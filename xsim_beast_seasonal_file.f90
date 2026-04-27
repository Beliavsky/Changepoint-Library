program xsim_beast_seasonal_file
use kind_mod, only: dp, long_int
use util_mod, only: print_wall_time
use compare_io_mod, only: get_data_file_arg, read_series_file
use rbeast_mod, only: solve_beast_harmonic_1d, top_cp_probabilities
implicit none

integer, parameter :: max_cp = 4           ! maximum number of trend changepoints
integer, parameter :: min_seg_len = 30     ! minimum segment length
integer, parameter :: ntop = 5             ! number of cp probabilities to print

character(len=256) :: data_file
integer :: n, ncp_best, top_idx(ntop), i
integer :: best_cps(max_cp)
integer(kind=long_int) :: start_count
real(kind=dp), parameter :: period = 12.0_dp ! seasonal period
real(kind=dp), allocatable :: y(:), cp_prob(:), fitted_mean(:), seasonal_fit(:), deseasonalized(:)
real(kind=dp) :: model_weights(max_cp + 1), top_prob(ntop), beta_season(2)

call get_data_file_arg('xbeast_seasonal_data.txt', data_file)
call read_series_file(trim(data_file), y)
n = size(y)
allocate(cp_prob(n), fitted_mean(n), seasonal_fit(n), deseasonalized(n))

call system_clock(start_count)
call solve_beast_harmonic_1d(y, period, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, fitted_mean, seasonal_fit, &
                             deseasonalized, model_weights, beta_season)
call top_cp_probabilities(cp_prob, top_idx, top_prob)

print '(a,a)', 'file = ', trim(data_file)
print '(a,i0)', 'n = ', n
print '(a)', 'season = harmonic'
print '(a,f0.1)', 'period = ', period
print '(a,i0)', 'max_cp = ', max_cp
print '(a,i0)', 'min_seg_len = ', min_seg_len
if (ncp_best > 0) then
    write(*, '(a)', advance='no') 'estimated changepoints ='
    print '(100(1x,i0))', best_cps(1:ncp_best)
else
    print '(a)', 'estimated changepoints ='
end if
print '(a,2(1x,f0.6))', 'seasonal beta =', beta_season
write(*, '(a)', advance='no') 'model weights ='
print '(100(1x,f0.6))', model_weights
write(*, '(a)', advance='no') 'top cp probabilities ='
do i = 1, ntop
    if (top_idx(i) <= 0) exit
    write(*, '(" (",i0,",",f0.6,")")', advance='no') top_idx(i), top_prob(i)
end do
print *
call print_wall_time(start_count)

deallocate(y, cp_prob, fitted_mean, seasonal_fit, deseasonalized)
end program xsim_beast_seasonal_file
