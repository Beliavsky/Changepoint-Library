program xsim_beast123_file
use kind_mod, only: dp, long_int
use util_mod, only: print_wall_time
use compare_io_mod, only: get_data_file_arg, read_matrix_file
use rbeast_mod, only: solve_beast123_regular_mv, top_cp_probabilities
implicit none

integer, parameter :: max_cp = 4           ! maximum number of trend changepoints
integer, parameter :: min_seg_len = 30     ! minimum segment length
integer, parameter :: ntop = 5             ! number of cp probabilities to print

character(len=256) :: data_file
integer :: n, p, i, j, top_idx(ntop)
integer(kind=long_int) :: start_count
real(kind=dp), allocatable :: x(:, :), cp_prob(:, :), mean_cp_prob(:), fitted_mean(:, :), model_weights(:, :)
real(kind=dp) :: top_prob(ntop)
integer, allocatable :: best_cps(:, :), ncp_best(:)

call get_data_file_arg('xbeast123_data.txt', data_file)
call read_matrix_file(trim(data_file), x)
n = size(x, 1)
p = size(x, 2)
allocate(best_cps(max_cp, p), ncp_best(p), cp_prob(n, p), mean_cp_prob(n), fitted_mean(n, p), model_weights(max_cp + 1, p))

call system_clock(start_count)
call solve_beast123_regular_mv(x, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, mean_cp_prob, fitted_mean, model_weights)
call top_cp_probabilities(mean_cp_prob, top_idx, top_prob)

print '(a,a)', 'file = ', trim(data_file)
print '(a,i0)', 'n = ', n
print '(a,i0)', 'p = ', p
print '(a)', 'season = none'
print '(a)', 'beast123 regular interface'
print '(a,i0)', 'max_cp = ', max_cp
print '(a,i0)', 'min_seg_len = ', min_seg_len
do j = 1, p
    write(*, '(a,i0,a)', advance='no') 'series ', j, ' changepoints ='
    if (ncp_best(j) > 0) then
        print '(100(1x,i0))', best_cps(1:ncp_best(j), j)
    else
        print *
    end if
end do
write(*, '(a)', advance='no') 'top mean cp probabilities ='
do i = 1, ntop
    if (top_idx(i) <= 0) exit
    write(*, '(" (",i0,",",f0.6,")")', advance='no') top_idx(i), top_prob(i)
end do
print *
call print_wall_time(start_count)

deallocate(x, best_cps, ncp_best, cp_prob, mean_cp_prob, fitted_mean, model_weights)
end program xsim_beast123_file
