program xsim_bocd_file
use kind_mod, only: dp, long_int
use util_mod, only: print_wall_time
use compare_io_mod, only: get_data_file_arg, read_series_file
use changepoint_mod, only: solve_bocd_student_t_1d
implicit none

real(kind=dp), parameter :: lam = 100.0_dp ! constant hazard scale

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: map_cps(:), map_cps_pkg(:), runlen_argmax(:)
integer(kind=long_int) :: start_count

call get_data_file_arg('xbocd_data.txt', data_file)
call read_series_file(trim(data_file), x)

call system_clock(start_count)
call solve_bocd_student_t_1d(x, lam, map_cps, runlen_argmax, mu0=0.0_dp, kappa0=1.0_dp, alpha0=1.0_dp, beta0=1.0_dp)
map_cps_pkg = pack(map_cps + 1, mask=map_cps > 0)

print '(a,a)', 'file = ', trim(data_file)
print '(a,i0)', 'n = ', size(x)
print '(a,f0.1)', 'lambda = ', lam
print '(a)', 'distribution = StudentT'
print '(a)', 'hazard = ConstantHazard'
write(*, '(a)', advance='no') 'estimated changepoints ='
if (size(map_cps_pkg) > 0) then
    print '(100(1x,i0))', map_cps_pkg
else
    print *
end if
print '(a,i0)', 'final MAP run length = ', runlen_argmax(size(runlen_argmax))
call print_wall_time(start_count)

deallocate(x, map_cps, map_cps_pkg, runlen_argmax)
end program xsim_bocd_file
