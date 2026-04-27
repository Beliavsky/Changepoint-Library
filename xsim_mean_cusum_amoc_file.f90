program xsim_mean_cusum_amoc_file
! read a univariate series from a text file and estimate one changepoint in the mean by AMOC CUSUM
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoint_mod, only: solve_amoc_mean_cusum_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: minseglen = 1 ! minimum segment length used by AMOC CUSUM
real(kind=dp), parameter :: penalty_value_in = 0.5_dp ! manual penalty passed to changepoint::cpt.mean

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: test_stat, penalty_value
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xcusum_data.txt", data_file)
call read_series_file(data_file, x)
call solve_amoc_mean_cusum_1d(x, cpt, test_stat, penalty_value, penalty="Manual", minseglen=minseglen, pen_value=penalty_value_in)

print *, "file               =", trim(data_file)
print *, "n                  =", size(x)
print *, "penalty            =", "Manual"
print *, "pen.value          =", penalty_value_in
print *, "minseglen          =", minseglen
print *, "estimated          =", cpt
print *, "test statistic     =", test_stat
print *, "effective penalty  =", penalty_value

deallocate(x)
call print_wall_time(t_start)

end program xsim_mean_cusum_amoc_file
