program xsim_mean_amoc_file
! read a univariate series from a text file and estimate one changepoint in the mean by AMOC
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoint_mod, only: solve_amoc_mean_norm_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: minseglen = 1 ! minimum segment length used by AMOC

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: conf_value, null_like, alt_like
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xamoc_data.txt", data_file)
call read_series_file(data_file, x)
call solve_amoc_mean_norm_1d(x, cpt, conf_value, null_like, alt_like, penalty="SIC", minseglen=minseglen)

print *, "file               =", trim(data_file)
print *, "n                  =", size(x)
print *, "penalty            =", "SIC"
print *, "minseglen          =", minseglen
print *, "estimated          =", cpt
print *, "conf.value         =", conf_value

deallocate(x)
call print_wall_time(t_start)

end program xsim_mean_amoc_file
