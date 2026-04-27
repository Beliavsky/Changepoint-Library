program xsim_window_file
! read a univariate series from a text file and estimate changepoints by the sliding-window l2 method
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_window_mean_shift_1d
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_bkps = 3
character(len=256) :: data_file
integer :: n, width
integer, allocatable :: bkps(:)
real(kind=dp), allocatable :: z(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xwindow_data.txt", data_file)
call read_series_file(data_file, z)

n = size(z)
width = max(20, n / 10)
call solve_window_mean_shift_1d(z, width, n_bkps, bkps)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "width       =", width
print *, "n_bkps      =", n_bkps
if (size(bkps) > 0) print *, "estimated   =", bkps

deallocate(z, bkps)
call print_wall_time(t_start)

end program xsim_window_file
