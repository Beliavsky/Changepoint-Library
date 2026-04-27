program xsim_clasp_alt_file
! read a second univariate series from a text file and segment it with fixed-K ClaSP binary segmentation
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use claspy_mod, only: solve_binary_clasp_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_segments = 5 ! target number of segments in the binary ClaSP segmentation
integer, parameter :: window_size = 10 ! subsequence window size used by ClaSP
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in ClaSP
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSP

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: cps(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xclasp_alt_data.txt", data_file)
call read_series_file(data_file, x)
call solve_binary_clasp_1d(x, window_size, n_segments, cps, k_neighbours, excl_radius)

print *, "file         =", trim(data_file)
print *, "n            =", size(x)
print *, "n_segments   =", n_segments
print *, "window_size  =", window_size
print *, "k_neighbours =", k_neighbours
print *, "excl_radius  =", excl_radius
if (size(cps) > 0) print *, "estimated    =", cps

deallocate(x, cps)
call print_wall_time(t_start)

end program xsim_clasp_alt_file
