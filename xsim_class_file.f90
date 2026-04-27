program xsim_class_file
! read a univariate series from a text file and score a single ClaSS split
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use claspy_mod, only: solve_class_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: window_size = 12 ! subsequence window size used by ClaSS
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in ClaSS
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSS
real(kind=dp), parameter :: threshold = 0.50_dp ! score-threshold validation cutoff

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: score
integer :: cp
logical :: found
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xclass_data.txt", data_file)
call read_series_file(data_file, x)
call solve_class_1d(x, window_size, cp, score, found, k_neighbours, excl_radius, threshold)

print *, "file         =", trim(data_file)
print *, "n            =", size(x)
print *, "window_size  =", window_size
print *, "k_neighbours =", k_neighbours
print *, "excl_radius  =", excl_radius
print *, "threshold    =", threshold
if (found) then
    print *, "estimated    =", cp
else
    print *, "estimated    = none"
end if
print *, "score        =", score

deallocate(x)
call print_wall_time(t_start)
end program xsim_class_file
