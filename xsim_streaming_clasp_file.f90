program xsim_streaming_clasp_file
! read a univariate series from a text file and segment it with a recompute-based streaming ClaSP analog
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use claspy_mod, only: solve_streaming_class_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_timepoints = 200 ! streaming sliding-window length
integer, parameter :: n_warmup = 100 ! warmup observations before streaming decisions
integer, parameter :: window_size = 12 ! subsequence window size used by ClaSS
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in ClaSS
integer, parameter :: jump = 5 ! streaming step between candidate evaluations
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSS
real(kind=dp), parameter :: threshold = 0.50_dp ! score-threshold validation cutoff

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: cps(:)
integer :: last_cp
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xstreaming_clasp_data.txt", data_file)
call read_series_file(data_file, x)
call solve_streaming_class_1d(x, n_timepoints, n_warmup, window_size, jump, cps, last_cp, &
    k_neighbours, excl_radius, threshold)

print *, "file         =", trim(data_file)
print *, "n            =", size(x)
print *, "n_timepoints =", n_timepoints
print *, "n_warmup     =", n_warmup
print *, "window_size  =", window_size
print *, "k_neighbours =", k_neighbours
print *, "jump         =", jump
print *, "excl_radius  =", excl_radius
print *, "threshold    =", threshold
if (size(cps) > 0) write(*,'(a,*(i0,1x))') " estimated    = ", cps
print *, "last_cp      =", last_cp

deallocate(x, cps)
call print_wall_time(t_start)
end program xsim_streaming_clasp_file
