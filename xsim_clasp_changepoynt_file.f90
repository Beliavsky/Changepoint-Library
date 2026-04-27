program xsim_clasp_changepoynt_file
! read a univariate series from a text file and segment it with changepoynt-style BinaryClaSPSegmentation defaults
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use claspy_mod, only: solve_binary_clasp_changepoynt_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_estimators = 10 ! number of ClaSP ensemble members
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in ClaSP
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSP
real(kind=dp), parameter :: threshold = 1.0e-15_dp ! significance-test threshold used by changepoynt-style CLASP

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: cps(:)
integer :: window_size_used, n_segments_used
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xclasp_changepoynt_data.txt", data_file)
call read_series_file(data_file, x)
call solve_binary_clasp_changepoynt_1d(x, cps, window_size_used, n_segments_used, &
    k_neighbours, excl_radius, n_estimators, threshold, 2357)

print *, "file             =", trim(data_file)
print *, "n                =", size(x)
print *, "n_segments       =", "learn"
print *, "window_size      =", "suss"
print *, "n_estimators     =", n_estimators
print *, "k_neighbours     =", k_neighbours
print *, "excl_radius      =", excl_radius
print *, "validation       =", "significance_test"
print *, "threshold        =", threshold
if (size(cps) > 0) print *, "estimated        =", cps
print *, "window_size used =", window_size_used
print *, "n_segments used  =", n_segments_used

deallocate(x, cps)
call print_wall_time(t_start)

end program xsim_clasp_changepoynt_file
