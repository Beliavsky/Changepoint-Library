program xsim_clasp_ensemble_file
! read a univariate series from a text file and score a deterministic ClaSPEnsemble-style split
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use claspy_mod, only: solve_clasp_ensemble_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_estimators = 5 ! number of temporal constraints in the deterministic ensemble
integer, parameter :: window_size = 12 ! subsequence window size used by each ClaSP member
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in each ClaSP member
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSP

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: score
integer :: cp
logical :: found
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xclasp_ensemble_data.txt", data_file)
call read_series_file(data_file, x)
call solve_clasp_ensemble_1d(x, window_size, cp, score, found, n_estimators, k_neighbours, excl_radius)

print *, "file              =", trim(data_file)
print *, "n                 =", size(x)
print *, "n_estimators      =", n_estimators
print *, "window_size       =", window_size
print *, "k_neighbours      =", k_neighbours
print *, "excl_radius       =", excl_radius
if (found) then
    print *, "estimated         =", cp
else
    print *, "estimated         = none"
end if
print *, "max profile score =", score

deallocate(x)
call print_wall_time(t_start)
end program xsim_clasp_ensemble_file
