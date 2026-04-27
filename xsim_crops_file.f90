program xsim_crops_file
! read a univariate series from a text file and estimate changepoints by CROPS
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoint_mod, only: solve_crops_mean_shift_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: min_penalty = 5.0_dp ! lower end of the CROPS penalty interval
real(kind=dp), parameter :: max_penalty = 60.0_dp ! upper end of the CROPS penalty interval
integer, parameter :: min_segment_length = 10 ! minimum allowed segment length in PELT

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: bkps(:)
integer(kind=long_int) :: t_start
real(kind=dp) :: best_penalty
integer :: n_path_solutions

call system_clock(t_start)
call get_data_file_arg("xcrops_data.txt", data_file)
call read_series_file(data_file, x)
call solve_crops_mean_shift_1d(x, min_penalty, max_penalty, bkps, best_penalty, n_path_solutions, min_segment_length)

print *, "file           =", trim(data_file)
print *, "n              =", size(x)
print *, "optimal penalty=", best_penalty
print *, "path solutions =", n_path_solutions
if (size(bkps) > 0) print *, "estimated      =", bkps

deallocate(x, bkps)
call print_wall_time(t_start)

end program xsim_crops_file
