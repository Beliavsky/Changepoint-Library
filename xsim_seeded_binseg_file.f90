program xsim_seeded_binseg_file
! read a univariate series from a text file and estimate changepoints by seeded binary segmentation
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoint_mod, only: solve_seeded_binseg_cusum_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: growth_factor = 1.5_dp ! seeded interval growth factor
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: bkps(:)
integer(kind=long_int) :: t_start
real(kind=dp) :: penalty
integer :: max_interval_length

call system_clock(t_start)
call get_data_file_arg("xseeded_binseg_data.txt", data_file)
call read_series_file(data_file, x)

penalty = 2.0_dp * log(real(size(x), dp))
max_interval_length = max(200, size(x) / 8)
call solve_seeded_binseg_cusum_1d(x, penalty, bkps, max_interval_length, growth_factor, "greedy")

print *, "file        =", trim(data_file)
print *, "n           =", size(x)
print *, "penalty     =", penalty
print *, "max_int_len =", max_interval_length
if (size(bkps) > 0) print *, "estimated   =", bkps

deallocate(x, bkps)
call print_wall_time(t_start)

end program xsim_seeded_binseg_file
