program xsim_mean_bottomup_file
! read a univariate series from a text file and estimate changepoints by bottom-up segmentation
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_bottomup_mean_shift_1d, refine_bkps_mean_shift_1d
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_bkps = 4
character(len=256) :: data_file
integer :: n, min_seg_len, jump
integer, allocatable :: bkps_initial(:), bkps(:)
real(kind=dp), allocatable :: z(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xbottomup_data.txt", data_file)
call read_series_file(data_file, z)

n = size(z)
min_seg_len = max(5, n / 25)
jump = 5

call solve_bottomup_mean_shift_1d(z, n_bkps, bkps_initial, min_seg_len, jump)
call refine_bkps_mean_shift_1d(z, bkps_initial, bkps, min_seg_len)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "n_bkps      =", n_bkps
print *, "min_seg_len =", min_seg_len
print *, "jump        =", jump
if (size(bkps_initial) > 0) print *, "initial     =", bkps_initial
if (size(bkps) > 0) print *, "estimated   =", bkps

deallocate(z, bkps_initial, bkps)
call print_wall_time(t_start)

end program xsim_mean_bottomup_file
