program xsim_mean_pelt_file
! read a univariate series from a text file and estimate changepoints by PELT with l2 cost
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_pelt_mean_shift_1d
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: pen = 10.0_dp
integer, parameter :: jump = 5
character(len=256) :: data_file
integer :: n, min_seg_len, n_bkps_found
integer, allocatable :: bkps(:)
real(kind=dp), allocatable :: z(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xmean_shift_data.txt", data_file)
call read_series_file(data_file, z)

n = size(z)
min_seg_len = max(5, n / 20)

call solve_pelt_mean_shift_1d(z, pen, bkps, n_bkps_found, min_seg_len=min_seg_len, jump=jump)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "pen         =", pen
print *, "min_seg_len =", min_seg_len
print *, "jump        =", jump
print *, "n_bkps      =", n_bkps_found
if (size(bkps) > 0) print *, "estimated   =", bkps

deallocate(z, bkps)
call print_wall_time(t_start)

end program xsim_mean_pelt_file
