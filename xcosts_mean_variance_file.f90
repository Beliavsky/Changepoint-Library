program xcosts_mean_variance_file
! read a univariate series from a text file and compare several dynp costs
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_dynp_cost_1d
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_bkps = 3  ! number of changepoints to estimate
integer, parameter :: max_m = n_bkps + 1  ! number of fitted segments
integer, parameter :: min_seg_len = 30  ! minimum allowed segment length
character(len=256) :: data_file
real(kind=dp), allocatable :: z(:)
integer :: n
integer :: seg_ends_l1(max_m), seg_ends_l2(max_m), seg_ends_normal(max_m), seg_ends_rbf(max_m)
integer(kind=long_int) :: t_start, t0, t1, rate

call system_clock(t_start, rate)
call get_data_file_arg("xcosts_mean_variance_data.txt", data_file)
call read_series_file(data_file, z)
n = size(z)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "n_bkps      =", n_bkps
print *, "min_seg_len =", min_seg_len
print *

call system_clock(t0)
call solve_dynp_cost_1d(z, max_m, "l1", seg_ends_l1, min_seg_len)
call system_clock(t1)
print *, "model       = l1"
print *, "note        = robust location shifts"
print *, "estimated   =", seg_ends_l1
write(*,'(a,f8.3)') "elapsed seconds = ", real(t1 - t0, dp) / real(rate, dp)
print *

call system_clock(t0)
call solve_dynp_cost_1d(z, max_m, "l2", seg_ends_l2, min_seg_len)
call system_clock(t1)
print *, "model       = l2"
print *, "note        = mean shifts"
print *, "estimated   =", seg_ends_l2
write(*,'(a,f8.3)') "elapsed seconds = ", real(t1 - t0, dp) / real(rate, dp)
print *

call system_clock(t0)
call solve_dynp_cost_1d(z, max_m, "normal", seg_ends_normal, min_seg_len)
call system_clock(t1)
print *, "model       = normal"
print *, "note        = mean and variance shifts"
print *, "estimated   =", seg_ends_normal
write(*,'(a,f8.3)') "elapsed seconds = ", real(t1 - t0, dp) / real(rate, dp)
print *

call system_clock(t0)
call solve_dynp_cost_1d(z, max_m, "rbf", seg_ends_rbf, min_seg_len)
call system_clock(t1)
print *, "model       = rbf"
print *, "note        = general distribution shifts"
print *, "estimated   =", seg_ends_rbf
write(*,'(a,f8.3)') "elapsed seconds = ", real(t1 - t0, dp) / real(rate, dp)
print *

deallocate(z)
call print_wall_time(t_start)

end program xcosts_mean_variance_file
