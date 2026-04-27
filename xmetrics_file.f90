program xmetrics_file
! read a univariate series from a text file and evaluate changepoint estimates with metrics
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_bottomup_mean_shift_1d, solve_window_mean_shift_1d, &
                           solve_dynp_cost_1d
use changepoint_binseg_mod, only: solve_binseg_mean_shift
use compare_io_mod, only: read_series_file, read_true_bkps_file, get_data_file_arg
use changepoint_metrics_mod, only: print_metric_block
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_bkps = 3  ! number of changepoints to estimate
integer, parameter :: max_m = n_bkps + 1  ! number of fitted segments
integer, parameter :: min_seg_len = 40  ! minimum allowed segment length
integer, parameter :: margin = 40  ! tolerance used in precision/recall
character(len=256) :: data_file
real(kind=dp), allocatable :: z(:)
integer, allocatable :: true_bkps(:), bkps_bottomup_init(:), bkps_window(:), bkps_binseg(:), bkps_binseg_out(:)
integer :: n, width, ncp_found
integer :: bkps_dynp(max_m)
integer(kind=long_int) :: t_start, t0, t1, rate

call system_clock(t_start, rate)
call get_data_file_arg("xmetrics_data.txt", data_file)
call read_series_file(data_file, z)
call read_true_bkps_file(data_file, true_bkps)

n = size(z)
width = max(80, n / 10)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "true bkps   =", true_bkps
print *, "n_bkps      =", n_bkps
print *, "min_seg_len =", min_seg_len
print *, "margin      =", margin
print *

allocate(bkps_binseg(n_bkps))
call system_clock(t0)
call solve_binseg_mean_shift(z, n_bkps, bkps_binseg, ncp_found, min_seg_len, 0.0_dp)
call system_clock(t1)
allocate(bkps_binseg_out(ncp_found + 1))
if (ncp_found > 0) bkps_binseg_out(1:ncp_found) = bkps_binseg(1:ncp_found)
bkps_binseg_out(ncp_found + 1) = n
call print_metric_block("binseg", bkps_binseg_out, true_bkps, margin, t0, t1, rate)
print *

call system_clock(t0)
call solve_bottomup_mean_shift_1d(z, n_bkps, bkps_bottomup_init, min_seg_len, 5)
call system_clock(t1)
call print_metric_block("bottomup", bkps_bottomup_init, true_bkps, margin, t0, t1, rate)
print *

call system_clock(t0)
call solve_window_mean_shift_1d(z, width, n_bkps, bkps_window)
call system_clock(t1)
call print_metric_block("window", bkps_window, true_bkps, margin, t0, t1, rate)
print *

call system_clock(t0)
call solve_dynp_cost_1d(z, max_m, "l2", bkps_dynp, min_seg_len)
call system_clock(t1)
call print_metric_block("dynp", bkps_dynp, true_bkps, margin, t0, t1, rate)

deallocate(z, true_bkps, bkps_bottomup_init, bkps_window, bkps_binseg, bkps_binseg_out)
call print_wall_time(t_start)

end program xmetrics_file
