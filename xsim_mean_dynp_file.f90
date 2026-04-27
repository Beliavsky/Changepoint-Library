program xsim_mean_dynp_file
! read a univariate series from a text file and estimate changepoints by dynamic programming
use kind_mod, only: dp, long_int
use changepoint_mod, only: mean_shift_cost_matrix, solve_changepoints, segment_ends
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_bkps = 3
integer, parameter :: max_m = n_bkps + 1
character(len=256) :: data_file
integer :: n, min_seg_len
integer, allocatable :: parent(:,:)
integer :: seg_ends(max_m)
real(kind=dp), allocatable :: z(:), cost(:,:), dp_table(:,:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xdynp_data.txt", data_file)
call read_series_file(data_file, z)

n = size(z)
min_seg_len = max(5, n / 20)

allocate(cost(n, n), dp_table(n, max_m), parent(n, max_m))
cost = mean_shift_cost_matrix(z, min_seg_len=min_seg_len)
call solve_changepoints(max_m, cost, dp_table, parent)
seg_ends = segment_ends(parent, max_m)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "n_bkps      =", n_bkps
print *, "min_seg_len =", min_seg_len
print *, "estimated   =", seg_ends

deallocate(z, cost, dp_table, parent)
call print_wall_time(t_start)

end program xsim_mean_dynp_file
