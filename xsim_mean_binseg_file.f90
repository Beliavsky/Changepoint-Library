program xsim_mean_binseg_file
! read a univariate series from a text file and estimate changepoints by binary segmentation
use kind_mod, only: dp, long_int
use changepoint_binseg_mod, only: solve_binseg_mean_shift
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

integer, parameter :: max_cp = 4
character(len=256) :: data_file
integer :: n, ncp_found, min_seg_len
integer, allocatable :: cp(:)
real(kind=dp), allocatable :: z(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)

call get_data_file_arg("xbinseg_data.txt", data_file)
call read_series_file(data_file, z)

n = size(z)
min_seg_len = max(5, n / 25)

allocate(cp(max_cp))
call solve_binseg_mean_shift(z, max_cp, cp, ncp_found, min_seg_len, 0.0_dp)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "max_cp      =", max_cp
print *, "min_seg_len =", min_seg_len
print *, "ncp_found   =", ncp_found
if (ncp_found > 0) print *, "estimated   =", cp(1:ncp_found)

deallocate(z, cp)
call print_wall_time(t_start)

end program xsim_mean_binseg_file
