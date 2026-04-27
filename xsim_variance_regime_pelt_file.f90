program xsim_variance_regime_pelt_file
! read squared returns from a text file and estimate variance-regime changepoints by RBF PELT
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_pelt_rbf_mv
use util_mod, only: print_wall_time, read_matrix
implicit none

real(kind=dp), parameter :: pen = 5.0_dp
integer, parameter :: jump = 5
character(len=256) :: data_file
integer :: n, p, min_seg_len, n_bkps_found
integer, allocatable :: bkps(:)
real(kind=dp), allocatable :: x(:,:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file(data_file)
call read_matrix(trim(data_file), x)

n = size(x, 1)
p = size(x, 2)
min_seg_len = max(5, n / 20)

call solve_pelt_rbf_mv(x, pen, bkps, n_bkps_found, min_seg_len=min_seg_len, jump=jump)

print *, "file        =", trim(data_file)
print *, "n           =", n
print *, "p           =", p
print *, "pen         =", pen
print *, "min_seg_len =", min_seg_len
print *, "jump        =", jump
print *, "n_bkps      =", n_bkps_found
if (size(bkps) > 0) print *, "estimated   =", bkps

deallocate(x, bkps)
call print_wall_time(t_start)

contains

subroutine get_data_file(filename)
character(len=*), intent(out) :: filename
integer :: status

filename = "xvariance_regime_data.txt"
call get_command_argument(1, filename, status=status)
if (status /= 0 .or. len_trim(filename) == 0) filename = "xvariance_regime_data.txt"

end subroutine get_data_file

end program xsim_variance_regime_pelt_file
