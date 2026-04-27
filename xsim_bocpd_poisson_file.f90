program xsim_bocpd_poisson_file
! read a univariate count series from a text file and estimate changepoints by BOCPD
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_bocpd_poisson_gamma_1d
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: lam = 120.0_dp ! characteristic run length in the geometric hazard
character(len=256) :: data_file
real(kind=dp), allocatable :: x_raw(:)
integer, allocatable :: x(:), map_cps(:), runlen_argmax(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xbocpd_poisson_data.txt", data_file)
call read_series_file(data_file, x_raw)

allocate(x(size(x_raw)))
x = nint(x_raw)

call solve_bocpd_poisson_gamma_1d(x, lam, map_cps, runlen_argmax)

print *, "file        =", trim(data_file)
print *, "n           =", size(x)
print *, "lam         =", lam
if (size(map_cps) > 0) print *, "estimated   =", map_cps

deallocate(x_raw, x, map_cps, runlen_argmax)
call print_wall_time(t_start)

end program xsim_bocpd_poisson_file
