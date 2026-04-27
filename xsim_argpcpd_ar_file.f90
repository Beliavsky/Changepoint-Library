program xsim_argpcpd_ar_file
! read a univariate series from a text file and estimate changepoints by an AR-based online detector
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_bocpd_ar_1d
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: logistic_hazard_h = -5.0_dp ! hazard scale in logit units
real(kind=dp), parameter :: logistic_hazard_a = 1.0_dp ! logistic hazard slope parameter
real(kind=dp), parameter :: logistic_hazard_b = 1.0_dp ! logistic hazard offset parameter
real(kind=dp), parameter :: epsilon = 1.0e-10_dp ! tail-mass threshold for truncating active run lengths
integer, parameter :: max_lag = 12 ! maximum autoregressive lag used in each local AR fit

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: map_cps(:), runlen_argmax(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xargpcpd_data.txt", data_file)
call read_series_file(data_file, x)

call solve_bocpd_ar_1d(x, logistic_hazard_h, map_cps, runlen_argmax, max_lag, logistic_hazard_a, logistic_hazard_b, epsilon)

print *, "file        =", trim(data_file)
print *, "n           =", size(x)
print *, "max_lag     =", max_lag
print *, "hazard_h    =", logistic_hazard_h
if (size(map_cps) > 0) print *, "estimated   =", map_cps

deallocate(x, map_cps, runlen_argmax)
call print_wall_time(t_start)

end program xsim_argpcpd_ar_file
