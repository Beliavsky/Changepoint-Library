program xsim_roerich_classifier_file
! read a univariate series from a text file and score a deterministic roerich classifier-style detector
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use roerich_mod, only: solve_cpdc_qda_klsym_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: window_size = 40 ! reference/test window length
integer, parameter :: periods = 1 ! autoregressive embedding length
integer, parameter :: step = 1 ! score evaluation stride

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: best_score
integer :: cp
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xroerich_classifier_data.txt", data_file)
call read_series_file(data_file, x)
call solve_cpdc_qda_klsym_1d(x, window_size, cp, best_score, periods, step)

print *, "file               =", trim(data_file)
print *, "n                  =", size(x)
print *, "window_size        =", window_size
print *, "periods            =", periods
print *, "step               =", step
print *, "base_classifier    = qda"
print *, "metric             = klsym"
print *, "estimated          =", cp
print *, "max raw score      =", best_score

deallocate(x)
call print_wall_time(t_start)
end program xsim_roerich_classifier_file
