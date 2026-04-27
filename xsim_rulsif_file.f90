program xsim_rulsif_file
! read a univariate series from a text file and score it with a deterministic Gaussian-kernel RuLSIF detector
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoynt_mod, only: solve_rulsif_gaussian_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: window_length = 10 ! Hankel subsequence length
integer, parameter :: n_windows = 50 ! number of reference/test subsequences
integer, parameter :: lag = 50 ! gap between reference and test windows
integer, parameter :: scoring_step = 1 ! score evaluation stride
real(kind=dp), parameter :: alpha = 0.1_dp ! relative density-ratio smoothing parameter
real(kind=dp), parameter :: sigma = 2.0_dp ! Gaussian kernel width
real(kind=dp), parameter :: lambda_reg = 1.0e-2_dp ! ridge regularization

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:), score(:)
real(kind=dp) :: best_score
integer :: cp
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xrulsif_data.txt", data_file)
call read_series_file(data_file, x)
allocate(score(size(x)))
call solve_rulsif_gaussian_1d(x, window_length, n_windows, lag, alpha, sigma, lambda_reg, score, cp, best_score, scoring_step)

print *, "file               =", trim(data_file)
print *, "n                  =", size(x)
print *, "window_length      =", window_length
print *, "n_windows          =", n_windows
print *, "lag                =", lag
print *, "alpha              =", alpha
print *, "sigma              =", sigma
print *, "lambda             =", lambda_reg
print *, "scoring_step       =", scoring_step
print *, "estimated          =", cp
print *, "max raw score      =", best_score

deallocate(score, x)
call print_wall_time(t_start)
end program xsim_rulsif_file
