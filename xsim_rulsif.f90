program xsim_rulsif
! simulate a univariate series and score it with a deterministic Gaussian-kernel RuLSIF detector
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoynt_mod, only: solve_rulsif_gaussian_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 400 ! series length
integer, parameter :: true_cp = 200 ! true change location
integer, parameter :: window_length = 10 ! Hankel subsequence length
integer, parameter :: n_windows = 50 ! number of reference/test subsequences
integer, parameter :: lag = 50 ! gap between reference and test windows
integer, parameter :: scoring_step = 1 ! score evaluation stride

integer, parameter :: regime_starts(2) = [0, true_cp] ! starting indices of the two simulated regimes
real(kind=dp), parameter :: means(2) = [0.0_dp, 3.0_dp] ! regime means
real(kind=dp), parameter :: sds(2) = [1.0_dp, 1.0_dp] ! regime standard deviations
real(kind=dp), parameter :: alpha = 0.1_dp ! relative density-ratio smoothing parameter
real(kind=dp), parameter :: sigma = 2.0_dp ! Gaussian kernel width
real(kind=dp), parameter :: lambda_reg = 1.0e-2_dp ! ridge regularization

real(kind=dp) :: x(n), score(n), best_score
integer :: cp
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(101)
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)
call solve_rulsif_gaussian_1d(x, window_length, n_windows, lag, alpha, sigma, lambda_reg, score, cp, best_score, scoring_step)

print *, "n                  =", n
print *, "true changepoint   =", true_cp
print *, "window_length      =", window_length
print *, "n_windows          =", n_windows
print *, "lag                =", lag
print *, "alpha              =", alpha
print *, "sigma              =", sigma
print *, "lambda             =", lambda_reg
print *, "scoring_step       =", scoring_step
print *, "estimated          =", cp
print *, "max raw score      =", best_score

call print_wall_time(t_start)
end program xsim_rulsif
