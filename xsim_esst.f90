program xsim_esst
! simulate a frequency-change signal and score it with a deterministic exact-SVD ESST analog
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_sine_1d
use changepoynt_mod, only: solve_esst_exact_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n_per_segment = 320 ! length of each sinusoidal regime
integer, parameter :: n = 2 * n_per_segment ! total series length
integer, parameter :: true_cp = n_per_segment ! true changepoint location
integer, parameter :: window_length = 48 ! Hankel window length
integer, parameter :: n_windows = 24 ! number of past/future Hankel columns
integer, parameter :: lag = 24 ! comparison lag between past and future windows
integer, parameter :: rank = 2 ! retained singular-vector count
integer, parameter :: scoring_step = 1 ! score evaluation stride

integer, parameter :: regime_starts(2) = [0, true_cp] ! starting indices of the two sinusoidal regimes
real(kind=dp), parameter :: periods(2) = [48.0_dp, 14.0_dp] ! sinusoidal periods before and after the change
real(kind=dp), parameter :: sigma = 0.02_dp ! additive noise standard deviation

real(kind=dp) :: x(n), score(n), best_score
integer :: cp
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(1234)
call simulate_piecewise_sine_1d(n, regime_starts, periods, sigma, x)
call solve_esst_exact_1d(x, window_length, n_windows, lag, rank, score, cp, best_score, scoring_step)

print *, "n                  =", n
print *, "true changepoint   =", true_cp
print *, "window_length      =", window_length
print *, "n_windows          =", n_windows
print *, "lag                =", lag
print *, "rank               =", rank
print *, "scoring_step       =", scoring_step
print *, "method             = exact"
print *, "estimated          =", cp
print *, "max raw score      =", best_score

call print_wall_time(t_start)
end program xsim_esst
