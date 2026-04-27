program xsim_roerich_rulsif
! simulate a univariate series and score a deterministic roerich RuLSIF-style detector
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use roerich_mod, only: solve_rulsif_linear_pesym_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 400 ! series length
integer, parameter :: true_cp = 200 ! true change location in the simulated series
integer, parameter :: window_size = 40 ! reference/test window length
integer, parameter :: periods = 1 ! autoregressive embedding length
integer, parameter :: step = 1 ! score evaluation stride
real(kind=dp), parameter :: alpha = 0.05_dp ! RuLSIF alpha parameter
real(kind=dp), parameter :: l2 = 1.0e-3_dp ! ridge regularization for the linear RuLSIF fit

integer, parameter :: regime_starts(2) = [0, true_cp] ! starting indices of the two simulated regimes
real(kind=dp), parameter :: means(2) = [0.0_dp, 3.0_dp] ! regime means
real(kind=dp), parameter :: sds(2) = [1.0_dp, 1.0_dp] ! regime standard deviations

real(kind=dp) :: x(n), best_score
integer :: cp
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(101)
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)
call solve_rulsif_linear_pesym_1d(x, window_size, cp, best_score, periods, step, 2357, alpha, l2)

print *, "n                  =", n
print *, "true changepoint   =", true_cp
print *, "window_size        =", window_size
print *, "periods            =", periods
print *, "step               =", step
print *, "base_regressor     = linear_rulsif"
print *, "metric             = pesym"
print *, "alpha              =", alpha
print *, "l2                 =", l2
print *, "estimated          =", cp
print *, "max raw score      =", best_score

call print_wall_time(t_start)
end program xsim_roerich_rulsif
