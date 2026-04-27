program xsim_bocpd_changepoynt
! simulate a univariate Gaussian series and score it with changepoynt BOCPD
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoynt_mod, only: solve_bocpd_gaussian_mean_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 600 ! series length
integer, parameter :: run_length = 120 ! constant-hazard reciprocal and sliding-window fit length
integer, parameter :: regime_starts(3) = [0, 200, 400] ! starting indices of the simulated regimes
real(kind=dp), parameter :: means(3) = [0.0_dp, 3.5_dp, -1.5_dp] ! regime means
real(kind=dp), parameter :: sds(3) = [1.0_dp, 1.0_dp, 1.6_dp] ! regime standard deviations

real(kind=dp) :: x(n), score(n), best_score, prior_mean, prior_var, signal_var
integer :: cp_raw, cp_post
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(21)
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)
call solve_bocpd_gaussian_mean_1d(x, run_length, score, cp_raw, cp_post, best_score, prior_mean, prior_var, signal_var)

print *, "n                  =", n
print *, "run_length         =", run_length
print *, "raw argmax         =", cp_raw
print *, "post_warmup argmax =", cp_post
print *, "max score          =", best_score
print *, "prior_mean         =", prior_mean
print *, "prior_var          =", prior_var
print *, "signal_var         =", signal_var

call print_wall_time(t_start)
end program xsim_bocpd_changepoynt
