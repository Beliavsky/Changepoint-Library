program xsim_bocpd_beta_bernoulli
! simulate a piecewise-Bernoulli binary series and estimate changepoints by BOCPD
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_bernoulli_1d
use changepoint_mod, only: solve_bocpd_beta_bernoulli_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 51 ! RNG seed used for the simulated binary series
integer, parameter :: n = 600 ! number of observations in the simulated series
real(kind=dp), parameter :: lam = 120.0_dp ! characteristic run length in the geometric hazard
integer, parameter :: nreg = 3 ! number of piecewise-constant regimes
integer, parameter :: true_cps(nreg) = [0, 200, 400] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: probs(nreg) = [0.15_dp, 0.80_dp, 0.35_dp] ! segment Bernoulli success probabilities

logical :: x(n)
integer, allocatable :: map_cps(:), runlen_argmax(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(seed)
call simulate_piecewise_bernoulli_1d(n, true_cps, probs, x)
call solve_bocpd_beta_bernoulli_1d(x, lam, map_cps, runlen_argmax)

print *, "n           =", n
print *, "lam         =", lam
print *, "true cps    =", true_cps
if (size(map_cps) > 0) print *, "estimated   =", map_cps

deallocate(map_cps, runlen_argmax)
call print_wall_time(t_start)

end program xsim_bocpd_beta_bernoulli
