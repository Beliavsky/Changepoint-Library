program xsim_argpcpd_gp
! simulate a univariate Gaussian series and estimate changepoints by ArgpCpd
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_argpcpd_gp_rust_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 31 ! RNG seed used for the simulated series
integer, parameter :: n = 600 ! number of observations in the simulated series
integer, parameter :: nreg = 3 ! number of piecewise-constant regimes
integer, parameter :: true_cps(nreg) = [0, 200, 400] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: means(nreg) = [0.0_dp, 4.0_dp, -3.0_dp] ! segment means
real(kind=dp), parameter :: sds(nreg) = [0.4_dp, 0.4_dp, 0.5_dp] ! segment standard deviations
real(kind=dp), parameter :: scale = 3.0_dp ! scale of the constant kernel factor
real(kind=dp), parameter :: length_scale = 10.0_dp ! length scale of the RBF kernel
real(kind=dp), parameter :: noise_level = 0.01_dp ! white-noise kernel level
integer, parameter :: max_lag = 12 ! maximum autoregressive lag in the lag-feature representation
real(kind=dp), parameter :: alpha0 = 2.0_dp ! scale-gamma alpha parameter
real(kind=dp), parameter :: beta0 = 1.0_dp ! scale-gamma beta parameter
real(kind=dp), parameter :: logistic_hazard_h = -5.0_dp ! hazard scale in logit units
real(kind=dp), parameter :: logistic_hazard_a = 1.0_dp ! logistic hazard slope parameter
real(kind=dp), parameter :: logistic_hazard_b = 1.0_dp ! logistic hazard offset parameter
real(kind=dp), parameter :: epsilon = 1.0e-10_dp ! tail-probability threshold for truncating run-length states

real(kind=dp) :: x(n)
integer, allocatable :: map_cps(:), runlen_argmax(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, true_cps, means, sds, x)
call solve_argpcpd_gp_rust_1d(x, map_cps, runlen_argmax, scale, length_scale, noise_level, max_lag, alpha0, beta0, &
                              logistic_hazard_h, logistic_hazard_a, logistic_hazard_b, epsilon)

print *, "n           =", n
print *, "true cps    =", true_cps
print *, "max_lag     =", max_lag
print *, "hazard_h    =", logistic_hazard_h
if (size(map_cps) > 0) print *, "estimated   =", map_cps

deallocate(map_cps, runlen_argmax)
call print_wall_time(t_start)

end program xsim_argpcpd_gp
