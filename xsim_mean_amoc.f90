program xsim_mean_amoc
! simulate a univariate Gaussian series and estimate one changepoint in the mean by AMOC
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_amoc_mean_norm_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 1 ! RNG seed used for the simulated series
integer, parameter :: n = 200 ! number of observations in the simulated series
integer, parameter :: minseglen = 1 ! minimum segment length used by AMOC
integer, parameter :: nreg = 2 ! number of mean regimes in the simulated data
integer, parameter :: true_cps(nreg) = [0, 100] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: means(nreg) = [0.0_dp, 5.0_dp] ! segment means
real(kind=dp), parameter :: sds(nreg) = [1.0_dp, 1.0_dp] ! segment standard deviations

real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: conf_value, null_like, alt_like
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, true_cps, means, sds, x)
call solve_amoc_mean_norm_1d(x, cpt, conf_value, null_like, alt_like, penalty="SIC", minseglen=minseglen)

print *, "n                  =", n
print *, "penalty            =", "SIC"
print *, "minseglen          =", minseglen
print *, "true changepoint   =", true_cps(2)
print *, "estimated          =", cpt
print *, "conf.value         =", conf_value

deallocate(x)
call print_wall_time(t_start)

end program xsim_mean_amoc
