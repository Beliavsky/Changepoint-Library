program xsim_var_css_amoc
! simulate a univariate Gaussian series and estimate one changepoint in variance by CSS AMOC
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_amoc_var_css_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 1 ! RNG seed used for the simulated series
integer, parameter :: n = 200 ! number of observations in the simulated series
integer, parameter :: minseglen = 2 ! minimum segment length used by CSS AMOC
integer, parameter :: nreg = 2 ! number of variance regimes in the simulated data
integer, parameter :: true_cps(nreg) = [0, 100] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: means(nreg) = [0.0_dp, 0.0_dp] ! segment means
real(kind=dp), parameter :: sds(nreg) = [1.0_dp, 10.0_dp] ! segment standard deviations
real(kind=dp), parameter :: pen_value = log(2.0_dp * log(real(n, dp))) ! manual CSS penalty used in the R example

real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: test_stat, penalty_value
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, true_cps, means, sds, x)
call solve_amoc_var_css_1d(x, cpt, test_stat, penalty_value, penalty="Manual", minseglen=minseglen, pen_value=pen_value)

print *, "n                  =", n
print *, "penalty            =", "Manual"
print *, "pen.value          =", pen_value
print *, "minseglen          =", minseglen
print *, "true changepoint   =", true_cps(2)
print *, "estimated          =", cpt
print *, "test statistic     =", test_stat

deallocate(x)
call print_wall_time(t_start)

end program xsim_var_css_amoc
