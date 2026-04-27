program xsim_seeded_binseg
! simulate a univariate Gaussian series and estimate changepoints by seeded binary segmentation
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_seeded_binseg_cusum_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 61 ! RNG seed used for the simulated series
integer, parameter :: n = 20000 ! number of observations in the simulated series
integer, parameter :: nreg = 5 ! number of piecewise-constant regimes
integer, parameter :: true_cps(nreg) = [0, 4000, 8000, 12000, 16000] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: means(nreg) = [0.0_dp, 4.0_dp, -2.0_dp, 3.0_dp, -1.5_dp] ! segment means
real(kind=dp), parameter :: sds(nreg) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp] ! segment standard deviations
integer, parameter :: max_interval_length = n / 8 ! maximum seeded interval length
real(kind=dp), parameter :: growth_factor = 1.5_dp ! seeded interval growth factor

real(kind=dp), allocatable :: x(:)
integer, allocatable :: bkps(:)
integer(kind=long_int) :: t_start
real(kind=dp) :: penalty

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, true_cps, means, sds, x)
penalty = 2.0_dp * log(real(n, dp))
call solve_seeded_binseg_cusum_1d(x, penalty, bkps, max_interval_length, growth_factor, "greedy")

print *, "n           =", n
print *, "penalty     =", penalty
print *, "true cps    =", true_cps
if (size(bkps) > 0) print *, "estimated   =", bkps

deallocate(x, bkps)
call print_wall_time(t_start)

end program xsim_seeded_binseg
