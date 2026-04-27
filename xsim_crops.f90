program xsim_crops
! simulate a univariate Gaussian series and estimate changepoints by CROPS
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_crops_mean_shift_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 71 ! RNG seed used for the simulated series
integer, parameter :: n = 4000 ! number of observations in the simulated series
integer, parameter :: nreg = 5 ! number of piecewise-constant regimes
integer, parameter :: true_cps(nreg) = [0, 800, 1600, 2400, 3200] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: means(nreg) = [0.0_dp, 4.0_dp, -2.0_dp, 3.5_dp, -1.5_dp] ! segment means
real(kind=dp), parameter :: sds(nreg) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp] ! segment standard deviations
real(kind=dp), parameter :: min_penalty = 5.0_dp ! lower end of the CROPS penalty interval
real(kind=dp), parameter :: max_penalty = 60.0_dp ! upper end of the CROPS penalty interval
integer, parameter :: min_segment_length = 10 ! minimum allowed segment length in PELT

real(kind=dp), allocatable :: x(:)
integer, allocatable :: bkps(:)
integer(kind=long_int) :: t_start
real(kind=dp) :: best_penalty
integer :: n_path_solutions

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, true_cps, means, sds, x)
call solve_crops_mean_shift_1d(x, min_penalty, max_penalty, bkps, best_penalty, n_path_solutions, min_segment_length)

print *, "n              =", n
print *, "true cps       =", true_cps
print *, "optimal penalty=", best_penalty
print *, "path solutions =", n_path_solutions
if (size(bkps) > 0) print *, "estimated      =", bkps

deallocate(x, bkps)
call print_wall_time(t_start)

end program xsim_crops
