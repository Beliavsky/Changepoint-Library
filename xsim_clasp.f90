program xsim_clasp
! simulate a univariate Gaussian series and segment it with a fixed-K ClaSP binary segmentation
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use claspy_mod, only: solve_binary_clasp_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 101 ! RNG seed used for the simulated series
integer, parameter :: n = 500 ! number of observations in the simulated series
integer, parameter :: n_segments = 5 ! target number of segments in the binary ClaSP segmentation
integer, parameter :: window_size = 10 ! subsequence window size used by ClaSP
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in ClaSP
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSP
integer, parameter :: nreg = 5 ! number of piecewise-constant regimes in the simulated data
integer, parameter :: true_cps(nreg) = [0, 100, 200, 300, 400] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: means(nreg) = [0.0_dp, 4.0_dp, -3.0_dp, 3.5_dp, -2.5_dp] ! segment means
real(kind=dp), parameter :: sds(nreg) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp] ! segment standard deviations

real(kind=dp), allocatable :: x(:)
integer, allocatable :: cps(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, true_cps, means, sds, x)
call solve_binary_clasp_1d(x, window_size, n_segments, cps, k_neighbours, excl_radius)

print *, "n            =", n
print *, "n_segments   =", n_segments
print *, "window_size  =", window_size
print *, "k_neighbours =", k_neighbours
print *, "excl_radius  =", excl_radius
print *, "true cps     =", true_cps
if (size(cps) > 0) print *, "estimated    =", cps

deallocate(x, cps)
call print_wall_time(t_start)

end program xsim_clasp
