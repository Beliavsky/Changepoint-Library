program xsim_clasp_changepoynt
! simulate a piecewise-sinusoidal series and segment it with changepoynt-style BinaryClaSPSegmentation defaults
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_sine_1d
use claspy_mod, only: solve_binary_clasp_changepoynt_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 123 ! RNG seed used for the simulated series
integer, parameter :: n = 880 ! number of observations in the simulated series
integer, parameter :: n_estimators = 10 ! number of ClaSP ensemble members
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in ClaSP
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSP
integer, parameter :: nreg = 4 ! number of sinusoidal regimes in the simulated data
integer, parameter :: true_cps(nreg) = [0, 220, 440, 660] ! regime start indices in 0-based changepoint notation
real(kind=dp), parameter :: periods(nreg) = [12.0_dp, 30.0_dp, 18.0_dp, 42.0_dp] ! sinusoid periods for each regime
real(kind=dp), parameter :: sigma = 0.2_dp ! innovation standard deviation in each regime
real(kind=dp), parameter :: threshold = 1.0e-15_dp ! significance-test threshold used by changepoynt-style CLASP

real(kind=dp), allocatable :: x(:)
integer, allocatable :: cps(:)
integer :: window_size_used, n_segments_used
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(seed)
call simulate_piecewise_sine_1d(n, true_cps, periods, sigma, x)
call solve_binary_clasp_changepoynt_1d(x, cps, window_size_used, n_segments_used, &
    k_neighbours, excl_radius, n_estimators, threshold, 2357)

print *, "n                =", n
print *, "n_segments       =", "learn"
print *, "window_size      =", "suss"
print *, "n_estimators     =", n_estimators
print *, "k_neighbours     =", k_neighbours
print *, "excl_radius      =", excl_radius
print *, "validation       =", "significance_test"
print *, "threshold        =", threshold
print *, "true cps         =", true_cps
if (size(cps) > 0) print *, "estimated        =", cps
print *, "window_size used =", window_size_used
print *, "n_segments used  =", n_segments_used

deallocate(x, cps)
call print_wall_time(t_start)

end program xsim_clasp_changepoynt
