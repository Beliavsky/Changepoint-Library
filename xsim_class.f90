program xsim_class
! simulate a univariate series and score a single ClaSS split
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use claspy_mod, only: solve_class_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 240 ! series length
integer, parameter :: true_cp = 120 ! true change location in the simulated series
integer, parameter :: window_size = 12 ! subsequence window size used by ClaSS
integer, parameter :: k_neighbours = 3 ! number of nearest neighbours in ClaSS
integer, parameter :: excl_radius = 5 ! minimum-segment multiplier used by ClaSS
real(kind=dp), parameter :: threshold = 0.50_dp ! score-threshold validation cutoff

integer, parameter :: regime_starts(2) = [0, true_cp] ! starting indices of the two simulated regimes
real(kind=dp), parameter :: means(2) = [0.0_dp, 4.0_dp] ! regime means
real(kind=dp), parameter :: sds(2) = [1.0_dp, 1.0_dp] ! regime standard deviations

real(kind=dp) :: x(n), score
integer :: cp
logical :: found
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(101)
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)
call solve_class_1d(x, window_size, cp, score, found, k_neighbours, excl_radius, threshold)

print *, "n                 =", n
print *, "window_size       =", window_size
print *, "k_neighbours      =", k_neighbours
print *, "excl_radius       =", excl_radius
print *, "threshold         =", threshold
print *, "true changepoint  =", true_cp
if (found) then
    print *, "estimated         =", cp
else
    print *, "estimated         = none"
end if
print *, "score             =", score

call print_wall_time(t_start)
end program xsim_class
