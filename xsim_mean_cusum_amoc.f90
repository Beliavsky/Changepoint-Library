program xsim_mean_cusum_amoc
! simulate a univariate mean-shift series and estimate one changepoint by AMOC CUSUM
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_amoc_mean_cusum_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 200 ! number of observations
integer, parameter :: true_cp = 100 ! true changepoint after this observation
integer, parameter :: minseglen = 1 ! minimum segment length used by AMOC CUSUM
real(kind=dp), parameter :: penalty_value_in = 0.5_dp ! manual penalty passed to changepoint::cpt.mean

integer, dimension(2) :: regime_starts
real(kind=dp), dimension(2) :: means, sds
real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: test_stat, penalty_value
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(1)

regime_starts = [1, true_cp + 1]
means = [0.0_dp, 5.0_dp]
sds = [1.0_dp, 1.0_dp]
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)
call solve_amoc_mean_cusum_1d(x, cpt, test_stat, penalty_value, penalty="Manual", minseglen=minseglen, pen_value=penalty_value_in)

print *, "n                  =", n
print *, "penalty            =", "Manual"
print *, "pen.value          =", penalty_value_in
print *, "minseglen          =", minseglen
print *, "true_cp            =", true_cp
print *, "estimated          =", cpt
print *, "test statistic     =", test_stat
print *, "effective penalty  =", penalty_value

deallocate(x)
call print_wall_time(t_start)

end program xsim_mean_cusum_amoc
