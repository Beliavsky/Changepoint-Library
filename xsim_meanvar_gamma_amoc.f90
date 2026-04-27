program xsim_meanvar_gamma_amoc
! simulate a univariate Gamma series and estimate one changepoint by AMOC
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_gamma_1d
use changepoint_mod, only: solve_amoc_meanvar_gamma_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 200 ! number of observations
integer, parameter :: true_cp = 100 ! true changepoint after this observation
integer, parameter :: minseglen = 2 ! minimum segment length used by Gamma AMOC
real(kind=dp), parameter :: shape = 2.0_dp ! fixed Gamma shape parameter

integer, dimension(2) :: regime_starts
real(kind=dp), dimension(2) :: means
real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: penalty_value, null_like, alt_like, alt_mbic
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(1)

regime_starts = [1, true_cp + 1]
means = [5.0_dp, 12.0_dp]
call simulate_piecewise_gamma_1d(n, regime_starts, shape, means, x)
call solve_amoc_meanvar_gamma_1d(x, shape, cpt, penalty_value, null_like, alt_like, penalty="MBIC", minseglen=minseglen)
if (cpt < n) then
    alt_mbic = alt_like + log(real(cpt, dp)) + log(real(n - cpt + 1, dp))
else
    alt_mbic = alt_like
end if

print *, "n                  =", n
print *, "shape              =", shape
print *, "penalty            =", "MBIC"
print *, "minseglen          =", minseglen
print *, "true_cp            =", true_cp
print *, "estimated          =", cpt
print *, "null               =", null_like
print *, "alt                =", alt_like
print *, "alt.mbic           =", alt_mbic
print *, "pen.value          =", penalty_value

deallocate(x)
call print_wall_time(t_start)

end program xsim_meanvar_gamma_amoc
