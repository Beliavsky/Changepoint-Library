program xsim_meanvar_exp_amoc
! simulate a univariate Exponential series and estimate one changepoint by AMOC
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_exponential_1d
use changepoint_mod, only: solve_amoc_meanvar_exp_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 200 ! number of observations
integer, parameter :: true_cp = 100 ! true changepoint after this observation
integer, parameter :: minseglen = 1 ! minimum segment length used by Exponential AMOC

integer, dimension(2) :: regime_starts
real(kind=dp), dimension(2) :: means
real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: p_value, null_like, alt_like, decision_penalty, alt_mbic
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(1)

regime_starts = [1, true_cp + 1]
means = [5.0_dp, 12.0_dp]
call simulate_piecewise_exponential_1d(n, regime_starts, means, x)
call solve_amoc_meanvar_exp_1d(x, cpt, p_value, null_like, alt_like, decision_penalty, penalty="MBIC", minseglen=minseglen)
if (cpt < n) then
    alt_mbic = alt_like + log(real(cpt, dp)) + log(real(n - cpt + 1, dp))
else
    alt_mbic = alt_like
end if

print *, "n                  =", n
print *, "penalty            =", "MBIC"
print *, "minseglen          =", minseglen
print *, "true_cp            =", true_cp
print *, "estimated          =", cpt
print *, "p.value            =", p_value
print *, "null               =", null_like
print *, "alt                =", alt_like
print *, "alt.mbic           =", alt_mbic
print *, "decision pen.value =", decision_penalty

deallocate(x)
call print_wall_time(t_start)

end program xsim_meanvar_exp_amoc
