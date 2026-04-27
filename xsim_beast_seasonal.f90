program xsim_beast_seasonal
use kind_mod, only: dp, long_int
use util_mod, only: print_wall_time
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_linear_harmonic_1d
use rbeast_mod, only: solve_beast_harmonic_1d, top_cp_probabilities
implicit none

integer, parameter :: n = 240              ! number of observations
integer, parameter :: max_cp = 4           ! maximum number of trend changepoints
integer, parameter :: min_seg_len = 30     ! minimum segment length
integer, parameter :: nreg = 3             ! number of linear regimes
integer, parameter :: ntop = 5             ! number of cp probabilities to print
integer, parameter :: seed = 20260422      ! RNG seed for reproducibility

integer :: regime_starts(nreg), true_cps(nreg - 1), best_cps(max_cp), ncp_best
integer :: top_idx(ntop), i
integer(kind=long_int) :: start_count
real(kind=dp) :: intercepts(nreg), slopes(nreg), sigma, period, amp_sin, amp_cos
real(kind=dp) :: y(n), cp_prob(n), fitted_mean(n), seasonal_fit(n), deseasonalized(n), model_weights(max_cp + 1), top_prob(ntop)
real(kind=dp) :: beta_season(2)

regime_starts = [0, 80, 160]
true_cps = [80, 160]
intercepts = [1.0_dp, 4.6_dp, -1.2_dp]
slopes = [0.03_dp, -0.02_dp, 0.015_dp]
sigma = 0.35_dp
period = 12.0_dp
amp_sin = 1.4_dp
amp_cos = -0.8_dp

call seed_rng_fixed(seed)
call simulate_piecewise_linear_harmonic_1d(n, regime_starts, intercepts, slopes, sigma, period, amp_sin, amp_cos, y)

call system_clock(start_count)
call solve_beast_harmonic_1d(y, period, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, fitted_mean, seasonal_fit, &
                             deseasonalized, model_weights, beta_season)
call top_cp_probabilities(cp_prob, top_idx, top_prob)

print '(a,i0)', 'n = ', n
print '(a)', 'season = harmonic'
print '(a,f0.1)', 'period = ', period
print '(a,i0)', 'max_cp = ', max_cp
print '(a,i0)', 'min_seg_len = ', min_seg_len
print '(a,2(1x,i0))', 'true changepoints =', true_cps
if (ncp_best > 0) then
    write(*, '(a)', advance='no') 'estimated changepoints ='
    print '(100(1x,i0))', best_cps(1:ncp_best)
else
    print '(a)', 'estimated changepoints ='
end if
print '(a,2(1x,f0.6))', 'seasonal beta =', beta_season
write(*, '(a)', advance='no') 'model weights ='
print '(100(1x,f0.6))', model_weights
write(*, '(a)', advance='no') 'top cp probabilities ='
do i = 1, ntop
    if (top_idx(i) <= 0) exit
    write(*, '(" (",i0,",",f0.6,")")', advance='no') top_idx(i), top_prob(i)
end do
print *
call print_wall_time(start_count)

end program xsim_beast_seasonal
