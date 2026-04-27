program xsim_beast_irreg
use kind_mod, only: dp, long_int
use util_mod, only: print_wall_time
use compare_sim_mod, only: seed_rng_fixed, simulate_irregular_time_grid, simulate_piecewise_linear_irregular_1d
use rbeast_mod, only: solve_beast_irreg_trend_only_1d, top_cp_probabilities
implicit none

integer, parameter :: n = 220              ! number of observations
integer, parameter :: max_cp = 4           ! maximum number of trend changepoints
integer, parameter :: min_seg_len = 28     ! minimum segment length
integer, parameter :: nreg = 3             ! number of linear regimes
integer, parameter :: ntop = 5             ! number of cp probabilities to print
integer, parameter :: seed_time = 20260423 ! RNG seed for irregular time grid
integer, parameter :: seed_y = 20260424    ! RNG seed for response values

integer :: regime_starts(nreg), true_cps(nreg - 1), best_cps(max_cp), ncp_best
integer :: top_idx(ntop), i
integer(kind=long_int) :: start_count
real(kind=dp) :: intercepts(nreg), slopes(nreg), sigma
real(kind=dp) :: t(n), y(n), cp_prob(n), fitted_mean(n), model_weights(max_cp + 1), top_prob(ntop)

regime_starts = [0, 70, 145]
true_cps = [70, 145]
intercepts = [1.0_dp, 6.0_dp, -2.0_dp]
slopes = [0.05_dp, -0.01_dp, 0.03_dp]
sigma = 0.30_dp

call seed_rng_fixed(seed_time)
call simulate_irregular_time_grid(n, 0.0_dp, 1.0_dp, 0.35_dp, t)
call seed_rng_fixed(seed_y)
call simulate_piecewise_linear_irregular_1d(t, regime_starts, intercepts, slopes, sigma, y)

call system_clock(start_count)
call solve_beast_irreg_trend_only_1d(t, y, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, fitted_mean, model_weights)
call top_cp_probabilities(cp_prob, top_idx, top_prob)

print '(a,i0)', 'n = ', n
print '(a)', 'season = none'
print '(a)', 'irregular = True'
print '(a,i0)', 'max_cp = ', max_cp
print '(a,i0)', 'min_seg_len = ', min_seg_len
print '(a,2(1x,i0))', 'true changepoints =', true_cps
if (ncp_best > 0) then
    write(*, '(a)', advance='no') 'estimated changepoints ='
    print '(100(1x,i0))', best_cps(1:ncp_best)
else
    print '(a)', 'estimated changepoints ='
end if
write(*, '(a)', advance='no') 'model weights ='
print '(100(1x,f0.6))', model_weights
write(*, '(a)', advance='no') 'top cp probabilities ='
do i = 1, ntop
    if (top_idx(i) <= 0) exit
    write(*, '(" (",i0,",",f0.6,")")', advance='no') top_idx(i), top_prob(i)
end do
print *
call print_wall_time(start_count)

end program xsim_beast_irreg
