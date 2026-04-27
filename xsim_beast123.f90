program xsim_beast123
use kind_mod, only: dp, long_int
use util_mod, only: print_wall_time
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_linear_mv_1d
use rbeast_mod, only: solve_beast123_regular_mv, top_cp_probabilities
implicit none

integer, parameter :: n = 240              ! number of observations
integer, parameter :: p = 3                ! number of series
integer, parameter :: max_cp = 4           ! maximum number of trend changepoints
integer, parameter :: min_seg_len = 30     ! minimum segment length
integer, parameter :: nreg = 3             ! number of linear regimes
integer, parameter :: ntop = 5             ! number of cp probabilities to print
integer, parameter :: seed = 20260430      ! RNG seed for reproducibility

integer :: regime_starts(nreg), true_cps(nreg - 1), best_cps(max_cp, p), ncp_best(p)
integer :: top_idx(ntop), i, j
integer(kind=long_int) :: start_count
real(kind=dp) :: sigma, intercepts(nreg, p), slopes(nreg, p)
real(kind=dp) :: x(n, p), cp_prob(n, p), mean_cp_prob(n), fitted_mean(n, p), model_weights(max_cp + 1, p), top_prob(ntop)

regime_starts = [0, 80, 160]
true_cps = [80, 160]
sigma = 0.35_dp
intercepts = reshape([1.0_dp, 4.6_dp, -1.2_dp, -1.0_dp, 2.0_dp, 3.2_dp, 2.5_dp, -0.5_dp, 1.0_dp], [nreg, p])
slopes = reshape([0.03_dp, -0.02_dp, 0.015_dp, 0.015_dp, -0.01_dp, 0.02_dp, -0.02_dp, 0.025_dp, -0.015_dp], [nreg, p])

call seed_rng_fixed(seed)
call simulate_piecewise_linear_mv_1d(n, p, regime_starts, intercepts, slopes, sigma, x)

call system_clock(start_count)
call solve_beast123_regular_mv(x, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, mean_cp_prob, fitted_mean, model_weights)
call top_cp_probabilities(mean_cp_prob, top_idx, top_prob)

print '(a,i0)', 'n = ', n
print '(a,i0)', 'p = ', p
print '(a)', 'season = none'
print '(a)', 'beast123 regular interface'
print '(a,i0)', 'max_cp = ', max_cp
print '(a,i0)', 'min_seg_len = ', min_seg_len
print '(a,2(1x,i0))', 'true changepoints =', true_cps
do j = 1, p
    write(*, '(a,i0,a)', advance='no') 'series ', j, ' changepoints ='
    if (ncp_best(j) > 0) then
        print '(100(1x,i0))', best_cps(1:ncp_best(j), j)
    else
        print *
    end if
end do
write(*, '(a)', advance='no') 'top mean cp probabilities ='
do i = 1, ntop
    if (top_idx(i) <= 0) exit
    write(*, '(" (",i0,",",f0.6,")")', advance='no') top_idx(i), top_prob(i)
end do
print *
call print_wall_time(start_count)

end program xsim_beast123
