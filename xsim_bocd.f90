program xsim_bocd
use kind_mod, only: dp, long_int
use util_mod, only: print_wall_time
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_bocd_student_t_1d
implicit none

integer, parameter :: n = 400              ! number of observations
real(kind=dp), parameter :: lam = 100.0_dp ! constant hazard scale
integer, parameter :: nreg = 3             ! number of Gaussian regimes
integer, parameter :: seed = 20260422      ! RNG seed for reproducibility

integer :: regime_starts(nreg), true_cps(nreg - 1)
integer, allocatable :: map_cps(:), map_cps_pkg(:), runlen_argmax(:)
integer(kind=long_int) :: start_count
real(kind=dp) :: means(nreg), sds(nreg), x(n)

regime_starts = [0, 150, 300]
true_cps = [150, 300]
means = [0.0_dp, 3.0_dp, -2.0_dp]
sds = [1.0_dp, 1.0_dp, 1.0_dp]

call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)

call system_clock(start_count)
call solve_bocd_student_t_1d(x, lam, map_cps, runlen_argmax, mu0=0.0_dp, kappa0=1.0_dp, alpha0=1.0_dp, beta0=1.0_dp)
map_cps_pkg = pack(map_cps + 1, mask=map_cps > 0)

print '(a,i0)', 'n = ', n
print '(a,f0.1)', 'lambda = ', lam
print '(a,2(1x,i0))', 'true changepoints =', true_cps
write(*, '(a)', advance='no') 'estimated changepoints ='
if (size(map_cps_pkg) > 0) then
    print '(100(1x,i0))', map_cps_pkg
else
    print *
end if
print '(a,i0)', 'final MAP run length = ', runlen_argmax(size(runlen_argmax))
call print_wall_time(start_count)

deallocate(map_cps, map_cps_pkg, runlen_argmax)
end program xsim_bocd
