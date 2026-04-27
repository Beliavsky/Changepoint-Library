program xsim_wbsip_cov
! simulate a multivariate series and an independent copy, then estimate covariance changepoints by WBSIP.cov
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_mv
use changepoints_pkg_mod, only: solve_wbsip_cov_1d, threshold_wbs_tree_1d
use util_mod, only: sort_int, print_wall_time
implicit none

integer, parameter :: n = 180 ! number of time points
integer, parameter :: p = 4 ! number of variables
integer, parameter :: n_regime = 2 ! number of covariance regimes
integer, parameter :: regime_starts(n_regime) = [0, 80] ! zero-based regime starts
integer, parameter :: m = 200 ! number of WBS intervals
integer, parameter :: delta = 5 ! minimum spacing parameter
real(kind=dp), parameter :: tau = 8.0_dp ! trimming threshold
real(kind=dp), parameter :: means(n_regime, p) = 0.0_dp ! segment means
real(kind=dp), parameter :: sds_x(n_regime, p) = reshape([ &
    1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
    2.0_dp, 0.6_dp, 0.8_dp, 1.7_dp ], [n_regime, p]) ! diagonal scales for X
real(kind=dp), parameter :: sds_xprime(n_regime, p) = reshape([ &
    1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
    2.0_dp, 0.6_dp, 0.8_dp, 1.7_dp ], [n_regime, p]) ! diagonal scales for X'

real(kind=dp) :: x_time_major(n, p), xprime_time_major(n, p)
real(kind=dp), allocatable :: x(:, :), x_prime(:, :), dval_nodes(:), dval_hat(:)
integer, allocatable :: alpha(:), beta(:), s_nodes(:), level_nodes(:), parent_nodes(:, :), parent_idx(:), cpt_hat(:)
integer :: n_nodes, n_keep, j
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(1234)
call simulate_piecewise_normal_mv(n, p, regime_starts, means, sds_x, x_time_major)
call seed_rng_fixed(4321)
call simulate_piecewise_normal_mv(n, p, regime_starts, means, sds_xprime, xprime_time_major)
allocate(x(p, n), x_prime(p, n), alpha(m), beta(m))
do j = 1, p
    x(j, :) = x_time_major(:, j)
    x_prime(j, :) = xprime_time_major(:, j)
end do
call generate_intervals(m, n, 202, alpha, beta)

call solve_wbsip_cov_1d(x, x_prime, alpha, beta, 1, n, delta, n_nodes, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx)
call threshold_wbs_tree_1d(n_nodes, s_nodes, dval_nodes, parent_idx, tau, n_keep, cpt_hat, dval_hat)
if (n_keep > 0) call sort_int(cpt_hat)

print *, "n                  =", n
print *, "p                  =", p
print *, "M                  =", m
print *, "delta              =", delta
print *, "tau                =", tau
if (n_keep > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
print *, "raw nodes          =", n_nodes

deallocate(x, x_prime, alpha, beta, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx, cpt_hat, dval_hat)
call print_wall_time(t_start)

contains

    subroutine generate_intervals(m_int, n_int, seed, alpha_out, beta_out)
        integer, intent(in) :: m_int, n_int, seed
        integer, intent(out) :: alpha_out(m_int), beta_out(m_int)
        integer :: i
        real(kind=dp) :: u1, u2

        call seed_rng_fixed(seed)
        do i = 1, m_int
            call random_number(u1)
            call random_number(u2)
            alpha_out(i) = 1 + int(u1 * real(n_int, dp))
            beta_out(i) = 1 + int(u2 * real(n_int, dp))
            alpha_out(i) = min(alpha_out(i), n_int)
            beta_out(i) = min(beta_out(i), n_int)
            if (alpha_out(i) > beta_out(i)) call swap_int(alpha_out(i), beta_out(i))
        end do
    end subroutine generate_intervals

    subroutine swap_int(a, b)
        integer, intent(inout) :: a, b
        integer :: tmp
        tmp = a
        a = b
        b = tmp
    end subroutine swap_int

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_wbsip_cov
