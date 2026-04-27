program xsim_bs_cov
! simulate a multivariate series with covariance changes and estimate changepoints by BS.cov
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_mv
use changepoints_pkg_mod, only: solve_bs_cov_1d, threshold_wbs_tree_1d
use util_mod, only: sort_int, print_wall_time
implicit none

integer, parameter :: n = 180 ! number of time points
integer, parameter :: p = 4 ! number of variables
integer, parameter :: n_regime = 3 ! number of covariance regimes
integer, parameter :: regime_starts(n_regime) = [0, 60, 120] ! zero-based regime starts
real(kind=dp), parameter :: tau = 7.0_dp ! threshold for trimming the BS tree
real(kind=dp), parameter :: means(n_regime, p) = reshape([ &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ], [n_regime, p]) ! segment means
real(kind=dp), parameter :: sds(n_regime, p) = reshape([ &
    1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
    2.4_dp, 0.5_dp, 1.8_dp, 0.6_dp, &
    0.6_dp, 2.1_dp, 0.5_dp, 1.9_dp ], [n_regime, p]) ! diagonal segment scales

real(kind=dp) :: x_time_major(n, p)
real(kind=dp), allocatable :: x(:, :), dval_nodes(:), dval_hat(:)
integer, allocatable :: s_nodes(:), level_nodes(:), parent_nodes(:, :), parent_idx(:), cpt_hat(:)
integer :: n_nodes, n_keep, j
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(1234)
call simulate_piecewise_normal_mv(n, p, regime_starts, means, sds, x_time_major)
allocate(x(p, n))
do j = 1, p
    x(j, :) = x_time_major(:, j)
end do

call solve_bs_cov_1d(x, 1, n, n_nodes, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx)
call threshold_wbs_tree_1d(n_nodes, s_nodes, dval_nodes, parent_idx, tau, n_keep, cpt_hat, dval_hat)
if (n_keep > 0) call sort_int(cpt_hat)

print *, "n                  =", n
print *, "p                  =", p
print *, "tau                =", tau
if (n_keep > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
print *, "raw nodes          =", n_nodes

deallocate(x, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx, cpt_hat, dval_hat)
call print_wall_time(t_start)

contains

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_bs_cov
