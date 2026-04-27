program xsim_wbs_univar
! simulate a univariate mean-change series and estimate changepoints by WBS
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoints_pkg_mod, only: solve_wbs_univar_1d, threshold_wbs_tree_1d
use util_mod, only: sort_int, print_wall_time
implicit none

integer, parameter :: n = 300 ! number of observations
integer, parameter :: delta = 5 ! minimum spacing used by WBS
integer, parameter :: m = 300 ! number of random intervals
real(kind=dp), parameter :: tau = 4.0_dp ! threshold for trimming the WBS tree

integer, dimension(4) :: regime_starts
real(kind=dp), dimension(4) :: means, sds
real(kind=dp), allocatable :: y(:)
integer, allocatable :: alpha(:), beta(:), s_nodes(:), level_nodes(:), parent_nodes(:, :), parent_idx(:), cpt_hat(:)
real(kind=dp), allocatable :: dval_nodes(:), dval_hat(:)
integer :: n_nodes, n_keep
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(y(n), alpha(m), beta(m))

call seed_rng_fixed(0)
regime_starts = [0, 20, 50, 170]
means = [0.0_dp, 2.0_dp, 0.0_dp, -2.0_dp]
sds = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, y)

call seed_rng_fixed(1)
call generate_wbs_intervals(m, 1, n, alpha, beta)

call solve_wbs_univar_1d(y, alpha, beta, 1, n, delta, n_nodes, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx)
call threshold_wbs_tree_1d(n_nodes, s_nodes, dval_nodes, parent_idx, tau, n_keep, cpt_hat, dval_hat)
if (n_keep > 0) call sort_int(cpt_hat)

print *, "n                  =", n
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

deallocate(y, alpha, beta, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx, cpt_hat, dval_hat)
call print_wall_time(t_start)

contains

    subroutine generate_wbs_intervals(m, lower, upper, alpha, beta)
        integer, intent(in) :: m, lower, upper
        integer, intent(out) :: alpha(m), beta(m)
        integer :: j, a, b
        real(kind=dp) :: u
        do j = 1, m
            call random_number(u)
            a = lower + int((upper - lower + 1) * u)
            if (a > upper) a = upper
            call random_number(u)
            b = lower + int((upper - lower + 1) * u)
            if (b > upper) b = upper
            alpha(j) = min(a, b)
            beta(j) = max(a, b)
        end do
    end subroutine generate_wbs_intervals

    subroutine print_int_list(x)
        integer, intent(in) :: x(:)
        integer :: i
        do i = 1, size(x)
            write (*,'(1X,I0)', advance='no') x(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_wbs_univar
