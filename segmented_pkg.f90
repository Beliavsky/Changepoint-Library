module segmented_pkg_mod
use kind_mod, only: dp
implicit none
private
public :: eval_segmented_lm1_draws, eval_segmented_lm2_draws, eval_stepmented_lm1_draws

contains

subroutine eval_segmented_lm1_draws(x, draws, fitted_mean, param_means)
real(kind=dp), intent(in) :: x(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), param_means(:)

integer :: n, ndraw, i, j
real(kind=dp) :: xv, b0, b1, du, psi

n = size(x)
ndraw = size(draws, 1)

allocate(fitted_mean(n), param_means(size(draws, 2)))
fitted_mean = 0.0_dp

do i = 1, ndraw
    psi = draws(i, 1)
    b0 = draws(i, 2)
    b1 = draws(i, 3)
    du = draws(i, 4)
    do j = 1, n
        xv = x(j)
        fitted_mean(j) = fitted_mean(j) + b0 + b1 * xv + du * max(0.0_dp, xv - psi)
    end do
end do

fitted_mean = fitted_mean / real(ndraw, dp)
param_means = sum(draws, dim = 1) / real(ndraw, dp)
end subroutine eval_segmented_lm1_draws

subroutine eval_segmented_lm2_draws(x, draws, fitted_mean, param_means)
real(kind=dp), intent(in) :: x(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), param_means(:)

integer :: n, ndraw, i, j
real(kind=dp) :: xv, psi1, psi2, b0, b1, du1, du2

n = size(x)
ndraw = size(draws, 1)

allocate(fitted_mean(n), param_means(size(draws, 2)))
fitted_mean = 0.0_dp

do i = 1, ndraw
    psi1 = draws(i, 1)
    psi2 = draws(i, 2)
    b0 = draws(i, 3)
    b1 = draws(i, 4)
    du1 = draws(i, 5)
    du2 = draws(i, 6)
    do j = 1, n
        xv = x(j)
        fitted_mean(j) = fitted_mean(j) + b0 + b1 * xv + du1 * max(0.0_dp, xv - psi1) + du2 * max(0.0_dp, xv - psi2)
    end do
end do

fitted_mean = fitted_mean / real(ndraw, dp)
param_means = sum(draws, dim = 1) / real(ndraw, dp)
end subroutine eval_segmented_lm2_draws

subroutine eval_stepmented_lm1_draws(x, draws, fitted_mean, param_means)
real(kind=dp), intent(in) :: x(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), param_means(:)

integer :: n, ndraw, i, j
real(kind=dp) :: xv, psi, b0, jump

n = size(x)
ndraw = size(draws, 1)

allocate(fitted_mean(n), param_means(size(draws, 2)))
fitted_mean = 0.0_dp

do i = 1, ndraw
    psi = draws(i, 1)
    b0 = draws(i, 2)
    jump = draws(i, 3)
    do j = 1, n
        xv = x(j)
        fitted_mean(j) = fitted_mean(j) + b0 + merge(jump, 0.0_dp, xv >= psi)
    end do
end do

fitted_mean = fitted_mean / real(ndraw, dp)
param_means = sum(draws, dim = 1) / real(ndraw, dp)
end subroutine eval_stepmented_lm1_draws

end module segmented_pkg_mod
