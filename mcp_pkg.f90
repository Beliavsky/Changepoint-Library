module mcp_pkg_mod
use kind_mod, only: dp
implicit none
private
public :: eval_mcp_demo_draws, eval_mcp_sigma_draws, eval_mcp_ar1_draws, eval_mcp_ar_change_draws, eval_mcp_arsigma_draws

contains

subroutine eval_mcp_demo_draws(time, draws, fitted_mean, cp_means, param_means)
real(kind=dp), intent(in) :: time(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), cp_means(:), param_means(:)

integer :: n, ndraw, i, j
real(kind=dp) :: t, cp1, cp2, int1, int3, time2, time3

n = size(time)
ndraw = size(draws, 1)

allocate(fitted_mean(n), cp_means(2), param_means(size(draws, 2)))
fitted_mean = 0.0_dp

do i = 1, ndraw
    cp1 = draws(i, 1)
    cp2 = draws(i, 2)
    int1 = draws(i, 3)
    int3 = draws(i, 4)
    time2 = draws(i, 6)
    time3 = draws(i, 7)
    do j = 1, n
        t = time(j)
        if (t < cp1) then
            fitted_mean(j) = fitted_mean(j) + int1
        else if (t < cp2) then
            fitted_mean(j) = fitted_mean(j) + int1 + time2 * (t - cp1)
        else
            fitted_mean(j) = fitted_mean(j) + int3 + time3 * (t - cp2)
        end if
    end do
end do

fitted_mean = fitted_mean / real(ndraw, dp)
cp_means = [sum(draws(:, 1)) / real(ndraw, dp), sum(draws(:, 2)) / real(ndraw, dp)]
param_means = sum(draws, dim = 1) / real(ndraw, dp)
end subroutine eval_mcp_demo_draws

subroutine eval_mcp_sigma_draws(x, draws, fitted_mean, sigma_mean, cp_means, param_means)
real(kind=dp), intent(in) :: x(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), sigma_mean(:), cp_means(:), param_means(:)

integer :: n, ndraw, i, j
real(kind=dp) :: xv, cp1, cp2, int1, sigma1, sigma2, sigma_x2, x3, sig

n = size(x)
ndraw = size(draws, 1)

allocate(fitted_mean(n), sigma_mean(n), cp_means(2), param_means(size(draws, 2)))
fitted_mean = 0.0_dp
sigma_mean = 0.0_dp

do i = 1, ndraw
    cp1 = draws(i, 1)
    cp2 = draws(i, 2)
    int1 = draws(i, 3)
    sigma1 = draws(i, 4)
    sigma2 = draws(i, 5)
    sigma_x2 = draws(i, 6)
    x3 = draws(i, 7)
    do j = 1, n
        xv = x(j)
        if (xv < cp1) then
            fitted_mean(j) = fitted_mean(j) + int1
            sigma_mean(j) = sigma_mean(j) + sigma1
        else
            if (xv < cp2) then
                fitted_mean(j) = fitted_mean(j) + int1
            else
                fitted_mean(j) = fitted_mean(j) + int1 + x3 * (xv - cp2)
            end if
            sig = sigma2 + sigma_x2 * (xv - cp1)
            if (sig < 0.0_dp) sig = 0.0_dp
            sigma_mean(j) = sigma_mean(j) + sig
        end if
    end do
end do

fitted_mean = fitted_mean / real(ndraw, dp)
sigma_mean = sigma_mean / real(ndraw, dp)
cp_means = [sum(draws(:, 1)) / real(ndraw, dp), sum(draws(:, 2)) / real(ndraw, dp)]
param_means = sum(draws, dim = 1) / real(ndraw, dp)
end subroutine eval_mcp_sigma_draws

subroutine eval_mcp_ar1_draws(time, draws, fitted_mean, param_means)
real(kind=dp), intent(in) :: time(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), param_means(:)

integer :: n, ndraw
real(kind=dp) :: int1_mean

n = size(time)
ndraw = size(draws, 1)

allocate(fitted_mean(n), param_means(size(draws, 2)))
param_means = sum(draws, dim = 1) / real(ndraw, dp)
int1_mean = param_means(2)
fitted_mean = int1_mean
end subroutine eval_mcp_ar1_draws

subroutine eval_mcp_ar_change_draws(x, draws, fitted_mean, cp_means, param_means)
real(kind=dp), intent(in) :: x(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), cp_means(:), param_means(:)

integer :: n, ndraw, i, j
real(kind=dp) :: xv, cp1, cp2, int1, x2, plateau

n = size(x)
ndraw = size(draws, 1)

allocate(fitted_mean(n), cp_means(2), param_means(size(draws, 2)))
fitted_mean = 0.0_dp

do i = 1, ndraw
    cp1 = draws(i, 5)
    cp2 = draws(i, 6)
    int1 = draws(i, 7)
    x2 = draws(i, 9)
    plateau = int1 + x2 * (cp2 - cp1)
    do j = 1, n
        xv = x(j)
        if (xv < cp1) then
            fitted_mean(j) = fitted_mean(j) + int1
        else if (xv < cp2) then
            fitted_mean(j) = fitted_mean(j) + int1 + x2 * (xv - cp1)
        else
            fitted_mean(j) = fitted_mean(j) + plateau
        end if
    end do
end do

fitted_mean = fitted_mean / real(ndraw, dp)
cp_means = [sum(draws(:, 5)) / real(ndraw, dp), sum(draws(:, 6)) / real(ndraw, dp)]
param_means = sum(draws, dim = 1) / real(ndraw, dp)
end subroutine eval_mcp_ar_change_draws

subroutine eval_mcp_arsigma_draws(x, draws, fitted_mean, sigma_mean, cp_mean, param_means)
real(kind=dp), intent(in) :: x(:), draws(:, :)
real(kind=dp), allocatable, intent(out) :: fitted_mean(:), sigma_mean(:), param_means(:)
real(kind=dp), intent(out) :: cp_mean

integer :: n, ndraw, i, j
real(kind=dp) :: xv, cp1, int1, sigma1, sigma2, sigma_x1, x2, sig

n = size(x)
ndraw = size(draws, 1)

allocate(fitted_mean(n), sigma_mean(n), param_means(size(draws, 2)))
fitted_mean = 0.0_dp
sigma_mean = 0.0_dp

do i = 1, ndraw
    cp1 = draws(i, 6)
    int1 = draws(i, 7)
    sigma1 = draws(i, 8)
    sigma2 = draws(i, 9)
    sigma_x1 = draws(i, 10)
    x2 = draws(i, 11)
    do j = 1, n
        xv = x(j)
        if (xv < cp1) then
            fitted_mean(j) = fitted_mean(j) + int1
            sig = sigma1 + sigma_x1 * xv
            if (sig < 0.0_dp) sig = 0.0_dp
            sigma_mean(j) = sigma_mean(j) + sig
        else
            fitted_mean(j) = fitted_mean(j) + int1 + x2 * (xv - cp1)
            sigma_mean(j) = sigma_mean(j) + sigma2
        end if
    end do
end do

fitted_mean = fitted_mean / real(ndraw, dp)
sigma_mean = sigma_mean / real(ndraw, dp)
cp_mean = sum(draws(:, 6)) / real(ndraw, dp)
param_means = sum(draws, dim = 1) / real(ndraw, dp)
end subroutine eval_mcp_arsigma_draws

end module mcp_pkg_mod
