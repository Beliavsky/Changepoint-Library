program xsim_strucchange_fstats
use kind_mod, only: dp
use strucchange_pkg_mod, only: solve_fstats_regression_1d
implicit none

real(kind=dp), parameter :: from = 0.2_dp, to = 0.8_dp
integer, parameter :: n = 120
real(kind=dp) :: x(n), y(n), beta0(n), beta1(n), u1, u2, z
real(kind=dp), allocatable :: xx(:, :), fstats(:)
integer, allocatable :: point_idx(:)
integer :: i, breakpoint
real(kind=dp) :: min_rss

call random_seed()
do i = 1, n
    x(i) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
end do
beta0(1:40) = -0.5_dp
beta0(41:80) = 1.0_dp
beta0(81:n) = -1.2_dp
beta1(1:40) = 1.2_dp
beta1(41:80) = -0.8_dp
beta1(81:n) = 0.9_dp
do i = 1, n
    call random_number(u1)
    call random_number(u2)
    z = sqrt(-2.0_dp * log(max(u1, 1.0e-12_dp))) * cos(2.0_dp * acos(-1.0_dp) * u2)
    y(i) = beta0(i) + beta1(i) * x(i) + 0.18_dp * z
end do

allocate(xx(n, 1))
xx(:, 1) = x
call solve_fstats_regression_1d(y, xx, from, to, fstats, point_idx, breakpoint, min_rss)

print *, "breakpoint         =", breakpoint
print *, "max F              =", maxval(fstats)
print *, "min RSS            =", min_rss

deallocate(xx, fstats, point_idx)
end program xsim_strucchange_fstats
