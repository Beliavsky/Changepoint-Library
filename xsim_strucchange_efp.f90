program xsim_strucchange_efp
use kind_mod, only: dp
use strucchange_pkg_mod, only: solve_efp_ols_cusum_1d
implicit none

integer, parameter :: n = 120
real(kind=dp) :: x(n), y(n), beta0(n), beta1(n), u1, u2, z
real(kind=dp), allocatable :: xx(:, :), process(:)
integer :: i, imax(1)

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
call solve_efp_ols_cusum_1d(y, xx, process)
imax = maxloc(abs(process))

print *, "argmax abs process =", imax(1) - 1
print *, "max abs process    =", maxval(abs(process))

deallocate(xx, process)
end program xsim_strucchange_efp
