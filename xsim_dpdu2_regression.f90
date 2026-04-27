program xsim_dpdu2_regression
! simulate regression changepoints and estimate them by DPDU2.regression
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_dpdu_mod, only: solve_dpdu2_regression_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 120 ! number of observations
integer, parameter :: p = 3 ! number of predictors
real(kind=dp), parameter :: lambda = 0.5_dp ! lasso penalty parameter
integer, parameter :: zeta = 20 ! minimum segment length / l0 penalty parameter
real(kind=dp), parameter :: eps = 0.001_dp ! lasso convergence tolerance

real(kind=dp) :: x(n, p), y(n)
real(kind=dp), parameter :: beta1(p) = [1.5_dp, 0.0_dp, -1.0_dp]
real(kind=dp), parameter :: beta2(p) = [0.2_dp, 1.8_dp, -0.2_dp]
real(kind=dp), parameter :: beta3(p) = [-1.2_dp, 0.5_dp, 1.1_dp]
real(kind=dp), allocatable :: beta_mat(:, :)
integer, allocatable :: partition(:), cpt_hat(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_normal_matrix(x)
call fill_response(y, x)
call solve_dpdu2_regression_1d(y, x, lambda, zeta, eps, partition, cpt_hat, beta_mat)

print *, "n                  = ", n
print *, "p                  = ", p
print *, "lambda             = ", lambda
print *, "zeta               = ", zeta
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          = "
end if
print *, "beta checksum      = ", sum(beta_mat)

deallocate(partition, cpt_hat, beta_mat)
call print_wall_time(t_start)

contains

    subroutine fill_normal_matrix(a)
        real(kind=dp), intent(out) :: a(:, :)
        integer :: i, j
        real(kind=dp) :: u1, u2, r, theta, z1, z2

        i = 1
        j = 1
        do while (i <= size(a, 1))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, 1.0e-12_dp)
            r = sqrt(-2.0_dp * log(u1))
            theta = 2.0_dp * acos(-1.0_dp) * u2
            z1 = r * cos(theta)
            z2 = r * sin(theta)
            a(i, j) = z1
            j = j + 1
            if (j > size(a, 2)) then
                j = 1
                i = i + 1
                if (i > size(a, 1)) exit
            end if
            a(i, j) = z2
            j = j + 1
            if (j > size(a, 2)) then
                j = 1
                i = i + 1
            end if
        end do
    end subroutine fill_normal_matrix

    subroutine fill_response(y_out, x_in)
        real(kind=dp), intent(out) :: y_out(:)
        real(kind=dp), intent(in) :: x_in(:, :)
        integer :: i

        do i = 1, 40
            y_out(i) = dot_product(x_in(i, :), beta1) + 0.5_dp * sample_standard_normal()
        end do
        do i = 41, 80
            y_out(i) = dot_product(x_in(i, :), beta2) + 0.5_dp * sample_standard_normal()
        end do
        do i = 81, 120
            y_out(i) = dot_product(x_in(i, :), beta3) + 0.5_dp * sample_standard_normal()
        end do
    end subroutine fill_response

    real(kind=dp) function sample_standard_normal() result(z)
        real(kind=dp) :: u1, u2
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, 1.0e-12_dp)
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
    end function sample_standard_normal

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_dpdu2_regression
