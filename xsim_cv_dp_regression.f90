program xsim_cv_dp_regression
! simulate regression changepoints and run cross-validated grid search for DP.regression
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_mod, only: solve_cv_dp_regression_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 120 ! number of observations
integer, parameter :: p = 3 ! number of predictors
integer, parameter :: delta = 5 ! minimum spacing
real(kind=dp), parameter :: eps = 0.001_dp ! lasso convergence tolerance
integer, parameter :: n_gamma = 3 ! number of candidate l0 penalties
integer, parameter :: n_lambda = 3 ! number of candidate lasso penalties
real(kind=dp), parameter :: gamma_set(n_gamma) = [1.0_dp, 2.0_dp, 3.0_dp] ! candidate l0 penalties
real(kind=dp), parameter :: lambda_set(n_lambda) = [0.2_dp, 0.5_dp, 1.0_dp] ! candidate lasso penalties

real(kind=dp) :: x(n, p), y(n)
real(kind=dp), parameter :: beta1(p) = [1.5_dp, 0.0_dp, -1.0_dp]
real(kind=dp), parameter :: beta2(p) = [0.2_dp, 1.8_dp, -0.2_dp]
real(kind=dp), parameter :: beta3(p) = [-1.2_dp, 0.5_dp, 1.1_dp]
integer, allocatable :: cpt_hat(:, :), k_hat(:)
real(kind=dp), allocatable :: test_error(:), train_error(:)
integer(kind=long_int) :: t_start
integer :: j, i, min_idx(2)

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_normal_matrix(x)
call fill_response(y, x)
call evaluate_grid(y, x, gamma_set, lambda_set, delta, eps, cpt_hat, k_hat, test_error, train_error)

min_idx = minloc(reshape(test_error, [n_gamma, n_lambda]))
print *, "n                  =", n
print *, "p                  =", p
print *, "gamma_set          =", gamma_set
print *, "lambda_set         =", lambda_set
print *, "delta              =", delta
print *, "eps                =", eps
do i = 1, n_lambda
    do j = 1, n_gamma
        print *, "gamma              =", gamma_set(j)
        print *, "lambda             =", lambda_set(i)
        if (k_hat((i - 1) * n_gamma + j) > 0) then
            write (*,'(A)', advance='no') " estimated          ="
            call print_int_list(cpt_hat(1:k_hat((i - 1) * n_gamma + j), (i - 1) * n_gamma + j))
        else
            print *, "estimated          ="
        end if
        print *, "n changepoints     =", k_hat((i - 1) * n_gamma + j)
        print *, "test error         =", test_error((i - 1) * n_gamma + j)
        print *, "train error        =", train_error((i - 1) * n_gamma + j)
    end do
end do
print *, "best gamma         =", gamma_set(min_idx(1))
print *, "best lambda        =", lambda_set(min_idx(2))
print *, "best index         =", min_idx
if (k_hat((min_idx(2) - 1) * n_gamma + min_idx(1)) > 0) then
    write (*,'(A)', advance='no') " best changepoints  ="
    call print_int_list(cpt_hat(1:k_hat((min_idx(2) - 1) * n_gamma + min_idx(1)), (min_idx(2) - 1) * n_gamma + min_idx(1)))
else
    print *, "best changepoints  ="
end if
print *, "best test error    =", test_error((min_idx(2) - 1) * n_gamma + min_idx(1))

deallocate(cpt_hat, k_hat, test_error, train_error)
call print_wall_time(t_start)

contains

    subroutine evaluate_grid(y, x, gamma_set, lambda_set, delta, eps, cpt_hat, k_hat, test_error, train_error)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: gamma_set(:), lambda_set(:), eps
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: cpt_hat(:, :), k_hat(:)
        real(kind=dp), allocatable, intent(out) :: test_error(:), train_error(:)

        integer :: i, j, idx, max_k_local
        integer, allocatable :: cpt_one(:)

        allocate(k_hat(size(gamma_set) * size(lambda_set)), test_error(size(gamma_set) * size(lambda_set)), &
                 train_error(size(gamma_set) * size(lambda_set)))
        max_k_local = 0
        idx = 0
        do i = 1, size(lambda_set)
            do j = 1, size(gamma_set)
                idx = idx + 1
                call solve_cv_dp_regression_1d(y, x, gamma_set(j), lambda_set(i), delta, eps, cpt_one, k_hat(idx), test_error(idx), train_error(idx))
                max_k_local = max(max_k_local, k_hat(idx))
                deallocate(cpt_one)
            end do
        end do

        allocate(cpt_hat(max_k_local, size(gamma_set) * size(lambda_set)))
        cpt_hat = 0
        idx = 0
        do i = 1, size(lambda_set)
            do j = 1, size(gamma_set)
                idx = idx + 1
                call solve_cv_dp_regression_1d(y, x, gamma_set(j), lambda_set(i), delta, eps, cpt_one, k_hat(idx), test_error(idx), train_error(idx))
                if (k_hat(idx) > 0) cpt_hat(1:k_hat(idx), idx) = cpt_one
                deallocate(cpt_one)
            end do
        end do
    end subroutine evaluate_grid

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

end program xsim_cv_dp_regression
