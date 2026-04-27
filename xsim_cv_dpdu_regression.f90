program xsim_cv_dpdu_regression
! simulate regression changepoints and run cross-validated grid search for DPDU.regression
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_dpdu_mod, only: solve_cv_dpdu_regression_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 120 ! number of observations
integer, parameter :: p = 3 ! number of predictors
real(kind=dp), parameter :: eps = 0.001_dp ! lasso convergence tolerance
integer, parameter :: n_lambda = 3 ! number of candidate lasso penalties
integer, parameter :: n_zeta = 3 ! number of candidate zeta values
real(kind=dp), parameter :: lambda_set(n_lambda) = [0.2_dp, 0.5_dp, 1.0_dp] ! candidate lasso penalties
integer, parameter :: zeta_set(n_zeta) = [10, 20, 30] ! candidate minimum segment lengths

real(kind=dp) :: x(n, p), y(n), test_error(n_lambda * n_zeta), train_error(n_lambda * n_zeta), test_error_mat(n_zeta, n_lambda)
real(kind=dp), parameter :: beta1(p) = [1.5_dp, 0.0_dp, -1.0_dp]
real(kind=dp), parameter :: beta2(p) = [0.2_dp, 1.8_dp, -0.2_dp]
real(kind=dp), parameter :: beta3(p) = [-1.2_dp, 0.5_dp, 1.1_dp]
integer, allocatable :: cpt_hat(:, :), k_hat(:)
integer(kind=long_int) :: t_start
integer :: i, j, idx, max_k_local
integer, allocatable :: cpt_one(:)
integer :: min_idx(2)

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_normal_matrix(x)
call fill_response(y, x)

allocate(k_hat(n_lambda * n_zeta))
max_k_local = 0
idx = 0
do i = 1, n_lambda
    do j = 1, n_zeta
        idx = idx + 1
        call solve_cv_dpdu_regression_1d(y, x, lambda_set(i), zeta_set(j), eps, cpt_one, k_hat(idx), test_error(idx), train_error(idx))
        max_k_local = max(max_k_local, k_hat(idx))
        deallocate(cpt_one)
    end do
end do

allocate(cpt_hat(max_k_local, n_lambda * n_zeta))
cpt_hat = 0
idx = 0
do i = 1, n_lambda
    do j = 1, n_zeta
        idx = idx + 1
        call solve_cv_dpdu_regression_1d(y, x, lambda_set(i), zeta_set(j), eps, cpt_one, k_hat(idx), test_error(idx), train_error(idx))
        if (k_hat(idx) > 0) cpt_hat(1:k_hat(idx), idx) = cpt_one
        deallocate(cpt_one)
        test_error_mat(j, i) = test_error(idx)
    end do
end do
min_idx = minloc(test_error_mat)

print *, "n                  = ", n
print *, "p                  = ", p
print *, "lambda_set         = ", lambda_set
print *, "zeta_set           = ", zeta_set
print *, "eps                = ", eps
idx = 0
do i = 1, n_lambda
    do j = 1, n_zeta
        idx = idx + 1
        print *, "lambda             = ", lambda_set(i)
        print *, "zeta               = ", zeta_set(j)
        if (k_hat(idx) > 0) then
            write (*,'(A)', advance='no') " estimated          ="
            call print_int_list(cpt_hat(1:k_hat(idx), idx))
        else
            print *, "estimated          = "
        end if
        print *, "n changepoints     = ", k_hat(idx)
        print *, "test error         = ", test_error(idx)
        print *, "train error        = ", train_error(idx)
    end do
end do
print *, "best lambda        = ", lambda_set(min_idx(2))
print *, "best zeta          = ", zeta_set(min_idx(1))
print *, "best index         = ", min_idx
idx = (min_idx(2) - 1) * n_zeta + min_idx(1)
if (k_hat(idx) > 0) then
    write (*,'(A)', advance='no') " best changepoints  ="
    call print_int_list(cpt_hat(1:k_hat(idx), idx))
else
    print *, "best changepoints  = "
end if
print *, "best test error    = ", test_error(idx)

deallocate(cpt_hat, k_hat)
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

end program xsim_cv_dpdu_regression
