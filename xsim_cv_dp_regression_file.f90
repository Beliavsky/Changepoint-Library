program xsim_cv_dp_regression_file
! read regression data from a text file and run cross-validated grid search for DP.regression
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_cv_dp_regression_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: delta = 5 ! minimum spacing
real(kind=dp), parameter :: eps = 0.001_dp ! lasso convergence tolerance
integer, parameter :: n_gamma = 3 ! number of candidate l0 penalties
integer, parameter :: n_lambda = 3 ! number of candidate lasso penalties
real(kind=dp), parameter :: gamma_set(n_gamma) = [1.0_dp, 2.0_dp, 3.0_dp] ! candidate l0 penalties
real(kind=dp), parameter :: lambda_set(n_lambda) = [0.2_dp, 0.5_dp, 1.0_dp] ! candidate lasso penalties

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), test_error(:), train_error(:), test_error_mat(:, :)
integer, allocatable :: cpt_hat(:, :), k_hat(:)
integer(kind=long_int) :: t_start
integer :: i, j, idx, min_idx(2)

call system_clock(t_start)
call get_data_file_arg("xcv_dp_regression_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))
call evaluate_grid(y, x, gamma_set, lambda_set, delta, eps, cpt_hat, k_hat, test_error, train_error)

allocate(test_error_mat(n_gamma, n_lambda))
idx = 0
do i = 1, n_lambda
    do j = 1, n_gamma
        idx = idx + 1
        test_error_mat(j, i) = test_error(idx)
    end do
end do
min_idx = minloc(test_error_mat)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "p                  =", size(x, 2)
print *, "gamma_set          =", gamma_set
print *, "lambda_set         =", lambda_set
print *, "delta              =", delta
print *, "eps                =", eps
idx = 0
do i = 1, n_lambda
    do j = 1, n_gamma
        idx = idx + 1
        print *, "gamma              =", gamma_set(j)
        print *, "lambda             =", lambda_set(i)
        if (k_hat(idx) > 0) then
            write (*,'(A)', advance='no') " estimated          ="
            call print_int_list(cpt_hat(1:k_hat(idx), idx))
        else
            print *, "estimated          ="
        end if
        print *, "n changepoints     =", k_hat(idx)
        print *, "test error         =", test_error(idx)
        print *, "train error        =", train_error(idx)
    end do
end do
print *, "best gamma         =", gamma_set(min_idx(1))
print *, "best lambda        =", lambda_set(min_idx(2))
print *, "best index         =", min_idx
idx = (min_idx(2) - 1) * n_gamma + min_idx(1)
if (k_hat(idx) > 0) then
    write (*,'(A)', advance='no') " best changepoints  ="
    call print_int_list(cpt_hat(1:k_hat(idx), idx))
else
    print *, "best changepoints  ="
end if
print *, "best test error    =", test_error(idx)

deallocate(dat, y, x, cpt_hat, k_hat, test_error, train_error, test_error_mat)
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

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_cv_dp_regression_file
