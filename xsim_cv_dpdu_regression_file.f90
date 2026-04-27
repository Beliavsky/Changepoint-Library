program xsim_cv_dpdu_regression_file
! read regression data from a text file and run cross-validated grid search for DPDU.regression
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_dpdu_mod, only: solve_cv_dpdu_regression_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: eps = 0.001_dp ! lasso convergence tolerance
integer, parameter :: n_lambda = 3 ! number of candidate lasso penalties
integer, parameter :: n_zeta = 3 ! number of candidate zeta values
real(kind=dp), parameter :: lambda_set(n_lambda) = [0.2_dp, 0.5_dp, 1.0_dp] ! candidate lasso penalties
integer, parameter :: zeta_set(n_zeta) = [10, 20, 30] ! candidate minimum segment lengths

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), test_error(:), train_error(:), test_error_mat(:, :)
integer, allocatable :: cpt_hat(:, :), k_hat(:)
integer(kind=long_int) :: t_start
integer :: i, j, idx, max_k_local
integer, allocatable :: cpt_one(:)
integer :: min_idx(2)

call system_clock(t_start)
call get_data_file_arg("xcv_dpdu_regression_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))

allocate(k_hat(n_lambda * n_zeta), test_error(n_lambda * n_zeta), train_error(n_lambda * n_zeta))
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

allocate(cpt_hat(max_k_local, n_lambda * n_zeta), test_error_mat(n_zeta, n_lambda))
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

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(y)
print *, "p                  = ", size(x, 2)
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

deallocate(dat, y, x, cpt_hat, k_hat, test_error, train_error, test_error_mat)
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

end program xsim_cv_dpdu_regression_file
