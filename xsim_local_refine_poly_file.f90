program xsim_local_refine_poly_file
! read univariate data, cross-validate DP.poly, and refine changepoints locally
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_cv_dp_poly_1d, solve_local_refine_poly_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: r = 2 ! polynomial order
integer, parameter :: n_gamma = 4 ! number of candidate penalties
real(kind=dp), parameter :: gamma_set(n_gamma) = [0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp] ! candidate l0 penalties
integer, parameter :: delta = 5 ! minimum spacing
integer, parameter :: delta_lr = 5 ! local refinement spacing

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), test_error(:), train_error(:)
integer, allocatable :: cpt_hat(:, :), k_hat(:), cpt_one(:), cpt_refined(:)
integer(kind=long_int) :: t_start
integer :: j, max_k_local, min_idx(1)

call system_clock(t_start)
call get_data_file_arg("xcv_dp_poly_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)))
y = dat(:, 1)

allocate(k_hat(n_gamma), test_error(n_gamma), train_error(n_gamma))
max_k_local = 0
do j = 1, n_gamma
    call solve_cv_dp_poly_1d(y, r, gamma_set(j), delta, cpt_one, k_hat(j), test_error(j), train_error(j))
    max_k_local = max(max_k_local, k_hat(j))
    deallocate(cpt_one)
end do
allocate(cpt_hat(max_k_local, n_gamma))
cpt_hat = 0
do j = 1, n_gamma
    call solve_cv_dp_poly_1d(y, r, gamma_set(j), delta, cpt_one, k_hat(j), test_error(j), train_error(j))
    if (k_hat(j) > 0) cpt_hat(1:k_hat(j), j) = cpt_one
    deallocate(cpt_one)
end do
min_idx = minloc(test_error)
call solve_local_refine_poly_1d(cpt_hat(1:k_hat(min_idx(1)), min_idx(1)), y, r, delta_lr, cpt_refined)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(y)
print *, "r                  = ", r
print *, "gamma_set          = ", gamma_set
print *, "delta              = ", delta
print *, "delta_lr           = ", delta_lr
print *, "best gamma         = ", gamma_set(min_idx(1))
if (k_hat(min_idx(1)) > 0) then
    write (*,'(A)', advance='no') " initial            ="
    call print_int_list(cpt_hat(1:k_hat(min_idx(1)), min_idx(1)))
else
    print *, "initial            = "
end if
if (size(cpt_refined) > 0) then
    write (*,'(A)', advance='no') " refined            ="
    call print_int_list(cpt_refined)
else
    print *, "refined            = "
end if

deallocate(dat, y, cpt_hat, k_hat, cpt_refined)
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

end program xsim_local_refine_poly_file
