program xsim_cv_dp_poly_file
! read univariate data from a text file and run cross-validated grid search for DP.poly
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_cv_dp_poly_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: r = 2 ! polynomial order
integer, parameter :: n_gamma = 4 ! number of candidate penalties
real(kind=dp), parameter :: gamma_set(n_gamma) = [0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp] ! candidate l0 penalties
integer, parameter :: delta = 5 ! minimum spacing

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), test_error(:), train_error(:)
integer, allocatable :: cpt_hat(:, :), k_hat(:), cpt_one(:)
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

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(y)
print *, "r                  = ", r
print *, "gamma_set          = ", gamma_set
print *, "delta              = ", delta
do j = 1, n_gamma
    print *, "gamma              = ", gamma_set(j)
    if (k_hat(j) > 0) then
        write (*,'(A)', advance='no') " estimated          ="
        call print_int_list(cpt_hat(1:k_hat(j), j))
    else
        print *, "estimated          = "
    end if
    print *, "n changepoints     = ", k_hat(j)
    print *, "test error         = ", test_error(j)
    print *, "train error        = ", train_error(j)
end do
print *, "best gamma         = ", gamma_set(min_idx(1))
print *, "best index         = ", min_idx(1)
if (k_hat(min_idx(1)) > 0) then
    write (*,'(A)', advance='no') " best changepoints  ="
    call print_int_list(cpt_hat(1:k_hat(min_idx(1)), min_idx(1)))
else
    print *, "best changepoints  = "
end if
print *, "best test error    = ", test_error(min_idx(1))

deallocate(dat, y, cpt_hat, k_hat, test_error, train_error)
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

end program xsim_cv_dp_poly_file
