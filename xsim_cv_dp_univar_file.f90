program xsim_cv_dp_univar_file
! read a univariate series from a text file and run cross-validated dynamic programming over a gamma grid
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_cv_dp_univar_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: delta = 5 ! minimum spacing parameter in CV.search.DP.univar
integer, parameter :: n_gamma = 8 ! number of candidate l0 penalties
integer, parameter :: gamma_set(n_gamma) = [1, 2, 3, 4, 5, 6, 7, 8] ! candidate l0 penalties

character(len=256) :: data_file
real(kind=dp), allocatable :: y(:), test_error(:), train_error(:)
integer, allocatable :: cpt_hat(:, :), k_hat(:)
integer(kind=long_int) :: t_start
integer :: j, min_idx

call system_clock(t_start)
call get_data_file_arg("xcv_dp_univar_data.txt", data_file)
call read_series_file(data_file, y)
call evaluate_gamma_grid(y, gamma_set, delta, cpt_hat, k_hat, test_error, train_error)

min_idx = minloc(test_error, dim=1)
print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "gamma_set          =", gamma_set
print *, "delta              =", delta
do j = 1, n_gamma
    print *, "gamma              =", gamma_set(j)
    if (k_hat(j) > 0) then
        write (*,'(A)', advance='no') " estimated          ="
        call print_int_list(cpt_hat(1:k_hat(j), j))
    else
        print *, "estimated          ="
    end if
    print *, "n changepoints     =", k_hat(j)
    print *, "test error         =", test_error(j)
    print *, "train error        =", train_error(j)
end do
print *, "best gamma         =", gamma_set(min_idx)
print *, "best index         =", min_idx
if (k_hat(min_idx) > 0) then
    write (*,'(A)', advance='no') " best changepoints  ="
    call print_int_list(cpt_hat(1:k_hat(min_idx), min_idx))
else
    print *, "best changepoints  ="
end if
print *, "best test error    =", test_error(min_idx)

deallocate(y, cpt_hat, k_hat, test_error, train_error)
call print_wall_time(t_start)

contains

    subroutine evaluate_gamma_grid(y, gamma_set, delta, cpt_hat, k_hat, test_error, train_error)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: gamma_set(:), delta
        integer, allocatable, intent(out) :: cpt_hat(:, :), k_hat(:)
        real(kind=dp), allocatable, intent(out) :: test_error(:), train_error(:)

        integer :: j, max_k_local
        integer, allocatable :: cpt_one(:)

        allocate(k_hat(size(gamma_set)), test_error(size(gamma_set)), train_error(size(gamma_set)))
        max_k_local = 0
        do j = 1, size(gamma_set)
            call solve_cv_dp_univar_1d(y, real(gamma_set(j), dp), delta, cpt_one, k_hat(j), test_error(j), train_error(j))
            max_k_local = max(max_k_local, k_hat(j))
            deallocate(cpt_one)
        end do

        allocate(cpt_hat(max_k_local, size(gamma_set)))
        cpt_hat = 0
        do j = 1, size(gamma_set)
            call solve_cv_dp_univar_1d(y, real(gamma_set(j), dp), delta, cpt_one, k_hat(j), test_error(j), train_error(j))
            if (k_hat(j) > 0) cpt_hat(1:k_hat(j), j) = cpt_one
            deallocate(cpt_one)
        end do
    end subroutine evaluate_gamma_grid

    subroutine print_int_list(x)
        integer, intent(in) :: x(:)
        integer :: i
        do i = 1, size(x)
            write (*,'(1X,I0)', advance='no') x(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_cv_dp_univar_file
