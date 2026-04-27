program xsim_local_refine_regression_file
! read regression data from a text file, estimate changepoints by DP.regression, and refine them locally
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_dp_regression_1d, solve_local_refine_regression_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: gamma = 2.0_dp ! l0 penalty parameter
real(kind=dp), parameter :: lambda = 0.5_dp ! lasso penalty parameter
integer, parameter :: delta = 5 ! minimum spacing
real(kind=dp), parameter :: eps = 0.001_dp ! lasso convergence tolerance
real(kind=dp), parameter :: zeta = 0.5_dp ! group-lasso refinement parameter

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :)
integer, allocatable :: partition(:), cpt_init(:), cpt_refined(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xdp_regression_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))
call solve_dp_regression_1d(y, x, gamma, lambda, delta, eps, partition, cpt_init)
call solve_local_refine_regression_1d(cpt_init, y, x, zeta, cpt_refined)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "p                  =", size(x, 2)
print *, "gamma              =", gamma
print *, "lambda             =", lambda
print *, "delta              =", delta
print *, "zeta               =", zeta
if (size(cpt_init) > 0) then
    write (*,'(A)', advance='no') " initial            ="
    call print_int_list(cpt_init)
else
    print *, "initial            ="
end if
if (size(cpt_refined) > 0) then
    write (*,'(A)', advance='no') " refined            ="
    call print_int_list(cpt_refined)
else
    print *, "refined            ="
end if

deallocate(dat, y, x, partition, cpt_init, cpt_refined)
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

end program xsim_local_refine_regression_file
