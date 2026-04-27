program xsim_local_refine_dpdu_regression_file
! read regression data from a text file, estimate changepoints by DPDU.regression, and refine them locally
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_dpdu_mod, only: solve_dpdu_regression_1d, solve_local_refine_dpdu_regression_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: lambda = 0.5_dp ! lasso penalty parameter
integer, parameter :: zeta = 20 ! minimum segment length / l0 penalty parameter
real(kind=dp), parameter :: eps = 0.001_dp ! lasso convergence tolerance
real(kind=dp), parameter :: w = 0.9_dp ! local refinement window weight

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), beta_mat(:, :), beta_seg(:, :)
integer, allocatable :: partition(:), cpt_init(:), cpt_refined(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xdpdu_regression_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))
call solve_dpdu_regression_1d(y, x, lambda, zeta, eps, partition, cpt_init, beta_mat)
call extract_segment_betas(beta_mat, cpt_init, beta_seg)
call solve_local_refine_dpdu_regression_1d(cpt_init, beta_seg, y, x, w, cpt_refined)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(y)
print *, "p                  = ", size(x, 2)
print *, "lambda             = ", lambda
print *, "zeta               = ", zeta
print *, "w                  = ", w
if (size(cpt_init) > 0) then
    write (*,'(A)', advance='no') " initial            ="
    call print_int_list(cpt_init)
else
    print *, "initial            = "
end if
if (size(cpt_refined) > 0) then
    write (*,'(A)', advance='no') " refined            ="
    call print_int_list(cpt_refined)
else
    print *, "refined            = "
end if

deallocate(dat, y, x, partition, cpt_init, cpt_refined, beta_mat, beta_seg)
call print_wall_time(t_start)

contains

    subroutine extract_segment_betas(beta_full, cps, beta_seg)
        real(kind=dp), intent(in) :: beta_full(:, :)
        integer, intent(in) :: cps(:)
        real(kind=dp), allocatable, intent(out) :: beta_seg(:, :)
        integer :: nseg

        nseg = size(cps) + 1
        allocate(beta_seg(size(beta_full, 1), nseg))
        if (size(cps) > 0) beta_seg(:, 1:size(cps)) = beta_full(:, cps)
        beta_seg(:, nseg) = beta_full(:, size(beta_full, 2))
    end subroutine extract_segment_betas

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_local_refine_dpdu_regression_file
