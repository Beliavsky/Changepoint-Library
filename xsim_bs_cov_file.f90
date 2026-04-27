program xsim_bs_cov_file
! read a multivariate series from a text file and estimate covariance changepoints by BS.cov
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_bs_cov_1d, threshold_wbs_tree_1d
use util_mod, only: sort_int, print_wall_time
implicit none

real(kind=dp), parameter :: tau = 7.0_dp ! threshold for trimming the BS tree

character(len=256) :: data_file
real(kind=dp), allocatable :: x_raw(:, :), x(:, :), dval_nodes(:), dval_hat(:)
integer, allocatable :: s_nodes(:), level_nodes(:), parent_nodes(:, :), parent_idx(:), cpt_hat(:)
integer :: n_nodes, n_keep, i
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xbs_cov_data.txt", data_file)
call read_matrix_file(data_file, x_raw)
allocate(x(size(x_raw, 2), size(x_raw, 1)))
do i = 1, size(x_raw, 1)
    x(:, i) = x_raw(i, :)
end do
call solve_bs_cov_1d(x, 1, size(x, 2), n_nodes, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx)
call threshold_wbs_tree_1d(n_nodes, s_nodes, dval_nodes, parent_idx, tau, n_keep, cpt_hat, dval_hat)
if (n_keep > 0) call sort_int(cpt_hat)

print *, "file               =", trim(data_file)
print *, "n                  =", size(x, 2)
print *, "p                  =", size(x, 1)
print *, "tau                =", tau
if (n_keep > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
print *, "raw nodes          =", n_nodes

deallocate(x_raw, x, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx, cpt_hat, dval_hat)
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

end program xsim_bs_cov_file
