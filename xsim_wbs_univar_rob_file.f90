program xsim_wbs_univar_rob_file
! read a contaminated univariate series and WBS intervals from text files and estimate changepoints by robust WBS
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_wbs_univar_rob_1d, threshold_wbs_tree_1d
use util_mod, only: sort_int, print_wall_time
implicit none

integer, parameter :: delta = 5 ! minimum spacing used by WBS.uni.rob
real(kind=dp), parameter :: k_huber = 1.345_dp ! Huber robustification parameter
real(kind=dp), parameter :: tau = 8.0_dp ! threshold for trimming the WBS tree

character(len=256) :: data_file, interval_file
real(kind=dp), allocatable :: y(:), interval_mat(:, :), dval_nodes(:), dval_hat(:)
integer, allocatable :: alpha(:), beta(:), s_nodes(:), level_nodes(:), parent_nodes(:, :), parent_idx(:), cpt_hat(:)
integer :: n_nodes, n_keep
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xwbs_rob_data.txt", data_file)
interval_file = derive_interval_file(data_file)
call read_series_file(data_file, y)
call read_matrix_file(interval_file, interval_mat)
allocate(alpha(size(interval_mat, 1)), beta(size(interval_mat, 1)))
alpha = nint(interval_mat(:, 1))
beta = nint(interval_mat(:, 2))

call solve_wbs_univar_rob_1d(y, alpha, beta, 1, size(y), k_huber, delta, n_nodes, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx)
call threshold_wbs_tree_1d(n_nodes, s_nodes, dval_nodes, parent_idx, tau, n_keep, cpt_hat, dval_hat)
if (n_keep > 0) call sort_int(cpt_hat)

print *, "file               =", trim(data_file)
print *, "interval_file      =", trim(interval_file)
print *, "n                  =", size(y)
print *, "M                  =", size(alpha)
print *, "K                  =", k_huber
print *, "delta              =", delta
print *, "tau                =", tau
if (n_keep > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
print *, "raw nodes          =", n_nodes

deallocate(y, interval_mat, alpha, beta, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx, cpt_hat, dval_hat)
call print_wall_time(t_start)

contains

    function derive_interval_file(data_file) result(interval_file)
        character(len=*), intent(in) :: data_file
        character(len=256) :: interval_file
        integer :: n
        n = len_trim(data_file)
        interval_file = trim(data_file)
        if (n >= 4 .and. data_file(n-3:n) == ".txt") then
            interval_file = data_file(:n-4) // "_intervals.txt"
        else
            interval_file = trim(data_file) // "_intervals.txt"
        end if
    end function derive_interval_file

    subroutine print_int_list(x)
        integer, intent(in) :: x(:)
        integer :: i
        do i = 1, size(x)
            write (*,'(1X,I0)', advance='no') x(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_wbs_univar_rob_file
