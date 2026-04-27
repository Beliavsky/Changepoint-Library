program xsim_wbsip_cov_file
! read a multivariate series, its independent copy, and WBS intervals from text files, then estimate covariance changepoints by WBSIP.cov
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_wbsip_cov_1d, threshold_wbs_tree_1d
use util_mod, only: sort_int, print_wall_time
implicit none

integer, parameter :: delta = 5 ! minimum spacing parameter

character(len=256) :: data_file, prime_file, interval_file
real(kind=dp), allocatable :: x_raw(:, :), xprime_raw(:, :), interval_mat(:, :), x(:, :), x_prime(:, :), dval_nodes(:), dval_hat(:)
integer, allocatable :: alpha(:), beta(:), s_nodes(:), level_nodes(:), parent_nodes(:, :), parent_idx(:), cpt_hat(:)
integer :: n_nodes, n_keep, i
real(kind=dp), parameter :: tau = 8.0_dp ! trimming threshold
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xwbsip_cov_data.txt", data_file)
prime_file = derive_prime_file(data_file)
interval_file = derive_interval_file(data_file)
call read_matrix_file(data_file, x_raw)
call read_matrix_file(prime_file, xprime_raw)
call read_matrix_file(interval_file, interval_mat)
allocate(x(size(x_raw, 2), size(x_raw, 1)), x_prime(size(xprime_raw, 2), size(xprime_raw, 1)))
do i = 1, size(x_raw, 1)
    x(:, i) = x_raw(i, :)
    x_prime(:, i) = xprime_raw(i, :)
end do
allocate(alpha(size(interval_mat, 1)), beta(size(interval_mat, 1)))
alpha = nint(interval_mat(:, 1))
beta = nint(interval_mat(:, 2))
call solve_wbsip_cov_1d(x, x_prime, alpha, beta, 1, size(x, 2), delta, n_nodes, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx)
call threshold_wbs_tree_1d(n_nodes, s_nodes, dval_nodes, parent_idx, tau, n_keep, cpt_hat, dval_hat)
if (n_keep > 0) call sort_int(cpt_hat)

print *, "file               =", trim(data_file)
print *, "prime_file         =", trim(prime_file)
print *, "interval_file      =", trim(interval_file)
print *, "n                  =", size(x, 2)
print *, "p                  =", size(x, 1)
print *, "M                  =", size(alpha)
print *, "delta              =", delta
print *, "tau                =", tau
if (n_keep > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
print *, "raw nodes          =", n_nodes

deallocate(x_raw, xprime_raw, interval_mat, x, x_prime, alpha, beta, s_nodes, dval_nodes, level_nodes, parent_nodes, parent_idx, cpt_hat, dval_hat)
call print_wall_time(t_start)

contains

    function derive_prime_file(data_file_in) result(prime_file_out)
        character(len=*), intent(in) :: data_file_in
        character(len=256) :: prime_file_out
        integer :: ntrim
        ntrim = len_trim(data_file_in)
        prime_file_out = trim(data_file_in)
        if (ntrim >= 4 .and. data_file_in(ntrim-3:ntrim) == ".txt") then
            prime_file_out = data_file_in(:ntrim-4) // "_prime.txt"
        else
            prime_file_out = trim(data_file_in) // "_prime.txt"
        end if
    end function derive_prime_file

    function derive_interval_file(data_file_in) result(interval_file_out)
        character(len=*), intent(in) :: data_file_in
        character(len=256) :: interval_file_out
        integer :: ntrim
        ntrim = len_trim(data_file_in)
        interval_file_out = trim(data_file_in)
        if (ntrim >= 4 .and. data_file_in(ntrim-3:ntrim) == ".txt") then
            interval_file_out = data_file_in(:ntrim-4) // "_intervals.txt"
        else
            interval_file_out = trim(data_file_in) // "_intervals.txt"
        end if
    end function derive_interval_file

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: j
        do j = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(j)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_wbsip_cov_file
