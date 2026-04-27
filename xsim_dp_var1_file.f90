program xsim_dp_var1_file
! read multivariate VAR(1) data and estimate changepoints by DP.VAR1 with lambda = 0
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_var1_mod, only: solve_dp_var1_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: gamma = 0.2_dp ! l0 penalty parameter
real(kind=dp), parameter :: lambda = 0.0_dp ! lasso penalty parameter
integer, parameter :: delta = 5 ! minimum spacing

character(len=256) :: data_file
real(kind=dp), allocatable :: data(:, :), data_t(:, :), x_curr(:, :), x_futu(:, :)
integer, allocatable :: partition(:), cpt_hat(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xdp_var1_data.txt", data_file)
call read_matrix_file(data_file, data)
allocate(data_t(size(data, 2), size(data, 1)))
data_t = transpose(data)
allocate(x_curr(size(data_t, 1), size(data_t, 2) - 1), x_futu(size(data_t, 1), size(data_t, 2) - 1))
x_curr = data_t(:, 1:size(data_t, 2) - 1)
x_futu = data_t(:, 2:size(data_t, 2))
call solve_dp_var1_1d(x_futu, x_curr, gamma, lambda, delta, partition, cpt_hat)

print *, "file               = ", trim(data_file)
print *, "p                  = ", size(data_t, 1)
print *, "n transitions      = ", size(x_curr, 2)
print *, "gamma              = ", gamma
print *, "lambda             = ", lambda
print *, "delta              = ", delta
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          = "
end if
print *, "partition checksum = ", sum(partition)

deallocate(data, data_t, x_curr, x_futu, partition, cpt_hat)
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

end program xsim_dp_var1_file
