program xsim_dp_poly_file
! read univariate data from a text file and estimate polynomial changepoints by DP.poly
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_dp_poly_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: r = 2 ! polynomial order
real(kind=dp), parameter :: gamma = 1.0_dp ! l0 penalty parameter
integer, parameter :: delta = 5 ! minimum spacing

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), yhat(:)
integer, allocatable :: partition(:), cpt_hat(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xdp_poly_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)))
y = dat(:, 1)
call solve_dp_poly_1d(y, r, gamma, delta, partition, yhat, cpt_hat)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(y)
print *, "r                  = ", r
print *, "gamma              = ", gamma
print *, "delta              = ", delta
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          = "
end if
print *, "yhat checksum      = ", sum(yhat)

deallocate(dat, y, partition, cpt_hat, yhat)
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

end program xsim_dp_poly_file
