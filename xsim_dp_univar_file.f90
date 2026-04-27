program xsim_dp_univar_file
! read a univariate series from a text file and estimate changepoints by l0-penalized dynamic programming
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_dp_univar_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: delta = 5 ! minimum spacing parameter in DP.univar
real(kind=dp), parameter :: gamma = 5.0_dp ! l0 penalty used by DP.univar

character(len=256) :: data_file
real(kind=dp), allocatable :: y(:), yhat(:)
integer, allocatable :: partition(:), cpt_hat(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xdp_univar_data.txt", data_file)
call read_series_file(data_file, y)
call solve_dp_univar_1d(y, gamma, delta, partition, yhat, cpt_hat)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "gamma              =", gamma
print *, "delta              =", delta
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
print *, "n changepoints     =", size(cpt_hat)
print *, "partition last     =", partition(size(partition))
print *, "yhat checksum      =", sum(yhat)

deallocate(y, partition, yhat, cpt_hat)
call print_wall_time(t_start)

contains

    subroutine print_int_list(x)
        integer, intent(in) :: x(:)
        integer :: i
        do i = 1, size(x)
            write (*,'(1X,I0)', advance='no') x(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_dp_univar_file
