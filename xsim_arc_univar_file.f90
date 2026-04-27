program xsim_arc_univar_file
! read a contaminated univariate series from a text file and estimate changepoints by ARC
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_arc_univar_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: h = 120 ! block size used by ARC
integer, parameter :: block_num = 1 ! number of blocks used in local maximiser search
real(kind=dp), parameter :: epsilon = 0.1_dp ! contamination proportion
logical, parameter :: gaussian = .true. ! use Gaussian threshold calibration

character(len=256) :: data_file
real(kind=dp), allocatable :: y(:), cusum(:)
integer, allocatable :: cpt_hat(:)
real(kind=dp) :: lambda
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xarc_data.txt", data_file)
call read_series_file(data_file, y)
call solve_arc_univar_1d(y, h, block_num, epsilon, gaussian, cpt_hat, lambda, cusum)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "h                  =", h
print *, "block_num          =", block_num
print *, "epsilon            =", epsilon
print *, "gaussian           =", gaussian
print *, "lambda             =", lambda
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if

deallocate(y, cusum, cpt_hat)
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

end program xsim_arc_univar_file
