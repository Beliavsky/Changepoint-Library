program xsim_aarc_univar_file
! read a contaminated univariate series from a text file and estimate changepoints by aARC
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_aarc_univar_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: t_dat_n = 200 ! number of training observations used by aARC
integer, parameter :: h = 120 ! block size used by aARC
integer, parameter :: block_num = 1 ! number of blocks used in local maximiser search
real(kind=dp), parameter :: guess_true = 0.05_dp ! preferred contamination root

character(len=256) :: data_file
real(kind=dp), allocatable :: y(:)
integer, allocatable :: cpt_hat(:)
real(kind=dp) :: eps_hat, arc_epsilon, lambda
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xaarc_data.txt", data_file)
call read_series_file(data_file, y)
call solve_aarc_univar_1d(y, y(1:t_dat_n), guess_true, h, block_num, cpt_hat, eps_hat, arc_epsilon, lambda)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "t_dat_n            =", t_dat_n
print *, "h                  =", h
print *, "block_num          =", block_num
print *, "guess_true         =", guess_true
print *, "eps_hat            =", eps_hat
print *, "arc epsilon        =", arc_epsilon
print *, "lambda             =", lambda
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if

deallocate(y, cpt_hat)
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

end program xsim_aarc_univar_file
