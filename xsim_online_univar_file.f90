program xsim_online_univar_file
! read a univariate series and threshold vector from text files and detect the changepoint by online.univar
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoints_pkg_mod, only: solve_online_univar_1d
use util_mod, only: print_wall_time
implicit none

character(len=256) :: data_file, threshold_file
real(kind=dp), allocatable :: y(:), b_vec(:)
integer(kind=long_int) :: t_start
integer :: cpt_hat

call system_clock(t_start)
call get_data_file_arg("xonline_univar_data.txt", data_file)
threshold_file = derive_threshold_file(data_file)
call read_series_file(data_file, y)
call read_series_file(threshold_file, b_vec)
call solve_online_univar_1d(y, b_vec, cpt_hat)

print *, "file               =", trim(data_file)
print *, "threshold_file     =", trim(threshold_file)
print *, "n                  =", size(y)
print *, "estimated          =", cpt_hat

deallocate(y, b_vec)
call print_wall_time(t_start)

contains

    function derive_threshold_file(data_file) result(threshold_file)
        character(len=*), intent(in) :: data_file
        character(len=256) :: threshold_file
        integer :: n
        n = len_trim(data_file)
        threshold_file = trim(data_file)
        if (n >= 4 .and. data_file(n-3:n) == ".txt") then
            threshold_file = data_file(:n-4) // "_thresholds.txt"
        else
            threshold_file = trim(data_file) // "_thresholds.txt"
        end if
    end function derive_threshold_file

end program xsim_online_univar_file
