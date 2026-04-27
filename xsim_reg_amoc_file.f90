program xsim_reg_amoc_file
! read a regression dataset from a text file and estimate one changepoint by changepoint::cpt.reg(method="AMOC")
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoint_mod, only: solve_amoc_reg_norm_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: minseglen = 5 ! minimum segment length used by cpt.reg AMOC

character(len=256) :: data_file
real(kind=dp), allocatable :: data(:, :)
real(kind=dp) :: penalty_value, null_rss, alt_rss
integer :: cpt
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xreg_data.txt", data_file)
call read_matrix_file(data_file, data)
call solve_amoc_reg_norm_1d(data, cpt, penalty_value, null_rss, alt_rss, penalty="SIC", minseglen=minseglen)

print *, "file               =", trim(data_file)
print *, "n                  =", size(data, 1)
print *, "penalty            =", "SIC"
print *, "minseglen          =", minseglen
print *, "estimated          =", cpt
print *, "pen.value          =", penalty_value

deallocate(data)
call print_wall_time(t_start)

end program xsim_reg_amoc_file
