program xsim_meanvar_exp_amoc_file
! read a univariate Exponential series from a text file and estimate one changepoint by AMOC
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoint_mod, only: solve_amoc_meanvar_exp_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: minseglen = 1 ! minimum segment length used by Exponential AMOC

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer :: cpt
real(kind=dp) :: p_value, null_like, alt_like, decision_penalty, alt_mbic
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xexp_data.txt", data_file)
call read_series_file(data_file, x)
call solve_amoc_meanvar_exp_1d(x, cpt, p_value, null_like, alt_like, decision_penalty, penalty="MBIC", minseglen=minseglen)
if (cpt < size(x)) then
    alt_mbic = alt_like + log(real(cpt, dp)) + log(real(size(x) - cpt + 1, dp))
else
    alt_mbic = alt_like
end if

print *, "file               =", trim(data_file)
print *, "n                  =", size(x)
print *, "penalty            =", "MBIC"
print *, "minseglen          =", minseglen
print *, "estimated          =", cpt
print *, "p.value            =", p_value
print *, "null               =", null_like
print *, "alt                =", alt_like
print *, "alt.mbic           =", alt_mbic
print *, "decision pen.value =", decision_penalty

deallocate(x)
call print_wall_time(t_start)

end program xsim_meanvar_exp_amoc_file
