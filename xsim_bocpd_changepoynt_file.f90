program xsim_bocpd_changepoynt_file
! read a univariate series from a text file and score it with changepoynt BOCPD
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoynt_mod, only: solve_bocpd_gaussian_mean_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: run_length = 120 ! constant-hazard reciprocal and sliding-window fit length

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:), score(:)
real(kind=dp) :: best_score, prior_mean, prior_var, signal_var
integer :: cp_raw, cp_post
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xbocpd_changepoynt_data.txt", data_file)
call read_series_file(data_file, x)
allocate(score(size(x)))
call solve_bocpd_gaussian_mean_1d(x, run_length, score, cp_raw, cp_post, best_score, prior_mean, prior_var, signal_var)

print *, "file               =", trim(data_file)
print *, "n                  =", size(x)
print *, "run_length         =", run_length
print *, "raw argmax         =", cp_raw
print *, "post_warmup argmax =", cp_post
print *, "max score          =", best_score
print *, "prior_mean         =", prior_mean
print *, "prior_var          =", prior_var
print *, "signal_var         =", signal_var

deallocate(score, x)
call print_wall_time(t_start)
end program xsim_bocpd_changepoynt_file
