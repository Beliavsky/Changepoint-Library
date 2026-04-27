program xsim_sst_file
! read a univariate series from a text file and score it with changepoynt SST using the slow naive method
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoynt_mod, only: solve_sst_naive_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: window_length = 48 ! Hankel window length
integer, parameter :: n_windows = 48 ! number of Hankel columns
integer, parameter :: lag = 16 ! comparison lag between past and future windows
integer, parameter :: rank = 2 ! retained subspace rank
integer, parameter :: scoring_step = 1 ! score evaluation stride

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:), score(:)
real(kind=dp) :: best_score
integer :: cp
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xsst_data.txt", data_file)
call read_series_file(data_file, x)
allocate(score(size(x)))
call solve_sst_naive_1d(x, window_length, n_windows, lag, rank, score, cp, best_score, scoring_step)

print *, "file               =", trim(data_file)
print *, "n                  =", size(x)
print *, "window_length      =", window_length
print *, "n_windows          =", n_windows
print *, "lag                =", lag
print *, "rank               =", rank
print *, "scoring_step       =", scoring_step
print *, "method             = naive"
print *, "estimated          =", cp
print *, "max raw score      =", best_score

deallocate(score, x)
call print_wall_time(t_start)
end program xsim_sst_file
