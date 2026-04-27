program xsim_clap_file
! read a labeled univariate series from a text file and score a CLaP-style centroid classifier
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_labels_file, get_data_file_arg
use claspy_mod, only: solve_clap_centroid_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: window_size = 20 ! subsequence window size used for state classification
integer, parameter :: n_splits = 5 ! number of cross-validation folds
integer, parameter :: sample_cap = 1000 ! maximum retained windows per label

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: state_labels(:), y_true(:), y_pred(:)
real(kind=dp) :: score
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xclap_data.txt", data_file)
call read_series_labels_file(data_file, x, state_labels)
call solve_clap_centroid_1d(x, state_labels, window_size, n_splits, y_true, y_pred, score, sample_cap)

print *, "file           =", trim(data_file)
print *, "n              =", size(x)
print *, "window_size    =", window_size
print *, "n_splits       =", n_splits
print *, "sample_cap     =", sample_cap
print '(A,*(1X,I0))', "true labels    =", y_true
print '(A,*(1X,I0))', "pred labels    =", y_pred
print *, "macro_f1       =", score

deallocate(x, state_labels, y_true, y_pred)
call print_wall_time(t_start)
end program xsim_clap_file
