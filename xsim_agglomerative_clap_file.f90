program xsim_agglomerative_clap_file
! read a repeated-state series from a file and run AgglomerativeCLaP-style merging
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, read_true_bkps_file, get_data_file_arg
use claspy_mod, only: solve_agglomerative_clap_1d, collapse_segment_process_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: window_size = 20 ! subsequence window size used by the deterministic CLaP classifier
integer, parameter :: n_splits = 5 ! number of cross-validation folds
integer, parameter :: sample_cap = 1000 ! maximum retained windows per label

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: gain
integer, allocatable :: true_bkps(:), segment_labels(:), cps_out(:), labels_out(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xagglomerative_clap_data.txt", data_file)
call read_series_file(data_file, x)
call read_true_bkps_file(data_file, true_bkps)
call solve_agglomerative_clap_1d(x, true_bkps, window_size, n_splits, segment_labels, gain, sample_cap)
call collapse_segment_process_1d(true_bkps, segment_labels, cps_out, labels_out)

print *, "file                =", trim(data_file)
print *, "n                   =", size(x)
print *, "true_bkps           =", true_bkps
print *, "window_size         =", window_size
print *, "n_splits            =", n_splits
print *, "sample_cap          =", sample_cap
print '(A,*(1X,I0))', "segment labels      =", segment_labels
print '(A,*(1X,I0))', "collapsed bkps      =", cps_out
print '(A,*(1X,I0))', "collapsed labels    =", labels_out
print *, "classification gain =", gain

deallocate(x, true_bkps, segment_labels, cps_out, labels_out)
call print_wall_time(t_start)
end program xsim_agglomerative_clap_file
