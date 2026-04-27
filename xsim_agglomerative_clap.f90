program xsim_agglomerative_clap
! simulate a repeated-state series and run AgglomerativeCLaP-style merging
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use claspy_mod, only: solve_agglomerative_clap_1d, collapse_segment_process_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 480 ! series length
integer, parameter :: window_size = 20 ! subsequence window size used by the deterministic CLaP classifier
integer, parameter :: n_splits = 5 ! number of cross-validation folds
integer, parameter :: sample_cap = 1000 ! maximum retained windows per label

integer, parameter :: regime_starts(6) = [0, 80, 160, 240, 320, 400] ! starting indices of the simulated regimes
integer, parameter :: true_bkps(5) = [80, 160, 240, 320, 400] ! initial segment boundaries passed to the agglomerative merger
integer, parameter :: true_states(6) = [1, 2, 3, 1, 2, 3] ! true recurring-state identities
real(kind=dp), parameter :: means(6) = [0.0_dp, 2.0_dp, -2.0_dp, 0.0_dp, 2.0_dp, -2.0_dp] ! regime means
real(kind=dp), parameter :: sds(6) = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp] ! regime standard deviations

real(kind=dp) :: x(n), gain
integer, allocatable :: segment_labels(:), cps_out(:), labels_out(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(101)
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)
call solve_agglomerative_clap_1d(x, true_bkps, window_size, n_splits, segment_labels, gain, sample_cap)
call collapse_segment_process_1d(true_bkps, segment_labels, cps_out, labels_out)

print *, "n                   =", n
print *, "true_bkps           =", true_bkps
print *, "true states         =", true_states
print *, "window_size         =", window_size
print *, "n_splits            =", n_splits
print *, "sample_cap          =", sample_cap
print '(A,*(1X,I0))', "segment labels      =", segment_labels
print '(A,*(1X,I0))', "collapsed bkps      =", cps_out
print '(A,*(1X,I0))', "collapsed labels    =", labels_out
print *, "classification gain =", gain

deallocate(segment_labels, cps_out, labels_out)
call print_wall_time(t_start)
end program xsim_agglomerative_clap
