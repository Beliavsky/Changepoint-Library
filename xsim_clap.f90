program xsim_clap
! simulate a labeled univariate series and score a CLaP-style centroid classifier
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_ar1_1d
use claspy_mod, only: solve_clap_centroid_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 480 ! series length
integer, parameter :: window_size = 20 ! subsequence window size used for state classification
integer, parameter :: n_splits = 5 ! number of cross-validation folds
integer, parameter :: sample_cap = 1000 ! maximum retained windows per label

integer, parameter :: regime_starts(6) = [0, 80, 160, 240, 320, 400] ! starting indices of the simulated regimes
real(kind=dp), parameter :: phis(6) = [0.8_dp, -0.6_dp, 0.2_dp, 0.8_dp, -0.6_dp, 0.2_dp] ! AR(1) coefficients for the repeated states
real(kind=dp), parameter :: sigmas(6) = [0.8_dp, 0.8_dp, 0.8_dp, 0.8_dp, 0.8_dp, 0.8_dp] ! innovation standard deviations
integer, parameter :: regime_labels(6) = [1, 2, 3, 1, 2, 3] ! repeated state labels assigned to the regimes

real(kind=dp) :: x(n), score
integer :: i, i1, i2
integer :: state_labels(n)
integer, allocatable :: y_true(:), y_pred(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(211)
call simulate_piecewise_ar1_1d(n, regime_starts, phis, sigmas, x)

do i = 1, size(regime_starts)
    i1 = regime_starts(i) + 1
    if (i < size(regime_starts)) then
        i2 = regime_starts(i + 1)
    else
        i2 = n
    end if
    state_labels(i1:i2) = regime_labels(i)
end do

call solve_clap_centroid_1d(x, state_labels, window_size, n_splits, y_true, y_pred, score, sample_cap)

print *, "n                =", n
print *, "window_size      =", window_size
print *, "n_splits         =", n_splits
print *, "sample_cap       =", sample_cap
print *, "regime_starts    =", regime_starts
print *, "regime_labels    =", regime_labels
print *, "phis             =", phis
print *, "sigmas           =", sigmas
print '(A,*(1X,I0))', "true labels      =", y_true
print '(A,*(1X,I0))', "pred labels      =", y_pred
print *, "macro_f1         =", score

deallocate(y_true, y_pred)
call print_wall_time(t_start)
end program xsim_clap
