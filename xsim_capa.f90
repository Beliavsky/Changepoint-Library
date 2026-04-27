program xsim_capa
! simulate a zero-baseline Gaussian series and detect collective and point anomalies by CAPA
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoint_mod, only: solve_capa_l2_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 81 ! RNG seed used for the simulated series
integer, parameter :: n = 5000 ! number of observations in the simulated series
integer, parameter :: nbase = 1 ! number of baseline Gaussian regimes
integer, parameter :: nseg_true = 3 ! number of injected collective anomalies
integer, parameter :: npoint_true = 1 ! number of injected point anomalies
integer, parameter :: baseline_starts(nbase) = [0] ! start index of the baseline regime in 0-based notation
real(kind=dp), parameter :: baseline_means(nbase) = [0.0_dp] ! baseline mean level
real(kind=dp), parameter :: baseline_sds(nbase) = [0.7_dp] ! baseline standard deviation
integer, parameter :: true_seg_starts(nseg_true) = [1000, 2800, 4000] ! 0-based collective anomaly start indices
integer, parameter :: true_seg_ends(nseg_true) = [1100, 2950, 4200] ! 0-based collective anomaly end indices
real(kind=dp), parameter :: true_seg_shifts(nseg_true) = [5.0_dp, -4.5_dp, 6.0_dp] ! mean shifts added on anomalous segments
integer, parameter :: true_point_locs(npoint_true) = [2000] ! 0-based point anomaly locations
real(kind=dp), parameter :: true_point_shifts(npoint_true) = [9.0_dp] ! additive shifts at point anomalies
real(kind=dp), parameter :: segment_penalty = 24.0_dp ! CAPA penalty for collective anomalies
real(kind=dp), parameter :: point_penalty = 14.0_dp ! CAPA penalty for point anomalies
integer, parameter :: min_segment_length = 10 ! minimum allowed collective anomaly length
integer, parameter :: max_segment_length = 500 ! maximum allowed collective anomaly length

real(kind=dp), allocatable :: x(:)
integer, allocatable :: segment_starts(:), segment_ends(:), point_locs(:)
integer(kind=long_int) :: t_start
integer :: i

call system_clock(t_start)
allocate(x(n))
call seed_rng_fixed(seed)
call simulate_piecewise_normal_1d(n, baseline_starts, baseline_means, baseline_sds, x)
do i = 1, nseg_true
    x(true_seg_starts(i) + 1:true_seg_ends(i)) = x(true_seg_starts(i) + 1:true_seg_ends(i)) + true_seg_shifts(i)
end do
do i = 1, npoint_true
    x(true_point_locs(i) + 1) = x(true_point_locs(i) + 1) + true_point_shifts(i)
end do

call solve_capa_l2_1d(x, segment_penalty, point_penalty, segment_starts, segment_ends, point_locs, min_segment_length, max_segment_length)

print *, "n                =", n
print *, "segment penalty  =", segment_penalty
print *, "point penalty    =", point_penalty
write(*, '(A)', advance='no') 'true segments    ='
do i = 1, nseg_true
    write(*, '(" [", I0, ",", I0, ")")', advance='no') true_seg_starts(i), true_seg_ends(i)
end do
print *
write(*, '(A)', advance='no') 'true points      ='
do i = 1, npoint_true
    write(*, '(" [", I0, ",", I0, ")")', advance='no') true_point_locs(i), true_point_locs(i) + 1
end do
print *
write(*, '(A)', advance='no') 'segment anomalies='
if (size(segment_starts) > 0) then
    do i = 1, size(segment_starts)
        write(*, '(" [", I0, ",", I0, ")")', advance='no') segment_starts(i), segment_ends(i)
    end do
else
    write(*, '(A)', advance='no') ' none'
end if
print *
write(*, '(A)', advance='no') 'point anomalies  ='
if (size(point_locs) > 0) then
    do i = 1, size(point_locs)
        write(*, '(" [", I0, ",", I0, ")")', advance='no') point_locs(i), point_locs(i) + 1
    end do
else
    write(*, '(A)', advance='no') ' none'
end if
print *

deallocate(x, segment_starts, segment_ends, point_locs)
call print_wall_time(t_start)

end program xsim_capa
