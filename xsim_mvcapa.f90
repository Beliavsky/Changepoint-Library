program xsim_mvcapa
! simulate multivariate zero-baseline data and detect subset anomalies by CAPA
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_mv
use changepoint_mod, only: solve_mvcapa_l2_mv
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 91 ! RNG seed used for the simulated data
integer, parameter :: n = 1200 ! number of observations in the simulated series
integer, parameter :: p = 3 ! number of variables
integer, parameter :: nbase = 1 ! number of baseline regimes
integer, parameter :: nseg_true = 2 ! number of collective anomalies
integer, parameter :: npoint_true = 1 ! number of point anomalies
integer, parameter :: baseline_starts(nbase) = [0] ! start indices of baseline regimes in 0-based notation
real(kind=dp), parameter :: baseline_means(nbase, p) = reshape([0.0_dp, 0.0_dp, 0.0_dp], [nbase, p]) ! baseline means by regime and variable
real(kind=dp), parameter :: baseline_sds(nbase, p) = reshape([0.7_dp, 0.7_dp, 0.7_dp], [nbase, p]) ! baseline standard deviations by regime and variable
integer, parameter :: true_seg_starts(nseg_true) = [200, 700] ! 0-based start indices of collective anomalies
integer, parameter :: true_seg_ends(nseg_true) = [260, 780] ! 0-based end indices of collective anomalies
integer, parameter :: true_seg_ncomponents(nseg_true) = [2, 1] ! number of affected variables in each collective anomaly
integer, parameter :: true_seg_components(nseg_true, p) = reshape([0, 1, 2, -1, -1, -1], [nseg_true, p]) ! affected variable indices for each collective anomaly
real(kind=dp), parameter :: true_seg_shifts(nseg_true) = [5.0_dp, -6.0_dp] ! additive mean shifts for collective anomalies
integer, parameter :: true_point_locs(npoint_true) = [950] ! 0-based point anomaly locations
integer, parameter :: true_point_ncomponents(npoint_true) = [2] ! number of affected variables in each point anomaly
integer, parameter :: true_point_components(npoint_true, p) = reshape([0, 1, -1], [npoint_true, p]) ! affected variable indices for each point anomaly
real(kind=dp), parameter :: true_point_shifts(npoint_true) = [8.0_dp] ! additive shifts for point anomalies
real(kind=dp), parameter :: segment_penalty(p) = [18.0_dp, 26.0_dp, 34.0_dp] ! CAPA penalty array for collective anomalies
real(kind=dp), parameter :: point_penalty(p) = [14.0_dp, 20.0_dp, 26.0_dp] ! CAPA penalty array for point anomalies
integer, parameter :: min_segment_length = 5 ! minimum allowed collective anomaly length
integer, parameter :: max_segment_length = 120 ! maximum allowed collective anomaly length

real(kind=dp), allocatable :: x(:, :)
integer, allocatable :: segment_starts(:), segment_ends(:), segment_components(:, :), segment_ncomponents(:)
integer, allocatable :: point_locs(:), point_components(:, :), point_ncomponents(:)
integer(kind=long_int) :: t_start
integer :: i, j

call system_clock(t_start)
allocate(x(n, p))
call seed_rng_fixed(seed)
call simulate_piecewise_normal_mv(n, p, baseline_starts, baseline_means, baseline_sds, x)

do i = 1, nseg_true
    do j = 1, true_seg_ncomponents(i)
        x(true_seg_starts(i) + 1:true_seg_ends(i), true_seg_components(i, j) + 1) = &
            x(true_seg_starts(i) + 1:true_seg_ends(i), true_seg_components(i, j) + 1) + true_seg_shifts(i)
    end do
end do
do i = 1, npoint_true
    do j = 1, true_point_ncomponents(i)
        x(true_point_locs(i) + 1, true_point_components(i, j) + 1) = &
            x(true_point_locs(i) + 1, true_point_components(i, j) + 1) + true_point_shifts(i)
    end do
end do

call solve_mvcapa_l2_mv(x, segment_penalty, point_penalty, segment_starts, segment_ends, segment_components, segment_ncomponents, &
                        point_locs, point_components, point_ncomponents, min_segment_length, max_segment_length)

print *, "n                =", n
print *, "p                =", p
print *, "segment penalty  =", segment_penalty
print *, "point penalty    =", point_penalty
write(*, '(A)', advance='no') 'true segments    ='
do i = 1, nseg_true
    write(*, '(" [", I0, ",", I0, "):")', advance='no') true_seg_starts(i), true_seg_ends(i)
    do j = 1, true_seg_ncomponents(i)
        write(*, '(" ", I0)', advance='no') true_seg_components(i, j)
    end do
end do
print *
write(*, '(A)', advance='no') 'true points      ='
do i = 1, npoint_true
    write(*, '(" [", I0, ",", I0, "):")', advance='no') true_point_locs(i), true_point_locs(i) + 1
    do j = 1, true_point_ncomponents(i)
        write(*, '(" ", I0)', advance='no') true_point_components(i, j)
    end do
end do
print *
write(*, '(A)', advance='no') 'segment anomalies='
if (size(segment_starts) > 0) then
    do i = 1, size(segment_starts)
        write(*, '(" [", I0, ",", I0, "):")', advance='no') segment_starts(i), segment_ends(i)
        do j = 1, segment_ncomponents(i)
            write(*, '(" ", I0)', advance='no') segment_components(i, j)
        end do
    end do
else
    write(*, '(A)', advance='no') ' none'
end if
print *
write(*, '(A)', advance='no') 'point anomalies  ='
if (size(point_locs) > 0) then
    do i = 1, size(point_locs)
        write(*, '(" [", I0, ",", I0, "):")', advance='no') point_locs(i), point_locs(i) + 1
        do j = 1, point_ncomponents(i)
            write(*, '(" ", I0)', advance='no') point_components(i, j)
        end do
    end do
else
    write(*, '(A)', advance='no') ' none'
end if
print *

deallocate(x, segment_starts, segment_ends, segment_components, segment_ncomponents, point_locs, point_components, point_ncomponents)
call print_wall_time(t_start)

end program xsim_mvcapa
