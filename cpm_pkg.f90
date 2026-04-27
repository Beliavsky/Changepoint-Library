module cpm_pkg_mod
use kind_mod, only: dp
implicit none
private

integer, parameter :: cpm_student_kind = 1
integer, parameter :: cpm_bartlett_kind = 2
integer, parameter :: cpm_joint_kind = 3
integer, parameter :: cpm_exponential_kind = 4
integer, parameter :: cpm_poisson_kind = 5
integer, parameter :: cpm_mw_kind = 6
integer, parameter :: cpm_mood_kind = 7
integer, parameter :: cpm_lepage_kind = 8
integer, parameter :: cpm_fet_kind = 9
integer, parameter :: cpm_ks_kind = 10
integer, parameter :: cpm_cvm_kind = 11
integer, parameter :: cpm_joint_adjusted_kind = 12
integer, parameter :: cpm_exponential_adjusted_kind = 13
integer, parameter :: cpm_joint_hawkins_kind = 14

public :: solve_cpm_detect_student_1d, solve_cpm_process_student_1d, solve_cpm_batch_student_1d
public :: solve_cpm_detect_bartlett_1d, solve_cpm_process_bartlett_1d, solve_cpm_batch_bartlett_1d
public :: solve_cpm_detect_joint_1d, solve_cpm_process_joint_1d, solve_cpm_batch_joint_1d
public :: solve_cpm_detect_exponential_1d, solve_cpm_process_exponential_1d, solve_cpm_batch_exponential_1d
public :: solve_cpm_detect_joint_adjusted_1d, solve_cpm_process_joint_adjusted_1d, solve_cpm_batch_joint_adjusted_1d
public :: solve_cpm_detect_exponential_adjusted_1d, solve_cpm_process_exponential_adjusted_1d, solve_cpm_batch_exponential_adjusted_1d
public :: solve_cpm_detect_joint_hawkins_1d, solve_cpm_batch_joint_hawkins_1d
public :: solve_cpm_detect_poisson_1d, solve_cpm_process_poisson_1d, solve_cpm_batch_poisson_1d
public :: solve_cpm_detect_mw_1d, solve_cpm_process_mw_1d, solve_cpm_batch_mw_1d
public :: solve_cpm_detect_mood_1d, solve_cpm_process_mood_1d, solve_cpm_batch_mood_1d
public :: solve_cpm_detect_lepage_1d, solve_cpm_process_lepage_1d, solve_cpm_batch_lepage_1d
public :: solve_cpm_detect_fet_1d, solve_cpm_process_fet_1d, solve_cpm_batch_fet_1d
public :: solve_cpm_detect_ks_1d, solve_cpm_process_ks_1d, solve_cpm_batch_ks_1d
public :: solve_cpm_detect_cvm_1d, solve_cpm_process_cvm_1d, solve_cpm_batch_cvm_1d

contains

subroutine solve_cpm_detect_student_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_student_kind, ds, cp, dt)
end subroutine solve_cpm_detect_student_1d

subroutine solve_cpm_process_student_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_student_kind, cps, dts)
end subroutine solve_cpm_process_student_1d

subroutine solve_cpm_batch_student_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_student_kind, ds, cp)
end subroutine solve_cpm_batch_student_1d

subroutine solve_cpm_detect_bartlett_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_bartlett_kind, ds, cp, dt)
end subroutine solve_cpm_detect_bartlett_1d

subroutine solve_cpm_process_bartlett_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_bartlett_kind, cps, dts)
end subroutine solve_cpm_process_bartlett_1d

subroutine solve_cpm_batch_bartlett_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_bartlett_kind, ds, cp)
end subroutine solve_cpm_batch_bartlett_1d

subroutine solve_cpm_detect_joint_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_joint_kind, ds, cp, dt)
end subroutine solve_cpm_detect_joint_1d

subroutine solve_cpm_process_joint_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_joint_kind, cps, dts)
end subroutine solve_cpm_process_joint_1d

subroutine solve_cpm_batch_joint_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_joint_kind, ds, cp)
end subroutine solve_cpm_batch_joint_1d

subroutine solve_cpm_detect_joint_adjusted_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_joint_adjusted_kind, ds, cp, dt)
end subroutine solve_cpm_detect_joint_adjusted_1d

subroutine solve_cpm_process_joint_adjusted_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_joint_adjusted_kind, cps, dts)
end subroutine solve_cpm_process_joint_adjusted_1d

subroutine solve_cpm_batch_joint_adjusted_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_joint_adjusted_kind, ds, cp)
end subroutine solve_cpm_batch_joint_adjusted_1d

subroutine solve_cpm_detect_exponential_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_exponential_kind, ds, cp, dt)
end subroutine solve_cpm_detect_exponential_1d

subroutine solve_cpm_process_exponential_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_exponential_kind, cps, dts)
end subroutine solve_cpm_process_exponential_1d

subroutine solve_cpm_batch_exponential_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_exponential_kind, ds, cp)
end subroutine solve_cpm_batch_exponential_1d

subroutine solve_cpm_detect_exponential_adjusted_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_exponential_adjusted_kind, ds, cp, dt)
end subroutine solve_cpm_detect_exponential_adjusted_1d

subroutine solve_cpm_process_exponential_adjusted_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_exponential_adjusted_kind, cps, dts)
end subroutine solve_cpm_process_exponential_adjusted_1d

subroutine solve_cpm_batch_exponential_adjusted_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_exponential_adjusted_kind, ds, cp)
end subroutine solve_cpm_batch_exponential_adjusted_1d

subroutine solve_cpm_detect_joint_hawkins_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_joint_hawkins_kind, ds, cp, dt)
end subroutine solve_cpm_detect_joint_hawkins_1d

subroutine solve_cpm_batch_joint_hawkins_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_joint_hawkins_kind, ds, cp)
end subroutine solve_cpm_batch_joint_hawkins_1d

subroutine solve_cpm_detect_poisson_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_poisson_kind, ds, cp, dt)
end subroutine solve_cpm_detect_poisson_1d

subroutine solve_cpm_process_poisson_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_poisson_kind, cps, dts)
end subroutine solve_cpm_process_poisson_1d

subroutine solve_cpm_batch_poisson_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_poisson_kind, ds, cp)
end subroutine solve_cpm_batch_poisson_1d

subroutine solve_cpm_detect_mw_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_mw_kind, ds, cp, dt)
end subroutine solve_cpm_detect_mw_1d

subroutine solve_cpm_process_mw_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_mw_kind, cps, dts)
end subroutine solve_cpm_process_mw_1d

subroutine solve_cpm_batch_mw_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_mw_kind, ds, cp)
end subroutine solve_cpm_batch_mw_1d

subroutine solve_cpm_detect_mood_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_mood_kind, ds, cp, dt)
end subroutine solve_cpm_detect_mood_1d

subroutine solve_cpm_process_mood_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_mood_kind, cps, dts)
end subroutine solve_cpm_process_mood_1d

subroutine solve_cpm_batch_mood_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_mood_kind, ds, cp)
end subroutine solve_cpm_batch_mood_1d

subroutine solve_cpm_detect_lepage_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_lepage_kind, ds, cp, dt)
end subroutine solve_cpm_detect_lepage_1d

subroutine solve_cpm_process_lepage_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_lepage_kind, cps, dts)
end subroutine solve_cpm_process_lepage_1d

subroutine solve_cpm_batch_lepage_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_lepage_kind, ds, cp)
end subroutine solve_cpm_batch_lepage_1d

subroutine solve_cpm_detect_fet_1d(x, thresholds, startup, lambda, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:), lambda
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_fet_kind, ds, cp, dt, lambda)
end subroutine solve_cpm_detect_fet_1d

subroutine solve_cpm_process_fet_1d(x, thresholds, startup, lambda, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:), lambda
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_fet_kind, cps, dts, lambda)
end subroutine solve_cpm_process_fet_1d

subroutine solve_cpm_batch_fet_1d(x, threshold, lambda, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold, lambda
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_fet_kind, ds, cp, lambda)
end subroutine solve_cpm_batch_fet_1d

subroutine solve_cpm_detect_ks_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_ks_kind, ds, cp, dt)
end subroutine solve_cpm_detect_ks_1d

subroutine solve_cpm_process_ks_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_ks_kind, cps, dts)
end subroutine solve_cpm_process_ks_1d

subroutine solve_cpm_batch_ks_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_ks_kind, ds, cp)
end subroutine solve_cpm_batch_ks_1d

subroutine solve_cpm_detect_cvm_1d(x, thresholds, startup, ds, cp, dt)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt
call solve_cpm_detect_1d(x, thresholds, startup, cpm_cvm_kind, ds, cp, dt)
end subroutine solve_cpm_detect_cvm_1d

subroutine solve_cpm_process_cvm_1d(x, thresholds, startup, cps, dts)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup
integer, allocatable, intent(out) :: cps(:), dts(:)
call solve_cpm_process_1d(x, thresholds, startup, cpm_cvm_kind, cps, dts)
end subroutine solve_cpm_process_cvm_1d

subroutine solve_cpm_batch_cvm_1d(x, threshold, ds, cp)
real(kind=dp), intent(in) :: x(:), threshold
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp
call solve_cpm_batch_1d(x, threshold, cpm_cvm_kind, ds, cp)
end subroutine solve_cpm_batch_cvm_1d

subroutine solve_cpm_detect_1d(x, thresholds, startup, stat_kind, ds, cp, dt, lambda)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup, stat_kind
real(kind=dp), intent(in), optional :: lambda
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp, dt

real(kind=dp), allocatable :: s(:), w(:)
real(kind=dp) :: u, threshold
integer :: n_state, i, maxind

allocate(ds(size(x)))
allocate(s(size(x)), w(size(x)))
ds = 0.0_dp
call reset_state_1d(n_state)

do i = 1, size(x)
    call process_point_1d(x(i), stat_kind, s, w, n_state)
    if (n_state >= startup) then
        call cpm_mle_max_1d(stat_kind, s(1:n_state), w(1:n_state), n_state, u, maxind, lambda)
        ds(i) = u
        threshold = get_threshold_1d(thresholds, n_state, 9999999.0_dp)
        if (u > threshold) then
            cp = maxind + 1
            dt = i
            return
        end if
    end if
end do

cp = 0
dt = 0
end subroutine solve_cpm_detect_1d

subroutine solve_cpm_process_1d(x, thresholds, startup, stat_kind, cps, dts, lambda)
real(kind=dp), intent(in) :: x(:), thresholds(:)
integer, intent(in) :: startup, stat_kind
real(kind=dp), intent(in), optional :: lambda
integer, allocatable, intent(out) :: cps(:), dts(:)

real(kind=dp), allocatable :: s(:), w(:)
real(kind=dp) :: u, threshold
integer, allocatable :: cps_work(:), dts_work(:)
integer :: n_state, i, count_changes, maxind, cp, dt, last_change

allocate(s(size(x)), w(size(x)))
allocate(cps_work(size(x)), dts_work(size(x)))
call reset_state_1d(n_state)
count_changes = 0
last_change = 0
i = 1

do while (i <= size(x))
    call process_point_1d(x(i), stat_kind, s, w, n_state)
    if (n_state >= startup) then
        call cpm_mle_max_1d(stat_kind, s(1:n_state), w(1:n_state), n_state, u, maxind, lambda)
        threshold = get_threshold_1d(thresholds, n_state, 99999.0_dp)
        if (u > threshold) then
            dt = i
            cp = last_change + maxind + 1
            count_changes = count_changes + 1
            cps_work(count_changes) = cp
            dts_work(count_changes) = dt
            last_change = cp
            i = cp + 1
            call reset_state_1d(n_state)
            cycle
        end if
    end if
    i = i + 1
end do

if (count_changes > 0) then
    allocate(cps(count_changes), dts(count_changes))
    cps = cps_work(1:count_changes)
    dts = dts_work(1:count_changes)
else
    allocate(cps(0), dts(0))
end if
end subroutine solve_cpm_process_1d

subroutine solve_cpm_batch_1d(x, threshold, stat_kind, ds, cp, lambda)
real(kind=dp), intent(in) :: x(:), threshold
integer, intent(in) :: stat_kind
real(kind=dp), intent(in), optional :: lambda
real(kind=dp), allocatable, intent(out) :: ds(:)
integer, intent(out) :: cp

real(kind=dp), allocatable :: s(:), w(:)
integer :: n_state, i, imax(1)

allocate(s(size(x)), w(size(x)))
call reset_state_1d(n_state)
do i = 1, size(x)
    call process_point_1d(x(i), stat_kind, s, w, n_state)
end do

 call cpm_mle_vector_1d(stat_kind, s(1:n_state), w(1:n_state), n_state, ds, lambda)
imax = maxloc(ds)
cp = imax(1)
if (maxval(ds) <= threshold) cp = 0
end subroutine solve_cpm_batch_1d

subroutine process_point_t_1d(obs, s, w, n_state)
real(kind=dp), intent(in) :: obs
real(kind=dp), intent(inout) :: s(:), w(:)
integer, intent(inout) :: n_state
real(kind=dp) :: temp

n_state = n_state + 1
if (n_state == 1) then
    s(1) = obs
    w(1) = 0.0_dp
else
    s(n_state) = obs + s(n_state - 1)
    temp = real(n_state - 1, dp) * obs - s(n_state - 1)
    w(n_state) = w(n_state - 1) + temp * temp / (real(n_state, dp) * real(n_state - 1, dp))
end if
end subroutine process_point_t_1d

subroutine process_point_sum_1d(obs, s, w, n_state)
real(kind=dp), intent(in) :: obs
real(kind=dp), intent(inout) :: s(:), w(:)
integer, intent(inout) :: n_state
n_state = n_state + 1
if (n_state == 1) then
    s(1) = obs
else
    s(n_state) = s(n_state - 1) + obs
end if
w(n_state) = 0.0_dp
end subroutine process_point_sum_1d

subroutine process_point_rank_1d(obs, s, w, n_state)
real(kind=dp), intent(in) :: obs
real(kind=dp), intent(inout) :: s(:), w(:)
integer, intent(inout) :: n_state
integer :: i, tie_count
real(kind=dp) :: rank
integer, allocatable :: ties(:)

n_state = n_state + 1
if (n_state == 1) then
    s(1) = obs
    w(1) = 1.0_dp
    return
end if

allocate(ties(n_state - 1))
rank = 1.0_dp
tie_count = 0
do i = 1, n_state - 1
    if (s(i) > obs) then
        w(i) = w(i) + 1.0_dp
    else if (obs > s(i)) then
        rank = rank + 1.0_dp
    else
        tie_count = tie_count + 1
        ties(tie_count) = i
    end if
end do
if (tie_count > 0) then
    rank = (rank + rank + real(tie_count, dp)) / 2.0_dp
    do i = 1, tie_count
        w(ties(i)) = rank
    end do
end if
s(n_state) = obs
w(n_state) = rank
deallocate(ties)
end subroutine process_point_rank_1d

subroutine process_point_1d(obs, stat_kind, s, w, n_state)
real(kind=dp), intent(in) :: obs
integer, intent(in) :: stat_kind
real(kind=dp), intent(inout) :: s(:), w(:)
integer, intent(inout) :: n_state
select case (stat_kind)
case (cpm_student_kind, cpm_bartlett_kind, cpm_joint_kind, cpm_joint_adjusted_kind, cpm_joint_hawkins_kind)
    call process_point_t_1d(obs, s, w, n_state)
case (cpm_exponential_kind, cpm_poisson_kind, cpm_fet_kind, cpm_exponential_adjusted_kind)
    call process_point_sum_1d(obs, s, w, n_state)
case (cpm_mw_kind, cpm_mood_kind, cpm_lepage_kind)
    call process_point_rank_1d(obs, s, w, n_state)
case (cpm_ks_kind, cpm_cvm_kind)
    call process_point_order_1d(obs, s, w, n_state)
end select
end subroutine process_point_1d

subroutine process_point_order_1d(obs, s, w, n_state)
real(kind=dp), intent(in) :: obs
real(kind=dp), intent(inout) :: s(:), w(:)
integer, intent(inout) :: n_state
integer :: sz, rank, i, j

sz = n_state
n_state = n_state + 1
if (sz == 0) then
    s(1) = obs
    w(1) = 1.0_dp
    return
end if

rank = 0
do i = 1, sz
    if (obs > s(i)) rank = rank + 1
end do

if (rank == sz) then
    w(sz + 1) = real(sz + 1, dp)
else
    do j = sz + 1, rank + 2, -1
        w(j) = w(j - 1)
    end do
    w(rank + 1) = real(sz + 1, dp)
end if
s(sz + 1) = obs
end subroutine process_point_order_1d

subroutine reset_state_1d(n_state)
integer, intent(out) :: n_state
n_state = 0
end subroutine reset_state_1d

real(kind=dp) function get_threshold_1d(thresholds, n_state, default_value)
real(kind=dp), intent(in) :: thresholds(:), default_value
integer, intent(in) :: n_state
if (size(thresholds) <= 0) then
    get_threshold_1d = default_value
else if (n_state >= size(thresholds)) then
    get_threshold_1d = thresholds(size(thresholds))
else
    get_threshold_1d = thresholds(n_state)
end if
end function get_threshold_1d

subroutine cpm_mle_max_1d(stat_kind, s, w, n_state, maxvalue, maxindex, lambda)
integer, intent(in) :: stat_kind, n_state
real(kind=dp), intent(in) :: s(:), w(:)
real(kind=dp), intent(in), optional :: lambda
real(kind=dp), intent(out) :: maxvalue
integer, intent(out) :: maxindex
real(kind=dp), allocatable :: us(:)
integer :: i

call cpm_mle_vector_1d(stat_kind, s, w, n_state, us, lambda)
maxvalue = 0.0_dp
maxindex = 0
do i = 2, size(us) - 2
    if (us(i) > maxvalue) then
        maxvalue = us(i)
        maxindex = i - 1
    end if
end do
end subroutine cpm_mle_max_1d

subroutine cpm_mle_vector_1d(stat_kind, s, w, n_state, us, lambda)
integer, intent(in) :: stat_kind, n_state
real(kind=dp), intent(in) :: s(:), w(:)
real(kind=dp), intent(in), optional :: lambda
real(kind=dp), allocatable, intent(out) :: us(:)

integer :: i, n1, n2
real(kind=dp) :: j, nn, sigma, temp, e, denom
real(kind=dp) :: mu1, mu2, v1, v2, sigma1, sigma2, sigma_all, c, g
real(kind=dp) :: sok, son, skn

allocate(us(n_state))
us = 0.0_dp
nn = real(n_state, dp)

select case (stat_kind)
case (cpm_student_kind)
    sigma = sqrt((nn - 2.0_dp) / (nn - 4.0_dp))
    do i = 2, n_state - 2
        j = real(i, dp)
        temp = nn * s(i) - j * s(n_state)
        e = temp * temp / (nn * j * (nn - j))
        denom = w(n_state) - e
        if (denom > 0.0_dp .and. e > 0.0_dp) us(i) = sqrt((nn - 2.0_dp) * e / denom) / sigma
    end do
case (cpm_bartlett_kind)
    do i = 2, n_state - 2
        n1 = i
        n2 = n_state - n1
        mu1 = s(n1) / real(n1, dp)
        mu2 = (s(n_state) - s(n1)) / real(n2, dp)
        v1 = w(n1)
        v2 = w(n_state) - w(n1) - real(n1 * n2, dp) * (mu1 - mu2) * (mu1 - mu2) / real(n_state, dp)
        sigma1 = v1 / real(n1 - 1, dp)
        sigma2 = v2 / real(n2 - 1, dp)
        sigma_all = (v1 + v2) / real(n_state - 2, dp)
        c = 1.0_dp + (1.0_dp / real(n1 - 1, dp) + 1.0_dp / real(n2 - 1, dp) - 1.0_dp / real(n_state - 2, dp)) / 3.0_dp
        if (sigma1 > 0.0_dp .and. sigma2 > 0.0_dp .and. sigma_all > 0.0_dp) then
            g = (real(n1 - 1, dp) * log(sigma_all / sigma1) + real(n2 - 1, dp) * log(sigma_all / sigma2)) / c
            us(i) = g
        end if
    end do
case (cpm_joint_kind)
    do i = 2, n_state - 2
        n1 = i
        n2 = n_state - n1
        mu1 = s(n1) / real(n1, dp)
        mu2 = (s(n_state) - s(n1)) / real(n2, dp)
        sok = w(n1) / real(n1, dp)
        son = w(n_state) / real(n_state, dp)
        skn = w(n_state) - w(n1) - real(n1 * n2, dp) * (mu1 - mu2) * (mu1 - mu2) / real(n_state, dp)
        skn = skn / real(n2, dp)
        c = 1.0_dp + 11.0_dp / 12.0_dp * (1.0_dp / real(n1, dp) + 1.0_dp / real(n2, dp) - 1.0_dp / real(n_state, dp)) + &
            (1.0_dp / real(n1 * n1, dp) + 1.0_dp / real(n2 * n2, dp) - 1.0_dp / real(n_state * n_state, dp))
        if (sok > 0.0_dp .and. son > 0.0_dp .and. skn > 0.0_dp) then
            g = (real(n1, dp) * log(son / sok) + real(n2, dp) * log(son / skn)) / c
            us(i) = g
        end if
    end do
case (cpm_joint_adjusted_kind)
    do i = 2, n_state - 2
        n1 = i
        n2 = n_state - n1
        mu1 = s(n1) / real(n1, dp)
        mu2 = (s(n_state) - s(n1)) / real(n2, dp)
        sok = w(n1) / real(n1, dp)
        son = w(n_state) / real(n_state, dp)
        skn = w(n_state) - w(n1) - real(n1 * n2, dp) * (mu1 - mu2) * (mu1 - mu2) / real(n_state, dp)
        skn = skn / real(n2, dp)
        if (sok > 0.0_dp .and. son > 0.0_dp .and. skn > 0.0_dp) then
            c = real(n_state, dp) * (log(2.0_dp / real(n_state, dp)) + digamma_half_index_1d(n_state - 1)) - &
                real(n1, dp) * (log(2.0_dp / real(n1, dp)) + digamma_half_index_1d(n1 - 1)) - &
                real(n2, dp) * (log(2.0_dp / real(n2, dp)) + digamma_half_index_1d(n2 - 1))
            if (abs(c) > 1.0e-12_dp) then
                g = 2.0_dp * (real(n1, dp) * log(son / sok) + real(n2, dp) * log(son / skn)) / c
                us(i) = g
            end if
        end if
    end do
case (cpm_joint_hawkins_kind)
    do i = 2, n_state - 2
        n1 = i
        n2 = n_state - n1
        mu1 = s(n1) / real(n1, dp)
        mu2 = (s(n_state) - s(n1)) / real(n2, dp)
        sok = w(n1) / real(n1, dp)
        son = w(n_state) / real(n_state, dp)
        skn = w(n_state) - w(n1) - real(n1 * n2, dp) * (mu1 - mu2) * (mu1 - mu2) / real(n_state, dp)
        skn = skn / real(n2, dp)
        c = 1.0_dp + 11.0_dp / 12.0_dp * (1.0_dp / real(n1, dp) + 1.0_dp / real(n2, dp) - 1.0_dp / real(n_state, dp)) + &
            (1.0_dp / real(n1 * n1, dp) + 1.0_dp / real(n2 * n2, dp) - 1.0_dp / real(n_state * n_state, dp))
        if (sok > 0.0_dp .and. son > 0.0_dp .and. skn > 0.0_dp) then
            g = (real(n1, dp) * log(son / sok) + real(n2, dp) * log(son / skn)) / c
            us(i) = g
        end if
    end do
    call apply_hawkins_edge_fixups_1d(us)
case (cpm_exponential_kind)
    do i = 1, n_state - 1
        j = real(i, dp)
        e = nn - j
        if (s(i) > 0.0_dp .and. (s(n_state) - s(i)) > 0.0_dp) then
            temp = -nn * log(s(n_state)) + j * log(s(i)) + e * log(s(n_state) - s(i))
            g = nn * log(nn) - j * log(j) - e * log(e)
            us(i) = -2.0_dp * (temp + g)
        end if
    end do
    us(n_state) = 0.0_dp
case (cpm_exponential_adjusted_kind)
    do i = 1, n_state - 1
        j = real(i, dp)
        e = nn - j
        if (s(i) > 0.0_dp .and. (s(n_state) - s(i)) > 0.0_dp) then
            temp = nn * log(nn / s(n_state)) - j * log(j / s(i)) - e * log(e / (s(n_state) - s(i)))
            g = nn * log(nn) - j * log(j) - e * log(e)
            c = j * digamma_half_index_1d(i) + e * digamma_half_index_1d(n_state - i) - nn * digamma_half_index_1d(n_state) + g
            if (abs(c) > 1.0e-12_dp) us(i) = temp / c
        end if
    end do
    us(n_state) = 0.0_dp
case (cpm_poisson_kind)
    do i = 2, n_state - 2
        j = real(i, dp)
        e = nn - j
        if (s(i) <= 0.0_dp) then
            mu1 = 0.5_dp
        else
            mu1 = s(i)
        end if
        if ((s(n_state) - s(i)) <= 0.0_dp) then
            mu2 = 0.5_dp
        else
            mu2 = s(n_state) - s(i)
        end if
        temp = e / j
        us(i) = abs((log(mu2 / mu1) - log(temp)) / sqrt(1.0_dp / mu1 + 1.0_dp / mu2))
    end do
case (cpm_mw_kind)
    call cpm_mle_mw_1d(w(1:n_state), us)
case (cpm_mood_kind)
    call cpm_mle_mood_1d(w(1:n_state), us)
case (cpm_lepage_kind)
    call cpm_mle_lepage_1d(w(1:n_state), us)
case (cpm_fet_kind)
    call cpm_mle_fet_1d(s(1:n_state), us, lambda)
case (cpm_ks_kind)
    call cpm_mle_ks_1d(w(1:n_state), us)
case (cpm_cvm_kind)
    call cpm_mle_cvm_1d(w(1:n_state), us)
end select
end subroutine cpm_mle_vector_1d

subroutine apply_hawkins_edge_fixups_1d(us)
real(kind=dp), intent(inout) :: us(:)
real(kind=dp), parameter :: mean_adjustments(3) = [2.2989_dp, 2.0814_dp, 2.0335_dp]
real(kind=dp), parameter :: sd_adjustments(3) = [2.3151_dp, 2.0871_dp, 2.0368_dp]
integer :: n_state, k

n_state = size(us)
if (n_state < 10) return
do k = 1, 3
    us(k + 1) = ((us(k + 1) - mean_adjustments(k)) / sd_adjustments(k)) * 2.0_dp + 2.0_dp
    us(n_state - k - 1) = ((us(n_state - k - 1) - mean_adjustments(k)) / sd_adjustments(k)) * 2.0_dp + 2.0_dp
end do
end subroutine apply_hawkins_edge_fixups_1d

real(kind=dp) function digamma_half_index_1d(n) result(val)
integer, intent(in) :: n
real(kind=dp), parameter :: euler_gamma = 0.5772156649015328606_dp
integer :: m, j

if (n <= 0) then
    val = 0.0_dp
else if (mod(n, 2) == 0) then
    m = n / 2
    val = -euler_gamma
    do j = 1, m - 1
        val = val + 1.0_dp / real(j, dp)
    end do
else
    m = (n - 1) / 2
    val = -euler_gamma - 2.0_dp * log(2.0_dp)
    do j = 1, m
        val = val + 2.0_dp / real(2 * j - 1, dp)
    end do
end if
end function digamma_half_index_1d

subroutine cpm_mle_mw_1d(ranks, us)
real(kind=dp), intent(in) :: ranks(:)
real(kind=dp), intent(inout) :: us(:)
real(kind=dp), allocatable :: cumsums(:)
integer :: i, n_state
real(kind=dp) :: n0, n1, r1, u, mu, sd

n_state = size(ranks)
allocate(cumsums(n_state))
cumsums(1) = ranks(1)
do i = 2, n_state
    cumsums(i) = cumsums(i - 1) + ranks(i)
end do
do i = 2, n_state - 2
    n0 = real(i, dp)
    n1 = real(n_state, dp) - n0
    r1 = cumsums(i)
    u = r1 - n0 * (n0 + 1.0_dp) / 2.0_dp
    mu = n0 * n1 / 2.0_dp
    sd = sqrt(n0 * n1 * (n0 + n1 + 1.0_dp) / 12.0_dp)
    us(i) = abs((u - mu) / sd)
end do
deallocate(cumsums)
end subroutine cpm_mle_mw_1d

subroutine cpm_mle_mood_1d(ranks, us)
real(kind=dp), intent(in) :: ranks(:)
real(kind=dp), intent(inout) :: us(:)
real(kind=dp), allocatable :: cumsums(:)
integer :: i, n_state
real(kind=dp) :: n0, n1, m, mu, sd, med

n_state = size(ranks)
allocate(cumsums(n_state))
med = (real(n_state, dp) + 1.0_dp) / 2.0_dp
cumsums(1) = (ranks(1) - med) * (ranks(1) - med)
do i = 2, n_state
    cumsums(i) = cumsums(i - 1) + (ranks(i) - med) * (ranks(i) - med)
end do
do i = 2, n_state - 2
    n0 = real(i, dp)
    n1 = real(n_state, dp) - n0
    m = cumsums(i)
    mu = n0 * (real(n_state * n_state, dp) - 1.0_dp) / 12.0_dp
    sd = sqrt(n0 * n1 * (real(n_state, dp) + 1.0_dp) * (real(n_state * n_state, dp) - 4.0_dp) / 180.0_dp)
    us(i) = abs((m - mu) / sd)
end do
deallocate(cumsums)
end subroutine cpm_mle_mood_1d

subroutine cpm_mle_lepage_1d(ranks, us)
real(kind=dp), intent(in) :: ranks(:)
real(kind=dp), intent(inout) :: us(:)
real(kind=dp), allocatable :: us_mw(:), us_mood(:)
integer :: i

allocate(us_mw(size(us)), us_mood(size(us)))
us_mw = 0.0_dp
us_mood = 0.0_dp
call cpm_mle_mw_1d(ranks, us_mw)
call cpm_mle_mood_1d(ranks, us_mood)
do i = 2, size(us) - 2
    us(i) = us_mw(i) * us_mw(i) + us_mood(i) * us_mood(i)
end do
deallocate(us_mw, us_mood)
end subroutine cpm_mle_lepage_1d

subroutine cpm_mle_fet_1d(sums, us, lambda)
real(kind=dp), intent(in) :: sums(:)
real(kind=dp), intent(inout) :: us(:)
real(kind=dp), intent(in), optional :: lambda
integer :: i, n_state, n0, n1, s0, s1
real(kind=dp) :: lam

n_state = size(sums)
lam = 0.1_dp
if (present(lambda)) lam = lambda

do i = 2, n_state - 2
    n0 = i
    n1 = n_state - n0
    s0 = nint(sums(i))
    s1 = nint(sums(n_state) - sums(i))
    us(i) = hypergeom_upper_tail_gt_1d(s0, s0 + s1, n0 + n1 - s0 - s1, n0)
end do

if (n_state > 3 .and. lam < 1.0_dp) then
    do i = 3, n_state - 2
        us(i) = (1.0_dp - lam) * us(i - 1) + lam * us(i)
    end do
end if
end subroutine cpm_mle_fet_1d

subroutine cpm_mle_ks_1d(order_vec, us)
real(kind=dp), intent(in) :: order_vec(:)
real(kind=dp), intent(inout) :: us(:)
real(kind=dp), allocatable :: cumsums(:)
integer :: i, j, n_state
real(kind=dp) :: n0, n1, a, b, statistic, temp, z, correction, big_n, swapv

n_state = size(order_vec)
big_n = real(n_state, dp)
allocate(cumsums(n_state))

do i = 2, n_state - 2
    n0 = real(i, dp)
    n1 = big_n - real(i, dp)
    a = 1.0_dp / n0
    b = -1.0_dp / n1
    do j = 1, n_state
        if (order_vec(j) <= n0) then
            cumsums(j) = a
        else
            cumsums(j) = b
        end if
    end do
    do j = 2, n_state
        cumsums(j) = cumsums(j - 1) + cumsums(j)
    end do
    statistic = 0.0_dp
    do j = 1, n_state
        temp = abs(cumsums(j))
        if (temp > statistic) statistic = temp
    end do

    if (n1 > n0) then
        swapv = n1
        n1 = n0
        n0 = swapv
    end if
    if (n0 > 2.0_dp * n1) then
        correction = 1.0_dp / (2.0_dp * sqrt(n0))
    else
        if (mod(int(n0), int(n1)) == 0) then
            correction = 2.0_dp / (3.0_dp * sqrt(n0))
        else
            correction = 2.0_dp / (5.0_dp * sqrt(n0))
        end if
    end if
    z = statistic * sqrt((n0 * n1) / (n0 + n1)) + correction
    z = z * z
    us(i) = 1.0_dp - 2.0_dp * (exp(-2.0_dp * z) - exp(-8.0_dp * z))
end do
deallocate(cumsums)
end subroutine cpm_mle_ks_1d

subroutine cpm_mle_cvm_1d(order_vec, us)
real(kind=dp), intent(in) :: order_vec(:)
real(kind=dp), intent(inout) :: us(:)
real(kind=dp), allocatable :: cumsums(:)
integer :: i, j, n_state
real(kind=dp) :: n0, n1, a, b, statistic, mu, sigma, big_n, prod

n_state = size(order_vec)
big_n = real(n_state, dp)
allocate(cumsums(n_state))
mu = 1.0_dp / 6.0_dp + 1.0_dp / (6.0_dp * big_n)

do i = 2, n_state - 2
    n0 = real(i, dp)
    n1 = big_n - real(i + 1, dp) + 1.0_dp
    a = 1.0_dp / n0
    b = -1.0_dp / n1
    do j = 1, n_state
        if (order_vec(j) <= n0) then
            cumsums(j) = a
        else
            cumsums(j) = b
        end if
    end do
    do j = 2, n_state
        cumsums(j) = cumsums(j - 1) + cumsums(j)
    end do
    statistic = 0.0_dp
    do j = 1, n_state
        statistic = statistic + cumsums(j) * cumsums(j)
    end do
    prod = n0 * n1
    sigma = sqrt((1.0_dp / 45.0_dp) * (big_n + 1.0_dp) / (big_n * big_n) * &
        (4.0_dp * prod * big_n - 3.0_dp * (n1 * n1 + n0 * n0) - 2.0_dp * prod) / (4.0_dp * prod))
    us(i) = (statistic * prod / (big_n * big_n) - mu) / sigma
end do

deallocate(cumsums)
end subroutine cpm_mle_cvm_1d

real(kind=dp) function hypergeom_upper_tail_gt_1d(q, white, black, draws)
integer, intent(in) :: q, white, black, draws
integer :: x, xmax, xmin
real(kind=dp) :: logp, p

xmin = max(0, draws - black)
xmax = min(draws, white)
if (q >= xmax) then
    hypergeom_upper_tail_gt_1d = 0.0_dp
    return
end if
if (q < xmin) then
    hypergeom_upper_tail_gt_1d = 1.0_dp
    return
end if

hypergeom_upper_tail_gt_1d = 0.0_dp
do x = q + 1, xmax
    logp = log_binom_coeff_1d(white, x) + log_binom_coeff_1d(black, draws - x) - log_binom_coeff_1d(white + black, draws)
    p = exp(logp)
    hypergeom_upper_tail_gt_1d = hypergeom_upper_tail_gt_1d + p
end do
if (hypergeom_upper_tail_gt_1d < 0.0_dp) hypergeom_upper_tail_gt_1d = 0.0_dp
if (hypergeom_upper_tail_gt_1d > 1.0_dp) hypergeom_upper_tail_gt_1d = 1.0_dp
end function hypergeom_upper_tail_gt_1d

real(kind=dp) function log_binom_coeff_1d(n, k)
integer, intent(in) :: n, k
if (k < 0 .or. k > n) then
    log_binom_coeff_1d = -huge(1.0_dp)
else
    log_binom_coeff_1d = log_gamma(real(n + 1, dp)) - log_gamma(real(k + 1, dp)) - log_gamma(real(n - k + 1, dp))
end if
end function log_binom_coeff_1d

end module cpm_pkg_mod
