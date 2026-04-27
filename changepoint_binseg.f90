module changepoint_binseg_mod
use kind_mod, only: dp
implicit none
private
public :: solve_binseg_cost, solve_binseg_mean_shift

contains

subroutine solve_binseg_cost(cost, max_cp, cp, ncp_found, min_gain)
! perform binary segmentation on a precomputed segment cost matrix
real(kind=dp), intent(in) :: cost(:,:)
integer, intent(in) :: max_cp
integer, intent(out) :: cp(max_cp)
integer, intent(out) :: ncp_found
real(kind=dp), intent(in), optional :: min_gain

integer :: n
integer :: seg_start(max_cp+1), seg_end(max_cp+1)
integer :: nseg
integer :: iter
integer :: iseg
integer :: best_seg
integer :: s, e, k
integer :: best_k_seg
real(kind=dp) :: gain
real(kind=dp) :: best_gain
real(kind=dp) :: min_gain_
real(kind=dp), parameter :: huge_cost = 1.0e19_dp

n = size(cost, 1)
cp = 0
ncp_found = 0
min_gain_ = 0.0_dp
if (present(min_gain)) min_gain_ = min_gain

seg_start(1) = 1
seg_end(1) = n
nseg = 1

do iter = 1, max_cp

    best_gain = -1.0_dp
    best_seg = 0
    best_k_seg = 0

    do iseg = 1, nseg
        s = seg_start(iseg)
        e = seg_end(iseg)

        if (s >= e) cycle
        if (cost(s, e) >= huge_cost) cycle

        do k = s, e - 1
            if (cost(s, k) >= huge_cost) cycle
            if (cost(k+1, e) >= huge_cost) cycle
            gain = cost(s, e) - cost(s, k) - cost(k+1, e)
            if (gain > best_gain) then
                best_gain = gain
                best_seg = iseg
                best_k_seg = k
            end if
        end do
    end do

    if (best_seg == 0) exit
    if (best_gain <= min_gain_) exit

    ncp_found = ncp_found + 1
    cp(ncp_found) = best_k_seg

    call split_segment(seg_start, seg_end, nseg, best_seg, best_k_seg)

end do

if (ncp_found > 1) call sort_int_inplace(cp, ncp_found)

end subroutine solve_binseg_cost

subroutine split_segment(seg_start, seg_end, nseg, iseg, ksplit)
! replace one segment by its two child segments
integer, intent(inout) :: seg_start(:), seg_end(:)
integer, intent(inout) :: nseg
integer, intent(in) :: iseg, ksplit

integer :: s, e
integer :: j

s = seg_start(iseg)
e = seg_end(iseg)

do j = nseg, iseg + 1, -1
    seg_start(j+1) = seg_start(j)
    seg_end(j+1) = seg_end(j)
end do

seg_start(iseg) = s
seg_end(iseg) = ksplit
seg_start(iseg+1) = ksplit + 1
seg_end(iseg+1) = e
nseg = nseg + 1

end subroutine split_segment

subroutine sort_int_inplace(x, n)
! sort the first n elements of x in ascending order
integer, intent(inout) :: x(:)
integer, intent(in) :: n

integer :: i, j, tmp

do i = 1, n - 1
    do j = i + 1, n
        if (x(j) < x(i)) then
            tmp = x(i)
            x(i) = x(j)
            x(j) = tmp
        end if
    end do
end do

end subroutine sort_int_inplace

subroutine solve_binseg_mean_shift(z, max_cp, cp, ncp_found, min_seg_len, min_gain)
! perform binary segmentation for changepoints in the mean of a 1d series
real(kind=dp), intent(in) :: z(:)
integer, intent(in) :: max_cp
integer, intent(out) :: cp(max_cp)
integer, intent(out) :: ncp_found
integer, intent(in), optional :: min_seg_len
real(kind=dp), intent(in), optional :: min_gain

real(kind=dp), allocatable :: sz(:), szz(:)
integer :: n
integer :: min_seg_len_
real(kind=dp) :: min_gain_
real(kind=dp), parameter :: huge_cost = 1.0e19_dp
integer :: seg_start(max_cp+1), seg_end(max_cp+1)
integer :: nseg
integer :: iter
integer :: iseg
integer :: best_seg
integer :: s, e, k
integer :: best_k_seg
real(kind=dp) :: gain
real(kind=dp) :: best_gain
real(kind=dp) :: parent_cost, left_cost, right_cost

min_seg_len_ = 50
if (present(min_seg_len)) min_seg_len_ = min_seg_len

n = size(z)
cp = 0
ncp_found = 0
min_gain_ = 0.0_dp
if (present(min_gain)) min_gain_ = min_gain

allocate(sz(0:n), szz(0:n))
call build_mean_shift_prefix(z, sz, szz)

seg_start(1) = 1
seg_end(1) = n
nseg = 1

do iter = 1, max_cp

    best_gain = -1.0_dp
    best_seg = 0
    best_k_seg = 0

    do iseg = 1, nseg
        s = seg_start(iseg)
        e = seg_end(iseg)
        parent_cost = mean_shift_segment_cost(sz, szz, s, e, min_seg_len_, huge_cost)
        if (parent_cost >= huge_cost) cycle

        do k = s, e - 1
            left_cost = mean_shift_segment_cost(sz, szz, s, k, min_seg_len_, huge_cost)
            if (left_cost >= huge_cost) cycle

            right_cost = mean_shift_segment_cost(sz, szz, k + 1, e, min_seg_len_, huge_cost)
            if (right_cost >= huge_cost) cycle

            gain = parent_cost - left_cost - right_cost
            if (gain > best_gain) then
                best_gain = gain
                best_seg = iseg
                best_k_seg = k
            end if
        end do
    end do

    if (best_seg == 0) exit
    if (best_gain <= min_gain_) exit

    ncp_found = ncp_found + 1
    cp(ncp_found) = best_k_seg

    call split_segment(seg_start, seg_end, nseg, best_seg, best_k_seg)

end do

if (ncp_found > 1) call sort_int_inplace(cp, ncp_found)

deallocate(sz, szz)

end subroutine solve_binseg_mean_shift

subroutine build_mean_shift_prefix(z, sz, szz)
real(kind=dp), intent(in) :: z(:)
real(kind=dp), intent(out) :: sz(0:), szz(0:)
integer :: i, n

n = size(z)
sz(0) = 0.0_dp
szz(0) = 0.0_dp
do i = 1, n
    sz(i) = sz(i - 1) + z(i)
    szz(i) = szz(i - 1) + z(i) * z(i)
end do

end subroutine build_mean_shift_prefix

pure function mean_shift_segment_cost(sz, szz, s, e, min_seg_len, huge_cost) result(cost)
real(kind=dp), intent(in) :: sz(0:), szz(0:)
integer, intent(in) :: s, e, min_seg_len
real(kind=dp), intent(in) :: huge_cost
real(kind=dp) :: cost
real(kind=dp) :: m, sum_z, sum_zz, var_z

if (e < s) then
    cost = huge_cost
    return
end if

m = real(e - s + 1, dp)
if (m < real(min_seg_len, dp)) then
    cost = huge_cost
    return
end if

sum_z = sz(e) - sz(s - 1)
sum_zz = szz(e) - szz(s - 1)
var_z = sum_zz / m - (sum_z / m) ** 2
cost = 0.5_dp * m * log(max(var_z, 1.0e-20_dp))

end function mean_shift_segment_cost

end module changepoint_binseg_mod
