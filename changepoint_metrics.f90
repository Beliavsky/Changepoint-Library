module changepoint_metrics_mod
use kind_mod, only: dp, long_int
implicit none
private
public :: precision_recall_metric, hausdorff_metric, randindex_metric, print_metric_block

contains

subroutine precision_recall_metric(true_bkps, my_bkps, margin, precision, recall)
integer, intent(in) :: true_bkps(:), my_bkps(:), margin
real(kind=dp), intent(out) :: precision, recall
logical, allocatable :: used(:)
integer :: i, j, tp

if (size(my_bkps) == 1) then
    precision = 0.0_dp
    recall = 0.0_dp
    return
end if

allocate(used(size(my_bkps) - 1))
used = .false.
tp = 0
do i = 1, size(true_bkps) - 1
    do j = 1, size(my_bkps) - 1
        if (.not. used(j) .and. my_bkps(j) - margin < true_bkps(i) .and. true_bkps(i) < my_bkps(j) + margin) then
            used(j) = .true.
            tp = tp + 1
            exit
        end if
    end do
end do

precision = real(tp, dp) / real(size(my_bkps) - 1, dp)
recall = real(tp, dp) / real(size(true_bkps) - 1, dp)
deallocate(used)
end subroutine precision_recall_metric

real(kind=dp) function hausdorff_metric(bkps1, bkps2) result(hd)
integer, intent(in) :: bkps1(:), bkps2(:)
integer :: i, j
real(kind=dp) :: dmin, dmax1, dmax2

dmax1 = 0.0_dp
do i = 1, size(bkps1) - 1
    dmin = huge(1.0_dp)
    do j = 1, size(bkps2) - 1
        dmin = min(dmin, abs(real(bkps1(i) - bkps2(j), dp)))
    end do
    dmax1 = max(dmax1, dmin)
end do

dmax2 = 0.0_dp
do i = 1, size(bkps2) - 1
    dmin = huge(1.0_dp)
    do j = 1, size(bkps1) - 1
        dmin = min(dmin, abs(real(bkps2(i) - bkps1(j), dp)))
    end do
    dmax2 = max(dmax2, dmin)
end do

hd = max(dmax1, dmax2)
end function hausdorff_metric

real(kind=dp) function randindex_metric(bkps1, bkps2) result(ri)
integer, intent(in) :: bkps1(:), bkps2(:)
integer :: i, j, start1, end1, start2, end2, nij, beginj, n_samples, disagreement

n_samples = bkps1(size(bkps1))
disagreement = 0
beginj = 1
do i = 1, size(bkps1)
    if (i == 1) then
        start1 = 0
    else
        start1 = bkps1(i - 1)
    end if
    end1 = bkps1(i)
    do j = beginj, size(bkps2)
        if (j == 1) then
            start2 = 0
        else
            start2 = bkps2(j - 1)
        end if
        end2 = bkps2(j)
        nij = max(min(end1, end2) - max(start1, start2), 0)
        disagreement = disagreement + nij * abs(end1 - end2)
        if (end1 < end2) then
            exit
        else
            beginj = j + 1
        end if
    end do
end do
ri = 1.0_dp - real(disagreement, dp) / (real(n_samples, dp) * real(n_samples - 1, dp) / 2.0_dp)
end function randindex_metric

subroutine print_metric_block(name, bkps, truth, margin, t0, t1, rate)
character(len=*), intent(in) :: name
integer, intent(in) :: bkps(:), truth(:), margin
integer(kind=long_int), intent(in) :: t0, t1, rate
real(kind=dp) :: precision, recall, hd, ri

call precision_recall_metric(truth, bkps, margin, precision, recall)
hd = hausdorff_metric(truth, bkps)
ri = randindex_metric(truth, bkps)
print *, "method      =", trim(name)
print *, "estimated   =", bkps
write(*,'(a,f8.3)') "precision = ", precision
write(*,'(a,f8.3)') "recall = ", recall
write(*,'(a,f8.3)') "hausdorff = ", hd
write(*,'(a,f10.6)') "randindex = ", ri
write(*,'(a,f8.3)') "elapsed seconds = ", real(t1 - t0, dp) / real(rate, dp)
end subroutine print_metric_block

end module changepoint_metrics_mod
