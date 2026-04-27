module ecp_pkg_mod
use kind_mod, only: dp
implicit none
private

public :: solve_e_divisive_fixedk_2d, solve_e_divisive_permtest_2d, solve_e_agglo_default_2d, solve_e_cp3o_2d, solve_e_cp3o_delta_2d

contains

subroutine solve_e_divisive_fixedk_2d(x, k, min_size, alpha, estimates, cluster)
real(kind=dp), intent(in) :: x(:, :)
integer, intent(in) :: k, min_size
real(kind=dp), intent(in) :: alpha
integer, allocatable, intent(out) :: estimates(:), cluster(:)

real(kind=dp), allocatable :: d(:, :)
integer, allocatable :: changes(:)
integer :: n, step, best_seg_idx, best_split
real(kind=dp) :: best_stat

n = size(x, 1)
allocate(d(n, n))
call build_distance_matrix_2d(x, alpha, d)

allocate(changes(2))
changes = [1, n + 1]

do step = 1, k
    call best_global_split_1d(changes, d, min_size, best_seg_idx, best_split, best_stat)
    if (best_seg_idx <= 0 .or. best_split <= 0) exit
    call insert_change_1d(changes, best_split)
end do

allocate(estimates(size(changes)))
estimates = changes
call build_cluster_1d(changes, cluster)

deallocate(changes, d)
end subroutine solve_e_divisive_fixedk_2d

subroutine solve_e_divisive_permtest_2d(x, sig_lvl, r, min_size, alpha, perms_in, order_found, estimates, considered_last, p_values, permutations, cluster)
real(kind=dp), intent(in) :: x(:, :)
real(kind=dp), intent(in) :: sig_lvl, alpha
integer, intent(in) :: r, min_size
integer, intent(in) :: perms_in(:, :)
integer, allocatable, intent(out) :: order_found(:), estimates(:), permutations(:), cluster(:)
integer, intent(out) :: considered_last
real(kind=dp), allocatable, intent(out) :: p_values(:)

real(kind=dp), allocatable :: d(:, :), d1(:, :)
integer, allocatable :: changes(:), p_work(:), perm(:)
real(kind=dp), allocatable :: pval_work(:)
integer :: n, stage, stage_count_max, best_seg_idx, best_split, over, f, row_idx
integer :: perm_best_seg_idx, perm_best_split
real(kind=dp) :: best_stat, perm_best_stat, pval

n = size(x, 1)
allocate(d(n, n), d1(n, n), perm(n))
call build_distance_matrix_2d(x, alpha, d)

allocate(changes(2))
changes = [1, n + 1]
allocate(p_work(max(1, size(perms_in, 1) / max(1, r))))
allocate(pval_work(size(p_work)))
stage_count_max = size(p_work)
stage = 0
considered_last = -1

do while (stage < stage_count_max)
    call best_global_split_1d(changes, d, min_size, best_seg_idx, best_split, best_stat)
    if (best_seg_idx <= 0 .or. best_split <= 0) exit
    considered_last = best_split
    stage = stage + 1
    over = 0
    do f = 1, r
        row_idx = (stage - 1) * r + f
        if (row_idx > size(perms_in, 1)) exit
        perm = perms_in(row_idx, :)
        call permute_distance_1d(d, changes, perm, d1)
        call best_global_split_1d(changes, d1, min_size, perm_best_seg_idx, perm_best_split, perm_best_stat)
        if (perm_best_stat >= best_stat) over = over + 1
    end do
    pval = real(1 + over, dp) / real(r + 1, dp)
    p_work(stage) = r
    pval_work(stage) = pval
    if (pval > sig_lvl) exit
    call append_change_unsorted_1d(changes, best_split)
end do

allocate(order_found(size(changes)))
order_found = changes
call sorted_copy_1d(changes, estimates)
call build_cluster_1d(estimates, cluster)

if (stage > 0) then
    allocate(permutations(stage), p_values(stage))
    permutations = p_work(1:stage)
    p_values = pval_work(1:stage)
else
    allocate(permutations(0), p_values(0))
end if

deallocate(changes, d, d1, perm, p_work, pval_work)
end subroutine solve_e_divisive_permtest_2d

subroutine solve_e_agglo_default_2d(x, alpha, estimates, cluster, fit, progression, merged)
real(kind=dp), intent(in) :: x(:, :)
real(kind=dp), intent(in) :: alpha
integer, allocatable, intent(out) :: estimates(:), cluster(:), progression(:, :), merged(:, :)
real(kind=dp), allocatable, intent(out) :: fit(:)

real(kind=dp), allocatable :: dobs(:, :), distm(:, :)
real(kind=dp) :: best_val, xval
integer, allocatable :: sizes(:), left(:), right(:), lm(:), est_row(:)
logical, allocatable :: open(:)
integer :: n, k, i, j, best_i, best_j, row_idx, best_row

n = size(x, 1)
allocate(dobs(n, n))
call build_distance_matrix_2d(x, alpha, dobs)

allocate(distm(2 * n, 2 * n))
distm = huge(1.0_dp)
do i = 1, 2 * n
    distm(i, i) = 0.0_dp
end do
do i = 1, n
    do j = i + 1, n
        distm(i, j) = 2.0_dp * dobs(i, j)
        distm(j, i) = distm(i, j)
    end do
end do

allocate(sizes(2 * n))
sizes = 0
sizes(1:n) = 1

allocate(left(2 * n - 1), right(2 * n - 1), lm(2 * n - 1))
allocate(open(2 * n - 1))
left = 0
right = 0
open = .true.
lm = 0
lm(1:n) = [(i, i = 1, n)]

do i = 2, n - 1
    left(i) = i - 1
    right(i) = i + 1
end do
left(1) = n
right(1) = 2
left(n) = n - 1
right(n) = 1

allocate(merged(n - 1, 2))
merged = 0
allocate(progression(n, n + 1))
progression = 0
progression(1, :) = [(i, i = 1, n + 1)]
allocate(fit(n))
fit(1) = 0.0_dp
do i = 1, n
    fit(1) = fit(1) + distm(i, left(i)) + distm(i, right(i))
end do

do k = n, 2 * n - 2
    best_val = -huge(1.0_dp)
    best_i = 0
    best_j = 0
    do i = 1, k
        if (open(i)) then
            xval = gof_update_agglo_1d(i, distm, sizes, left, right, fit(k - n + 1))
            if (xval > best_val) then
                best_val = xval
                best_i = i
                best_j = right(i)
            end if
        end if
    end do
    row_idx = k - n + 2
    fit(row_idx) = best_val
    call update_distance_agglo_1d(best_i, best_j, k, n, distm, sizes, left, right, open, lm, progression, merged)
end do

best_row = maxloc(fit, dim = 1)
allocate(est_row(count(progression(best_row, :) > 0)))
est_row = pack(progression(best_row, :), progression(best_row, :) > 0)
if (est_row(1) /= 1) then
    allocate(estimates(size(est_row) - 1))
    estimates = est_row(1:size(est_row) - 1)
else
    allocate(estimates(size(est_row)))
    estimates = est_row
end if
call build_cluster_from_estimates_1d(estimates, n, cluster)

deallocate(dobs, distm, sizes, left, right, lm, open, est_row)
end subroutine solve_e_agglo_default_2d

subroutine solve_e_cp3o_2d(x, kmax, minsize, alpha, number, estimates, gofm, cploc, cploc_lens)
real(kind=dp), intent(in) :: x(:, :)
integer, intent(in) :: kmax, minsize
real(kind=dp), intent(in) :: alpha
integer, intent(out) :: number
integer, allocatable, intent(out) :: estimates(:), cploc(:, :), cploc_lens(:)
real(kind=dp), allocatable, intent(out) :: gofm(:)

real(kind=dp), allocatable :: d(:, :), pref(:, :), ff(:, :)
integer, allocatable :: a(:, :)
real(kind=dp), allocatable :: poi(:)
integer :: n, kfit, tf, ef, uf, best_idx, kk, lenk
real(kind=dp) :: stat, best_sse, sse

n = size(x, 1)
allocate(d(n, n), pref(0:n, 0:n))
call build_distance_matrix_2d(x, alpha, d)
call build_prefix_sum_2d(d, pref)

allocate(ff(kmax, n), a(kmax, n))
ff = -huge(1.0_dp)
a = 0

do tf = 2 * minsize, n
    do ef = minsize + 1, tf - minsize + 1
        stat = cp3o_seg_stat_1d(1, ef, tf, pref)
        if (stat > ff(1, tf)) then
            ff(1, tf) = stat
            a(1, tf) = ef
        end if
    end do
end do

do kfit = 2, kmax
    do tf = (kfit + 1) * minsize, n
        do ef = kfit * minsize + 1, tf - minsize + 1
            if (ff(kfit - 1, ef - 1) <= -huge(1.0_dp) / 2.0_dp) cycle
            uf = a(kfit - 1, ef - 1)
            if (uf <= 0) cycle
            stat = cp3o_seg_stat_1d(uf, ef, tf, pref) + ff(kfit - 1, ef - 1)
            if (stat > ff(kfit, tf)) then
                ff(kfit, tf) = stat
                a(kfit, tf) = ef
            end if
        end do
    end do
end do

allocate(gofm(kmax))
gofm = ff(:, n)

if (kmax == 1) then
    number = 1
else
    allocate(poi(kmax - 2))
    do kk = 1, kmax - 2
        call cp3o_poi_sse_1d(gofm, kk + 1, sse)
        poi(kk) = sse
    end do
    best_idx = 1
    best_sse = poi(1)
    do kk = 2, size(poi)
        if (poi(kk) < best_sse) then
            best_sse = poi(kk)
            best_idx = kk
        end if
    end do
    number = best_idx + 1
    deallocate(poi)
end if

allocate(cploc(kmax, kmax), cploc_lens(kmax))
cploc = 0
cploc_lens = 0
do kfit = 1, kmax
    call cp3o_backtrack_1d(a, n, kfit, cploc(kfit, :), lenk)
    cploc_lens(kfit) = lenk
end do

allocate(estimates(cploc_lens(number)))
estimates = cploc(number, 1:cploc_lens(number))

deallocate(d, pref, ff, a)
end subroutine solve_e_cp3o_2d

subroutine solve_e_cp3o_delta_2d(x, kmax, delta, alpha, number, estimates, gofm, cploc, cploc_lens)
real(kind=dp), intent(in) :: x(:, :)
integer, intent(in) :: kmax, delta
real(kind=dp), intent(in) :: alpha
integer, intent(out) :: number
integer, allocatable, intent(out) :: estimates(:), cploc(:, :), cploc_lens(:)
real(kind=dp), allocatable, intent(out) :: gofm(:)

real(kind=dp), allocatable :: d(:, :), ff(:, :), dll(:), drr(:), dlr(:), csum(:), gof_work(:)
integer, allocatable :: a(:, :)
logical, allocatable :: cand(:, :), keep(:)
real(kind=dp), allocatable :: poi(:)
integer :: n, kuse, minsize, kfit, t0, s0, best_idx, kk, lenk
integer :: t, u0, flag, prev_flag, a2, b2
real(kind=dp) :: stat, best_sse, sse, stat2

n = size(x, 1)
minsize = delta + 1
kuse = min(kmax, n / max(1, minsize))
if (kuse < 1) kuse = 1

allocate(d(n, n), dll(n), drr(n), dlr(n), csum(n), gof_work(kuse))
call build_distance_matrix_2d(x, alpha, d)
call build_delta_precompute_2d(d, delta, dll, drr, dlr, csum)

allocate(ff(2, n), a(kuse, n), cand(n, n), keep(n))
ff = -huge(1.0_dp)
a = -1
cand = .false.
gof_work = -huge(1.0_dp)

if (kuse == 1) then
    t0 = n - 1
    t = n
    do s0 = delta, t0 - minsize
        stat = cp3o_delta_seg_stat_1d(-1, s0, t0, delta, dll, drr, dlr, csum)
        if (stat > ff(1, t)) then
            ff(1, t) = stat
            a(1, t) = s0
        end if
    end do
    gof_work(1) = ff(1, n)
else
    do t0 = 2 * minsize - 1, n - 1
        t = t0 + 1
        keep = .false.
        do s0 = delta, t0 - minsize
            keep(s0 + 1) = .true.
            stat = cp3o_delta_seg_stat_1d(-1, s0, t0, delta, dll, drr, dlr, csum)
            if (stat > ff(1, t)) then
                ff(1, t) = stat
                a(1, t) = s0
            end if
        end do
        a2 = t0 - minsize
        stat2 = cp3o_delta_seg_stat_1d(a(1, a2 + 1), a2, t0, delta, dll, drr, dlr, csum)
        do s0 = delta, t0 - minsize
            stat = cp3o_delta_seg_stat_1d(a(1, s0 + 1), s0, t0, delta, dll, drr, dlr, csum)
            if (ff(1, s0 + 1) + stat < ff(1, a2 + 1) + stat2) keep(s0 + 1) = .false.
        end do
        cand(:, t) = keep
    end do
    gof_work(1) = ff(1, n)

    flag = 2
    prev_flag = 1
    do kfit = 2, kuse
        if (kfit >= kuse) then
            t0 = n - 1
            t = n
            ff(flag, t) = -huge(1.0_dp)
            a(kfit, t) = -1
            do s0 = delta, t0 - minsize
                if (.not. cand(s0 + 1, t)) cycle
                u0 = a(kfit - 1, s0 + 1)
                stat = cp3o_delta_seg_stat_1d(u0, s0, t0, delta, dll, drr, dlr, csum)
                if (u0 > 0) stat = stat + ff(prev_flag, s0 + 1)
                if (stat > ff(flag, t)) then
                    ff(flag, t) = stat
                    a(kfit, t) = s0
                end if
            end do
        else
            do t0 = 2 * minsize - 1, n - 1
                t = t0 + 1
                ff(flag, t) = -huge(1.0_dp)
                a(kfit, t) = -1
                do s0 = delta, t0 - minsize
                    if (.not. cand(s0 + 1, t)) cycle
                    u0 = a(kfit - 1, s0 + 1)
                    stat = cp3o_delta_seg_stat_1d(u0, s0, t0, delta, dll, drr, dlr, csum)
                    if (u0 > 0) stat = stat + ff(prev_flag, s0 + 1)
                    if (stat > ff(flag, t)) then
                        ff(flag, t) = stat
                        a(kfit, t) = s0
                    end if
                end do
                a2 = t0 - minsize
                b2 = a(kfit, a2 + 1)
                stat2 = cp3o_delta_seg_stat_1d(b2, a2, t0, delta, dll, drr, dlr, csum)
                keep = cand(:, t)
                do s0 = delta, t0 - minsize
                    if (.not. keep(s0 + 1)) cycle
                    stat = cp3o_delta_seg_stat_1d(a(kfit, s0 + 1), s0, t0, delta, dll, drr, dlr, csum)
                    if (ff(flag, s0 + 1) + stat < ff(flag, a2 + 1) + stat2) keep(s0 + 1) = .false.
                end do
                cand(:, t) = keep
            end do
        end if
        gof_work(kfit) = ff(flag, n)
        prev_flag = flag
        flag = 3 - flag
    end do
end if

allocate(gofm(kuse))
gofm = gof_work

if (kuse == 1) then
    number = 1
else
    allocate(poi(max(0, kuse - 2)))
    if (size(poi) == 0) then
        number = 1
    else
        do kk = 1, size(poi)
            call cp3o_poi_sse_1d(gofm, kk + 1, sse)
            poi(kk) = sse
        end do
        best_idx = 1
        best_sse = poi(1)
        do kk = 2, size(poi)
            if (poi(kk) < best_sse) then
                best_sse = poi(kk)
                best_idx = kk
            end if
        end do
        number = best_idx + 1
    end if
    deallocate(poi)
end if

allocate(cploc(kuse, kuse), cploc_lens(kuse))
cploc = 0
cploc_lens = 0
do kfit = 1, kuse
    call cp3o_delta_backtrack_1d(a, n, kfit, cploc(kfit, :), lenk)
    cploc_lens(kfit) = lenk
end do

allocate(estimates(cploc_lens(number)))
estimates = cploc(number, 1:cploc_lens(number))

deallocate(d, ff, a, dll, drr, dlr, csum, cand, keep, gof_work)
end subroutine solve_e_cp3o_delta_2d

subroutine build_delta_precompute_2d(d, delta, dll, drr, dlr, csum)
real(kind=dp), intent(in) :: d(:, :)
integer, intent(in) :: delta
real(kind=dp), intent(out) :: dll(:), drr(:), dlr(:), csum(:)

real(kind=dp), allocatable :: dll0(:), drr0(:), dlr0(:), csum0(:), left0(:, :), right0(:, :)
integer :: n, s, i, j, minsize
real(kind=dp) :: r1, r2, r3, a1, a2, a3

n = size(d, 1)
minsize = delta + 1
dll = 0.0_dp
drr = 0.0_dp
dlr = 0.0_dp
csum = 0.0_dp

allocate(dll0(0:n-1), drr0(0:n-1), dlr0(0:n-1), csum0(0:n-1), left0(0:n-1, 0:1), right0(0:n-1, 0:1))
dll0 = 0.0_dp
drr0 = 0.0_dp
dlr0 = 0.0_dp
csum0 = 0.0_dp
left0 = 0.0_dp
right0 = 0.0_dp

do s = delta, n - 1
    dll0(s) = delta_sum_1d(d, s - delta + 2, s + 1)
    drr0(s - delta) = dll0(s)
end do

do i = delta, n - delta - 1
    do j = i - delta, i - 1
        left0(i, 0) = left0(i, 0) + d(i + 1, j + 1)
    end do
    do j = i + 1, i + delta
        right0(i, 0) = right0(i, 0) + d(i + 1, j + 1)
    end do
    if (i >= 2 * delta - 1) then
        do j = i - 2 * delta + 1, i - delta
            left0(i, 1) = left0(i, 1) + d(i + 1, j + 1)
        end do
    end if
    if (i + 2 * delta - 1 < n) then
        do j = i + delta, i + 2 * delta - 1
            right0(i, 1) = right0(i, 1) + d(i + 1, j + 1)
        end do
    end if
end do

do i = 1, minsize - 1
    do j = minsize, minsize + delta - 1
        dlr0(minsize - 1) = dlr0(minsize - 1) + d(i + 1, j + 1)
    end do
end do

do s = minsize, n - delta - 1
    r1 = left0(s, 0)
    r2 = right0(s - delta, 1)
    r3 = d(s + 1, s + delta + 1)
    a1 = left0(s + delta, 1)
    a2 = right0(s, 0)
    a3 = d(s + 1, s - delta + 1)
    dlr0(s) = dlr0(s - 1) - r1 - r2 - r3 + a1 + a2 + a3
end do

do i = 1, n - 1
    csum0(i) = csum0(i - 1) + d(i + 1, i)
end do

dll = dll0
drr = drr0
dlr = dlr0
csum = csum0

deallocate(dll0, drr0, dlr0, csum0, left0, right0)
end subroutine build_delta_precompute_2d

real(kind=dp) function delta_sum_1d(d, a, b) result(s)
real(kind=dp), intent(in) :: d(:, :)
integer, intent(in) :: a, b
integer :: i, j

s = 0.0_dp
if (b <= a) return
do i = a, b - 1
    do j = i + 1, b
        s = s + d(i, j)
    end do
end do
end function delta_sum_1d

real(kind=dp) function cp3o_delta_seg_stat_1d(u0, s0, t0, delta, dll, drr, dlr, csum) result(stat)
integer, intent(in) :: u0, s0, t0, delta
real(kind=dp), intent(in) :: dll(:), drr(:), dlr(:), csum(:)
integer :: cll, crr, clr
real(kind=dp) :: dllv, drrv, dlrv, num, dnom

cll = delta * (delta - 1) / 2
crr = cll
clr = delta * delta
dllv = dll(s0 + 1)
drrv = drr(s0 + 1)
dlrv = dlr(s0 + 1)

dllv = dllv + csum(s0 - delta + 1)
if (u0 > 0) dllv = dllv - csum(u0)
drrv = drrv + (csum(t0 + 1) - csum(s0 + delta + 1))
cll = cll + (s0 - delta - u0)
crr = crr + (t0 - s0 - delta)

stat = 2.0_dp * dlrv / real(clr, dp) - drrv / real(crr, dp) - dllv / real(cll, dp)
num = real((s0 - u0) * (t0 - s0), dp)
dnom = real((t0 - u0) * (t0 - u0), dp)
stat = stat * (num / dnom)
end function cp3o_delta_seg_stat_1d

subroutine cp3o_delta_backtrack_1d(a, n, kfit, out_row, out_len)
integer, intent(in) :: a(:, :), n, kfit
integer, intent(out) :: out_row(:), out_len
integer, allocatable :: tmp(:)
integer :: cp0, kk, i, j, hold

out_row = 0
allocate(tmp(kfit))
cp0 = a(kfit, n)
kk = kfit
out_len = 0
do while (kk >= 1 .and. cp0 >= 0)
    out_len = out_len + 1
    tmp(out_len) = cp0 + 2
    kk = kk - 1
    if (kk >= 1) cp0 = a(kk, cp0 + 1)
end do
do i = 1, out_len - 1
    do j = i + 1, out_len
        if (tmp(j) < tmp(i)) then
            hold = tmp(i)
            tmp(i) = tmp(j)
            tmp(j) = hold
        end if
    end do
end do
out_row(1:out_len) = tmp(1:out_len)
deallocate(tmp)
end subroutine cp3o_delta_backtrack_1d

subroutine build_distance_matrix_2d(x, alpha, d)
real(kind=dp), intent(in) :: x(:, :)
real(kind=dp), intent(in) :: alpha
real(kind=dp), intent(out) :: d(:, :)
integer :: i, j
real(kind=dp) :: s

d = 0.0_dp
do i = 1, size(x, 1)
    do j = i + 1, size(x, 1)
        s = sqrt(sum((x(i, :) - x(j, :)) ** 2))
        if (alpha /= 1.0_dp) s = s ** alpha
        d(i, j) = s
        d(j, i) = s
    end do
end do
end subroutine build_distance_matrix_2d

subroutine build_prefix_sum_2d(d, pref)
real(kind=dp), intent(in) :: d(:, :)
real(kind=dp), intent(out) :: pref(0:, 0:)
integer :: i, j, n

n = size(d, 1)
pref = 0.0_dp
do i = 1, n
    do j = 1, n
        pref(i, j) = d(i, j) + pref(i - 1, j) + pref(i, j - 1) - pref(i - 1, j - 1)
    end do
end do
end subroutine build_prefix_sum_2d

real(kind=dp) function rect_sum_2d(pref, r1, r2, c1, c2) result(s)
real(kind=dp), intent(in) :: pref(0:, 0:)
integer, intent(in) :: r1, r2, c1, c2
s = pref(r2, c2) - pref(r1 - 1, c2) - pref(r2, c1 - 1) + pref(r1 - 1, c1 - 1)
end function rect_sum_2d

real(kind=dp) function within_sum_1d(pref, a, b) result(s)
real(kind=dp), intent(in) :: pref(0:, 0:)
integer, intent(in) :: a, b
if (b < a) then
    s = 0.0_dp
else
    s = 0.5_dp * rect_sum_2d(pref, a, b, a, b)
end if
end function within_sum_1d

real(kind=dp) function cp3o_seg_stat_1d(uf, ef, tf, pref) result(stat)
integer, intent(in) :: uf, ef, tf
real(kind=dp), intent(in) :: pref(0:, 0:)
integer :: nseg, mseg, tend
real(kind=dp) :: lr, ll, rr

tend = tf
nseg = ef - uf
mseg = tend - ef + 1
ll = within_sum_1d(pref, uf, ef - 1)
rr = within_sum_1d(pref, ef, tend)
lr = rect_sum_2d(pref, uf, ef - 1, ef, tend)
stat = 2.0_dp * lr / real(nseg * mseg, dp) - 2.0_dp * ll / real(nseg * (nseg - 1), dp) - 2.0_dp * rr / real(mseg * (mseg - 1), dp)
stat = stat * real((ef - uf) * (tend - ef + 1), dp) / real((tend - uf + 1) * (tend - uf + 1), dp)
end function cp3o_seg_stat_1d

subroutine cp3o_poi_sse_1d(gofm, split_idx, sse)
real(kind=dp), intent(in) :: gofm(:)
integer, intent(in) :: split_idx
real(kind=dp), intent(out) :: sse
integer :: j, n
real(kind=dp) :: x1_mean, y1_mean, beta1_num, beta1_den, beta1, alpha1
real(kind=dp) :: x2_mean, y2_mean, beta2_num, beta2_den, beta2, alpha2

n = size(gofm)
x1_mean = 0.0_dp
do j = 0, split_idx - 1
    x1_mean = x1_mean + real(j, dp)
end do
x1_mean = x1_mean / real(split_idx, dp)
y1_mean = sum(gofm(1:split_idx)) / real(split_idx, dp)
beta1_num = 0.0_dp
beta1_den = 0.0_dp
do j = 0, split_idx - 1
    beta1_num = beta1_num + (real(j, dp) - x1_mean) * (gofm(j + 1) - y1_mean)
    beta1_den = beta1_den + (real(j, dp) - x1_mean) ** 2
end do
beta1 = beta1_num / beta1_den
alpha1 = y1_mean - beta1 * x1_mean

x2_mean = 0.0_dp
do j = split_idx - 1, n - 1
    x2_mean = x2_mean + real(j, dp)
end do
x2_mean = x2_mean / real(n - split_idx + 1, dp)
y2_mean = sum(gofm(split_idx:n)) / real(n - split_idx + 1, dp)
beta2_num = 0.0_dp
beta2_den = 0.0_dp
do j = split_idx - 1, n - 1
    beta2_num = beta2_num + (real(j, dp) - x2_mean) * (gofm(j + 1) - y2_mean)
    beta2_den = beta2_den + (real(j, dp) - x2_mean) ** 2
end do
beta2 = beta2_num / beta2_den
alpha2 = y2_mean - beta2 * x2_mean

sse = 0.0_dp
do j = 0, split_idx - 1
    sse = sse + (alpha1 + beta1 * real(j, dp) - gofm(j + 1)) ** 2
end do
do j = split_idx - 1, n - 1
    sse = sse + (alpha2 + beta2 * real(j, dp) - gofm(j + 1)) ** 2
end do
end subroutine cp3o_poi_sse_1d

subroutine cp3o_backtrack_1d(a, n, kfit, out_row, out_len)
integer, intent(in) :: a(:, :), n, kfit
integer, intent(out) :: out_row(:), out_len
integer, allocatable :: tmp(:)
integer :: cp, kk, i, j, hold

out_row = 0
allocate(tmp(kfit))
cp = a(kfit, n)
kk = kfit
out_len = 0
do while (kk >= 1 .and. cp > 0)
    out_len = out_len + 1
    tmp(out_len) = cp
    kk = kk - 1
    if (kk >= 1) cp = a(kk, cp - 1)
end do
do i = 1, out_len - 1
    do j = i + 1, out_len
        if (tmp(j) < tmp(i)) then
            hold = tmp(i)
            tmp(i) = tmp(j)
            tmp(j) = hold
        end if
    end do
end do
out_row(1:out_len) = tmp(1:out_len)
deallocate(tmp)
end subroutine cp3o_backtrack_1d

real(kind=dp) function gof_update_agglo_1d(i, dmat, sizes, left, right, fit_curr) result(val)
integer, intent(in) :: i, sizes(:), left(:), right(:)
real(kind=dp), intent(in) :: dmat(:, :), fit_curr
integer :: j, ll, rr, n1, n2, n3
real(kind=dp) :: ktmp

j = right(i)
rr = right(j)
ll = left(i)
val = fit_curr - 2.0_dp * (dmat(i, j) + dmat(i, ll) + dmat(j, rr))
n1 = sizes(i)
n2 = sizes(j)
n3 = sizes(ll)
ktmp = (real(n1 + n3, dp) * dmat(i, ll) + real(n2 + n3, dp) * dmat(j, ll) - real(n3, dp) * dmat(i, j)) / real(n1 + n2 + n3, dp)
val = val + 2.0_dp * ktmp
n3 = sizes(rr)
ktmp = (real(n1 + n3, dp) * dmat(i, rr) + real(n2 + n3, dp) * dmat(j, rr) - real(n3, dp) * dmat(i, j)) / real(n1 + n2 + n3, dp)
val = val + 2.0_dp * ktmp
end function gof_update_agglo_1d

subroutine update_distance_agglo_1d(i, j, k, n, dmat, sizes, left, right, open, lm, progression, merged)
integer, intent(in) :: i, j, k, n
real(kind=dp), intent(inout) :: dmat(:, :)
integer, intent(inout) :: sizes(:), left(:), right(:), lm(:), progression(:, :), merged(:, :)
logical, intent(inout) :: open(:)
integer :: ll, rr, n1, n2, n3, kk, row_idx
real(kind=dp) :: denom

row_idx = k - n + 1
if (i <= n) then
    merged(row_idx, 1) = -i
else
    merged(row_idx, 1) = i - n
end if
if (j <= n) then
    merged(row_idx, 2) = -j
else
    merged(row_idx, 2) = j - n
end if

ll = left(i)
rr = right(j)
left(k + 1) = ll
right(k + 1) = rr
right(ll) = k + 1
left(rr) = k + 1
open(i) = .false.
open(j) = .false.
open(k + 1) = .true.

n1 = sizes(i)
n2 = sizes(j)
sizes(k + 1) = n1 + n2

progression(row_idx + 1, :) = progression(row_idx, :)
progression(row_idx + 1, lm(j)) = 0
lm(k + 1) = lm(i)

do kk = 1, k
    if (open(kk)) then
        n3 = sizes(kk)
        denom = real(n1 + n2 + n3, dp)
        dmat(k + 1, kk) = (real(n1 + n3, dp) * dmat(i, kk) + real(n2 + n3, dp) * dmat(j, kk) - real(n3, dp) * dmat(i, j)) / denom
        dmat(kk, k + 1) = dmat(k + 1, kk)
    end if
end do
end subroutine update_distance_agglo_1d

subroutine best_global_split_1d(changes, d, min_size, best_seg_idx, best_split, best_stat)
integer, intent(in) :: changes(:), min_size
real(kind=dp), intent(in) :: d(:, :)
integer, intent(out) :: best_seg_idx, best_split
real(kind=dp), intent(out) :: best_stat

integer, allocatable :: splits(:)
integer :: i, s, e, cand_split
real(kind=dp) :: cand_stat

best_seg_idx = -1
best_split = -1
best_stat = -huge(1.0_dp)

call sorted_copy_1d(changes, splits)

do i = 1, size(splits) - 1
    s = splits(i)
    e = splits(i + 1) - 1
    call best_segment_split_1d(s, e, d, min_size, cand_split, cand_stat)
    if (cand_stat > best_stat) then
        best_stat = cand_stat
        best_seg_idx = i
        best_split = cand_split
    end if
end do
deallocate(splits)
end subroutine best_global_split_1d

subroutine best_segment_split_1d(s, e, d, min_size, best_split, best_stat)
integer, intent(in) :: s, e, min_size
real(kind=dp), intent(in) :: d(:, :)
integer, intent(out) :: best_split
real(kind=dp), intent(out) :: best_stat

real(kind=dp), allocatable :: b(:), ab(:)
integer :: n, t1, t2, i, j
real(kind=dp) :: a_left, b1, ab1, add_a, add_b, tmp

best_split = -1
best_stat = -huge(1.0_dp)
n = e - s + 1
if (n < 2 * min_size) return

t1 = min_size
t2 = 2 * min_size

a_left = 0.0_dp
do i = s, s + t1 - 1
    do j = i + 1, s + t1 - 1
        a_left = a_left + d(i, j)
    end do
end do

b1 = 0.0_dp
do i = s + t1, s + t2 - 1
    do j = i + 1, s + t2 - 1
        b1 = b1 + d(i, j)
    end do
end do

ab1 = 0.0_dp
do i = s, s + t1 - 1
    do j = s + t1, s + t2 - 1
        ab1 = ab1 + d(i, j)
    end do
end do

tmp = 2.0_dp * ab1 / real((t2 - t1) * t1, dp) - 2.0_dp * b1 / real((t2 - t1 - 1) * (t2 - t1), dp) - 2.0_dp * a_left / real((t1 - 1) * t1, dp)
tmp = tmp * real(t1 * (t2 - t1), dp) / real(t2, dp)
if (tmp > best_stat) then
    best_split = s + t1
    best_stat = tmp
end if

allocate(b(0:n), ab(0:n))
b = b1
ab = ab1

do t2 = 2 * min_size + 1, n
    b(t2) = b(t2 - 1)
    do j = s + t1, s + t2 - 2
        b(t2) = b(t2) + d(s + t2 - 1, j)
    end do
    ab(t2) = ab(t2 - 1)
    do j = s, s + t1 - 1
        ab(t2) = ab(t2) + d(s + t2 - 1, j)
    end do
    tmp = 2.0_dp * ab(t2) / real((t2 - t1) * t1, dp) - 2.0_dp * b(t2) / real((t2 - t1 - 1) * (t2 - t1), dp) - 2.0_dp * a_left / real(t1 * (t1 - 1), dp)
    tmp = tmp * real(t1 * (t2 - t1), dp) / real(t2, dp)
    if (tmp > best_stat) then
        best_split = s + t1
        best_stat = tmp
    end if
end do

do t1 = min_size + 1, n
    t2 = t1 + min_size
    if (t2 > n) exit
    add_a = 0.0_dp
    do j = s, s + t1 - 2
        add_a = add_a + d(s + t1 - 1, j)
    end do
    a_left = a_left + add_a

    add_b = 0.0_dp
    do j = s + t1, s + t2 - 2
        add_b = add_b + d(s + t1 - 1, j)
    end do

    do while (t2 <= n)
        add_b = add_b + d(s + t1 - 1, s + t2 - 1)
        b(t2) = b(t2) - add_b
        ab(t2) = ab(t2) + (add_b - add_a)
        tmp = 2.0_dp * ab(t2) / real((t2 - t1) * t1, dp) - 2.0_dp * b(t2) / real((t2 - t1 - 1) * (t2 - t1), dp) - 2.0_dp * a_left / real((t1 - 1) * t1, dp)
        tmp = tmp * real(t1 * (t2 - t1), dp) / real(t2, dp)
        if (tmp > best_stat) then
            best_split = s + t1
            best_stat = tmp
        end if
        t2 = t2 + 1
    end do
end do

deallocate(b, ab)
end subroutine best_segment_split_1d

subroutine insert_change_1d(changes, val)
integer, allocatable, intent(inout) :: changes(:)
integer, intent(in) :: val
integer, allocatable :: tmp(:)
integer :: pos

allocate(tmp(size(changes) + 1))
pos = count(changes < val) + 1
if (pos > 1) tmp(1:pos-1) = changes(1:pos-1)
tmp(pos) = val
if (pos <= size(changes)) tmp(pos+1:) = changes(pos:)
call move_alloc(tmp, changes)
end subroutine insert_change_1d

subroutine append_change_unsorted_1d(changes, val)
integer, allocatable, intent(inout) :: changes(:)
integer, intent(in) :: val
integer, allocatable :: tmp(:)

allocate(tmp(size(changes) + 1))
tmp(1:size(changes)) = changes
tmp(size(tmp)) = val
call move_alloc(tmp, changes)
end subroutine append_change_unsorted_1d

subroutine sorted_copy_1d(vals, sorted_vals)
integer, intent(in) :: vals(:)
integer, allocatable, intent(out) :: sorted_vals(:)
integer :: i, j, tmp

allocate(sorted_vals(size(vals)))
sorted_vals = vals
do i = 1, size(sorted_vals) - 1
    do j = i + 1, size(sorted_vals)
        if (sorted_vals(j) < sorted_vals(i)) then
            tmp = sorted_vals(i)
            sorted_vals(i) = sorted_vals(j)
            sorted_vals(j) = tmp
        end if
    end do
end do
end subroutine sorted_copy_1d

subroutine permute_distance_1d(d, changes, perm, d1)
real(kind=dp), intent(in) :: d(:, :)
integer, intent(in) :: changes(:)
integer, intent(in) :: perm(:)
real(kind=dp), intent(out) :: d1(:, :)
integer, allocatable :: splits(:)
integer :: i, j, k, s, e

d1 = d
call sorted_copy_1d(changes, splits)
do k = 1, size(splits) - 1
    s = splits(k)
    e = splits(k + 1) - 1
    do i = s, e
        do j = s, e
            d1(i, j) = d(perm(i), perm(j))
        end do
    end do
end do
deallocate(splits)
end subroutine permute_distance_1d

subroutine build_cluster_1d(estimates, cluster)
integer, intent(in) :: estimates(:)
integer, allocatable, intent(out) :: cluster(:)
integer :: n, i, a, b, idx

n = estimates(size(estimates)) - 1
allocate(cluster(n))
idx = 1
do i = 1, size(estimates) - 1
    a = estimates(i)
    b = estimates(i + 1) - 1
    cluster(a:b) = idx
    idx = idx + 1
end do
end subroutine build_cluster_1d

subroutine build_cluster_from_estimates_1d(estimates, n, cluster)
integer, intent(in) :: estimates(:), n
integer, allocatable, intent(out) :: cluster(:)
integer, allocatable :: tmp(:)
integer :: k

if (estimates(1) == 1) then
    call build_cluster_1d(estimates, cluster)
else
    allocate(tmp(size(estimates) + 1))
    tmp(1) = 1
    tmp(2:) = estimates
    call build_cluster_1d(tmp, cluster)
    k = n - size(cluster)
    if (k > 0) cluster = [cluster, spread(1, 1, k)]
    deallocate(tmp)
end if
end subroutine build_cluster_from_estimates_1d

end module ecp_pkg_mod
