module rbeast_mod
use kind_mod, only: dp
implicit none
private
public :: solve_beast_trend_only_1d, solve_beast_harmonic_1d, solve_beast_irreg_trend_only_1d, &
          solve_beast123_regular_mv, top_cp_probabilities

contains

subroutine solve_beast_trend_only_1d(y, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, fitted_mean, model_weights)
    !> Trend-only BEAST-style model averaging over piecewise-linear segmentations.
    real(kind=dp), intent(in) :: y(:)
    integer, intent(in) :: max_cp
    integer, intent(in) :: min_seg_len
    integer, intent(out) :: best_cps(max_cp)
    integer, intent(out) :: ncp_best
    real(kind=dp), intent(out) :: cp_prob(size(y))
    real(kind=dp), intent(out) :: fitted_mean(size(y))
    real(kind=dp), intent(out) :: model_weights(max_cp + 1)

    integer :: n, max_m, m, i, k, best_model, ncp_model
    real(kind=dp), allocatable :: sx(:), sx2(:), sy(:), sy2(:), sxy(:)
    real(kind=dp), allocatable :: dp_cost(:, :), bic(:), fitted(:, :)
    integer, allocatable :: parent(:, :), cps(:, :)
    real(kind=dp) :: rss, bic_min, wsum, nreal

    n = size(y)
    max_m = max_cp + 1
    nreal = real(n, dp)
    best_cps = 0
    cp_prob = 0.0_dp
    fitted_mean = 0.0_dp
    model_weights = 0.0_dp
    ncp_best = 0

    call build_linear_prefix(y, sx, sx2, sy, sy2, sxy)
    allocate(dp_cost(n, max_m), parent(n, max_m), bic(max_m), fitted(n, max_m), cps(max_cp, max_m))
    dp_cost = huge(1.0_dp)
    parent = 0
    bic = huge(1.0_dp)
    fitted = 0.0_dp
    cps = 0

    do i = min_seg_len, n
        dp_cost(i, 1) = segment_linear_rss(1, i, sx, sx2, sy, sy2, sxy)
    end do

    do m = 2, max_m
        do i = m * min_seg_len, n
            do k = (m - 1) * min_seg_len, i - min_seg_len
                rss = dp_cost(k, m - 1) + segment_linear_rss(k + 1, i, sx, sx2, sy, sy2, sxy)
                if (rss < dp_cost(i, m)) then
                    dp_cost(i, m) = rss
                    parent(i, m) = k
                end if
            end do
        end do
    end do

    do m = 1, max_m
        if (.not. isfinite(dp_cost(n, m))) cycle
        ncp_model = m - 1
        if (ncp_model > 0) call backtrack_cps(parent, n, m, cps(1:ncp_model, m))
        call fit_piecewise_linear(y, cps(1:ncp_model, m), fitted(:, m))
        rss = max(sum((y - fitted(:, m))**2), 1.0e-12_dp)
        bic(m) = nreal * log(rss / nreal) + real(2 * m, dp) * log(nreal)
    end do

    bic_min = minval(bic)
    do m = 1, max_m
        if (isfinite(bic(m))) model_weights(m) = exp(-0.5_dp * (bic(m) - bic_min))
    end do
    wsum = sum(model_weights)
    if (wsum > 0.0_dp) model_weights = model_weights / wsum

    best_model = maxloc(model_weights, dim=1)
    ncp_best = best_model - 1
    if (ncp_best > 0) best_cps(1:ncp_best) = cps(1:ncp_best, best_model)

    do m = 1, max_m
        if (model_weights(m) <= 0.0_dp) cycle
        fitted_mean = fitted_mean + model_weights(m) * fitted(:, m)
        ncp_model = m - 1
        do k = 1, ncp_model
            cp_prob(cps(k, m)) = cp_prob(cps(k, m)) + model_weights(m)
        end do
    end do

    deallocate(sx, sx2, sy, sy2, sxy, dp_cost, parent, bic, fitted, cps)
end subroutine solve_beast_trend_only_1d

subroutine solve_beast_harmonic_1d(y, period, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, fitted_mean, seasonal_fit, &
                                   deseasonalized, model_weights, beta_season)
    !> Harmonic-seasonal plus piecewise-linear trend BEAST-style approximation.
    real(kind=dp), intent(in) :: y(:)
    real(kind=dp), intent(in) :: period
    integer, intent(in) :: max_cp
    integer, intent(in) :: min_seg_len
    integer, intent(out) :: best_cps(max_cp)
    integer, intent(out) :: ncp_best
    real(kind=dp), intent(out) :: cp_prob(size(y))
    real(kind=dp), intent(out) :: fitted_mean(size(y))
    real(kind=dp), intent(out) :: seasonal_fit(size(y))
    real(kind=dp), intent(out) :: deseasonalized(size(y))
    real(kind=dp), intent(out) :: model_weights(max_cp + 1)
    real(kind=dp), intent(out) :: beta_season(2)

    integer :: i, n
    real(kind=dp) :: omega

    n = size(y)
    omega = 2.0_dp * acos(-1.0_dp) / period
    call fit_harmonic_component(y, omega, beta_season)
    do i = 1, n
        seasonal_fit(i) = beta_season(1) * sin(omega * real(i, dp)) + beta_season(2) * cos(omega * real(i, dp))
        deseasonalized(i) = y(i) - seasonal_fit(i)
    end do
    call solve_beast_trend_only_1d(deseasonalized, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, fitted_mean, model_weights)
    fitted_mean = fitted_mean + seasonal_fit
end subroutine solve_beast_harmonic_1d

subroutine solve_beast_irreg_trend_only_1d(t, y, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, fitted_mean, model_weights)
    !> Irregular-time trend-only BEAST-style model averaging over piecewise-linear segmentations.
    real(kind=dp), intent(in) :: t(:), y(:)
    integer, intent(in) :: max_cp
    integer, intent(in) :: min_seg_len
    integer, intent(out) :: best_cps(max_cp)
    integer, intent(out) :: ncp_best
    real(kind=dp), intent(out) :: cp_prob(size(y))
    real(kind=dp), intent(out) :: fitted_mean(size(y))
    real(kind=dp), intent(out) :: model_weights(max_cp + 1)

    integer :: n, max_m, m, i, k, best_model, ncp_model
    real(kind=dp), allocatable :: sx(:), sx2(:), sy(:), sy2(:), sxy(:)
    real(kind=dp), allocatable :: dp_cost(:, :), bic(:), fitted(:, :)
    integer, allocatable :: parent(:, :), cps(:, :)
    real(kind=dp) :: rss, bic_min, wsum, nreal

    n = size(y)
    max_m = max_cp + 1
    nreal = real(n, dp)
    best_cps = 0
    cp_prob = 0.0_dp
    fitted_mean = 0.0_dp
    model_weights = 0.0_dp
    ncp_best = 0

    call build_linear_prefix_irreg(t, y, sx, sx2, sy, sy2, sxy)
    allocate(dp_cost(n, max_m), parent(n, max_m), bic(max_m), fitted(n, max_m), cps(max_cp, max_m))
    dp_cost = huge(1.0_dp)
    parent = 0
    bic = huge(1.0_dp)
    fitted = 0.0_dp
    cps = 0

    do i = min_seg_len, n
        dp_cost(i, 1) = segment_linear_rss_irreg(1, i, sx, sx2, sy, sy2, sxy)
    end do

    do m = 2, max_m
        do i = m * min_seg_len, n
            do k = (m - 1) * min_seg_len, i - min_seg_len
                rss = dp_cost(k, m - 1) + segment_linear_rss_irreg(k + 1, i, sx, sx2, sy, sy2, sxy)
                if (rss < dp_cost(i, m)) then
                    dp_cost(i, m) = rss
                    parent(i, m) = k
                end if
            end do
        end do
    end do

    do m = 1, max_m
        if (.not. isfinite(dp_cost(n, m))) cycle
        ncp_model = m - 1
        if (ncp_model > 0) call backtrack_cps(parent, n, m, cps(1:ncp_model, m))
        call fit_piecewise_linear_irreg(t, y, cps(1:ncp_model, m), fitted(:, m))
        rss = max(sum((y - fitted(:, m))**2), 1.0e-12_dp)
        bic(m) = nreal * log(rss / nreal) + real(2 * m, dp) * log(nreal)
    end do

    bic_min = minval(bic)
    do m = 1, max_m
        if (isfinite(bic(m))) model_weights(m) = exp(-0.5_dp * (bic(m) - bic_min))
    end do
    wsum = sum(model_weights)
    if (wsum > 0.0_dp) model_weights = model_weights / wsum

    best_model = maxloc(model_weights, dim=1)
    ncp_best = best_model - 1
    if (ncp_best > 0) best_cps(1:ncp_best) = cps(1:ncp_best, best_model)

    do m = 1, max_m
        if (model_weights(m) <= 0.0_dp) cycle
        fitted_mean = fitted_mean + model_weights(m) * fitted(:, m)
        ncp_model = m - 1
        do k = 1, ncp_model
            cp_prob(cps(k, m)) = cp_prob(cps(k, m)) + model_weights(m)
        end do
    end do

    deallocate(sx, sx2, sy, sy2, sxy, dp_cost, parent, bic, fitted, cps)
end subroutine solve_beast_irreg_trend_only_1d

subroutine solve_beast123_regular_mv(y, max_cp, min_seg_len, best_cps, ncp_best, cp_prob, mean_cp_prob, fitted_mean, model_weights)
    !> Apply the regular trend-only BEAST-style solver column-wise to multiple series.
    real(kind=dp), intent(in) :: y(:, :)
    integer, intent(in) :: max_cp
    integer, intent(in) :: min_seg_len
    integer, intent(out) :: best_cps(max_cp, size(y, 2))
    integer, intent(out) :: ncp_best(size(y, 2))
    real(kind=dp), intent(out) :: cp_prob(size(y, 1), size(y, 2))
    real(kind=dp), intent(out) :: mean_cp_prob(size(y, 1))
    real(kind=dp), intent(out) :: fitted_mean(size(y, 1), size(y, 2))
    real(kind=dp), intent(out) :: model_weights(max_cp + 1, size(y, 2))

    integer :: j

    best_cps = 0
    ncp_best = 0
    cp_prob = 0.0_dp
    fitted_mean = 0.0_dp
    model_weights = 0.0_dp
    do j = 1, size(y, 2)
        call solve_beast_trend_only_1d(y(:, j), max_cp, min_seg_len, best_cps(:, j), ncp_best(j), cp_prob(:, j), &
                                       fitted_mean(:, j), model_weights(:, j))
    end do
    mean_cp_prob = sum(cp_prob, dim=2) / real(size(y, 2), dp)
end subroutine solve_beast123_regular_mv

subroutine top_cp_probabilities(cp_prob, top_idx, top_prob)
    !> Extract changepoint locations with the largest occurrence probabilities.
    real(kind=dp), intent(in) :: cp_prob(:)
    integer, intent(out) :: top_idx(:)
    real(kind=dp), intent(out) :: top_prob(:)

    integer :: n, i, j, pos
    real(kind=dp), allocatable :: work(:)

    n = size(cp_prob)
    allocate(work(n))
    work = cp_prob
    top_idx = 0
    top_prob = 0.0_dp
    do i = 1, size(top_idx)
        pos = maxloc(work, dim=1)
        if (work(pos) <= 0.0_dp) exit
        top_idx(i) = pos
        top_prob(i) = work(pos)
        work(pos) = -1.0_dp
    end do

    do i = 1, size(top_idx) - 1
        do j = i + 1, size(top_idx)
            if (top_idx(j) > 0 .and. (top_idx(i) == 0 .or. top_idx(j) < top_idx(i))) then
                call swap_int(top_idx(i), top_idx(j))
                call swap_real(top_prob(i), top_prob(j))
            end if
        end do
    end do
    deallocate(work)
end subroutine top_cp_probabilities

subroutine build_linear_prefix(y, sx, sx2, sy, sy2, sxy)
    real(kind=dp), intent(in) :: y(:)
    real(kind=dp), allocatable, intent(out) :: sx(:), sx2(:), sy(:), sy2(:), sxy(:)

    integer :: n, i
    real(kind=dp) :: x

    n = size(y)
    allocate(sx(0:n), sx2(0:n), sy(0:n), sy2(0:n), sxy(0:n))
    sx = 0.0_dp
    sx2 = 0.0_dp
    sy = 0.0_dp
    sy2 = 0.0_dp
    sxy = 0.0_dp
    do i = 1, n
        x = real(i, dp)
        sx(i) = sx(i - 1) + x
        sx2(i) = sx2(i - 1) + x * x
        sy(i) = sy(i - 1) + y(i)
        sy2(i) = sy2(i - 1) + y(i) * y(i)
        sxy(i) = sxy(i - 1) + x * y(i)
    end do
end subroutine build_linear_prefix

subroutine build_linear_prefix_irreg(t, y, sx, sx2, sy, sy2, sxy)
    real(kind=dp), intent(in) :: t(:), y(:)
    real(kind=dp), allocatable, intent(out) :: sx(:), sx2(:), sy(:), sy2(:), sxy(:)

    integer :: n, i

    n = size(y)
    allocate(sx(0:n), sx2(0:n), sy(0:n), sy2(0:n), sxy(0:n))
    sx = 0.0_dp
    sx2 = 0.0_dp
    sy = 0.0_dp
    sy2 = 0.0_dp
    sxy = 0.0_dp
    do i = 1, n
        sx(i) = sx(i - 1) + t(i)
        sx2(i) = sx2(i - 1) + t(i) * t(i)
        sy(i) = sy(i - 1) + y(i)
        sy2(i) = sy2(i - 1) + y(i) * y(i)
        sxy(i) = sxy(i - 1) + t(i) * y(i)
    end do
end subroutine build_linear_prefix_irreg

real(kind=dp) function segment_linear_rss(a, b, sx, sx2, sy, sy2, sxy) result(rss)
    integer, intent(in) :: a, b
    real(kind=dp), intent(in) :: sx(0:), sx2(0:), sy(0:), sy2(0:), sxy(0:)

    real(kind=dp) :: nseg, sumx, sumx2, sumy, sumy2, sumxy, denom, beta0, beta1

    nseg = real(b - a + 1, dp)
    sumx = sx(b) - sx(a - 1)
    sumx2 = sx2(b) - sx2(a - 1)
    sumy = sy(b) - sy(a - 1)
    sumy2 = sy2(b) - sy2(a - 1)
    sumxy = sxy(b) - sxy(a - 1)
    denom = nseg * sumx2 - sumx * sumx

    if (abs(denom) <= 1.0e-12_dp) then
        beta1 = 0.0_dp
        beta0 = sumy / nseg
    else
        beta1 = (nseg * sumxy - sumx * sumy) / denom
        beta0 = (sumy - beta1 * sumx) / nseg
    end if
    rss = max(sumy2 - beta0 * sumy - beta1 * sumxy, 0.0_dp)
end function segment_linear_rss

real(kind=dp) function segment_linear_rss_irreg(a, b, sx, sx2, sy, sy2, sxy) result(rss)
    integer, intent(in) :: a, b
    real(kind=dp), intent(in) :: sx(0:), sx2(0:), sy(0:), sy2(0:), sxy(0:)

    real(kind=dp) :: nseg, sumx, sumx2, sumy, sumy2, sumxy, denom, beta0, beta1

    nseg = real(b - a + 1, dp)
    sumx = sx(b) - sx(a - 1)
    sumx2 = sx2(b) - sx2(a - 1)
    sumy = sy(b) - sy(a - 1)
    sumy2 = sy2(b) - sy2(a - 1)
    sumxy = sxy(b) - sxy(a - 1)
    denom = nseg * sumx2 - sumx * sumx

    if (abs(denom) <= 1.0e-12_dp) then
        beta1 = 0.0_dp
        beta0 = sumy / nseg
    else
        beta1 = (nseg * sumxy - sumx * sumy) / denom
        beta0 = (sumy - beta1 * sumx) / nseg
    end if
    rss = max(sumy2 - beta0 * sumy - beta1 * sumxy, 0.0_dp)
end function segment_linear_rss_irreg

subroutine backtrack_cps(parent, n, m, cps)
    integer, intent(in) :: parent(:, :)
    integer, intent(in) :: n, m
    integer, intent(out) :: cps(:)

    integer :: idx, pos, t

    cps = 0
    t = n
    do idx = m, 2, -1
        pos = parent(t, idx)
        cps(idx - 1) = pos
        t = pos
    end do
end subroutine backtrack_cps

subroutine fit_piecewise_linear(y, cps, fitted)
    real(kind=dp), intent(in) :: y(:)
    integer, intent(in) :: cps(:)
    real(kind=dp), intent(out) :: fitted(:)

    integer :: start_idx, end_idx, seg, ncp
    real(kind=dp), allocatable :: sx(:), sx2(:), sy(:), sy2(:), sxy(:)

    ncp = count(cps > 0)
    call build_linear_prefix(y, sx, sx2, sy, sy2, sxy)
    start_idx = 1
    do seg = 1, ncp + 1
        if (seg <= ncp) then
            end_idx = cps(seg)
        else
            end_idx = size(y)
        end if
        call fit_segment(start_idx, end_idx, sx, sx2, sy, sxy, fitted)
        start_idx = end_idx + 1
    end do
    deallocate(sx, sx2, sy, sy2, sxy)
end subroutine fit_piecewise_linear

subroutine fit_piecewise_linear_irreg(t, y, cps, fitted)
    real(kind=dp), intent(in) :: t(:), y(:)
    integer, intent(in) :: cps(:)
    real(kind=dp), intent(out) :: fitted(:)

    integer :: start_idx, end_idx, seg, ncp

    ncp = count(cps > 0)
    start_idx = 1
    do seg = 1, ncp + 1
        if (seg <= ncp) then
            end_idx = cps(seg)
        else
            end_idx = size(y)
        end if
        call fit_segment_irreg(start_idx, end_idx, t, y, fitted)
        start_idx = end_idx + 1
    end do
end subroutine fit_piecewise_linear_irreg

subroutine fit_segment(a, b, sx, sx2, sy, sxy, fitted)
    integer, intent(in) :: a, b
    real(kind=dp), intent(in) :: sx(0:), sx2(0:), sy(0:), sxy(0:)
    real(kind=dp), intent(inout) :: fitted(:)

    integer :: i
    real(kind=dp) :: nseg, sumx, sumx2, sumy, sumxy, denom, beta0, beta1

    nseg = real(b - a + 1, dp)
    sumx = sx(b) - sx(a - 1)
    sumx2 = sx2(b) - sx2(a - 1)
    sumy = sy(b) - sy(a - 1)
    sumxy = sxy(b) - sxy(a - 1)
    denom = nseg * sumx2 - sumx * sumx

    if (abs(denom) <= 1.0e-12_dp) then
        beta1 = 0.0_dp
        beta0 = sumy / nseg
    else
        beta1 = (nseg * sumxy - sumx * sumy) / denom
        beta0 = (sumy - beta1 * sumx) / nseg
    end if
    do i = a, b
        fitted(i) = beta0 + beta1 * real(i, dp)
    end do
end subroutine fit_segment

subroutine fit_segment_irreg(a, b, t, y, fitted)
    integer, intent(in) :: a, b
    real(kind=dp), intent(in) :: t(:), y(:)
    real(kind=dp), intent(inout) :: fitted(:)

    real(kind=dp) :: xm, ym, denom, beta0, beta1

    xm = sum(t(a:b)) / real(b - a + 1, dp)
    ym = sum(y(a:b)) / real(b - a + 1, dp)
    denom = sum((t(a:b) - xm)**2)
    if (denom <= 1.0e-12_dp) then
        beta1 = 0.0_dp
        beta0 = ym
    else
        beta1 = sum((t(a:b) - xm) * (y(a:b) - ym)) / denom
        beta0 = ym - beta1 * xm
    end if
    fitted(a:b) = beta0 + beta1 * t(a:b)
end subroutine fit_segment_irreg

logical function isfinite(x)
    real(kind=dp), intent(in) :: x
    isfinite = (x < huge(x) .and. x > -huge(x))
end function isfinite

subroutine swap_int(a, b)
    integer, intent(inout) :: a, b
    integer :: tmp
    tmp = a
    a = b
    b = tmp
end subroutine swap_int

subroutine swap_real(a, b)
    real(kind=dp), intent(inout) :: a, b
    real(kind=dp) :: tmp
    tmp = a
    a = b
    b = tmp
end subroutine swap_real

subroutine fit_harmonic_component(y, omega, beta)
    real(kind=dp), intent(in) :: y(:)
    real(kind=dp), intent(in) :: omega
    real(kind=dp), intent(out) :: beta(2)

    integer :: i
    real(kind=dp) :: s11, s12, s22, rhs1, rhs2, det, s, c

    s11 = 0.0_dp
    s12 = 0.0_dp
    s22 = 0.0_dp
    rhs1 = 0.0_dp
    rhs2 = 0.0_dp
    do i = 1, size(y)
        s = sin(omega * real(i, dp))
        c = cos(omega * real(i, dp))
        s11 = s11 + s * s
        s12 = s12 + s * c
        s22 = s22 + c * c
        rhs1 = rhs1 + s * y(i)
        rhs2 = rhs2 + c * y(i)
    end do
    det = s11 * s22 - s12 * s12
    if (abs(det) <= 1.0e-12_dp) then
        beta = 0.0_dp
    else
        beta(1) = (rhs1 * s22 - rhs2 * s12) / det
        beta(2) = (s11 * rhs2 - s12 * rhs1) / det
    end if
end subroutine fit_harmonic_component

end module rbeast_mod
