module strucchange_pkg_mod
    use kind_mod, only: dp
    use pca_jacobi_mod, only: jacobi_eigen_sym
    implicit none
    private
    public :: solve_breakpoints_regression_1d, solve_fstats_regression_1d, solve_efp_ols_cusum_1d, solve_efp_ols_mosum_1d, solve_efp_rec_cusum_1d, solve_efp_rec_mosum_1d, solve_efp_re_1d, solve_efp_me_1d, solve_mefp_ols_cusum_1d, solve_mefp_ols_mosum_1d, solve_mefp_re_1d, solve_mefp_me_1d

contains

    subroutine solve_breakpoints_regression_1d(y, x, h, cpt_hat, rss_vec, bic_vec, best_m)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        integer, intent(in) :: h
        integer, allocatable, intent(out) :: cpt_hat(:)
        real(kind=dp), allocatable, intent(out) :: rss_vec(:), bic_vec(:)
        integer, intent(out) :: best_m

        integer :: n, p, k, max_breaks, m, i, j, best_i, ncpt, t
        real(kind=dp) :: best_val, total_rss
        real(kind=dp), allocatable :: rss_seg(:, :), dp_rss(:, :)
        integer, allocatable :: prev_break(:, :), last_break(:), cps_work(:)

        n = size(y)
        p = size(x, 2)
        k = p + 1
        if (h <= k) error stop "solve_breakpoints_regression_1d: h must be greater than number of regressors"
        if (h > n / 2) error stop "solve_breakpoints_regression_1d: h must be smaller than half the sample size"

        max_breaks = ceiling_div_int_1d(n, h) - 2
        if (max_breaks < 0) max_breaks = 0

        allocate(rss_vec(0:max_breaks), bic_vec(0:max_breaks), rss_seg(n, n))
        rss_seg = huge(1.0_dp)
        do i = 1, n
            do j = i, n
                if (j - i + 1 >= h) rss_seg(i, j) = regression_segment_rss_1d(y, x, i, j)
            end do
        end do

        rss_vec(0) = rss_seg(1, n)
        bic_vec(0) = bic_value_1d(n, k, 0, rss_vec(0))

        if (max_breaks == 0) then
            best_m = 0
            allocate(cpt_hat(0))
            return
        end if

        allocate(dp_rss(max_breaks, n), prev_break(max_breaks, n), last_break(max_breaks), cps_work(max_breaks))
        dp_rss = huge(1.0_dp)
        prev_break = 0
        last_break = 0

        do i = h, n - h
            dp_rss(1, i) = rss_seg(1, i)
        end do

        do m = 2, max_breaks
            do i = m * h, n - h
                best_val = huge(1.0_dp)
                best_i = 0
                do j = (m - 1) * h, i - h
                    total_rss = dp_rss(m - 1, j) + rss_seg(j + 1, i)
                    if (total_rss < best_val) then
                        best_val = total_rss
                        best_i = j
                    end if
                end do
                dp_rss(m, i) = best_val
                prev_break(m, i) = best_i
            end do
        end do

        do m = 1, max_breaks
            best_val = huge(1.0_dp)
            best_i = 0
            do i = m * h, n - h
                total_rss = dp_rss(m, i) + rss_seg(i + 1, n)
                if (total_rss < best_val) then
                    best_val = total_rss
                    best_i = i
                end if
            end do
            rss_vec(m) = best_val
            bic_vec(m) = bic_value_1d(n, k, m, best_val)
            last_break(m) = best_i
        end do

        best_m = 0
        best_val = bic_vec(0)
        do m = 1, max_breaks
            if (bic_vec(m) < best_val) then
                best_val = bic_vec(m)
                best_m = m
            end if
        end do

        if (best_m == 0) then
            allocate(cpt_hat(0))
        else
            ncpt = best_m
            cps_work = 0
            cps_work(best_m) = last_break(best_m)
            do m = best_m, 2, -1
                cps_work(m - 1) = prev_break(m, cps_work(m))
            end do
            allocate(cpt_hat(ncpt))
            do t = 1, ncpt
                cpt_hat(t) = cps_work(t)
            end do
        end if

        deallocate(dp_rss, prev_break, last_break, cps_work)
    end subroutine solve_breakpoints_regression_1d

    subroutine solve_fstats_regression_1d(y, x, from_in, to_in, fstats, point_idx, breakpoint, min_rss)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: from_in
        real(kind=dp), intent(in), optional :: to_in
        real(kind=dp), allocatable, intent(out) :: fstats(:)
        integer, allocatable, intent(out) :: point_idx(:)
        integer, intent(out) :: breakpoint
        real(kind=dp), intent(out) :: min_rss

        integer :: n, p, k, from_idx, to_idx, np, i, point
        real(kind=dp) :: sume2, ess, sigma2, max_f, to_use
        real(kind=dp), allocatable :: resid_full(:), fitted(:)

        n = size(y)
        p = size(x, 2)
        k = p + 1
        to_use = merge(to_in, -1.0_dp, present(to_in))

        call ols_fit_resid_1d(y, x, 1, n, fitted, resid_full)
        sume2 = sum(resid_full ** 2)
        deallocate(fitted, resid_full)

        if (from_in < 1.0_dp) then
            from_idx = floor(from_in * real(n, dp))
            if (present(to_in)) then
                to_idx = floor(to_use * real(n, dp))
            else
                to_idx = n - from_idx
            end if
        else
            from_idx = int(from_in)
            if (present(to_in)) then
                to_idx = int(to_use)
            else
                to_idx = n - from_idx
            end if
        end if

        if (from_idx < (k + 1)) from_idx = k + 1
        if (to_idx > (n - k - 1)) to_idx = n - k - 1
        if (from_idx > to_idx) error stop "solve_fstats_regression_1d: inadmissible from/to"

        np = to_idx - from_idx + 1
        allocate(fstats(np), point_idx(np))
        do i = 1, np
            point = from_idx + i - 1
            point_idx(i) = point
            ess = regression_segment_rss_1d(y, x, 1, point) + regression_segment_rss_1d(y, x, point + 1, n)
            sigma2 = ess / real(n - 2 * k, dp)
            fstats(i) = (sume2 - ess) / sigma2
        end do

        max_f = maxval(fstats)
        breakpoint = point_idx(maxloc(fstats, dim=1))
        min_rss = sume2 / (1.0_dp + max_f / real(n - 2 * k, dp))
    end subroutine solve_fstats_regression_1d

    subroutine solve_efp_ols_cusum_1d(y, x, process)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), allocatable, intent(out) :: process(:)

        integer :: n, p, df_resid, i
        real(kind=dp) :: sigma
        real(kind=dp), allocatable :: resid(:), fitted(:)

        n = size(y)
        p = size(x, 2)
        df_resid = n - (p + 1)
        call ols_fit_resid_1d(y, x, 1, n, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(df_resid, dp))

        allocate(process(n + 1))
        process(1) = 0.0_dp
        do i = 1, n
            process(i + 1) = process(i) + resid(i)
        end do
        process = process / (sigma * sqrt(real(n, dp)))

        deallocate(resid, fitted)
    end subroutine solve_efp_ols_cusum_1d

    subroutine solve_efp_ols_mosum_1d(y, x, h, process)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: h
        real(kind=dp), allocatable, intent(out) :: process(:)

        integer :: n, p, df_resid, nh, i
        real(kind=dp) :: sigma
        real(kind=dp), allocatable :: resid(:), fitted(:), cs(:)

        n = size(y)
        p = size(x, 2)
        df_resid = n - (p + 1)
        nh = floor(real(n, dp) * h)
        if (nh < 1) error stop "solve_efp_ols_mosum_1d: h too small"
        if (nh > n) error stop "solve_efp_ols_mosum_1d: h too large"

        call ols_fit_resid_1d(y, x, 1, n, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(df_resid, dp))

        allocate(cs(n + 1))
        cs(1) = 0.0_dp
        do i = 1, n
            cs(i + 1) = cs(i) + resid(i)
        end do

        allocate(process(n - nh + 1))
        do i = 1, n - nh + 1
            process(i) = (cs(i + nh) - cs(i)) / (sigma * sqrt(real(n, dp)))
        end do

        deallocate(resid, fitted, cs)
    end subroutine solve_efp_ols_mosum_1d

    subroutine solve_efp_rec_cusum_1d(y, x, process)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), allocatable, intent(out) :: process(:)

        integer :: n, p, k, nw, i
        real(kind=dp) :: sigma
        real(kind=dp), allocatable :: w(:)

        n = size(y)
        p = size(x, 2)
        k = p + 1
        call recresid_ols_1d(y, x, w)
        nw = size(w)
        sigma = sample_sd_1d(w)

        allocate(process(nw + 1))
        process(1) = 0.0_dp
        do i = 1, nw
            process(i + 1) = process(i) + w(i)
        end do
        process = process / (sigma * sqrt(real(n - k, dp)))

        deallocate(w)
    end subroutine solve_efp_rec_cusum_1d

    subroutine solve_efp_rec_mosum_1d(y, x, h, process)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: h
        real(kind=dp), allocatable, intent(out) :: process(:)

        integer :: n, p, k, nw, nh, i
        real(kind=dp) :: sigma
        real(kind=dp), allocatable :: w(:)

        n = size(y)
        p = size(x, 2)
        k = p + 1
        call recresid_ols_1d(y, x, w)
        nw = size(w)
        nh = floor(real(nw, dp) * h)
        if (nh < 1) error stop "solve_efp_rec_mosum_1d: h too small"
        if (nh >= nw) error stop "solve_efp_rec_mosum_1d: h too large"

        sigma = sqrt(sum((w - sum(w) / real(nw, dp)) ** 2) / real(nw - k, dp))
        allocate(process(nw - nh + 1))
        do i = 1, nw - nh + 1
            process(i) = sum(w(i:i + nh - 1)) / (sigma * sqrt(real(nw, dp)))
        end do

        deallocate(w)
    end subroutine solve_efp_rec_mosum_1d

    subroutine solve_efp_re_1d(y, x, process)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), allocatable, intent(out) :: process(:, :)

        integer :: n, p, k, i, info
        real(kind=dp) :: sigma
        real(kind=dp), allocatable :: fitted(:), resid(:), beta_hat(:), beta_i(:), xtx(:, :), xty(:), q12(:, :), qi12(:, :), tmp(:), work(:, :)

        n = size(y)
        p = size(x, 2)
        k = p + 1

        call ols_fit_resid_1d(y, x, 1, n, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(n - k, dp))

        allocate(xtx(k, k), xty(k), beta_hat(k))
        call normal_equations_1d(y, x, 1, n, xtx, xty)
        call solve_linear_system_square_1d(xtx, xty, beta_hat, info)
        if (info /= 0) error stop "solve_efp_re_1d: singular full fit"
        call symmetric_root_1d(xtx, q12)
        q12 = q12 / (sigma * sqrt(real(n, dp)))

        allocate(work(k, n - k + 1))
        work = 0.0_dp
        allocate(beta_i(k), tmp(k))
        do i = k, n - 1
            call normal_equations_1d(y, x, 1, i, xtx, xty)
            call solve_linear_system_square_1d(xtx, xty, beta_i, info)
            if (info /= 0) error stop "solve_efp_re_1d: singular recursive fit"
            call symmetric_root_1d(xtx, qi12)
            qi12 = qi12 / (sigma * sqrt(real(i, dp)))
            tmp = matmul(qi12, beta_i - beta_hat)
            work(:, i - k + 1) = tmp
            deallocate(qi12)
        end do

        allocate(process(n, k))
        process = 0.0_dp
        do i = 2, n
            process(i, :) = work(:, i - 1)
        end do
        do i = 1, n
            process(i, :) = process(i, :) * real(i, dp) / sqrt(real(n, dp))
        end do

        deallocate(fitted, resid, beta_hat, beta_i, xtx, xty, q12, tmp, work)
    end subroutine solve_efp_re_1d

    subroutine solve_efp_me_1d(y, x, h, process)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: h
        real(kind=dp), allocatable, intent(out) :: process(:, :)

        integer :: n, p, k, nh, i, info
        real(kind=dp) :: sigma
        real(kind=dp), allocatable :: fitted(:), resid(:), beta_hat(:), beta_i(:), xtx(:, :), xty(:), qnh12(:, :), temp(:)

        n = size(y)
        p = size(x, 2)
        k = p + 1
        nh = floor(real(n, dp) * h)
        if (nh < k) error stop "solve_efp_me_1d: h too small"
        if (nh > n) error stop "solve_efp_me_1d: h too large"

        call ols_fit_resid_1d(y, x, 1, n, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(n - k, dp))

        allocate(xtx(k, k), xty(k), beta_hat(k))
        call normal_equations_1d(y, x, 1, n, xtx, xty)
        call solve_linear_system_square_1d(xtx, xty, beta_hat, info)
        if (info /= 0) error stop "solve_efp_me_1d: singular full fit"

        allocate(process(n - nh + 1, k), beta_i(k), temp(k))
        do i = 0, n - nh
            call normal_equations_1d(y, x, i + 1, i + nh, xtx, xty)
            call solve_linear_system_square_1d(xtx, xty, beta_i, info)
            if (info /= 0) error stop "solve_efp_me_1d: singular window fit"
            call symmetric_root_1d(xtx, qnh12)
            qnh12 = qnh12 / (sigma * sqrt(real(nh, dp)))
            temp = matmul(qnh12, beta_i - beta_hat)
            process(i + 1, :) = temp * real(nh, dp) / sqrt(real(n, dp))
            deallocate(qnh12)
        end do

        deallocate(fitted, resid, beta_hat, beta_i, xtx, xty, temp)
    end subroutine solve_efp_me_1d

    subroutine solve_mefp_ols_cusum_1d(y_hist, x_hist, y_full, x_full, alpha, process, boundary, breakpoint, critval)
        real(kind=dp), intent(in) :: y_hist(:), x_hist(:, :)
        real(kind=dp), intent(in) :: y_full(:), x_full(:, :)
        real(kind=dp), intent(in) :: alpha
        real(kind=dp), allocatable, intent(out) :: process(:), boundary(:)
        integer, intent(out) :: breakpoint
        real(kind=dp), intent(out) :: critval

        integer :: histsize, nfull, p, i
        real(kind=dp) :: sigma, xk
        real(kind=dp), allocatable :: fitted(:), resid(:), histcoef(:), full_resid(:), xtx(:, :), xty(:), cs(:)

        histsize = size(y_hist)
        nfull = size(y_full)
        p = size(x_hist, 2)
        if (size(x_full, 2) /= p) error stop "solve_mefp_ols_cusum_1d: x dimension mismatch"
        if (nfull <= histsize) error stop "solve_mefp_ols_cusum_1d: full sample must exceed history sample"

        call mefp_critval_ols_cusum_1d(alpha, critval)
        call ols_fit_resid_1d(y_hist, x_hist, 1, histsize, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(histsize - (p + 1), dp))

        allocate(xtx(p + 1, p + 1), xty(p + 1), histcoef(p + 1))
        call normal_equations_1d(y_hist, x_hist, 1, histsize, xtx, xty)
        call solve_linear_system_square_1d(xtx, xty, histcoef, i)
        if (i /= 0) error stop "solve_mefp_ols_cusum_1d: singular history fit"

        allocate(full_resid(nfull), cs(nfull))
        do i = 1, nfull
            full_resid(i) = y_full(i) - histcoef(1)
            if (p > 0) full_resid(i) = full_resid(i) - dot_product(x_full(i, :), histcoef(2:))
        end do

        allocate(process(nfull - histsize), boundary(nfull - histsize))
        cs(1) = full_resid(1)
        do i = 2, nfull
            cs(i) = cs(i - 1) + full_resid(i)
        end do
        do i = 1, nfull - histsize
            process(i) = cs(histsize + i) / (sigma * sqrt(real(histsize, dp)))
        end do

        breakpoint = -1
        do i = histsize + 1, nfull
            xk = real(i, dp) / real(histsize, dp)
            boundary(i - histsize) = sqrt(xk * (xk - 1.0_dp) * (critval ** 2 + log(xk / (xk - 1.0_dp))))
            if (breakpoint < 0) then
                if (abs(process(i - histsize)) > boundary(i - histsize)) breakpoint = i
            end if
        end do

        deallocate(fitted, resid, histcoef, full_resid, xtx, xty, cs)
    end subroutine solve_mefp_ols_cusum_1d

    subroutine solve_mefp_ols_mosum_1d(y_hist, x_hist, y_full, x_full, h, critval, process, boundary, breakpoint)
        real(kind=dp), intent(in) :: y_hist(:), x_hist(:, :)
        real(kind=dp), intent(in) :: y_full(:), x_full(:, :)
        real(kind=dp), intent(in) :: h, critval
        real(kind=dp), allocatable, intent(out) :: process(:), boundary(:)
        integer, intent(out) :: breakpoint

        integer :: histsize, nfull, p, i, k_win, hist_proc_len, full_proc_len
        real(kind=dp) :: sigma, xk
        real(kind=dp), allocatable :: fitted(:), resid(:), histcoef(:), full_resid(:), xtx(:, :), xty(:), cs(:), full_proc(:)

        histsize = size(y_hist)
        nfull = size(y_full)
        p = size(x_hist, 2)
        if (size(x_full, 2) /= p) error stop "solve_mefp_ols_mosum_1d: x dimension mismatch"
        if (nfull <= histsize) error stop "solve_mefp_ols_mosum_1d: full sample must exceed history sample"

        k_win = floor(real(histsize, dp) * h)
        if (k_win < 1) error stop "solve_mefp_ols_mosum_1d: h too small"
        if (k_win > histsize) error stop "solve_mefp_ols_mosum_1d: h too large"

        call ols_fit_resid_1d(y_hist, x_hist, 1, histsize, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(histsize - (p + 1), dp))

        allocate(xtx(p + 1, p + 1), xty(p + 1), histcoef(p + 1))
        call normal_equations_1d(y_hist, x_hist, 1, histsize, xtx, xty)
        call solve_linear_system_square_1d(xtx, xty, histcoef, i)
        if (i /= 0) error stop "solve_mefp_ols_mosum_1d: singular history fit"

        allocate(full_resid(nfull), cs(nfull + 1))
        do i = 1, nfull
            full_resid(i) = y_full(i) - histcoef(1)
            if (p > 0) full_resid(i) = full_resid(i) - dot_product(x_full(i, :), histcoef(2:))
        end do
        cs(1) = 0.0_dp
        do i = 1, nfull
            cs(i + 1) = cs(i) + full_resid(i)
        end do

        full_proc_len = nfull - k_win + 1
        hist_proc_len = histsize - k_win + 1
        allocate(full_proc(full_proc_len))
        do i = 1, full_proc_len
            full_proc(i) = (cs(i + k_win) - cs(i)) / (sigma * sqrt(real(histsize, dp)))
        end do

        allocate(process(full_proc_len - hist_proc_len), boundary(nfull - histsize))
        process = full_proc(hist_proc_len + 1:full_proc_len)

        breakpoint = -1
        do i = histsize + 1, nfull
            xk = real(i, dp) / real(histsize, dp)
            if (xk <= exp(1.0_dp)) then
                boundary(i - histsize) = critval * sqrt(2.0_dp)
            else
                boundary(i - histsize) = critval * sqrt(2.0_dp * log(xk))
            end if
            if (breakpoint < 0) then
                if (abs(process(i - histsize)) > boundary(i - histsize)) breakpoint = i
            end if
        end do

        deallocate(fitted, resid, histcoef, full_resid, xtx, xty, cs, full_proc)
    end subroutine solve_mefp_ols_mosum_1d

    subroutine solve_mefp_re_1d(y_hist, x_hist, y_full, x_full, alpha, process, statistic, boundary, breakpoint, critval)
        real(kind=dp), intent(in) :: y_hist(:), x_hist(:, :)
        real(kind=dp), intent(in) :: y_full(:), x_full(:, :)
        real(kind=dp), intent(in) :: alpha
        real(kind=dp), allocatable, intent(out) :: process(:, :), statistic(:), boundary(:)
        integer, intent(out) :: breakpoint
        real(kind=dp), intent(out) :: critval

        integer :: histsize, nfull, p, k, i, info
        real(kind=dp) :: sigma, xk
        real(kind=dp), allocatable :: fitted(:), resid(:), histcoef(:), xtx(:, :), xty(:), q12(:, :), beta_i(:), temp(:)

        histsize = size(y_hist)
        nfull = size(y_full)
        p = size(x_hist, 2)
        k = p + 1
        if (size(x_full, 2) /= p) error stop "solve_mefp_re_1d: x dimension mismatch"
        if (nfull <= histsize) error stop "solve_mefp_re_1d: full sample must exceed history sample"

        call mefp_critval_re_1d(alpha / real(k, dp), critval)
        call ols_fit_resid_1d(y_hist, x_hist, 1, histsize, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(histsize - k, dp))

        allocate(xtx(k, k), xty(k), histcoef(k))
        call normal_equations_1d(y_hist, x_hist, 1, histsize, xtx, xty)
        call solve_linear_system_square_1d(xtx, xty, histcoef, info)
        if (info /= 0) error stop "solve_mefp_re_1d: singular history fit"
        call symmetric_root_1d(xtx, q12)
        q12 = q12 / (sigma * sqrt(real(histsize, dp)))

        allocate(process(nfull - histsize, k), statistic(nfull - histsize), boundary(nfull - histsize), beta_i(k), temp(k))
        breakpoint = -1
        do i = histsize + 1, nfull
            call normal_equations_1d(y_full, x_full, 1, i, xtx, xty)
            call solve_linear_system_square_1d(xtx, xty, beta_i, info)
            if (info /= 0) error stop "solve_mefp_re_1d: singular recursive fit"
            temp = matmul(q12, beta_i - histcoef)
            process(i - histsize, :) = temp * real(i, dp) / (sigma * sqrt(real(histsize, dp)))
            if (i == histsize + 1) then
                statistic(i - histsize) = maxval(abs(process(i - histsize, :)))
            else
                statistic(i - histsize) = max(statistic(i - histsize - 1), maxval(abs(process(i - histsize, :))))
            end if

            xk = real(i, dp) / real(histsize, dp)
            boundary(i - histsize) = sqrt(xk * (xk - 1.0_dp) * (critval ** 2 + log(xk / (xk - 1.0_dp))))
            if (breakpoint < 0) then
                if (statistic(i - histsize) > boundary(i - histsize)) breakpoint = i
            end if
        end do

        deallocate(fitted, resid, histcoef, xtx, xty, q12, beta_i, temp)
    end subroutine solve_mefp_re_1d

    subroutine solve_mefp_me_1d(y_hist, x_hist, y_full, x_full, h, critval, process, statistic, boundary, breakpoint)
        real(kind=dp), intent(in) :: y_hist(:), x_hist(:, :)
        real(kind=dp), intent(in) :: y_full(:), x_full(:, :)
        real(kind=dp), intent(in) :: h, critval
        real(kind=dp), allocatable, intent(out) :: process(:, :), statistic(:), boundary(:)
        integer, intent(out) :: breakpoint

        integer :: histsize, nfull, p, k, k_win, i, info, row
        real(kind=dp) :: sigma, xk
        real(kind=dp), allocatable :: fitted(:), resid(:), histcoef(:), xtx(:, :), xty(:), beta_i(:), qwin12(:, :), temp(:)

        histsize = size(y_hist)
        nfull = size(y_full)
        p = size(x_hist, 2)
        k = p + 1
        if (size(x_full, 2) /= p) error stop "solve_mefp_me_1d: x dimension mismatch"
        if (nfull <= histsize) error stop "solve_mefp_me_1d: full sample must exceed history sample"

        k_win = floor(real(histsize, dp) * h)
        if (k_win < k) error stop "solve_mefp_me_1d: h too small"
        if (k_win > histsize) error stop "solve_mefp_me_1d: h too large"

        call ols_fit_resid_1d(y_hist, x_hist, 1, histsize, fitted, resid)
        sigma = sqrt(sum(resid ** 2) / real(histsize - k, dp))

        allocate(xtx(k, k), xty(k), histcoef(k))
        call normal_equations_1d(y_hist, x_hist, 1, histsize, xtx, xty)
        call solve_linear_system_square_1d(xtx, xty, histcoef, info)
        if (info /= 0) error stop "solve_mefp_me_1d: singular history fit"

        allocate(process(nfull - histsize, k), statistic(nfull - histsize), boundary(nfull - histsize), beta_i(k), temp(k))
        breakpoint = -1
        do i = histsize + 1, nfull
            row = i - histsize
            call normal_equations_1d(y_full, x_full, i - k_win + 1, i, xtx, xty)
            call solve_linear_system_square_1d(xtx, xty, beta_i, info)
            if (info /= 0) error stop "solve_mefp_me_1d: singular moving-window fit"
            call symmetric_root_1d(xtx, qwin12)
            qwin12 = qwin12 / sqrt(real(k_win, dp))
            temp = matmul(qwin12, beta_i - histcoef)
            process(row, :) = temp * real(k_win, dp) / (sigma * sqrt(real(histsize, dp)))
            if (row == 1) then
                statistic(row) = maxval(abs(process(row, :)))
            else
                statistic(row) = max(statistic(row - 1), maxval(abs(process(row, :))))
            end if

            xk = real(i, dp) / real(histsize, dp)
            if (xk <= exp(1.0_dp)) then
                boundary(row) = critval * sqrt(2.0_dp)
            else
                boundary(row) = critval * sqrt(2.0_dp * log(xk))
            end if
            if (breakpoint < 0) then
                if (statistic(row) > boundary(row)) breakpoint = i
            end if
            deallocate(qwin12)
        end do

        deallocate(fitted, resid, histcoef, xtx, xty, beta_i, temp)
    end subroutine solve_mefp_me_1d

    real(kind=dp) function bic_value_1d(n, k, breaks, rss) result(val)
        integer, intent(in) :: n, k, breaks
        real(kind=dp), intent(in) :: rss
        val = real(n, dp) * (log(rss) + 1.0_dp - log(real(n, dp)) + log(2.0_dp * acos(-1.0_dp))) + &
              log(real(n, dp)) * real((k + 1) * (breaks + 1), dp)
    end function bic_value_1d

    real(kind=dp) function regression_segment_rss_1d(y, x, s, e) result(rss)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        integer, intent(in) :: s, e

        integer :: nseg, p, k, i, j, jj, info
        real(kind=dp), allocatable :: xtx(:, :), xty(:), beta(:)
        real(kind=dp) :: yhat

        nseg = e - s + 1
        p = size(x, 2)
        k = p + 1
        allocate(xtx(k, k), xty(k), beta(k))
        xtx = 0.0_dp
        xty = 0.0_dp

        do i = s, e
            xtx(1, 1) = xtx(1, 1) + 1.0_dp
            xty(1) = xty(1) + y(i)
            do j = 1, p
                xtx(1, j + 1) = xtx(1, j + 1) + x(i, j)
                xtx(j + 1, 1) = xtx(j + 1, 1) + x(i, j)
                xty(j + 1) = xty(j + 1) + x(i, j) * y(i)
                do jj = 1, p
                    xtx(j + 1, jj + 1) = xtx(j + 1, jj + 1) + x(i, j) * x(i, jj)
                end do
            end do
        end do

        call solve_linear_system_square_1d(xtx, xty, beta, info)
        if (info /= 0) then
            rss = huge(1.0_dp)
            deallocate(xtx, xty, beta)
            return
        end if

        rss = 0.0_dp
        do i = s, e
            yhat = beta(1)
            if (p > 0) yhat = yhat + dot_product(x(i, :), beta(2:k))
            rss = rss + (y(i) - yhat) ** 2
        end do

        deallocate(xtx, xty, beta)
    end function regression_segment_rss_1d

    subroutine ols_fit_resid_1d(y, x, s, e, fitted, resid)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        integer, intent(in) :: s, e
        real(kind=dp), allocatable, intent(out) :: fitted(:), resid(:)

        integer :: nseg, p, k, i, j, jj, info
        real(kind=dp), allocatable :: xtx(:, :), xty(:), beta(:)

        nseg = e - s + 1
        p = size(x, 2)
        k = p + 1
        allocate(xtx(k, k), xty(k), beta(k), fitted(nseg), resid(nseg))
        xtx = 0.0_dp
        xty = 0.0_dp

        do i = s, e
            xtx(1, 1) = xtx(1, 1) + 1.0_dp
            xty(1) = xty(1) + y(i)
            do j = 1, p
                xtx(1, j + 1) = xtx(1, j + 1) + x(i, j)
                xtx(j + 1, 1) = xtx(j + 1, 1) + x(i, j)
                xty(j + 1) = xty(j + 1) + x(i, j) * y(i)
                do jj = 1, p
                    xtx(j + 1, jj + 1) = xtx(j + 1, jj + 1) + x(i, j) * x(i, jj)
                end do
            end do
        end do

        call solve_linear_system_square_1d(xtx, xty, beta, info)
        if (info /= 0) error stop "ols_fit_resid_1d: singular segment"

        do i = 1, nseg
            fitted(i) = beta(1)
            if (p > 0) fitted(i) = fitted(i) + dot_product(x(s + i - 1, :), beta(2:k))
            resid(i) = y(s + i - 1) - fitted(i)
        end do

        deallocate(xtx, xty, beta)
    end subroutine ols_fit_resid_1d

    subroutine recresid_ols_1d(y, x, rval)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), allocatable, intent(out) :: rval(:)

        integer :: n, p, k, q, r, i
        real(kind=dp) :: fr, yhat
        real(kind=dp), allocatable :: X1(:, :), xr(:), xtx(:, :), xty(:), beta(:), tmpv(:), outer_tmp(:, :)

        n = size(y)
        p = size(x, 2)
        k = p + 1
        q = k
        allocate(rval(n - q), X1(k, k), xr(k), xtx(k, k), xty(k), beta(k), tmpv(k), outer_tmp(k, k))

        call normal_equations_1d(y, x, 1, q, xtx, xty)
        call solve_linear_system_square_1d(xtx, xty, beta, i)
        if (i /= 0) error stop "recresid_ols_1d: singular initial fit"
        call invert_matrix_square_1d(xtx, X1, i)
        if (i /= 0) error stop "recresid_ols_1d: singular initial crossproduct"

        call design_row_1d(x, q + 1, xr)
        fr = 1.0_dp + dot_product(xr, matmul(X1, xr))
        yhat = dot_product(xr, beta)
        rval(1) = (y(q + 1) - yhat) / sqrt(fr)

        if ((q + 1) < n) then
            do r = q + 2, n
                tmpv = matmul(X1, xr)
                outer_tmp = spread(tmpv, 2, k) * spread(tmpv, 1, k)
                X1 = X1 - outer_tmp / fr
                tmpv = matmul(X1, xr)
                beta = beta + tmpv * rval(r - q - 1) * sqrt(fr)

                call design_row_1d(x, r, xr)
                fr = 1.0_dp + dot_product(xr, matmul(X1, xr))
                yhat = dot_product(xr, beta)
                rval(r - q) = (y(r) - yhat) / sqrt(fr)
            end do
        end if

        deallocate(X1, xr, xtx, xty, beta, tmpv, outer_tmp)
    end subroutine recresid_ols_1d

    subroutine normal_equations_1d(y, x, s, e, xtx, xty)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        integer, intent(in) :: s, e
        real(kind=dp), intent(out) :: xtx(:, :), xty(:)

        integer :: p, k, i, j, jj

        p = size(x, 2)
        k = p + 1
        xtx = 0.0_dp
        xty = 0.0_dp
        do i = s, e
            xtx(1, 1) = xtx(1, 1) + 1.0_dp
            xty(1) = xty(1) + y(i)
            do j = 1, p
                xtx(1, j + 1) = xtx(1, j + 1) + x(i, j)
                xtx(j + 1, 1) = xtx(j + 1, 1) + x(i, j)
                xty(j + 1) = xty(j + 1) + x(i, j) * y(i)
                do jj = 1, p
                    xtx(j + 1, jj + 1) = xtx(j + 1, jj + 1) + x(i, j) * x(i, jj)
                end do
            end do
        end do
    end subroutine normal_equations_1d

    subroutine symmetric_root_1d(a, root)
        real(kind=dp), intent(in) :: a(:, :)
        real(kind=dp), allocatable, intent(out) :: root(:, :)

        real(kind=dp), allocatable :: evals(:), evecs(:, :), d(:, :)
        integer :: n, i

        n = size(a, 1)
        call jacobi_eigen_sym(a, evals, evecs)
        allocate(d(n, n), root(n, n))
        d = 0.0_dp
        do i = 1, n
            d(i, i) = sqrt(max(evals(i), 0.0_dp))
        end do
        root = matmul(evecs, matmul(d, transpose(evecs)))
        deallocate(evals, evecs, d)
    end subroutine symmetric_root_1d

    subroutine design_row_1d(x, idx, xr)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: idx
        real(kind=dp), intent(out) :: xr(:)

        xr(1) = 1.0_dp
        if (size(x, 2) > 0) xr(2:) = x(idx, :)
    end subroutine design_row_1d

    subroutine invert_matrix_square_1d(a, ainv, info)
        real(kind=dp), intent(in) :: a(:, :)
        real(kind=dp), intent(out) :: ainv(:, :)
        integer, intent(out) :: info

        integer :: n, j
        real(kind=dp), allocatable :: rhs(:), sol(:)

        n = size(a, 1)
        info = 0
        allocate(rhs(n), sol(n))
        do j = 1, n
            rhs = 0.0_dp
            rhs(j) = 1.0_dp
            call solve_linear_system_square_1d(a, rhs, sol, info)
            if (info /= 0) then
                ainv = 0.0_dp
                deallocate(rhs, sol)
                return
            end if
            ainv(:, j) = sol
        end do
        deallocate(rhs, sol)
    end subroutine invert_matrix_square_1d

    real(kind=dp) function sample_sd_1d(x) result(s)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: mu
        integer :: n
        n = size(x)
        mu = sum(x) / real(n, dp)
        s = sqrt(sum((x - mu) ** 2) / real(n - 1, dp))
    end function sample_sd_1d

    subroutine mefp_critval_ols_cusum_1d(alpha, critval)
        real(kind=dp), intent(in) :: alpha
        real(kind=dp), intent(out) :: critval

        integer :: iter
        real(kind=dp) :: lo, hi, mid, fmid

        lo = 0.0_dp
        hi = 10.0_dp
        do iter = 1, 100
            mid = 0.5_dp * (lo + hi)
            fmid = 2.0_dp * (normal_cdf_1d(mid) - mid * normal_pdf_1d(mid)) + alpha - 2.0_dp
            if (fmid > 0.0_dp) then
                hi = mid
            else
                lo = mid
            end if
        end do
        critval = 0.5_dp * (lo + hi)
    end subroutine mefp_critval_ols_cusum_1d

    subroutine mefp_critval_re_1d(alpha_elem, critval)
        real(kind=dp), intent(in) :: alpha_elem
        real(kind=dp), intent(out) :: critval

        integer :: iter
        real(kind=dp) :: lo, hi, mid, fmid

        lo = 0.0_dp
        hi = 10.0_dp
        do iter = 1, 100
            mid = 0.5_dp * (lo + hi)
            fmid = 2.0_dp * (normal_cdf_1d(mid) - mid * normal_pdf_1d(mid)) + alpha_elem - 2.0_dp
            if (fmid > 0.0_dp) then
                hi = mid
            else
                lo = mid
            end if
        end do
        critval = 0.5_dp * (lo + hi)
    end subroutine mefp_critval_re_1d

    real(kind=dp) function normal_pdf_1d(x) result(v)
        real(kind=dp), intent(in) :: x
        v = exp(-0.5_dp * x * x) / sqrt(2.0_dp * acos(-1.0_dp))
    end function normal_pdf_1d

    real(kind=dp) function normal_cdf_1d(x) result(v)
        real(kind=dp), intent(in) :: x
        v = 0.5_dp * (1.0_dp + erf(x / sqrt(2.0_dp)))
    end function normal_cdf_1d

    subroutine solve_linear_system_square_1d(a, b, x, info)
        real(kind=dp), intent(in) :: a(:, :), b(:)
        real(kind=dp), intent(out) :: x(:)
        integer, intent(out) :: info

        integer :: n, i, j, k, pivot_row
        real(kind=dp) :: pivot_abs, factor
        real(kind=dp), allocatable :: aug(:, :), tmp_row(:)

        n = size(b)
        info = 0
        allocate(aug(n, n + 1), tmp_row(n + 1))
        aug(:, 1:n) = a
        aug(:, n + 1) = b

        do k = 1, n
            pivot_row = k
            pivot_abs = abs(aug(k, k))
            do i = k + 1, n
                if (abs(aug(i, k)) > pivot_abs) then
                    pivot_abs = abs(aug(i, k))
                    pivot_row = i
                end if
            end do
            if (pivot_abs <= 1.0e-12_dp) then
                info = 1
                x = 0.0_dp
                deallocate(aug, tmp_row)
                return
            end if
            if (pivot_row /= k) then
                tmp_row = aug(k, :)
                aug(k, :) = aug(pivot_row, :)
                aug(pivot_row, :) = tmp_row
            end if
            aug(k, k:n + 1) = aug(k, k:n + 1) / aug(k, k)
            do i = k + 1, n
                factor = aug(i, k)
                aug(i, k:n + 1) = aug(i, k:n + 1) - factor * aug(k, k:n + 1)
            end do
        end do

        x = 0.0_dp
        do i = n, 1, -1
            x(i) = aug(i, n + 1)
            do j = i + 1, n
                x(i) = x(i) - aug(i, j) * x(j)
            end do
        end do

        deallocate(aug, tmp_row)
    end subroutine solve_linear_system_square_1d

    integer function ceiling_div_int_1d(a, b) result(v)
        integer, intent(in) :: a, b
        v = (a + b - 1) / b
    end function ceiling_div_int_1d

end module strucchange_pkg_mod
