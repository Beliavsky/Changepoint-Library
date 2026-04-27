module changepoint_mod
    use kind_mod, only: dp, long_int
    use util_mod, only: sort_int
    use median_mod, only: median
    use arma_mod, only: ar_fit, ar_pred
    implicit none
    private
    public :: log_likelihood_corr, cost_matrix, cost_matrix_slow, mean_shift_cost_matrix, &
              rbf_kernel_cost_matrix_1d, solve_rbf_kernel_changepoints_1d, solve_pelt_mean_shift_1d, &
              solve_pelt_mean_shift_mv, solve_pelt_rbf_mv, solve_window_mean_shift_1d, &
              solve_bottomup_mean_shift_1d, refine_bkps_mean_shift_1d, solve_seeded_binseg_cusum_1d, solve_crops_mean_shift_1d, &
              solve_capa_l2_1d, solve_mvcapa_l2_mv, &
              solve_amoc_mean_norm_1d, solve_amoc_mean_cusum_1d, solve_amoc_meanvar_poisson_1d, solve_amoc_meanvar_exp_1d, solve_amoc_meanvar_gamma_1d, solve_amoc_reg_norm_1d, solve_amoc_var_css_1d, &
              solve_dynp_cost_1d, solve_bocpd_normal_gamma_1d, solve_bocd_student_t_1d, solve_bocpd_poisson_gamma_1d, &
              solve_bocpd_beta_bernoulli_1d, solve_bocpd_ar_1d, solve_argpcpd_gp_1d, solve_argpcpd_gp_rust_1d, &
              map_changepoints_bocpd, reduce_changepoints_normal_1d, &
              multivar_cost_matrix, solve_changepoints, segment_ends
    public :: solve_continuous_pwl_sse, fit_continuous_pwl_given_cps

    type :: quad_state
        integer :: nquad = 0
        real(kind=dp), allocatable :: a(:), b(:), c(:)
        integer, allocatable :: parent_t(:), parent_q(:)
    end type quad_state

contains

    !> Computes the negative log-likelihood of bivariate normal data segment.
    pure function log_likelihood_corr(x, y, min_len) result(ll)
        real(kind=dp), intent(in) :: x(:), y(:)
        integer, intent(in), optional :: min_len
        real(kind=dp) :: ll, r, ma, mb, va, vb
        integer :: n, min_len_

        min_len_ = 50
        if (present(min_len)) min_len_ = min_len

        n = size(x)
        if (n < min_len_) then
            ll = -1.0e20_dp
            return
        end if
        ma = sum(x)/n
        mb = sum(y)/n
        va = sum((x-ma)**2) / n
        vb = sum((y-mb)**2) / n
        r = sum((x-ma)*(y-mb)) / (sqrt(max(1e-20_dp, va * vb)) * n)

        ll = -0.5_dp * n * log(max(1.0e-10_dp, 1.0_dp - min(r**2, 1.0_dp - 1.0e-10_dp)))
    end function log_likelihood_corr

    !> Constructs a cost matrix for all possible segments using negative log-likelihood (O(N^2)).
    function cost_matrix(x, y, min_seg_len) result(cost)
        real(kind=dp), intent(in) :: x(:), y(:)
        integer, intent(in), optional :: min_seg_len
        real(kind=dp) :: cost(size(x), size(x))
        real(kind=dp) :: sx(size(x)), sy(size(x)), sxx(size(x)), syy(size(x)), sxy(size(x))
        integer :: n, i, j, min_len

        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        n = size(x)

        sx(1) = x(1); sy(1) = y(1); sxx(1) = x(1)**2; syy(1) = y(1)**2; sxy(1) = x(1)*y(1)
        do i = 2, n
            sx(i) = sx(i-1) + x(i)
            sy(i) = sy(i-1) + y(i)
            sxx(i) = sxx(i-1) + x(i)**2
            syy(i) = syy(i-1) + y(i)**2
            sxy(i) = sxy(i-1) + x(i)*y(i)
        end do

        do i = 1, n
            do j = i, n
                if (j - i + 1 < min_len) then
                    cost(i, j) = 1.0e20_dp
                else
                    call fast_ll(i, j, sx, sy, sxx, syy, sxy, cost(i, j))
                end if
            end do
        end do
        do i = 2, n
            do j = 1, i - 1
                cost(i, j) = 1.0e20_dp
            end do
        end do
    end function cost_matrix

    subroutine fast_ll(i, j, sx, sy, sxx, syy, sxy, cost_val)
        integer, intent(in) :: i, j
        real(kind=dp), intent(in) :: sx(:), sy(:), sxx(:), syy(:), sxy(:)
        real(kind=dp), intent(out) :: cost_val
        real(kind=dp) :: n, r, ma, mb, va, vb, sx_s, sy_s, sxx_s, syy_s, sxy_s

        n = real(j - i + 1, dp)
        if (i == 1) then
            sx_s = sx(j); sy_s = sy(j); sxx_s = sxx(j); syy_s = syy(j); sxy_s = sxy(j)
        else
            sx_s = sx(j) - sx(i-1)
            sy_s = sy(j) - sy(i-1)
            sxx_s = sxx(j) - sxx(i-1)
            syy_s = syy(j) - syy(i-1)
            sxy_s = sxy(j) - sxy(i-1)
        end if
        ma = sx_s / n
        mb = sy_s / n
        va = (sxx_s - n * ma**2) / n
        vb = (syy_s - n * mb**2) / n
        r = (sxy_s - n * ma * mb) / (sqrt(max(1e-20_dp, va * vb)) * n)

        cost_val = 0.5_dp * n * log(max(1.0e-10_dp, 1.0_dp - min(r**2, 1.0_dp - 1.0e-10_dp)))
    end subroutine fast_ll

    !> Cost matrix for a 1-D series z using the profile normal mean-shift log-likelihood.
    !! cost(i,j) = (j-i+1)/2 * log(max(sample_variance(z(i:j)), eps))
    !! Use z = x*y for covariance changepoints; z = x*x for variance changepoints.
    function mean_shift_cost_matrix(z, min_seg_len) result(cost)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in), optional :: min_seg_len
        real(kind=dp) :: cost(size(z), size(z))
        real(kind=dp) :: sz(size(z)), szz(size(z))
        real(kind=dp) :: m, sz_s, szz_s, s2
        integer :: n, i, j, min_len

        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        n = size(z)

        sz(1) = z(1);  szz(1) = z(1)**2
        do i = 2, n
            sz(i)  = sz(i-1)  + z(i)
            szz(i) = szz(i-1) + z(i)**2
        end do

        do i = 1, n
            do j = i, n
                if (j - i + 1 < min_len) then
                    cost(i, j) = 1.0e20_dp
                else
                    m     = real(j - i + 1, dp)
                    sz_s  = sz(j)  - merge(sz(i-1),  0.0_dp, i > 1)
                    szz_s = szz(j) - merge(szz(i-1), 0.0_dp, i > 1)
                    s2    = szz_s / m - (sz_s / m)**2
                    cost(i, j) = 0.5_dp * m * log(max(s2, 1.0e-20_dp))
                end if
            end do
            do j = 1, i - 1
                cost(i, j) = 1.0e20_dp
            end do
        end do
    end function mean_shift_cost_matrix

    subroutine solve_pelt_mean_shift_1d(z, pen, bkps, n_bkps_found, min_seg_len, jump)
        real(kind=dp), intent(in) :: z(:)
        real(kind=dp), intent(in) :: pen
        integer, allocatable, intent(out) :: bkps(:)
        integer, intent(out) :: n_bkps_found
        integer, intent(in), optional :: min_seg_len, jump

        real(kind=dp), allocatable :: sz(:), szz(:), best_cost(:), cand_cost(:)
        integer, allocatable :: parent(:), admissible(:), next_adm(:), tmp_bkps(:), ind(:)
        integer :: n, min_len, jump_, n_ind, i, b, t, new_adm_pt, n_adm, n_next, best_t
        real(kind=dp) :: best_val, seg_cost, cand_val
        real(kind=dp), parameter :: huge_cost = 1.0e20_dp

        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        jump_ = 5
        if (present(jump)) jump_ = jump

        if (n <= 0) then
            allocate(bkps(0))
            n_bkps_found = 0
            return
        end if

        allocate(sz(0:n), szz(0:n), best_cost(0:n), parent(0:n))
        call build_mean_shift_prefix_1d(z, sz, szz)

        best_cost = huge_cost
        parent = -1
        best_cost(0) = 0.0_dp
        parent(0) = 0

        n_ind = count([(i, i = 0, n - 1, jump_) ] >= min_len) + 1
        allocate(ind(n_ind))
        n_ind = 0
        do i = 0, n - 1, jump_
            if (i >= min_len) then
                n_ind = n_ind + 1
                ind(n_ind) = i
            end if
        end do
        n_ind = n_ind + 1
        ind(n_ind) = n

        allocate(admissible(n + 1), next_adm(n + 1), cand_cost(n + 1))
        n_adm = 0

        do i = 1, n_ind
            b = ind(i)
            new_adm_pt = ((b - min_len) / jump_) * jump_
            if (new_adm_pt >= 0) then
                if (n_adm == 0 .or. admissible(n_adm) /= new_adm_pt) then
                    n_adm = n_adm + 1
                    admissible(n_adm) = new_adm_pt
                end if
            end if

            best_val = huge_cost
            best_t = -1
            cand_cost(1:n_adm) = huge_cost

            do t = 1, n_adm
                if (best_cost(admissible(t)) >= huge_cost) cycle
                seg_cost = mean_shift_segment_sse_1d(sz, szz, admissible(t) + 1, b, min_len, huge_cost)
                if (seg_cost >= huge_cost) cycle
                cand_val = best_cost(admissible(t)) + seg_cost + pen
                cand_cost(t) = cand_val
                if (cand_val < best_val) then
                    best_val = cand_val
                    best_t = admissible(t)
                end if
            end do

            best_cost(b) = best_val
            parent(b) = best_t

            n_next = 0
            do t = 1, n_adm
                if (cand_cost(t) <= best_val + pen) then
                    n_next = n_next + 1
                    next_adm(n_next) = admissible(t)
                end if
            end do
            if (n_next > 0) admissible(1:n_next) = next_adm(1:n_next)
            n_adm = n_next
        end do

        allocate(tmp_bkps(n))
        n_bkps_found = 0
        b = n
        do while (b > 0 .and. parent(b) >= 0)
            n_bkps_found = n_bkps_found + 1
            tmp_bkps(n_bkps_found) = b
            b = parent(b)
        end do

        if (n_bkps_found > 0) then
            allocate(bkps(n_bkps_found))
            do i = 1, n_bkps_found
                bkps(i) = tmp_bkps(n_bkps_found - i + 1)
            end do
            n_bkps_found = n_bkps_found - 1
        else
            allocate(bkps(0))
        end if

        deallocate(sz, szz, best_cost, parent, admissible, next_adm, cand_cost, tmp_bkps, ind)
    end subroutine solve_pelt_mean_shift_1d

    subroutine solve_pelt_mean_shift_mv(x, pen, bkps, n_bkps_found, min_seg_len, jump)
        real(kind=dp), intent(in) :: x(:,:)
        real(kind=dp), intent(in) :: pen
        integer, allocatable, intent(out) :: bkps(:)
        integer, intent(out) :: n_bkps_found
        integer, intent(in), optional :: min_seg_len, jump

        real(kind=dp), allocatable :: sx(:,:), sxx(:,:), best_cost(:), cand_cost(:)
        integer, allocatable :: parent(:), admissible(:), next_adm(:), tmp_bkps(:), ind(:)
        integer :: n, p, min_len, jump_, n_ind, i, b, t, new_adm_pt, n_adm, n_next, best_t
        real(kind=dp) :: best_val, seg_cost, cand_val
        real(kind=dp), parameter :: huge_cost = 1.0e20_dp

        n = size(x, 1)
        p = size(x, 2)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        jump_ = 5
        if (present(jump)) jump_ = jump

        if (n <= 0 .or. p <= 0) then
            allocate(bkps(0))
            n_bkps_found = 0
            return
        end if

        allocate(sx(0:n, p), sxx(0:n, p), best_cost(0:n), parent(0:n))
        call build_mean_shift_prefix_mv(x, sx, sxx)

        best_cost = huge_cost
        parent = -1
        best_cost(0) = 0.0_dp
        parent(0) = 0

        n_ind = count([(i, i = 0, n - 1, jump_)] >= min_len) + 1
        allocate(ind(n_ind))
        n_ind = 0
        do i = 0, n - 1, jump_
            if (i >= min_len) then
                n_ind = n_ind + 1
                ind(n_ind) = i
            end if
        end do
        n_ind = n_ind + 1
        ind(n_ind) = n

        allocate(admissible(n + 1), next_adm(n + 1), cand_cost(n + 1))
        n_adm = 0

        do i = 1, n_ind
            b = ind(i)
            new_adm_pt = ((b - min_len) / jump_) * jump_
            if (new_adm_pt >= 0) then
                if (n_adm == 0 .or. admissible(n_adm) /= new_adm_pt) then
                    n_adm = n_adm + 1
                    admissible(n_adm) = new_adm_pt
                end if
            end if

            best_val = huge_cost
            best_t = -1
            cand_cost(1:n_adm) = huge_cost

            do t = 1, n_adm
                if (best_cost(admissible(t)) >= huge_cost) cycle
                seg_cost = mean_shift_segment_sse_mv(sx, sxx, admissible(t) + 1, b, min_len, huge_cost)
                if (seg_cost >= huge_cost) cycle
                cand_val = best_cost(admissible(t)) + seg_cost + pen
                cand_cost(t) = cand_val
                if (cand_val < best_val) then
                    best_val = cand_val
                    best_t = admissible(t)
                end if
            end do

            best_cost(b) = best_val
            parent(b) = best_t

            n_next = 0
            do t = 1, n_adm
                if (cand_cost(t) <= best_val + pen) then
                    n_next = n_next + 1
                    next_adm(n_next) = admissible(t)
                end if
            end do
            if (n_next > 0) admissible(1:n_next) = next_adm(1:n_next)
            n_adm = n_next
        end do

        allocate(tmp_bkps(n))
        n_bkps_found = 0
        b = n
        do while (b > 0 .and. parent(b) >= 0)
            n_bkps_found = n_bkps_found + 1
            tmp_bkps(n_bkps_found) = b
            b = parent(b)
        end do

        if (n_bkps_found > 0) then
            allocate(bkps(n_bkps_found))
            do i = 1, n_bkps_found
                bkps(i) = tmp_bkps(n_bkps_found - i + 1)
            end do
            n_bkps_found = n_bkps_found - 1
        else
            allocate(bkps(0))
        end if

        deallocate(sx, sxx, best_cost, parent, admissible, next_adm, cand_cost, tmp_bkps, ind)
    end subroutine solve_pelt_mean_shift_mv

    subroutine solve_pelt_rbf_mv(x, pen, bkps, n_bkps_found, min_seg_len, jump)
        real(kind=dp), intent(in) :: x(:,:)
        real(kind=dp), intent(in) :: pen
        integer, allocatable, intent(out) :: bkps(:)
        integer, intent(out) :: n_bkps_found
        integer, intent(in), optional :: min_seg_len, jump

        real(kind=dp), allocatable :: pref(:,:), dists(:), best_cost(:), cand_cost(:)
        integer, allocatable :: parent(:), admissible(:), next_adm(:), tmp_bkps(:), ind(:)
        integer :: n, p, min_len, jump_, n_ind, i, b, t, new_adm_pt, n_adm, n_next, best_t
        integer :: ndist
        real(kind=dp) :: best_val, seg_cost, cand_val, gamma
        real(kind=dp), parameter :: huge_cost = 1.0e20_dp

        n = size(x, 1)
        p = size(x, 2)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        jump_ = 5
        if (present(jump)) jump_ = jump

        if (n <= 0 .or. p <= 0) then
            allocate(bkps(0))
            n_bkps_found = 0
            return
        end if

        allocate(best_cost(0:n), parent(0:n))
        best_cost = huge_cost
        parent = -1
        best_cost(0) = 0.0_dp
        parent(0) = 0

        ndist = n * (n - 1) / 2
        gamma = 1.0_dp
        if (ndist > 0) then
            allocate(dists(ndist))
            call fill_pairwise_sqdist_mv(x, dists)
            gamma = gamma_from_pairwise_median(dists)
            deallocate(dists)
        end if

        allocate(pref(0:n, 0:n))
        call build_rbf_prefix_sum_mv(x, gamma, pref)

        n_ind = count([(i, i = 0, n - 1, jump_)] >= min_len) + 1
        allocate(ind(n_ind))
        n_ind = 0
        do i = 0, n - 1, jump_
            if (i >= min_len) then
                n_ind = n_ind + 1
                ind(n_ind) = i
            end if
        end do
        n_ind = n_ind + 1
        ind(n_ind) = n

        allocate(admissible(n + 1), next_adm(n + 1), cand_cost(n + 1))
        n_adm = 0

        do i = 1, n_ind
            b = ind(i)
            new_adm_pt = ((b - min_len) / jump_) * jump_
            if (new_adm_pt >= 0) then
                if (n_adm == 0 .or. admissible(n_adm) /= new_adm_pt) then
                    n_adm = n_adm + 1
                    admissible(n_adm) = new_adm_pt
                end if
            end if

            best_val = huge_cost
            best_t = -1
            cand_cost(1:n_adm) = huge_cost

            do t = 1, n_adm
                if (best_cost(admissible(t)) >= huge_cost) cycle
                seg_cost = rbf_kernel_segment_cost(pref, admissible(t) + 1, b, min_len, huge_cost)
                if (seg_cost >= huge_cost) cycle
                cand_val = best_cost(admissible(t)) + seg_cost + pen
                cand_cost(t) = cand_val
                if (cand_val < best_val) then
                    best_val = cand_val
                    best_t = admissible(t)
                end if
            end do

            best_cost(b) = best_val
            parent(b) = best_t

            n_next = 0
            do t = 1, n_adm
                if (cand_cost(t) <= best_val + pen) then
                    n_next = n_next + 1
                    next_adm(n_next) = admissible(t)
                end if
            end do
            if (n_next > 0) admissible(1:n_next) = next_adm(1:n_next)
            n_adm = n_next
        end do

        allocate(tmp_bkps(n))
        n_bkps_found = 0
        b = n
        do while (b > 0 .and. parent(b) >= 0)
            n_bkps_found = n_bkps_found + 1
            tmp_bkps(n_bkps_found) = b
            b = parent(b)
        end do

        if (n_bkps_found > 0) then
            allocate(bkps(n_bkps_found))
            do i = 1, n_bkps_found
                bkps(i) = tmp_bkps(n_bkps_found - i + 1)
            end do
            n_bkps_found = n_bkps_found - 1
        else
            allocate(bkps(0))
        end if

        deallocate(pref, best_cost, parent, admissible, next_adm, cand_cost, tmp_bkps, ind)
    end subroutine solve_pelt_rbf_mv

    subroutine solve_window_mean_shift_1d(z, width, n_bkps, bkps)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: width, n_bkps
        integer, allocatable, intent(out) :: bkps(:)

        real(kind=dp), allocatable :: sz(:), szz(:), score(:), peak_gain(:)
        integer, allocatable :: inds(:), peak_bkps(:)
        integer :: n, width_, jump, min_size, order
        integer :: i, k, start0, end0, n_inds, n_peaks, p
        real(kind=dp) :: gain, sse_full, sse_left, sse_right

        n = size(z)
        width_ = 2 * (width / 2)
        jump = 5
        min_size = 2

        if (n <= 0) then
            allocate(bkps(0))
            return
        end if

        call build_mean_shift_prefix_1d(z, sz, szz)

        allocate(inds((n + jump - 1) / jump))
        n_inds = 0
        do k = 0, n - 1, jump
            if (k >= width_ / 2 .and. k < n - width_ / 2) then
                n_inds = n_inds + 1
                inds(n_inds) = k
            end if
        end do

        if (n_inds == 0) then
            allocate(bkps(1))
            bkps(1) = n
            deallocate(sz, szz, inds)
            return
        end if

        allocate(score(n_inds))
        do i = 1, n_inds
            k = inds(i)
            start0 = k - width_ / 2
            end0 = k + width_ / 2
            sse_full = mean_shift_segment_sse_half_open_1d(sz, szz, start0, end0)
            sse_left = mean_shift_segment_sse_half_open_1d(sz, szz, start0, k)
            sse_right = mean_shift_segment_sse_half_open_1d(sz, szz, k, end0)
            gain = sse_full - sse_left - sse_right
            score(i) = gain
        end do

        order = max(max(width_, 2 * min_size) / (2 * jump), 1)
        call find_window_peaks(score, inds(1:n_inds), order, peak_bkps, peak_gain)
        n_peaks = size(peak_bkps)

        if (n_peaks == 0) then
            allocate(bkps(1))
            bkps(1) = n
            deallocate(sz, szz, inds, score)
            return
        end if

        call sort_peaks_by_gain(peak_bkps, peak_gain)

        allocate(bkps(min(n_bkps, n_peaks) + 1))
        do p = 1, min(n_bkps, n_peaks)
            bkps(p) = peak_bkps(n_peaks - p + 1)
        end do
        bkps(min(n_bkps, n_peaks) + 1) = n
        if (size(bkps) > 1) call sort_int(bkps(1:size(bkps) - 1))

        deallocate(sz, szz, inds, score, peak_bkps, peak_gain)
    end subroutine solve_window_mean_shift_1d

    subroutine solve_bottomup_mean_shift_1d(z, n_bkps, bkps, min_seg_len, jump)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: n_bkps
        integer, allocatable, intent(out) :: bkps(:)
        integer, intent(in), optional :: min_seg_len, jump

        real(kind=dp), allocatable :: sz(:), szz(:), seg_cost(:)
        integer, allocatable :: seg_start(:), seg_end(:)
        integer :: n, min_len, jump_, nseg, i, best_i
        real(kind=dp) :: best_gain, gain, merged_cost

        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        jump_ = 5
        if (present(jump)) jump_ = jump

        if (n <= 0) then
            allocate(bkps(0))
            return
        end if

        call build_mean_shift_prefix_1d(z, sz, szz)
        allocate(seg_start(n), seg_end(n), seg_cost(n))
        nseg = 0
        call grow_bottomup_leaves(0, n, min_len, jump_, seg_start, seg_end, nseg)

        do i = 1, nseg
            seg_cost(i) = mean_shift_segment_sse_half_open_1d(sz, szz, seg_start(i), seg_end(i))
        end do

        do while (nseg > n_bkps + 1)
            best_gain = huge(1.0_dp)
            best_i = -1
            merged_cost = 0.0_dp
            do i = 1, nseg - 1
                gain = mean_shift_segment_sse_half_open_1d(sz, szz, seg_start(i), seg_end(i + 1)) - &
                       seg_cost(i) - seg_cost(i + 1)
                if (gain < best_gain) then
                    best_gain = gain
                    best_i = i
                    merged_cost = mean_shift_segment_sse_half_open_1d(sz, szz, seg_start(i), seg_end(i + 1))
                end if
            end do
            if (best_i < 0) exit

            seg_end(best_i) = seg_end(best_i + 1)
            seg_cost(best_i) = merged_cost
            if (best_i + 1 <= nseg - 1) then
                seg_start(best_i + 1:nseg - 1) = seg_start(best_i + 2:nseg)
                seg_end(best_i + 1:nseg - 1) = seg_end(best_i + 2:nseg)
                seg_cost(best_i + 1:nseg - 1) = seg_cost(best_i + 2:nseg)
            end if
            nseg = nseg - 1
        end do

        allocate(bkps(nseg))
        bkps = seg_end(1:nseg)

        deallocate(sz, szz, seg_start, seg_end, seg_cost)
    end subroutine solve_bottomup_mean_shift_1d

    subroutine solve_seeded_binseg_cusum_1d(z, penalty, bkps, max_interval_length, growth_factor, selection_method)
        !> Seeded binary segmentation with univariate CUSUM scores and constant penalty.
        real(kind=dp), intent(in) :: z(:), penalty
        integer, allocatable, intent(out) :: bkps(:)
        integer, intent(in), optional :: max_interval_length
        real(kind=dp), intent(in), optional :: growth_factor
        character(len=*), intent(in), optional :: selection_method

        integer :: n, max_len, nintv, i, j, start0, end0, split0, best_split, ncp
        real(kind=dp) :: growth, best_score, score
        character(len=16) :: selection
        real(kind=dp), allocatable :: sz(:), szz(:), max_scores(:)
        integer, allocatable :: starts(:), ends(:), argmax_scores(:), tmp(:)
        logical, allocatable :: active(:)

        n = size(z)
        if (n < 2) then
            allocate(bkps(0))
            return
        end if

        max_len = 200
        if (present(max_interval_length)) max_len = max_interval_length
        growth = 1.5_dp
        if (present(growth_factor)) growth = growth_factor
        selection = "greedy"
        if (present(selection_method)) selection = trim(selection_method)

        call build_seeded_intervals(n, 2, max_len, growth, starts, ends)
        nintv = size(starts)
        allocate(sz(0:n), szz(0:n), max_scores(nintv), argmax_scores(nintv))
        call build_mean_shift_prefix_1d(z, sz, szz)

        do i = 1, nintv
            start0 = starts(i)
            end0 = ends(i)
            best_score = -huge(1.0_dp)
            best_split = start0 + 1
            do split0 = start0 + 1, end0 - 1
                score = abs(cusum_weighted_difference_1d(sz, start0, split0, end0)) - penalty
                if (score > best_score) then
                    best_score = score
                    best_split = split0
                end if
            end do
            max_scores(i) = best_score
            argmax_scores(i) = best_split
        end do

        allocate(tmp(nintv), active(nintv))
        active = (max_scores > 0.0_dp)
        ncp = 0
        do while (any(active))
            if (selection == "narrowest") then
                j = narrowest_seeded_interval(starts, ends, active)
            else
                j = maxloc(merge(max_scores, -huge(1.0_dp), active), dim=1)
            end if
            if (.not. active(j)) exit
            ncp = ncp + 1
            tmp(ncp) = argmax_scores(j)
            active = active .and. .not. ((tmp(ncp) >= starts) .and. (tmp(ncp) < ends))
        end do

        allocate(bkps(ncp))
        if (ncp > 0) then
            bkps = max(tmp(1:ncp) - 1, 0)
            if (ncp > 1) call sort_int(bkps)
        end if

        deallocate(sz, szz, max_scores, argmax_scores, starts, ends, tmp, active)
    end subroutine solve_seeded_binseg_cusum_1d

    subroutine solve_crops_mean_shift_1d(z, min_penalty, max_penalty, bkps, best_penalty, n_path_solutions, min_seg_len, middle_penalty_nudge)
        !> CROPS path search for univariate mean-shift PELT with BIC-based selection.
        real(kind=dp), intent(in) :: z(:), min_penalty, max_penalty
        integer, allocatable, intent(out) :: bkps(:)
        real(kind=dp), intent(out) :: best_penalty
        integer, intent(out) :: n_path_solutions
        integer, intent(in), optional :: min_seg_len
        real(kind=dp), intent(in), optional :: middle_penalty_nudge

        integer :: n, min_len, store_cap, queue_cap, n_store, n_queue, idx_low, idx_high, idx_mid
        integer :: num_low, num_high, num_mid, i, best_idx, j, best_ncp, n_unique
        integer, allocatable :: cps_store(:,:), numcps(:), queue_low(:), queue_high(:), unique_idx(:), order_idx(:)
        real(kind=dp), allocatable :: penalties(:), seg_costs(:), sz(:), szz(:), bic_values(:), effective_penalties(:)
        real(kind=dp) :: nudge, threshold_penalty, middle_penalty, bic
        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        nudge = 1.0e-5_dp
        if (present(middle_penalty_nudge)) nudge = middle_penalty_nudge

        if (n <= 0) then
            allocate(bkps(0))
            best_penalty = min_penalty
            n_path_solutions = 0
            return
        end if

        allocate(sz(0:n), szz(0:n))
        call build_mean_shift_prefix_1d(z, sz, szz)

        store_cap = 32
        queue_cap = 64
        allocate(cps_store(n, store_cap), numcps(store_cap), penalties(store_cap), seg_costs(store_cap))
        allocate(queue_low(queue_cap), queue_high(queue_cap))
        cps_store = 0
        numcps = 0
        penalties = 0.0_dp
        seg_costs = 0.0_dp
        n_store = 0
        n_queue = 0

        call add_crops_solution(z, sz, szz, min_penalty, min_len, cps_store, numcps, penalties, seg_costs, n_store, store_cap, idx_low)
        call add_crops_solution(z, sz, szz, max_penalty, min_len, cps_store, numcps, penalties, seg_costs, n_store, store_cap, idx_high)

        num_low = numcps(idx_low)
        num_high = numcps(idx_high)
        if (num_low > num_high + 1) then
            n_queue = 1
            queue_low(1) = idx_low
            queue_high(1) = idx_high
        end if

        do while (n_queue > 0)
            idx_low = queue_low(1)
            idx_high = queue_high(1)
            if (n_queue > 1) then
                queue_low(1:n_queue - 1) = queue_low(2:n_queue)
                queue_high(1:n_queue - 1) = queue_high(2:n_queue)
            end if
            n_queue = n_queue - 1

            num_low = numcps(idx_low)
            num_high = numcps(idx_high)
            if (num_low <= num_high + 1) cycle

            threshold_penalty = (seg_costs(idx_high) - seg_costs(idx_low)) / real(num_low - num_high, dp)
            middle_penalty = min(threshold_penalty + (penalties(idx_high) - threshold_penalty) * nudge, threshold_penalty * (1.0_dp + nudge))

            call add_crops_solution(z, sz, szz, middle_penalty, min_len, cps_store, numcps, penalties, seg_costs, n_store, store_cap, idx_mid)
            num_mid = numcps(idx_mid)

            if (num_mid == num_high) then
                cycle
            else if (num_mid == num_low) then
                cycle
            else
                if (n_queue + 2 > queue_cap) call grow_int_pair_storage(queue_low, queue_high, queue_cap)
                n_queue = n_queue + 1
                queue_low(n_queue) = idx_low
                queue_high(n_queue) = idx_mid
                n_queue = n_queue + 1
                queue_low(n_queue) = idx_mid
                queue_high(n_queue) = idx_high
            end if
        end do

        allocate(bic_values(n_store), unique_idx(n_store), effective_penalties(n_store), order_idx(n_store))
        best_idx = -1
        bic_values = huge(1.0_dp)
        effective_penalties = penalties
        best_ncp = -1
        do i = 1, n_store
            order_idx(i) = i
        end do

        do i = 1, n_store - 1
            do j = i + 1, n_store
                if (numcps(order_idx(j)) > numcps(order_idx(i)) .or. &
                    (numcps(order_idx(j)) == numcps(order_idx(i)) .and. penalties(order_idx(j)) < penalties(order_idx(i)))) then
                    idx_mid = order_idx(i)
                    order_idx(i) = order_idx(j)
                    order_idx(j) = idx_mid
                end if
            end do
        end do

        n_unique = 0
        do i = 1, n_store
            if (i == 1 .or. numcps(order_idx(i)) /= numcps(order_idx(i - 1))) then
                n_unique = n_unique + 1
                unique_idx(n_unique) = order_idx(i)
            end if
        end do

        do i = 1, n_unique - 1
            if (numcps(unique_idx(i)) == numcps(unique_idx(i + 1)) + 1) then
                effective_penalties(unique_idx(i + 1)) = (seg_costs(unique_idx(i + 1)) - seg_costs(unique_idx(i))) / &
                                                         real(numcps(unique_idx(i)) - numcps(unique_idx(i + 1)), dp)
            end if
        end do

        do i = 1, n_unique
            bic = seg_costs(unique_idx(i)) + real(numcps(unique_idx(i)) + 1, dp) * log(real(n, dp))
            bic_values(unique_idx(i)) = bic
            if (best_idx < 0 .or. bic < bic_values(best_idx)) then
                best_idx = unique_idx(i)
                best_ncp = numcps(unique_idx(i))
            end if
        end do

        allocate(bkps(numcps(best_idx)))
        if (numcps(best_idx) > 0) bkps = cps_store(1:numcps(best_idx), best_idx)
        best_penalty = effective_penalties(best_idx)
        n_path_solutions = n_unique

        deallocate(sz, szz, cps_store, numcps, penalties, seg_costs, queue_low, queue_high, bic_values, unique_idx, effective_penalties, order_idx)
    end subroutine solve_crops_mean_shift_1d

    subroutine solve_capa_l2_1d(z, segment_penalty, point_penalty, segment_starts, segment_ends, point_locs, min_seg_len, max_seg_len)
        !> CAPA with univariate L2 saving around a zero baseline mean.
        real(kind=dp), intent(in) :: z(:), segment_penalty, point_penalty
        integer, allocatable, intent(out) :: segment_starts(:), segment_ends(:), point_locs(:)
        integer, intent(in), optional :: min_seg_len, max_seg_len

        real(kind=dp), allocatable :: sz(:), opt_savings(:), candidate_savings(:)
        integer, allocatable :: opt_anomaly_starts(:), starts(:), keep_starts(:)
        integer :: n, min_len, max_len, t, nstarts, i, best_idx, best_start
        real(kind=dp) :: best_segment_saving, best_point_saving, best_total_saving
        real(kind=dp) :: max_segment_penalty

        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        max_len = min(1000, max(1, n))
        if (present(max_seg_len)) max_len = max_seg_len

        if (n <= 0) then
            allocate(segment_starts(0), segment_ends(0), point_locs(0))
            return
        end if

        allocate(sz(0:n))
        sz(0) = 0.0_dp
        do i = 1, n
            sz(i) = sz(i - 1) + z(i)
        end do
        allocate(opt_savings(0:n), opt_anomaly_starts(0:n - 1), starts(n), keep_starts(n), candidate_savings(n))
        opt_savings = 0.0_dp
        opt_anomaly_starts = -1
        nstarts = 0
        max_segment_penalty = segment_penalty

        do t = min_len - 1, n - 1
            nstarts = nstarts + 1
            starts(nstarts) = t - min_len + 1

            best_segment_saving = -huge(1.0_dp)
            best_idx = 1
            do i = 1, nstarts
                candidate_savings(i) = opt_savings(starts(i)) + l2_saving_half_open_1d(sz, starts(i), t + 1) - segment_penalty
                if (candidate_savings(i) > best_segment_saving) then
                    best_segment_saving = candidate_savings(i)
                    best_idx = i
                end if
            end do
            best_start = starts(best_idx)

            best_point_saving = opt_savings(t) + z(t + 1) ** 2 - point_penalty
            best_total_saving = opt_savings(t)
            opt_anomaly_starts(t) = -1

            if (best_segment_saving > best_total_saving) then
                best_total_saving = best_segment_saving
                opt_anomaly_starts(t) = best_start
            end if
            if (best_point_saving > best_total_saving) then
                best_total_saving = best_point_saving
                opt_anomaly_starts(t) = t
            end if
            opt_savings(t + 1) = best_total_saving

            best_idx = 0
            do i = 1, nstarts
                if (candidate_savings(i) + max_segment_penalty <= opt_savings(t + 1)) cycle
                if (starts(i) < t - max_len + 2) cycle
                best_idx = best_idx + 1
                keep_starts(best_idx) = starts(i)
            end do
            nstarts = best_idx
            if (nstarts > 0) starts(1:nstarts) = keep_starts(1:nstarts)
        end do

        call reconstruct_capa_anomalies(opt_anomaly_starts, segment_starts, segment_ends, point_locs)

        deallocate(sz, opt_savings, opt_anomaly_starts, starts, keep_starts, candidate_savings)
    end subroutine solve_capa_l2_1d

    pure real(kind=dp) function l2_saving_half_open_1d(sz, start0, end0) result(saving)
        !> Univariate L2 saving against a zero baseline mean on [start0, end0).
        real(kind=dp), intent(in) :: sz(0:)
        integer, intent(in) :: start0, end0
        real(kind=dp) :: seg_sum
        integer :: nseg

        nseg = end0 - start0
        if (nseg <= 0) then
            saving = 0.0_dp
            return
        end if

        seg_sum = sz(end0) - sz(start0)
        saving = seg_sum * seg_sum / real(nseg, dp)
    end function l2_saving_half_open_1d

    subroutine reconstruct_capa_anomalies(opt_anomaly_starts, segment_starts, segment_ends, point_locs)
        !> Recover CAPA segment and point anomalies from the stored optimal starts.
        integer, intent(in) :: opt_anomaly_starts(0:)
        integer, allocatable, intent(out) :: segment_starts(:), segment_ends(:), point_locs(:)

        integer :: i, start_i, anomaly_size, nseg, npoint, iseg, ipoint

        nseg = 0
        npoint = 0
        i = ubound(opt_anomaly_starts, 1)
        do while (i >= 0)
            start_i = opt_anomaly_starts(i)
            if (start_i >= 0) then
                anomaly_size = i - start_i + 1
                if (anomaly_size > 1) then
                    nseg = nseg + 1
                    i = start_i
                else
                    npoint = npoint + 1
                end if
            end if
            i = i - 1
        end do

        allocate(segment_starts(nseg), segment_ends(nseg), point_locs(npoint))
        iseg = nseg
        ipoint = npoint
        i = ubound(opt_anomaly_starts, 1)
        do while (i >= 0)
            start_i = opt_anomaly_starts(i)
            if (start_i >= 0) then
                anomaly_size = i - start_i + 1
                if (anomaly_size > 1) then
                    segment_starts(iseg) = start_i
                    segment_ends(iseg) = i + 1
                    iseg = iseg - 1
                    i = start_i
                else
                    point_locs(ipoint) = i
                    ipoint = ipoint - 1
                end if
            end if
            i = i - 1
        end do
    end subroutine reconstruct_capa_anomalies

    subroutine solve_mvcapa_l2_mv(x, segment_penalty, point_penalty, segment_starts, segment_ends, segment_components, segment_ncomponents, &
                                  point_locs, point_components, point_ncomponents, min_seg_len, max_seg_len)
        !> Multivariate CAPA with component-wise L2 savings against a zero baseline mean.
        real(kind=dp), intent(in) :: x(:, :)
        real(kind=dp), intent(in) :: segment_penalty(:), point_penalty(:)
        integer, allocatable, intent(out) :: segment_starts(:), segment_ends(:), segment_components(:, :), segment_ncomponents(:)
        integer, allocatable, intent(out) :: point_locs(:), point_components(:, :), point_ncomponents(:)
        integer, intent(in), optional :: min_seg_len, max_seg_len

        real(kind=dp), allocatable :: sx(:, :), opt_savings(:), candidate_savings(:), score(:)
        integer, allocatable :: starts(:), keep_starts(:), opt_anomaly_starts(:), opt_anomaly_kind(:), opt_ncomponents(:), opt_components(:, :)
        integer, allocatable :: candidate_components(:), seg_best_components(:), point_best_components(:)
        integer :: n, p, min_len, max_len, t, nstarts, i, nbest, nbest_point, seg_best_ncomp
        real(kind=dp) :: best_segment_saving, best_point_saving, best_total_saving, max_segment_penalty

        n = size(x, 1)
        p = size(x, 2)
        if (size(segment_penalty) /= p .or. size(point_penalty) /= p) error stop "penalty arrays must match the number of variables"

        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        max_len = min(1000, max(1, n))
        if (present(max_seg_len)) max_len = max_seg_len

        if (n <= 0 .or. p <= 0) then
            allocate(segment_starts(0), segment_ends(0), segment_components(0, 0), segment_ncomponents(0))
            allocate(point_locs(0), point_components(0, 0), point_ncomponents(0))
            return
        end if

        allocate(sx(0:n, p), opt_savings(0:n), candidate_savings(n), score(p))
        allocate(starts(n), keep_starts(n), opt_anomaly_starts(0:n - 1), opt_anomaly_kind(0:n - 1), opt_ncomponents(0:n - 1), opt_components(0:n - 1, p))
        allocate(candidate_components(p), seg_best_components(p), point_best_components(p))
        sx(0, :) = 0.0_dp
        do i = 1, n
            sx(i, :) = sx(i - 1, :) + x(i, :)
        end do

        opt_savings = 0.0_dp
        opt_anomaly_starts = -1
        opt_anomaly_kind = 0
        opt_ncomponents = 0
        opt_components = -1
        nstarts = 0
        max_segment_penalty = maxval(segment_penalty)

        do t = min_len - 1, n - 1
            nstarts = nstarts + 1
            starts(nstarts) = t - min_len + 1

            best_segment_saving = -huge(1.0_dp)
            nbest = 0
            seg_best_ncomp = 0
            seg_best_components = -1
            do i = 1, nstarts
                score = ((sx(t + 1, :) - sx(starts(i), :)) ** 2) / real(t + 1 - starts(i), dp)
                call best_penalised_subset(score, segment_penalty, candidate_savings(i), candidate_components, nbest)
                candidate_savings(i) = opt_savings(starts(i)) + candidate_savings(i)
                if (candidate_savings(i) > best_segment_saving) then
                    best_segment_saving = candidate_savings(i)
                    opt_anomaly_starts(t) = starts(i)
                    opt_anomaly_kind(t) = 1
                    seg_best_ncomp = nbest
                    seg_best_components = candidate_components
                end if
            end do

            score = x(t + 1, :) ** 2
            call best_penalised_subset(score, point_penalty, best_point_saving, point_best_components, nbest_point)
            best_point_saving = opt_savings(t) + best_point_saving

            best_total_saving = opt_savings(t)
            if (best_segment_saving > best_total_saving) then
                best_total_saving = best_segment_saving
                opt_ncomponents(t) = seg_best_ncomp
                opt_components(t, :) = seg_best_components
            else
                opt_anomaly_kind(t) = 0
                opt_anomaly_starts(t) = -1
                opt_ncomponents(t) = 0
                opt_components(t, :) = -1
            end if
            if (best_point_saving > best_total_saving) then
                best_total_saving = best_point_saving
                opt_anomaly_starts(t) = t
                opt_anomaly_kind(t) = 2
                opt_ncomponents(t) = nbest_point
                opt_components(t, :) = point_best_components
            end if
            opt_savings(t + 1) = best_total_saving

            nbest = 0
            do i = 1, nstarts
                if (candidate_savings(i) + max_segment_penalty <= opt_savings(t + 1)) cycle
                if (starts(i) < t - max_len + 2) cycle
                nbest = nbest + 1
                keep_starts(nbest) = starts(i)
            end do
            nstarts = nbest
            if (nstarts > 0) starts(1:nstarts) = keep_starts(1:nstarts)
        end do

        call reconstruct_mvcapa_anomalies(opt_anomaly_starts, opt_anomaly_kind, opt_ncomponents, opt_components, &
                                          segment_starts, segment_ends, segment_components, segment_ncomponents, &
                                          point_locs, point_components, point_ncomponents)

        deallocate(sx, opt_savings, candidate_savings, score, starts, keep_starts, opt_anomaly_starts, &
                   candidate_components, seg_best_components, point_best_components, &
                   opt_anomaly_kind, opt_ncomponents, opt_components)
    end subroutine solve_mvcapa_l2_mv

    subroutine best_penalised_subset(score, penalty, best_value, best_components, nbest)
        !> Best penalised subset for additive component-wise scores and a penalty array.
        real(kind=dp), intent(in) :: score(:), penalty(:)
        real(kind=dp), intent(out) :: best_value
        integer, intent(out) :: best_components(:), nbest

        real(kind=dp), allocatable :: sorted_score(:)
        integer, allocatable :: order(:)
        integer :: p, i, j, tmp_i
        real(kind=dp) :: tmp_r, csum, cand

        p = size(score)
        if (size(penalty) /= p .or. size(best_components) /= p) error stop "invalid subset dimensions"

        allocate(sorted_score(p), order(p))
        sorted_score = score
        do i = 1, p
            order(i) = i - 1
        end do

        do i = 1, p - 1
            do j = i + 1, p
                if (sorted_score(j) > sorted_score(i)) then
                    tmp_r = sorted_score(i)
                    sorted_score(i) = sorted_score(j)
                    sorted_score(j) = tmp_r
                    tmp_i = order(i)
                    order(i) = order(j)
                    order(j) = tmp_i
                end if
            end do
        end do

        best_value = -huge(1.0_dp)
        nbest = 0
        csum = 0.0_dp
        best_components = -1
        do i = 1, p
            csum = csum + sorted_score(i)
            cand = csum - penalty(i)
            if (cand > best_value) then
                best_value = cand
                nbest = i
                best_components(1:i) = order(1:i)
                if (i > 1) call sort_int(best_components(1:i))
                if (i < p) best_components(i + 1:p) = -1
            end if
        end do

        deallocate(sorted_score, order)
    end subroutine best_penalised_subset

    subroutine reconstruct_mvcapa_anomalies(opt_anomaly_starts, opt_anomaly_kind, opt_ncomponents, opt_components, &
                                            segment_starts, segment_ends, segment_components, segment_ncomponents, &
                                            point_locs, point_components, point_ncomponents)
        !> Recover multivariate CAPA anomalies and affected components from the stored dynamic program state.
        integer, intent(in) :: opt_anomaly_starts(0:), opt_anomaly_kind(0:), opt_ncomponents(0:), opt_components(0:, :)
        integer, allocatable, intent(out) :: segment_starts(:), segment_ends(:), segment_components(:, :), segment_ncomponents(:)
        integer, allocatable, intent(out) :: point_locs(:), point_components(:, :), point_ncomponents(:)

        integer :: i, start_i, anomaly_size, nseg, npoint, iseg, ipoint, p

        p = size(opt_components, 2)
        nseg = 0
        npoint = 0
        i = ubound(opt_anomaly_starts, 1)
        do while (i >= 0)
            if (opt_anomaly_kind(i) == 1) then
                nseg = nseg + 1
                i = opt_anomaly_starts(i)
            else if (opt_anomaly_kind(i) == 2) then
                npoint = npoint + 1
            end if
            i = i - 1
        end do

        allocate(segment_starts(nseg), segment_ends(nseg), segment_components(nseg, p), segment_ncomponents(nseg))
        allocate(point_locs(npoint), point_components(npoint, p), point_ncomponents(npoint))
        segment_components = -1
        point_components = -1

        iseg = nseg
        ipoint = npoint
        i = ubound(opt_anomaly_starts, 1)
        do while (i >= 0)
            if (opt_anomaly_kind(i) == 1) then
                start_i = opt_anomaly_starts(i)
                anomaly_size = i - start_i + 1
                if (anomaly_size > 1) then
                    segment_starts(iseg) = start_i
                    segment_ends(iseg) = i + 1
                    segment_ncomponents(iseg) = opt_ncomponents(i)
                    if (opt_ncomponents(i) > 0) segment_components(iseg, 1:opt_ncomponents(i)) = opt_components(i, 1:opt_ncomponents(i))
                    iseg = iseg - 1
                    i = start_i
                end if
            else if (opt_anomaly_kind(i) == 2) then
                point_locs(ipoint) = i
                point_ncomponents(ipoint) = opt_ncomponents(i)
                if (opt_ncomponents(i) > 0) point_components(ipoint, 1:opt_ncomponents(i)) = opt_components(i, 1:opt_ncomponents(i))
                ipoint = ipoint - 1
            end if
            i = i - 1
        end do
    end subroutine reconstruct_mvcapa_anomalies

    subroutine solve_dynp_cost_1d(z, max_m, model, seg_ends, min_seg_len)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: max_m
        character(len=*), intent(in) :: model
        integer, intent(out) :: seg_ends(max_m)
        integer, intent(in), optional :: min_seg_len

        real(kind=dp), allocatable :: sz(:), szz(:), dp_table(:,:)
        integer, allocatable :: parent(:,:)
        integer :: n, min_len, i, k, m
        real(kind=dp) :: seg_cost
        real(kind=dp), parameter :: huge_cost = 1.0e20_dp

        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len

        if (n <= 0) then
            seg_ends = 0
            return
        end if

        if (trim(model) == "rbf") then
            allocate(dp_table(n, max_m), parent(n, max_m))
            call solve_rbf_kernel_changepoints_1d(z, max_m, dp_table, parent, min_seg_len)
            seg_ends = segment_ends(parent, max_m)
            deallocate(dp_table, parent)
            return
        end if

        allocate(sz(0:n), szz(0:n), dp_table(n, max_m), parent(n, max_m))
        call build_mean_shift_prefix_1d(z, sz, szz)

        dp_table = huge_cost
        parent = 0
        do i = 1, n
            seg_cost = dynp_segment_cost_1d(model, z, sz, szz, 1, i, min_len, huge_cost)
            dp_table(i, 1) = seg_cost
        end do

        do m = 2, max_m
            do i = m, n
                if (i < m * min_len) cycle
                do k = m - 1, i - 1
                    if (dp_table(k, m - 1) >= huge_cost) cycle
                    seg_cost = dynp_segment_cost_1d(model, z, sz, szz, k + 1, i, min_len, huge_cost)
                    if (seg_cost >= huge_cost) cycle
                    if (dp_table(k, m - 1) + seg_cost < dp_table(i, m)) then
                        dp_table(i, m) = dp_table(k, m - 1) + seg_cost
                        parent(i, m) = k
                    end if
                end do
            end do
        end do

        seg_ends = segment_ends(parent, max_m)

        deallocate(sz, szz, dp_table, parent)
    end subroutine solve_dynp_cost_1d

    subroutine refine_bkps_mean_shift_1d(z, bkps_in, bkps_out, min_seg_len, max_iter)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: bkps_in(:)
        integer, allocatable, intent(out) :: bkps_out(:)
        integer, intent(in), optional :: min_seg_len, max_iter

        real(kind=dp), allocatable :: sz(:), szz(:)
        integer :: n, min_len, max_iter_, iter, i, left0, right0, lo, hi, k, best_k
        real(kind=dp) :: best_cost, cand_cost
        logical :: changed

        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        max_iter_ = 20
        if (present(max_iter)) max_iter_ = max_iter

        allocate(bkps_out(size(bkps_in)))
        bkps_out = bkps_in
        if (size(bkps_out) == 0) return
        if (bkps_out(size(bkps_out)) /= n) error stop "refine_bkps_mean_shift_1d: final endpoint must be n"

        call build_mean_shift_prefix_1d(z, sz, szz)

        do iter = 1, max_iter_
            changed = .false.
            do i = 1, size(bkps_out) - 1
                if (i == 1) then
                    left0 = 0
                else
                    left0 = bkps_out(i - 1)
                end if
                right0 = bkps_out(i + 1)
                lo = left0 + min_len
                hi = right0 - min_len
                if (lo > hi) cycle

                best_k = bkps_out(i)
                best_cost = mean_shift_segment_sse_half_open_1d(sz, szz, left0, best_k) + &
                            mean_shift_segment_sse_half_open_1d(sz, szz, best_k, right0)
                do k = lo, hi
                    cand_cost = mean_shift_segment_sse_half_open_1d(sz, szz, left0, k) + &
                                mean_shift_segment_sse_half_open_1d(sz, szz, k, right0)
                    if (cand_cost < best_cost) then
                        best_cost = cand_cost
                        best_k = k
                    end if
                end do

                if (best_k /= bkps_out(i)) then
                    bkps_out(i) = best_k
                    changed = .true.
                end if
            end do
            if (.not. changed) exit
        end do

        deallocate(sz, szz)
    end subroutine refine_bkps_mean_shift_1d

    recursive subroutine grow_bottomup_leaves(start0, end0, min_seg_len, jump, seg_start, seg_end, nseg)
        integer, intent(in) :: start0, end0, min_seg_len, jump
        integer, intent(inout) :: seg_start(:), seg_end(:)
        integer, intent(inout) :: nseg

        integer :: b, best_bkp
        real(kind=dp) :: mid, best_dist, dist

        mid = 0.5_dp * real(start0 + end0, dp)
        best_bkp = -1
        best_dist = huge(1.0_dp)
        do b = start0, end0
            if (mod(b, jump) /= 0) cycle
            if (b - start0 < min_seg_len) cycle
            if (end0 - b < min_seg_len) cycle
            dist = abs(real(b, dp) - mid)
            if (dist < best_dist) then
                best_dist = dist
                best_bkp = b
            end if
        end do

        if (best_bkp < 0) then
            nseg = nseg + 1
            seg_start(nseg) = start0
            seg_end(nseg) = end0
        else
            call grow_bottomup_leaves(start0, best_bkp, min_seg_len, jump, seg_start, seg_end, nseg)
            call grow_bottomup_leaves(best_bkp, end0, min_seg_len, jump, seg_start, seg_end, nseg)
        end if
    end subroutine grow_bottomup_leaves

    !> RBF-kernel segment cost matrix for a 1-D series z.
    !! Matches ruptures CostRbf with the median heuristic for gamma:
    !!   cost(i,j) = sum_{t=i}^j K_tt - (1/m) * sum_{a=i}^j sum_{b=i}^j K_ab
    !! where K_ab = exp(-clip(gamma * (z_a-z_b)^2, 1e-2, 1e2)) and K_tt = 1.
    function rbf_kernel_cost_matrix_1d(z, min_seg_len) result(cost)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in), optional :: min_seg_len
        real(kind=dp) :: cost(size(z), size(z))
        real(kind=dp), allocatable :: gram(:,:), pref(:,:), dists(:)
        real(kind=dp) :: gamma, m_r, sub_sum
        integer :: n, i, j, min_len, ndist

        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len

        cost = 1.0e20_dp
        if (n <= 0) return

        ndist = n * (n - 1) / 2
        gamma = 1.0_dp
        if (ndist > 0) then
            allocate(dists(ndist))
            call fill_pairwise_sqdist_1d(z, dists)
            gamma = gamma_from_pairwise_median(dists)
            deallocate(dists)
        end if

        allocate(gram(n, n), pref(0:n, 0:n))
        call build_rbf_gram_1d(z, gamma, gram)
        call build_matrix_prefix_sum(gram, pref)

        do i = 1, n
            do j = i + min_len - 1, n
                m_r = real(j - i + 1, dp)
                sub_sum = pref(j, j) - pref(i - 1, j) - pref(j, i - 1) + pref(i - 1, i - 1)
                cost(i, j) = m_r - sub_sum / m_r
            end do
        end do

        deallocate(gram, pref)
    end function rbf_kernel_cost_matrix_1d

    subroutine solve_rbf_kernel_changepoints_1d(z, max_m, dp_table, parent, min_seg_len)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: max_m
        real(kind=dp), intent(out) :: dp_table(size(z), max_m)
        integer, intent(out) :: parent(size(z), max_m)
        integer, intent(in), optional :: min_seg_len

        real(kind=dp), allocatable :: pref(:,:), dists(:)
        real(kind=dp) :: gamma, seg_cost
        integer :: n, min_len, ndist, i, k, m
        real(kind=dp), parameter :: huge_cost = 1.0e20_dp

        n = size(z)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len

        dp_table = huge_cost
        parent = 0
        if (n <= 0) return

        ndist = n * (n - 1) / 2
        gamma = 1.0_dp
        if (ndist > 0) then
            allocate(dists(ndist))
            call fill_pairwise_sqdist_1d(z, dists)
            gamma = gamma_from_pairwise_median(dists)
            deallocate(dists)
        end if

        allocate(pref(0:n, 0:n))
        call build_rbf_prefix_sum_1d(z, gamma, pref)

        do i = 1, n
            seg_cost = rbf_kernel_segment_cost(pref, 1, i, min_len, huge_cost)
            dp_table(i, 1) = seg_cost
        end do

        do m = 2, max_m
            do i = m, n
                if (i < m * min_len) cycle
                do k = m - 1, i - 1
                    if (dp_table(k, m - 1) >= huge_cost) cycle
                    seg_cost = rbf_kernel_segment_cost(pref, k + 1, i, min_len, huge_cost)
                    if (seg_cost >= huge_cost) cycle
                    if (dp_table(k, m - 1) + seg_cost < dp_table(i, m)) then
                        dp_table(i, m) = dp_table(k, m - 1) + seg_cost
                        parent(i, m) = k
                    end if
                end do
            end do
        end do

        deallocate(pref)
    end subroutine solve_rbf_kernel_changepoints_1d

    subroutine fill_pairwise_sqdist_1d(z, dists)
        real(kind=dp), intent(in) :: z(:)
        real(kind=dp), intent(out) :: dists(:)
        integer :: n, i, j, k
        real(kind=dp) :: dz

        n = size(z)
        k = 0
        do i = 1, n - 1
            do j = i + 1, n
                k = k + 1
                dz = z(i) - z(j)
                dists(k) = dz * dz
            end do
        end do
    end subroutine fill_pairwise_sqdist_1d

    subroutine fill_pairwise_sqdist_mv(x, dists)
        real(kind=dp), intent(in) :: x(:,:)
        real(kind=dp), intent(out) :: dists(:)
        integer :: n, i, j, k

        n = size(x, 1)
        k = 0
        do i = 1, n - 1
            do j = i + 1, n
                k = k + 1
                dists(k) = sum((x(i, :) - x(j, :)) ** 2)
            end do
        end do
    end subroutine fill_pairwise_sqdist_mv

    subroutine build_mean_shift_prefix_1d(z, sz, szz)
        real(kind=dp), intent(in) :: z(:)
        real(kind=dp), allocatable, intent(out) :: sz(:), szz(:)
        integer :: i, n

        n = size(z)
        allocate(sz(0:n), szz(0:n))
        sz(0) = 0.0_dp
        szz(0) = 0.0_dp
        do i = 1, n
            sz(i) = sz(i - 1) + z(i)
            szz(i) = szz(i - 1) + z(i) * z(i)
        end do
    end subroutine build_mean_shift_prefix_1d

    subroutine build_mean_shift_prefix_mv(x, sx, sxx)
        real(kind=dp), intent(in) :: x(:,:)
        real(kind=dp), intent(out) :: sx(0:, :), sxx(0:, :)
        integer :: i, n

        n = size(x, 1)
        sx(0, :) = 0.0_dp
        sxx(0, :) = 0.0_dp
        do i = 1, n
            sx(i, :) = sx(i - 1, :) + x(i, :)
            sxx(i, :) = sxx(i - 1, :) + x(i, :) * x(i, :)
        end do
    end subroutine build_mean_shift_prefix_mv

    pure function mean_shift_segment_sse_1d(sz, szz, s, e, min_seg_len, huge_cost) result(cost)
        real(kind=dp), intent(in) :: sz(0:), szz(0:)
        integer, intent(in) :: s, e, min_seg_len
        real(kind=dp), intent(in) :: huge_cost
        real(kind=dp) :: cost
        real(kind=dp) :: m, sum_z, sum_zz

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
        cost = sum_zz - (sum_z * sum_z) / m
    end function mean_shift_segment_sse_1d

    pure function mean_shift_segment_sse_half_open_1d(sz, szz, start0, end0) result(cost)
        real(kind=dp), intent(in) :: sz(0:), szz(0:)
        integer, intent(in) :: start0, end0
        real(kind=dp) :: cost
        real(kind=dp) :: m, sum_z, sum_zz
        integer :: s, e

        if (end0 <= start0) then
            cost = 0.0_dp
            return
        end if

        s = start0 + 1
        e = end0
        m = real(end0 - start0, dp)
        sum_z = sz(e) - sz(s - 1)
        sum_zz = szz(e) - szz(s - 1)
        cost = sum_zz - (sum_z * sum_z) / m
    end function mean_shift_segment_sse_half_open_1d

    pure function mean_shift_segment_sse_mv(sx, sxx, s, e, min_seg_len, huge_cost) result(cost)
        real(kind=dp), intent(in) :: sx(0:, :), sxx(0:, :)
        integer, intent(in) :: s, e, min_seg_len
        real(kind=dp), intent(in) :: huge_cost
        real(kind=dp) :: cost
        real(kind=dp) :: m

        if (e < s) then
            cost = huge_cost
            return
        end if

        m = real(e - s + 1, dp)
        if (m < real(min_seg_len, dp)) then
            cost = huge_cost
            return
        end if

        cost = sum((sxx(e, :) - sxx(s - 1, :)) - ((sx(e, :) - sx(s - 1, :)) ** 2) / m)
    end function mean_shift_segment_sse_mv

    pure function normal_segment_cost_1d(sz, szz, s, e, min_seg_len, huge_cost) result(cost)
        real(kind=dp), intent(in) :: sz(0:), szz(0:)
        integer, intent(in) :: s, e, min_seg_len
        real(kind=dp), intent(in) :: huge_cost
        real(kind=dp) :: cost
        real(kind=dp) :: m, sum_z, sum_zz, var_hat
        real(kind=dp), parameter :: small_diag = 1.0e-6_dp

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
        var_hat = max(sum_zz / m - (sum_z / m) ** 2, 0.0_dp)
        cost = m * log(var_hat + small_diag)
    end function normal_segment_cost_1d

    subroutine solve_amoc_mean_norm_1d(x, cpt, conf_value, null_like, alt_like, penalty, minseglen, pen_value)
        !> Single changepoint in the mean for Normal data, matching changepoint::single.mean.norm.
        real(kind=dp), intent(in) :: x(:)
        integer, intent(out) :: cpt
        real(kind=dp), intent(out) :: conf_value, null_like, alt_like
        character(len=*), intent(in), optional :: penalty
        integer, intent(in), optional :: minseglen
        real(kind=dp), intent(in), optional :: pen_value

        integer :: n, minseglen_use, tau, best_tau
        real(kind=dp) :: sum_x, sum_xx, left_sum, left_ss, right_sum, right_ss
        real(kind=dp) :: alt_raw, alt_use, best_alt_raw, penalty_use, test_stat
        real(kind=dp), allocatable :: cs(:), css(:)
        character(len=32) :: penalty_use_name

        n = size(x)
        minseglen_use = 1
        if (present(minseglen)) minseglen_use = minseglen
        if (n < 2) error stop "solve_amoc_mean_norm_1d requires at least 2 observations"
        if (n < 2 * minseglen_use) error stop "minseglen too large for solve_amoc_mean_norm_1d"

        penalty_use_name = "SIC"
        if (present(penalty)) penalty_use_name = trim(penalty)

        allocate(cs(0:n), css(0:n))
        cs(0) = 0.0_dp
        css(0) = 0.0_dp
        do tau = 1, n
            cs(tau) = cs(tau - 1) + x(tau)
            css(tau) = css(tau - 1) + x(tau) * x(tau)
        end do

        sum_x = cs(n)
        sum_xx = css(n)
        null_like = sum_xx - sum_x * sum_x / real(n, dp)

        best_tau = minseglen_use
        best_alt_raw = huge(1.0_dp)
        do tau = minseglen_use, n - minseglen_use + 1
            left_sum = cs(tau)
            left_ss = css(tau)
            if (n - tau <= 0) then
                alt_raw = huge(1.0_dp)
            else
                right_sum = sum_x - left_sum
                right_ss = sum_xx - left_ss
                alt_raw = left_ss - left_sum * left_sum / real(tau, dp) + &
                          right_ss - right_sum * right_sum / real(n - tau, dp)
            end if
            if (alt_raw < best_alt_raw) then
                best_alt_raw = alt_raw
                best_tau = tau
            end if
        end do

        alt_like = best_alt_raw
        alt_use = best_alt_raw
        if (penalty_use_name == "MBIC") then
            alt_use = alt_use + log(real(best_tau, dp)) + log(real(n - best_tau + 1, dp))
        end if

        penalty_use = amoc_penalty_value(trim(penalty_use_name), n, 1, pen_value)
        test_stat = null_like - alt_use
        if (test_stat >= penalty_use) then
            cpt = best_tau
        else
            cpt = n
        end if

        conf_value = amoc_mean_norm_conf_value(n, null_like, best_alt_raw)

        deallocate(cs, css)
    end subroutine solve_amoc_mean_norm_1d

    subroutine solve_amoc_mean_cusum_1d(x, cpt, test_stat, penalty_value, penalty, minseglen, pen_value)
        !> Single changepoint in the mean using the CUSUM statistic, matching changepoint::single.mean.cusum.
        real(kind=dp), intent(in) :: x(:)
        integer, intent(out) :: cpt
        real(kind=dp), intent(out) :: test_stat, penalty_value
        character(len=*), intent(in), optional :: penalty
        integer, intent(in), optional :: minseglen
        real(kind=dp), intent(in), optional :: pen_value

        integer :: n, minseglen_use, tau, best_tau
        real(kind=dp) :: xbar, scale_n, cusum_abs, max_cusum_abs
        real(kind=dp), allocatable :: cs(:)
        character(len=32) :: penalty_use_name

        n = size(x)
        minseglen_use = 1
        if (present(minseglen)) minseglen_use = minseglen
        if (n < 2) error stop "solve_amoc_mean_cusum_1d requires at least 2 observations"
        if (n < 2 * minseglen_use) error stop "minseglen too large for solve_amoc_mean_cusum_1d"

        penalty_use_name = "Manual"
        if (present(penalty)) penalty_use_name = trim(penalty)

        allocate(cs(0:n))
        xbar = sum(x) / real(n, dp)
        cs(0) = 0.0_dp
        do tau = 1, n
            cs(tau) = cs(tau - 1) + (x(tau) - xbar)
        end do

        scale_n = real(n, dp)
        max_cusum_abs = -1.0_dp
        best_tau = minseglen_use
        do tau = minseglen_use, n - minseglen_use + 1
            ! R stores y = c(0, cumsum(x - mean(x))) / n and searches y[minseglen:(n - minseglen + 1)].
            cusum_abs = abs(cs(tau - 1) / scale_n)
            if (cusum_abs > max_cusum_abs) then
                max_cusum_abs = cusum_abs
                best_tau = tau
            end if
        end do

        test_stat = max_cusum_abs
        penalty_value = cusum_penalty_value(trim(penalty_use_name), n, 1, pen_value)
        if (test_stat >= penalty_value) then
            cpt = best_tau
        else
            cpt = n
        end if

        deallocate(cs)
    end subroutine solve_amoc_mean_cusum_1d

    real(kind=dp) function amoc_penalty_value(penalty, n, diffparam, pen_value) result(penalty_out)
        !> Penalty calculation matching changepoint::penalty_decision for common AMOC penalties.
        character(len=*), intent(in) :: penalty
        integer, intent(in) :: n, diffparam
        real(kind=dp), intent(in), optional :: pen_value

        select case (trim(penalty))
        case ("SIC0", "BIC0")
            penalty_out = real(diffparam, dp) * log(real(n, dp))
        case ("SIC", "BIC")
            penalty_out = real(diffparam + 1, dp) * log(real(n, dp))
        case ("MBIC")
            penalty_out = real(diffparam + 2, dp) * log(real(n, dp))
        case ("AIC0")
            penalty_out = 2.0_dp * real(diffparam, dp)
        case ("AIC")
            penalty_out = 2.0_dp * real(diffparam + 1, dp)
        case ("Hannan-Quinn0")
            penalty_out = 2.0_dp * real(diffparam, dp) * log(log(real(n, dp)))
        case ("Hannan-Quinn")
            penalty_out = 2.0_dp * real(diffparam + 1, dp) * log(log(real(n, dp)))
        case ("None")
            penalty_out = 0.0_dp
        case ("Manual")
            if (.not. present(pen_value)) error stop "Manual AMOC penalty requires pen_value"
            penalty_out = pen_value
        case default
            error stop "Unsupported AMOC penalty in solve_amoc_mean_norm_1d"
        end select

        if (penalty_out < 0.0_dp) error stop "Negative AMOC penalty is invalid"
    end function amoc_penalty_value

    real(kind=dp) function cusum_penalty_value(penalty, n, diffparam, pen_value) result(penalty_out)
        !> Penalty calculation matching changepoint::penalty_decision for mean CUSUM.
        character(len=*), intent(in) :: penalty
        integer, intent(in) :: n, diffparam
        real(kind=dp), intent(in), optional :: pen_value

        select case (trim(penalty))
        case ("Asymptotic")
            error stop "Asymptotic penalties have not been implemented yet for CUSUM"
        case ("MBIC")
            error stop "MBIC penalty is not valid for nonparametric CUSUM"
        case default
            penalty_out = amoc_penalty_value(trim(penalty), n, diffparam, pen_value)
        end select
    end function cusum_penalty_value

    real(kind=dp) function amoc_mean_norm_conf_value(n, null_like, alt_like_raw) result(conf_value)
        !> Confidence value matching changepoint::single.mean.norm for class=FALSE output.
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: null_like, alt_like_raw

        real(kind=dp) :: alogn, blogn, term_a, term_b, pi_val

        pi_val = acos(-1.0_dp)
        alogn = (2.0_dp * log(log(real(n, dp)))) ** (-0.5_dp)
        blogn = 1.0_dp / alogn + 0.5_dp * alogn * log(log(log(real(n, dp))))
        term_a = exp(-alogn * sqrt(abs(null_like - alt_like_raw)) + blogn / alogn)
        term_b = exp(blogn / alogn)
        conf_value = exp(-2.0_dp * sqrt(pi_val) * term_a) - exp(-2.0_dp * sqrt(pi_val) * term_b)
    end function amoc_mean_norm_conf_value

    subroutine solve_amoc_meanvar_poisson_1d(x, cpt, penalty_value, null_like, alt_like, penalty, minseglen, pen_value)
        !> Single changepoint in Poisson mean/variance, matching changepoint::single.meanvar.poisson.
        real(kind=dp), intent(in) :: x(:)
        integer, intent(out) :: cpt
        real(kind=dp), intent(out) :: penalty_value, null_like, alt_like
        character(len=*), intent(in), optional :: penalty
        integer, intent(in), optional :: minseglen
        real(kind=dp), intent(in), optional :: pen_value

        integer :: n, minseglen_use, tau, best_tau
        real(kind=dp) :: total_sum, left_sum, right_sum, alt_use, best_alt
        real(kind=dp), allocatable :: cs(:)
        character(len=32) :: penalty_use_name

        n = size(x)
        minseglen_use = 1
        if (present(minseglen)) minseglen_use = minseglen
        if (n < 4) error stop "solve_amoc_meanvar_poisson_1d requires at least 4 observations"
        if (n < 2 * minseglen_use) error stop "minseglen too large for solve_amoc_meanvar_poisson_1d"
        if (any(x < 0.0_dp)) error stop "Poisson data must be nonnegative"

        penalty_use_name = "MBIC"
        if (present(penalty)) penalty_use_name = trim(penalty)

        allocate(cs(0:n))
        cs(0) = 0.0_dp
        do tau = 1, n
            cs(tau) = cs(tau - 1) + x(tau)
        end do

        total_sum = cs(n)
        if (total_sum <= 0.0_dp) then
            null_like = huge(1.0_dp)
        else
            null_like = 2.0_dp * total_sum * log(real(n, dp)) - 2.0_dp * total_sum * log(total_sum)
        end if

        best_alt = huge(1.0_dp)
        best_tau = minseglen_use
        do tau = minseglen_use, n - minseglen_use
            left_sum = cs(tau)
            right_sum = total_sum - left_sum
            if (left_sum <= 0.0_dp .or. right_sum <= 0.0_dp) then
                alt_use = huge(1.0_dp)
            else
                alt_use = 2.0_dp * log(real(tau, dp)) * left_sum - 2.0_dp * left_sum * log(left_sum) + &
                          2.0_dp * log(real(n - tau, dp)) * right_sum - 2.0_dp * right_sum * log(right_sum)
            end if
            if (alt_use < best_alt) then
                best_alt = alt_use
                best_tau = tau
            end if
        end do

        alt_like = best_alt
        alt_use = best_alt
        if (trim(penalty_use_name) == "MBIC") then
            alt_use = alt_use + log(real(best_tau, dp)) + log(real(n - best_tau + 1, dp))
        end if

        penalty_value = amoc_penalty_value(trim(penalty_use_name), n, 1, pen_value)
        if (null_like - alt_use >= penalty_value) then
            cpt = best_tau
        else
            cpt = n
        end if

        deallocate(cs)
    end subroutine solve_amoc_meanvar_poisson_1d

    subroutine solve_amoc_meanvar_exp_1d(x, cpt, p_value, null_like, alt_like, decision_penalty, penalty, minseglen, pen_value)
        !> Single changepoint in Exponential mean/variance, matching changepoint::single.meanvar.exp.
        real(kind=dp), intent(in) :: x(:)
        integer, intent(out) :: cpt
        real(kind=dp), intent(out) :: p_value, null_like, alt_like, decision_penalty
        character(len=*), intent(in), optional :: penalty
        integer, intent(in), optional :: minseglen
        real(kind=dp), intent(in), optional :: pen_value

        integer :: n, minseglen_use, tau, best_tau
        real(kind=dp) :: total_sum, left_sum, right_sum, alt_use
        real(kind=dp) :: best_alt, an, bn, pi_val
        real(kind=dp), allocatable :: cs(:)
        character(len=32) :: penalty_use_name

        n = size(x)
        minseglen_use = 1
        if (present(minseglen)) minseglen_use = minseglen
        if (n < 4) error stop "solve_amoc_meanvar_exp_1d requires at least 4 observations"
        if (n < 2 * minseglen_use) error stop "minseglen too large for solve_amoc_meanvar_exp_1d"
        if (any(x < 0.0_dp)) error stop "Exponential data must be nonnegative"

        penalty_use_name = "MBIC"
        if (present(penalty)) penalty_use_name = trim(penalty)

        allocate(cs(0:n))
        cs(0) = 0.0_dp
        do tau = 1, n
            cs(tau) = cs(tau - 1) + x(tau)
        end do

        total_sum = cs(n)
        if (total_sum <= 0.0_dp) then
            null_like = huge(1.0_dp)
        else
            null_like = 2.0_dp * real(n, dp) * log(total_sum) - 2.0_dp * real(n, dp) * log(real(n, dp))
        end if

        best_alt = huge(1.0_dp)
        best_tau = minseglen_use
        do tau = minseglen_use, n - minseglen_use
            left_sum = cs(tau)
            right_sum = total_sum - left_sum
            if (left_sum <= 0.0_dp .or. right_sum <= 0.0_dp) then
                alt_use = huge(1.0_dp)
            else
                alt_use = 2.0_dp * real(tau, dp) * log(left_sum) - 2.0_dp * real(tau, dp) * log(real(tau, dp)) + &
                          2.0_dp * real(n - tau, dp) * log(right_sum) - 2.0_dp * real(n - tau, dp) * log(real(n - tau, dp))
            end if
            if (alt_use < best_alt) then
                best_alt = alt_use
                best_tau = tau
            end if
        end do

        alt_like = best_alt
        alt_use = best_alt
        if (trim(penalty_use_name) == "MBIC") then
            alt_use = alt_use + log(real(best_tau, dp)) + log(real(n - best_tau + 1, dp))
        end if

        ! Mirror the package's actual scalar implementation: decision() receives the raw pen.value input.
        decision_penalty = 0.0_dp
        if (present(pen_value)) decision_penalty = pen_value
        if (null_like - alt_use >= decision_penalty) then
            cpt = best_tau
        else
            cpt = n
        end if

        pi_val = acos(-1.0_dp)
        an = sqrt(2.0_dp * log(log(real(n, dp))))
        bn = 2.0_dp * log(log(real(n, dp))) + 0.5_dp * log(log(log(real(n, dp)))) - 0.5_dp * log(pi_val)
        p_value = exp(-2.0_dp * exp(-an * sqrt(abs(null_like - alt_use)) + an * bn)) - exp(-2.0_dp * exp(an * bn))

        deallocate(cs)
    end subroutine solve_amoc_meanvar_exp_1d

    subroutine solve_amoc_meanvar_gamma_1d(x, shape, cpt, penalty_value, null_like, alt_like, penalty, minseglen, pen_value)
        !> Single changepoint in Gamma mean/variance with fixed shape, matching changepoint::single.meanvar.gamma.
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(in) :: shape
        integer, intent(out) :: cpt
        real(kind=dp), intent(out) :: penalty_value, null_like, alt_like
        character(len=*), intent(in), optional :: penalty
        integer, intent(in), optional :: minseglen
        real(kind=dp), intent(in), optional :: pen_value

        integer :: n, minseglen_use, tau, best_tau
        real(kind=dp) :: total_sum, left_sum, right_sum, alt_use, best_alt
        real(kind=dp), allocatable :: cs(:)
        character(len=32) :: penalty_use_name

        n = size(x)
        minseglen_use = 1
        if (present(minseglen)) minseglen_use = minseglen
        if (n < 4) error stop "solve_amoc_meanvar_gamma_1d requires at least 4 observations"
        if (n < 2 * minseglen_use) error stop "minseglen too large for solve_amoc_meanvar_gamma_1d"
        if (any(x <= 0.0_dp)) error stop "Gamma data must be strictly positive"
        if (shape <= 0.0_dp) error stop "Gamma shape must be positive"

        penalty_use_name = "MBIC"
        if (present(penalty)) penalty_use_name = trim(penalty)

        allocate(cs(0:n))
        cs(0) = 0.0_dp
        do tau = 1, n
            cs(tau) = cs(tau - 1) + x(tau)
        end do

        total_sum = cs(n)
        null_like = 2.0_dp * real(n, dp) * shape * log(total_sum) - 2.0_dp * real(n, dp) * shape * log(real(n, dp) * shape)

        best_alt = huge(1.0_dp)
        best_tau = minseglen_use
        do tau = minseglen_use, n - minseglen_use
            left_sum = cs(tau)
            right_sum = total_sum - left_sum
            if (left_sum <= 0.0_dp .or. right_sum <= 0.0_dp) then
                alt_use = huge(1.0_dp)
            else
                alt_use = 2.0_dp * real(tau, dp) * shape * log(left_sum) - 2.0_dp * real(tau, dp) * shape * log(real(tau, dp) * shape) + &
                          2.0_dp * real(n - tau, dp) * shape * log(right_sum) - 2.0_dp * real(n - tau, dp) * shape * log(real(n - tau, dp) * shape)
            end if
            if (alt_use < best_alt) then
                best_alt = alt_use
                best_tau = tau
            end if
        end do

        alt_like = best_alt
        alt_use = best_alt
        if (trim(penalty_use_name) == "MBIC") then
            alt_use = alt_use + log(real(best_tau, dp)) + log(real(n - best_tau + 1, dp))
        end if

        penalty_value = amoc_penalty_value(trim(penalty_use_name), n, 1, pen_value)
        if (null_like - alt_use >= penalty_value) then
            cpt = best_tau
        else
            cpt = n
        end if

        deallocate(cs)
    end subroutine solve_amoc_meanvar_gamma_1d

    subroutine solve_amoc_reg_norm_1d(data, cpt, penalty_value, null_rss, alt_rss, penalty, minseglen, pen_value)
        !> Single changepoint in a Normal linear regression model, matching changepoint::cpt.reg(method="AMOC").
        real(kind=dp), intent(in) :: data(:, :)
        integer, intent(out) :: cpt
        real(kind=dp), intent(out) :: penalty_value, null_rss, alt_rss
        character(len=*), intent(in), optional :: penalty
        integer, intent(in), optional :: minseglen
        real(kind=dp), intent(in), optional :: pen_value

        integer :: n, p, minseglen_use, tau, best_tau, diffparam
        real(kind=dp) :: alt_use, best_alt, test_stat
        character(len=32) :: penalty_use_name

        n = size(data, 1)
        p = size(data, 2) - 1
        if (p <= 0) error stop "solve_amoc_reg_norm_1d requires response plus at least one regressor"

        minseglen_use = 3
        if (present(minseglen)) minseglen_use = minseglen
        if (minseglen_use < size(data, 2) - 1) minseglen_use = size(data, 2)
        if (n < 2 * minseglen_use) error stop "minseglen too large for solve_amoc_reg_norm_1d"

        penalty_use_name = "SIC"
        if (present(penalty)) penalty_use_name = trim(penalty)
        diffparam = size(data, 2)
        penalty_value = amoc_penalty_value(trim(penalty_use_name), n, diffparam, pen_value)

        null_rss = regression_segment_rss(data, 1, n)
        best_tau = minseglen_use
        best_alt = huge(1.0_dp)
        do tau = minseglen_use, n - minseglen_use
            alt_use = regression_segment_rss(data, 1, tau) + regression_segment_rss(data, tau + 1, n)
            if (alt_use < best_alt) then
                best_alt = alt_use
                best_tau = tau
            end if
        end do

        alt_rss = best_alt
        test_stat = null_rss - alt_rss
        if (test_stat >= penalty_value) then
            cpt = best_tau
        else
            cpt = n
        end if
    end subroutine solve_amoc_reg_norm_1d

    subroutine solve_amoc_var_css_1d(x, cpt, test_stat, penalty_value, penalty, minseglen, pen_value)
        !> Single changepoint in variance using the CSS statistic, matching changepoint::single.var.css.
        real(kind=dp), intent(in) :: x(:)
        integer, intent(out) :: cpt
        real(kind=dp), intent(out) :: test_stat, penalty_value
        character(len=*), intent(in), optional :: penalty
        integer, intent(in), optional :: minseglen
        real(kind=dp), intent(in), optional :: pen_value

        integer :: n, minseglen_use, tau, best_tau, tau_idx, best_tau_idx
        real(kind=dp) :: total_ss, css_val, d_max, diff
        real(kind=dp), allocatable :: cs2(:)
        character(len=32) :: penalty_use_name

        n = size(x)
        minseglen_use = 2
        if (present(minseglen)) minseglen_use = minseglen
        if (n < 4) error stop "solve_amoc_var_css_1d requires at least 4 observations"
        if (n < 2 * minseglen_use) error stop "minseglen too large for solve_amoc_var_css_1d"

        penalty_use_name = "Manual"
        if (present(penalty)) penalty_use_name = trim(penalty)

        allocate(cs2(0:n))
        cs2(0) = 0.0_dp
        do tau = 1, n
            cs2(tau) = cs2(tau - 1) + x(tau) * x(tau)
        end do

        total_ss = cs2(n)
        d_max = -1.0_dp
        best_tau = minseglen_use
        best_tau_idx = 1
        tau_idx = 0
        do tau = minseglen_use, n - minseglen_use + 1
            tau_idx = tau_idx + 1
            css_val = cs2(tau) / total_ss - real(tau, dp) / real(n, dp)
            diff = abs(css_val)
            if (diff > d_max) then
                d_max = diff
                best_tau = tau
                best_tau_idx = tau_idx
            end if
        end do

        test_stat = sqrt(real(n, dp) / 2.0_dp) * d_max
        penalty_value = css_penalty_value(trim(penalty_use_name), n, pen_value)
        if (test_stat >= penalty_value) then
            cpt = best_tau_idx
        else
            cpt = n
        end if

        deallocate(cs2)
    end subroutine solve_amoc_var_css_1d

    real(kind=dp) function regression_segment_rss(data, s, e) result(rss)
        !> Residual sum of squares for an OLS regression on rows s:e with response in column 1.
        real(kind=dp), intent(in) :: data(:, :)
        integer, intent(in) :: s, e

        integer :: nseg, p, i, j, k, info
        real(kind=dp), allocatable :: xtx(:,:), xty(:), beta(:), resid(:)

        p = size(data, 2) - 1
        nseg = e - s + 1
        if (nseg <= 0) then
            rss = huge(1.0_dp)
            return
        end if

        allocate(xtx(p, p), xty(p), beta(p), resid(nseg))
        xtx = 0.0_dp
        xty = 0.0_dp
        do i = s, e
            do j = 1, p
                xty(j) = xty(j) + data(i, j + 1) * data(i, 1)
                do k = 1, p
                    xtx(j, k) = xtx(j, k) + data(i, j + 1) * data(i, k + 1)
                end do
            end do
        end do

        call solve_linear_system_square(xtx, xty, beta, info)
        if (info /= 0) then
            rss = huge(1.0_dp)
            deallocate(xtx, xty, beta, resid)
            return
        end if

        rss = 0.0_dp
        do i = s, e
            resid(i - s + 1) = data(i, 1) - dot_product(data(i, 2:p + 1), beta)
            rss = rss + resid(i - s + 1) * resid(i - s + 1)
        end do

        deallocate(xtx, xty, beta, resid)
    end function regression_segment_rss

    subroutine solve_linear_system_square(a, b, x, info)
        !> Dense Gaussian elimination with partial pivoting for small square systems.
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
    end subroutine solve_linear_system_square

    real(kind=dp) function css_penalty_value(penalty, n, pen_value) result(penalty_out)
        !> Penalty calculation for CSS variance AMOC, following changepoint::penalty_decision where applicable.
        character(len=*), intent(in) :: penalty
        integer, intent(in) :: n
        real(kind=dp), intent(in), optional :: pen_value

        select case (trim(penalty))
        case ("None")
            penalty_out = 0.0_dp
        case ("Manual")
            if (.not. present(pen_value)) error stop "Manual CSS penalty requires pen_value"
            penalty_out = pen_value
        case ("Asymptotic")
            if (.not. present(pen_value)) error stop "Asymptotic CSS penalty requires pen_value"
            select case (nint(100.0_dp * pen_value))
            case (1)
                penalty_out = 1.628_dp
            case (5)
                penalty_out = 1.358_dp
            case (10)
                penalty_out = 1.224_dp
            case (25)
                penalty_out = 1.019_dp
            case (50)
                penalty_out = 0.828_dp
            case (75)
                penalty_out = 0.677_dp
            case (90)
                penalty_out = 0.571_dp
            case (95)
                penalty_out = 0.520_dp
            case default
                error stop "Unsupported asymptotic CSS alpha"
            end select
        case default
            error stop "Unsupported CSS AMOC penalty"
        end select
    end function css_penalty_value

    function l1_segment_cost_1d(z, s, e, min_seg_len, huge_cost) result(cost)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: s, e, min_seg_len
        real(kind=dp), intent(in) :: huge_cost
        real(kind=dp) :: cost
        real(kind=dp), allocatable :: work(:)
        real(kind=dp) :: med
        integer :: m

        if (e < s) then
            cost = huge_cost
            return
        end if

        m = e - s + 1
        if (m < min_seg_len) then
            cost = huge_cost
            return
        end if

        allocate(work(m))
        work = z(s:e)
        med = median(work)
        cost = sum(abs(z(s:e) - med))
        deallocate(work)
    end function l1_segment_cost_1d

    function dynp_segment_cost_1d(model, z, sz, szz, s, e, min_seg_len, huge_cost) result(cost)
        character(len=*), intent(in) :: model
        real(kind=dp), intent(in) :: z(:)
        real(kind=dp), intent(in) :: sz(0:), szz(0:)
        integer, intent(in) :: s, e, min_seg_len
        real(kind=dp), intent(in) :: huge_cost
        real(kind=dp) :: cost

        select case (trim(model))
        case ("l2")
            cost = mean_shift_segment_sse_1d(sz, szz, s, e, min_seg_len, huge_cost)
        case ("normal")
            cost = normal_segment_cost_1d(sz, szz, s, e, min_seg_len, huge_cost)
        case ("l1")
            cost = l1_segment_cost_1d(z, s, e, min_seg_len, huge_cost)
        case default
            error stop "dynp_segment_cost_1d: unsupported model"
        end select
    end function dynp_segment_cost_1d

    subroutine find_window_peaks(score, inds, order, peak_bkps, peak_gain)
        real(kind=dp), intent(in) :: score(:)
        integer, intent(in) :: inds(:), order
        integer, allocatable, intent(out) :: peak_bkps(:)
        real(kind=dp), allocatable, intent(out) :: peak_gain(:)

        integer :: n, i, j, idx, n_peaks
        logical :: is_peak

        n = size(score)
        n_peaks = 0
        do i = 1, n
            is_peak = .true.
            do j = -order, order
                if (j == 0) cycle
                idx = 1 + modulo(i - 1 + j, n)
                if (score(i) <= score(idx)) then
                    is_peak = .false.
                    exit
                end if
            end do
            if (is_peak) n_peaks = n_peaks + 1
        end do

        allocate(peak_bkps(n_peaks), peak_gain(n_peaks))
        n_peaks = 0
        do i = 1, n
            is_peak = .true.
            do j = -order, order
                if (j == 0) cycle
                idx = 1 + modulo(i - 1 + j, n)
                if (score(i) <= score(idx)) then
                    is_peak = .false.
                    exit
                end if
            end do
            if (is_peak) then
                n_peaks = n_peaks + 1
                peak_bkps(n_peaks) = inds(i)
                peak_gain(n_peaks) = score(i)
            end if
        end do
    end subroutine find_window_peaks

    subroutine sort_peaks_by_gain(peak_bkps, peak_gain)
        integer, intent(inout) :: peak_bkps(:)
        real(kind=dp), intent(inout) :: peak_gain(:)
        integer :: i, j, tmp_i
        real(kind=dp) :: tmp_g

        do i = 2, size(peak_gain)
            tmp_g = peak_gain(i)
            tmp_i = peak_bkps(i)
            j = i - 1
            do while (j >= 1 .and. peak_gain(j) > tmp_g)
                peak_gain(j + 1) = peak_gain(j)
                peak_bkps(j + 1) = peak_bkps(j)
                j = j - 1
            end do
            peak_gain(j + 1) = tmp_g
            peak_bkps(j + 1) = tmp_i
        end do
    end subroutine sort_peaks_by_gain

    function gamma_from_pairwise_median(dists) result(gamma)
        real(kind=dp), intent(in) :: dists(:)
        real(kind=dp) :: gamma
        real(kind=dp), allocatable :: work(:)
        integer :: n
        real(kind=dp) :: med_lo, med_hi

        n = size(dists)
        gamma = 1.0_dp
        if (n <= 0) return

        allocate(work(n))
        work = dists

        if (mod(n, 2) == 1) then
            call select_kth_real(work, (n + 1) / 2, med_lo)
            med_hi = med_lo
        else
            call select_kth_real(work, n / 2, med_lo)
            call select_kth_real(work, n / 2 + 1, med_hi)
        end if
        if (mod(n, 2) == 1) then
            if (med_lo /= 0.0_dp) gamma = 1.0_dp / med_lo
        else
            if (0.5_dp * (med_lo + med_hi) /= 0.0_dp) gamma = 1.0_dp / (0.5_dp * (med_lo + med_hi))
        end if

        deallocate(work)
    end function gamma_from_pairwise_median

    subroutine select_kth_real(x, k, kth)
        real(kind=dp), intent(inout) :: x(:)
        integer, intent(in) :: k
        real(kind=dp), intent(out) :: kth

        integer :: left, right, pivot_index, new_pivot_index

        left = 1
        right = size(x)

        do
            if (left == right) then
                kth = x(left)
                return
            end if

            pivot_index = (left + right) / 2
            call partition_real(x, left, right, pivot_index, new_pivot_index)
            pivot_index = new_pivot_index

            if (k == pivot_index) then
                kth = x(k)
                return
            else if (k < pivot_index) then
                right = pivot_index - 1
            else
                left = pivot_index + 1
            end if
        end do
    end subroutine select_kth_real

    subroutine partition_real(x, left, right, pivot_index_in, pivot_index_out)
        real(kind=dp), intent(inout) :: x(:)
        integer, intent(in) :: left, right, pivot_index_in
        integer, intent(out) :: pivot_index_out

        real(kind=dp) :: pivot_value
        integer :: store_index, i

        pivot_value = x(pivot_index_in)
        call swap_real(x(pivot_index_in), x(right))
        store_index = left

        do i = left, right - 1
            if (x(i) < pivot_value) then
                call swap_real(x(store_index), x(i))
                store_index = store_index + 1
            end if
        end do

        call swap_real(x(right), x(store_index))
        pivot_index_out = store_index
    end subroutine partition_real

    elemental subroutine swap_real(a, b)
        real(kind=dp), intent(inout) :: a, b
        real(kind=dp) :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine swap_real

    subroutine build_rbf_gram_1d(z, gamma, gram)
        real(kind=dp), intent(in) :: z(:)
        real(kind=dp), intent(in) :: gamma
        real(kind=dp), intent(out) :: gram(:,:)
        integer :: n, i, j
        real(kind=dp) :: dz, scaled_d2

        n = size(z)
        do i = 1, n
            gram(i, i) = 1.0_dp
            do j = i + 1, n
                dz = z(i) - z(j)
                scaled_d2 = gamma * dz * dz
                scaled_d2 = max(1.0e-2_dp, min(1.0e2_dp, scaled_d2))
                gram(i, j) = exp(-scaled_d2)
                gram(j, i) = gram(i, j)
            end do
        end do
    end subroutine build_rbf_gram_1d

    subroutine build_rbf_prefix_sum_1d(z, gamma, pref)
        real(kind=dp), intent(in) :: z(:)
        real(kind=dp), intent(in) :: gamma
        real(kind=dp), intent(out) :: pref(0:, 0:)
        integer :: n, i, j
        real(kind=dp) :: dz, scaled_d2, kij

        n = size(z)
        pref(0, :) = 0.0_dp
        pref(:, 0) = 0.0_dp
        do i = 1, n
            do j = 1, n
                if (i == j) then
                    kij = 1.0_dp
                else
                    dz = z(i) - z(j)
                    scaled_d2 = gamma * dz * dz
                    scaled_d2 = max(1.0e-2_dp, min(1.0e2_dp, scaled_d2))
                    kij = exp(-scaled_d2)
                end if
                pref(i, j) = kij + pref(i - 1, j) + pref(i, j - 1) - pref(i - 1, j - 1)
            end do
        end do
    end subroutine build_rbf_prefix_sum_1d

    subroutine build_rbf_prefix_sum_mv(x, gamma, pref)
        real(kind=dp), intent(in) :: x(:,:)
        real(kind=dp), intent(in) :: gamma
        real(kind=dp), intent(out) :: pref(0:, 0:)
        integer :: n, i, j
        real(kind=dp) :: scaled_d2, kij

        n = size(x, 1)
        pref(0, :) = 0.0_dp
        pref(:, 0) = 0.0_dp
        do i = 1, n
            do j = 1, n
                if (i == j) then
                    kij = 1.0_dp
                else
                    scaled_d2 = gamma * sum((x(i, :) - x(j, :)) ** 2)
                    scaled_d2 = max(1.0e-2_dp, min(1.0e2_dp, scaled_d2))
                    kij = exp(-scaled_d2)
                end if
                pref(i, j) = kij + pref(i - 1, j) + pref(i, j - 1) - pref(i - 1, j - 1)
            end do
        end do
    end subroutine build_rbf_prefix_sum_mv

    pure function rbf_kernel_segment_cost(pref, s, e, min_seg_len, huge_cost) result(cost)
        real(kind=dp), intent(in) :: pref(0:, 0:)
        integer, intent(in) :: s, e, min_seg_len
        real(kind=dp), intent(in) :: huge_cost
        real(kind=dp) :: cost
        real(kind=dp) :: m_r, sub_sum

        if (e < s) then
            cost = huge_cost
            return
        end if

        m_r = real(e - s + 1, dp)
        if (m_r < real(min_seg_len, dp)) then
            cost = huge_cost
            return
        end if

        sub_sum = pref(e, e) - pref(s - 1, e) - pref(e, s - 1) + pref(s - 1, s - 1)
        cost = m_r - sub_sum / m_r
    end function rbf_kernel_segment_cost

    subroutine build_matrix_prefix_sum(x, pref)
        real(kind=dp), intent(in) :: x(:,:)
        real(kind=dp), intent(out) :: pref(0:, 0:)
        integer :: n, i, j

        n = size(x, 1)
        pref(0, :) = 0.0_dp
        pref(:, 0) = 0.0_dp
        do i = 1, n
            do j = 1, n
                pref(i, j) = x(i, j) + pref(i - 1, j) + pref(i, j - 1) - pref(i - 1, j - 1)
            end do
        end do
    end subroutine build_matrix_prefix_sum

    !> Log-determinant of a symmetric positive-definite matrix via Cholesky.
    !! Returns -1e30 if the matrix is not positive definite.
    pure function log_det_chol(A) result(ld)
        real(kind=dp), intent(in) :: A(:,:)
        real(kind=dp) :: ld
        integer :: p, i, j
        real(kind=dp) :: s, L(size(A,1), size(A,1))
        p = size(A, 1)
        L = 0.0_dp
        do j = 1, p
            s = A(j,j) - sum(L(j, 1:j-1)**2)
            if (s <= 0.0_dp) then
                ld = -1.0e30_dp
                return
            end if
            L(j,j) = sqrt(s)
            do i = j+1, p
                L(i,j) = (A(i,j) - sum(L(i,1:j-1)*L(j,1:j-1))) / L(j,j)
            end do
        end do
        ld = 2.0_dp * sum([(log(L(i,i)), i=1,p)])
    end function log_det_chol

    !> Joint covariance-matrix changepoint cost matrix.
    !! cost(i,j) = m/2 * log|Sigma_hat(i:j)|  where Sigma_hat is the biased sample
    !! covariance matrix of the p-column return matrix R over rows i..j.
    !! Uses O(n*p^2) prefix sums; each cell costs O(p^2) to assemble + O(p^3) Cholesky.
    function multivar_cost_matrix(R, min_seg_len) result(cost)
        real(kind=dp), intent(in) :: R(:,:)      ! n × p
        integer, intent(in), optional :: min_seg_len
        real(kind=dp) :: cost(size(R,1), size(R,1))
        real(kind=dp), allocatable :: SR(:,:), SP(:,:,:)
        real(kind=dp) :: S(size(R,2), size(R,2)), sr_s(size(R,2)), m_r, ld
        integer :: n, p, i, j, a, b, min_len

        n = size(R, 1);  p = size(R, 2)
        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        allocate(SR(0:n, p), SP(0:n, p, p))
        cost = 1.0e20_dp

        SR(0,:) = 0.0_dp;  SP(0,:,:) = 0.0_dp
        do i = 1, n
            SR(i,:) = SR(i-1,:) + R(i,:)
            do a = 1, p
                SP(i,a,:) = SP(i-1,a,:) + R(i,a) * R(i,:)
            end do
        end do

        do j = min_len, n
            do i = 1, j - min_len + 1
                m_r  = real(j - i + 1, dp)
                sr_s = SR(j,:) - SR(i-1,:)
                do a = 1, p
                    do b = a, p
                        S(a,b) = (SP(j,a,b) - SP(i-1,a,b)) / m_r &
                               - (sr_s(a)/m_r) * (sr_s(b)/m_r)
                        S(b,a) = S(a,b)
                    end do
                end do
                ld = log_det_chol(S)
                if (ld > -1.0e19_dp) cost(i,j) = 0.5_dp * m_r * ld
            end do
        end do
    end function multivar_cost_matrix

    !> Constructs a cost matrix for all possible segments using negative log-likelihood (slow O(N^3) version).
    pure function cost_matrix_slow(x, y) result(cost)
        real(kind=dp), intent(in) :: x(:), y(:)
        real(kind=dp) :: cost(size(x), size(x))
        integer :: n, i, j
        n = size(x)
        do i = 1, n
            do j = i, n
                cost(i, j) = -log_likelihood_corr(x(i:j), y(i:j))
            end do
        end do
    end function cost_matrix_slow

    subroutine solve_changepoints(max_m, cost, dp_table, parent)
        !> Solves the changepoint problem using dynamic programming to minimize total cost.
        integer, intent(in) :: max_m
        real(kind=dp), intent(in) :: cost(:, :)
        real(kind=dp), intent(out) :: dp_table(size(cost, 1), max_m)
        integer, intent(out) :: parent(size(cost, 1), max_m)
        integer :: i, m, k, n

        n = size(cost, 1)

        dp_table = 1.0e20_dp
        parent = 0
        do i = 1, n
            dp_table(i, 1) = cost(1, i)
        end do

        do m = 2, max_m
            do i = m, n
                do k = m-1, i-1
                    if (dp_table(k, m-1) + cost(k+1, i) < dp_table(i, m)) then
                        dp_table(i, m) = dp_table(k, m-1) + cost(k+1, i)
                        parent(i, m) = k
                    end if
                end do
            end do
        end do
    end subroutine solve_changepoints

    !> Solve the exact continuous piecewise linear least-squares problem for y.
    !!
    !! The fitted curve is continuous and piecewise linear, with observations 1:n.
    !! For m segments there are m-1 internal changepoints. Segment 1 fits y(1:t1).
    !! Later segments fit y(t_prev+1:t_curr), so each observation is used exactly once.
    !!
    !! The dynamic programming state for each (m,t) is represented as a lower envelope
    !! of quadratics in the endpoint value at time t. This gives an exact optimizer for
    !! the least-squares continuous PWL objective, not the greedy knot-insertion heuristic.
    subroutine solve_continuous_pwl_sse(y, max_m, min_seg_len, sse, cp_store, n_models_fitted)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: max_m
        integer, intent(in), optional :: min_seg_len
        real(kind=dp), intent(out) :: sse(max_m)
        integer, intent(out) :: cp_store(max(1, max_m - 1), max_m)
        integer, intent(out) :: n_models_fitted

        type(quad_state), allocatable :: states(:, :)
        real(kind=dp), allocatable :: sy(:), siy(:), syy(:)
        real(kind=dp), allocatable :: raw_a(:), raw_b(:), raw_c(:)
        integer, allocatable :: raw_pt(:), raw_pq(:), keep(:)
        integer :: n, min_len, max_m_eff, m, t, s, j, nraw, nkeep, best_q
        real(kind=dp) :: a0, b0, c0, qa, qb, qc, qd, qe, qf
        real(kind=dp) :: best_cost, denom, lin

        min_len = 50
        if (present(min_seg_len)) min_len = min_seg_len

        n = size(y)
        if (max_m < 1) error stop 'solve_continuous_pwl_sse: max_m must be positive'
        if (n < min_len) error stop 'solve_continuous_pwl_sse: size(y) < min_seg_len'

        max_m_eff = min(max_m, n / min_len)
        sse = huge(1.0_dp)
        cp_store = 0
        n_models_fitted = 0

        call build_prefix_sums(y, sy, siy, syy)
        allocate(states(max_m_eff, n))

        do t = min_len, n
            call first_segment_quadratic(t, sy, siy, syy, qa, qb, qc)
            call set_state_single(states(1, t), qa, qb, qc)
        end do

        do m = 2, max_m_eff
            do t = m * min_len, n
                nraw = 0
                do s = (m - 1) * min_len, t - min_len
                    nraw = nraw + states(m - 1, s)%nquad
                end do
                if (nraw == 0) cycle

                allocate(raw_a(nraw), raw_b(nraw), raw_c(nraw), raw_pt(nraw), raw_pq(nraw))
                nraw = 0
                do s = (m - 1) * min_len, t - min_len
                    if (states(m - 1, s)%nquad == 0) cycle
                    call later_segment_coefficients(s, t, sy, siy, syy, qa, qb, qc, qd, qe, qf)
                    do j = 1, states(m - 1, s)%nquad
                        a0 = states(m - 1, s)%a(j)
                        b0 = states(m - 1, s)%b(j)
                        c0 = states(m - 1, s)%c(j)
                        denom = a0 + qa
                        lin = b0 + qd
                        nraw = nraw + 1
                        raw_a(nraw) = qc - (qb * qb) / denom
                        raw_b(nraw) = qe - qb * lin / denom
                        raw_c(nraw) = c0 + qf - 0.25_dp * lin * lin / denom
                        raw_pt(nraw) = s
                        raw_pq(nraw) = j
                    end do
                end do

                call keep_lower_envelope(raw_a, raw_b, raw_c, keep, nkeep)
                call set_state_subset(states(m, t), raw_a, raw_b, raw_c, raw_pt, raw_pq, keep, nkeep)
                deallocate(raw_a, raw_b, raw_c, raw_pt, raw_pq, keep)
            end do
        end do

        do m = 1, max_m_eff
            if (states(m, n)%nquad == 0) cycle
            call best_terminal_quadratic(states(m, n)%a, states(m, n)%b, states(m, n)%c, best_q, best_cost)
            sse(m) = best_cost
            call backtrack_pwl_cps(states, m, n, best_q, cp_store(1:m-1, m))
            n_models_fitted = m
        end do

        deallocate(states, sy, siy, syy)
    end subroutine solve_continuous_pwl_sse

    !> Fit a continuous piecewise linear regression to y for a fixed changepoint set.
    subroutine fit_continuous_pwl_given_cps(y, cps, fitted, rss)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: cps(:)
        real(kind=dp), intent(out) :: fitted(:)
        real(kind=dp), intent(out), optional :: rss

        integer :: n, p, i, k
        real(kind=dp), allocatable :: xtx(:, :), xty(:), beta(:)
        real(kind=dp) :: t

        n = size(y)
        if (size(fitted) /= n) error stop 'fit_continuous_pwl_given_cps: size(fitted) /= size(y)'

        p = size(cps) + 2
        allocate(xtx(p, p), xty(p), beta(p))
        xtx = 0.0_dp
        xty = 0.0_dp

        do i = 1, n
            t = real(i, dp)
            call accumulate_normal_equations(t, y(i), cps, xtx, xty)
        end do

        call solve_linear_system(xtx, xty, beta)

        do i = 1, n
            t = real(i, dp)
            fitted(i) = beta(1) + beta(2) * t
            do k = 1, size(cps)
                fitted(i) = fitted(i) + beta(k + 2) * max(0.0_dp, t - real(cps(k), dp))
            end do
        end do

        if (present(rss)) rss = sum((y - fitted)**2)

        deallocate(xtx, xty, beta)
    end subroutine fit_continuous_pwl_given_cps

    subroutine build_prefix_sums(y, sy, siy, syy)
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), allocatable, intent(out) :: sy(:), siy(:), syy(:)
        integer :: n, i

        n = size(y)
        allocate(sy(0:n), siy(0:n), syy(0:n))
        sy = 0.0_dp
        siy = 0.0_dp
        syy = 0.0_dp
        do i = 1, n
            sy(i) = sy(i-1) + y(i)
            siy(i) = siy(i-1) + real(i, dp) * y(i)
            syy(i) = syy(i-1) + y(i) * y(i)
        end do
    end subroutine build_prefix_sums

    subroutine first_segment_quadratic(t, sy, siy, syy, a, b, c)
        integer, intent(in) :: t
        real(kind=dp), intent(in) :: sy(0:), siy(0:), syy(0:)
        real(kind=dp), intent(out) :: a, b, c
        real(kind=dp) :: qa, qb, qc, qd, qe, qf, d, sumy, sumiy, sumjy, ya, yb, denom

        d = real(t - 1, dp)
        sumy = sy(t)
        sumiy = siy(t)
        sumjy = sumiy - sumy

        qa = (d + 1.0_dp) * (2.0_dp * d + 1.0_dp) / (6.0_dp * d)
        qb = (d * d - 1.0_dp) / (6.0_dp * d)
        qc = qa
        ya = sumy - sumjy / d
        yb = sumjy / d
        qd = -2.0_dp * ya
        qe = -2.0_dp * yb
        qf = syy(t)

        denom = qa
        a = qc - (qb * qb) / denom
        b = qe - qb * qd / denom
        c = qf - 0.25_dp * qd * qd / denom
    end subroutine first_segment_quadratic

    subroutine later_segment_coefficients(s, t, sy, siy, syy, a, b, c, dcoef, ecoef, fcoef)
        integer, intent(in) :: s, t
        real(kind=dp), intent(in) :: sy(0:), siy(0:), syy(0:)
        real(kind=dp), intent(out) :: a, b, c, dcoef, ecoef, fcoef
        real(kind=dp) :: len, sumy, sumiy, sumjy, ya, yb

        len = real(t - s, dp)
        sumy = sy(t) - sy(s)
        sumiy = siy(t) - siy(s)
        sumjy = sumiy - real(s, dp) * sumy

        a = (len - 1.0_dp) * (2.0_dp * len - 1.0_dp) / (6.0_dp * len)
        b = (len * len - 1.0_dp) / (6.0_dp * len)
        c = (len + 1.0_dp) * (2.0_dp * len + 1.0_dp) / (6.0_dp * len)
        ya = sumy - sumjy / len
        yb = sumjy / len
        dcoef = -2.0_dp * ya
        ecoef = -2.0_dp * yb
        fcoef = syy(t) - syy(s)
    end subroutine later_segment_coefficients

    subroutine set_state_single(state, a, b, c)
        type(quad_state), intent(inout) :: state
        real(kind=dp), intent(in) :: a, b, c

        state%nquad = 1
        allocate(state%a(1), state%b(1), state%c(1), state%parent_t(1), state%parent_q(1))
        state%a(1) = a
        state%b(1) = b
        state%c(1) = c
        state%parent_t(1) = 0
        state%parent_q(1) = 0
    end subroutine set_state_single

    subroutine set_state_subset(state, a, b, c, parent_t, parent_q, keep, nkeep)
        type(quad_state), intent(inout) :: state
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        integer, intent(in) :: parent_t(:), parent_q(:), keep(:), nkeep
        integer :: i, idx

        state%nquad = nkeep
        allocate(state%a(nkeep), state%b(nkeep), state%c(nkeep), state%parent_t(nkeep), state%parent_q(nkeep))
        do i = 1, nkeep
            idx = keep(i)
            state%a(i) = a(idx)
            state%b(i) = b(idx)
            state%c(i) = c(idx)
            state%parent_t(i) = parent_t(idx)
            state%parent_q(i) = parent_q(idx)
        end do
    end subroutine set_state_subset

    subroutine best_terminal_quadratic(a, b, c, best_q, best_cost)
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        integer, intent(out) :: best_q
        real(kind=dp), intent(out) :: best_cost
        integer :: j
        real(kind=dp) :: val

        best_q = 1
        best_cost = c(1) - 0.25_dp * b(1) * b(1) / a(1)
        do j = 2, size(a)
            val = c(j) - 0.25_dp * b(j) * b(j) / a(j)
            if (val < best_cost) then
                best_cost = val
                best_q = j
            end if
        end do
    end subroutine best_terminal_quadratic

    subroutine backtrack_pwl_cps(states, m, t, q, cps)
        type(quad_state), intent(in) :: states(:, :)
        integer, intent(in) :: m, t, q
        integer, intent(out) :: cps(:)
        integer :: seg, cur_t, cur_q

        if (size(cps) /= max(0, m - 1)) error stop 'backtrack_pwl_cps: invalid cps size'

        cur_t = t
        cur_q = q
        do seg = m, 2, -1
            cps(seg - 1) = states(seg, cur_t)%parent_t(cur_q)
            cur_q = states(seg, cur_t)%parent_q(cur_q)
            cur_t = cps(seg - 1)
        end do
    end subroutine backtrack_pwl_cps

    subroutine keep_lower_envelope(a, b, c, keep, nkeep)
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        integer, allocatable, intent(out) :: keep(:)
        integer, intent(out) :: nkeep

        integer :: k, i, j, nroot, nuniq, idx, nsurvive
        real(kind=dp), allocatable :: roots(:), roots_u(:), aa_s(:), bb_s(:), cc_s(:)
        logical, allocatable :: active(:), survive(:)
        real(kind=dp) :: aa, bb, cc, disc, r1, r2, x, delta

        k = size(a)
        if (k == 0) then
            allocate(keep(0))
            nkeep = 0
            return
        end if

        allocate(survive(k))
        survive = .true.
        do i = 1, k
            do j = 1, k
                if (i == j) cycle
                if (quadratic_dominated(a(i), b(i), c(i), a(j), b(j), c(j))) then
                    survive(i) = .false.
                    exit
                end if
            end do
        end do

        nsurvive = count(survive)
        allocate(aa_s(nsurvive), bb_s(nsurvive), cc_s(nsurvive), keep(nsurvive))
        idx = 0
        do i = 1, k
            if (survive(i)) then
                idx = idx + 1
                aa_s(idx) = a(i)
                bb_s(idx) = b(i)
                cc_s(idx) = c(i)
                keep(idx) = i
            end if
        end do
        deallocate(survive)

        allocate(active(nsurvive))
        active = .false.
        allocate(roots(max(1, nsurvive * (nsurvive - 1))))
        nroot = 0

        do i = 1, nsurvive - 1
            do j = i + 1, nsurvive
                aa = aa_s(i) - aa_s(j)
                bb = bb_s(i) - bb_s(j)
                cc = cc_s(i) - cc_s(j)
                if (abs(aa) <= 1.0e-12_dp) then
                    if (abs(bb) > 1.0e-12_dp) then
                        nroot = nroot + 1
                        roots(nroot) = -cc / bb
                    end if
                else
                    disc = bb * bb - 4.0_dp * aa * cc
                    if (disc >= -1.0e-10_dp) then
                        disc = max(disc, 0.0_dp)
                        r1 = (-bb - sqrt(disc)) / (2.0_dp * aa)
                        r2 = (-bb + sqrt(disc)) / (2.0_dp * aa)
                        nroot = nroot + 1
                        roots(nroot) = r1
                        if (abs(r2 - r1) > 1.0e-10_dp) then
                            nroot = nroot + 1
                            roots(nroot) = r2
                        end if
                    end if
                end if
            end do
        end do

        if (nroot > 0) then
            call sort_real(roots(1:nroot))
            allocate(roots_u(nroot))
            nuniq = 1
            roots_u(1) = roots(1)
            do i = 2, nroot
                if (abs(roots(i) - roots_u(nuniq)) > 1.0e-8_dp) then
                    nuniq = nuniq + 1
                    roots_u(nuniq) = roots(i)
                end if
            end do

            do i = 1, nuniq
                call mark_minimizers(roots_u(i), aa_s, bb_s, cc_s, active)
            end do

            delta = abs(roots_u(1)) + 1.0_dp
            call mark_minimizers(roots_u(1) - delta, aa_s, bb_s, cc_s, active)
            do i = 1, nuniq - 1
                x = 0.5_dp * (roots_u(i) + roots_u(i + 1))
                call mark_minimizers(x, aa_s, bb_s, cc_s, active)
            end do
            delta = abs(roots_u(nuniq)) + 1.0_dp
            call mark_minimizers(roots_u(nuniq) + delta, aa_s, bb_s, cc_s, active)
            deallocate(roots_u)
        else
            call mark_minimizers(0.0_dp, aa_s, bb_s, cc_s, active)
        end if

        nkeep = count(active)
        if (nkeep == 0) then
            call mark_minimizers(0.0_dp, aa_s, bb_s, cc_s, active)
            nkeep = count(active)
        end if

        keep(1:nkeep) = pack(keep, active)
        if (nkeep < size(keep)) keep = keep(1:nkeep)

        deallocate(active, roots, aa_s, bb_s, cc_s)
    end subroutine keep_lower_envelope

    logical pure function quadratic_dominated(a1, b1, c1, a2, b2, c2)
        real(kind=dp), intent(in) :: a1, b1, c1, a2, b2, c2
        real(kind=dp) :: aa, bb, cc, disc

        aa = a1 - a2
        bb = b1 - b2
        cc = c1 - c2
        if (aa < -1.0e-12_dp) then
            quadratic_dominated = .false.
        else if (abs(aa) <= 1.0e-12_dp) then
            if (abs(bb) <= 1.0e-12_dp) then
                quadratic_dominated = (cc >= -1.0e-12_dp)
            else
                quadratic_dominated = .false.
            end if
        else
            disc = bb * bb - 4.0_dp * aa * cc
            quadratic_dominated = (disc <= 1.0e-10_dp)
        end if
    end function quadratic_dominated

    subroutine mark_minimizers(x, a, b, c, active)
        real(kind=dp), intent(in) :: x
        real(kind=dp), intent(in) :: a(:), b(:), c(:)
        logical, intent(inout) :: active(:)
        integer :: i
        real(kind=dp) :: minv, v, tol

        minv = a(1) * x * x + b(1) * x + c(1)
        do i = 2, size(a)
            v = a(i) * x * x + b(i) * x + c(i)
            if (v < minv) minv = v
        end do
        tol = 1.0e-8_dp * max(1.0_dp, abs(minv))
        do i = 1, size(a)
            v = a(i) * x * x + b(i) * x + c(i)
            if (v <= minv + tol) active(i) = .true.
        end do
    end subroutine mark_minimizers

    subroutine sort_real(x)
        real(kind=dp), intent(inout) :: x(:)
        integer :: i, j
        real(kind=dp) :: tmp

        do i = 2, size(x)
            tmp = x(i)
            j = i - 1
            do while (j >= 1 .and. x(j) > tmp)
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = tmp
        end do
    end subroutine sort_real

    subroutine accumulate_normal_equations(t, yval, cps, xtx, xty)
        real(kind=dp), intent(in) :: t, yval
        integer, intent(in) :: cps(:)
        real(kind=dp), intent(inout) :: xtx(:, :), xty(:)

        integer :: p, i, j
        real(kind=dp), allocatable :: x(:)

        p = size(cps) + 2
        allocate(x(p))
        x(1) = 1.0_dp
        x(2) = t
        do i = 1, size(cps)
            x(i + 2) = max(0.0_dp, t - real(cps(i), dp))
        end do

        do i = 1, p
            xty(i) = xty(i) + x(i) * yval
            do j = 1, p
                xtx(i, j) = xtx(i, j) + x(i) * x(j)
            end do
        end do

        deallocate(x)
    end subroutine accumulate_normal_equations

    subroutine solve_linear_system(a, b, x)
        real(kind=dp), intent(in) :: a(:, :), b(:)
        real(kind=dp), intent(out) :: x(:)

        integer :: n, i, k, ipiv
        real(kind=dp), allocatable :: aa(:, :), bb(:), rowtmp(:)
        real(kind=dp) :: piv, factor, best

        n = size(b)
        allocate(aa(n, n), bb(n), rowtmp(n))
        aa = a
        bb = b

        do k = 1, n - 1
            ipiv = k
            best = abs(aa(k, k))
            do i = k + 1, n
                if (abs(aa(i, k)) > best) then
                    best = abs(aa(i, k))
                    ipiv = i
                end if
            end do
            if (best <= 1.0e-12_dp) error stop 'solve_linear_system: singular matrix'

            if (ipiv /= k) then
                rowtmp = aa(k, :)
                aa(k, :) = aa(ipiv, :)
                aa(ipiv, :) = rowtmp
                piv = bb(k)
                bb(k) = bb(ipiv)
                bb(ipiv) = piv
            end if

            do i = k + 1, n
                factor = aa(i, k) / aa(k, k)
                aa(i, k:n) = aa(i, k:n) - factor * aa(k, k:n)
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do
        if (abs(aa(n, n)) <= 1.0e-12_dp) error stop 'solve_linear_system: singular matrix'

        x(n) = bb(n) / aa(n, n)
        do i = n - 1, 1, -1
            x(i) = (bb(i) - sum(aa(i, i+1:n) * x(i+1:n))) / aa(i, i)
        end do

        deallocate(aa, bb, rowtmp)
    end subroutine solve_linear_system

    pure function segment_ends(parent, ms) result(seg_ends)
        !> Backtrack through the DP parent table to recover the ms segment end-points,
        !! then sort into ascending order.
        integer, intent(in) :: parent(:,:)
        integer, intent(in) :: ms
        integer :: seg_ends(ms)
        integer :: n, k, cp_idx
        n            = size(parent, 1)
        cp_idx       = n
        seg_ends(ms) = n
        do k = ms, 2, -1
            seg_ends(k-1) = parent(cp_idx, k)
            cp_idx        = seg_ends(k-1)
        end do
        if (ms > 1) call sort_int(seg_ends(1:ms-1))
    end function segment_ends

    subroutine solve_bocpd_normal_gamma_1d(x, lam, map_cps, runlen_argmax, m0, r0, s0, v0)
        !> Bayesian online changepoint detection for a univariate Gaussian series
        !! with a Normal-Gamma prior. This stores only the per-time MAP run length,
        !! which is sufficient to reconstruct changepoints via map_changepoints_bocpd.
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(in) :: lam
        integer, allocatable, intent(out) :: map_cps(:)
        integer, allocatable, intent(out) :: runlen_argmax(:)
        real(kind=dp), intent(in), optional :: m0, r0, s0, v0

        integer :: n, t, i, nold
        real(kind=dp) :: prior_m, prior_r, prior_s, prior_v, hazard
        real(kind=dp) :: r0_mass, r_sum, pp, r_seen
        real(kind=dp), parameter :: cdf_threshold = 1.0e-3_dp
        real(kind=dp), allocatable :: r_old(:), r_new(:)
        real(kind=dp), allocatable :: m_old(:), rr_old(:), s_old(:), v_old(:)
        real(kind=dp), allocatable :: m_new(:), rr_new(:), s_new(:), v_new(:)

        n = size(x)
        if (n <= 0) then
            allocate(map_cps(0), runlen_argmax(0))
            return
        end if

        prior_m = 0.0_dp
        prior_r = 1.0_dp
        prior_s = 1.0_dp
        prior_v = 1.0_dp
        if (present(m0)) prior_m = m0
        if (present(r0)) prior_r = r0
        if (present(s0)) prior_s = s0
        if (present(v0)) prior_v = v0
        hazard = 1.0_dp / lam

        allocate(runlen_argmax(n))
        allocate(r_old(n), r_new(n), m_old(n), rr_old(n), s_old(n), v_old(n))
        allocate(m_new(n), rr_new(n), s_new(n), v_new(n))

        r_old = 0.0_dp
        m_old = 0.0_dp
        rr_old = 0.0_dp
        s_old = 0.0_dp
        v_old = 0.0_dp

        runlen_argmax(1) = 0
        r_old(1) = 1.0_dp
        call normal_gamma_observe(prior_m, prior_r, prior_s, prior_v, x(1), m_old(1), rr_old(1), s_old(1), v_old(1))

        do t = 2, n
            nold = t - 1
            r_new(1:t) = 0.0_dp
            r0_mass = 0.0_dp
            r_sum = 0.0_dp
            r_seen = 0.0_dp

            do i = nold, 1, -1
                if (r_old(i) <= 0.0_dp) cycle
                pp = exp(normal_gamma_log_predictive(x(t), m_old(i), rr_old(i), s_old(i), v_old(i)))
                r_seen = r_seen + r_old(i)
                r_new(i + 1) = r_old(i) * pp * (1.0_dp - hazard)
                r0_mass = r0_mass + r_old(i) * pp * hazard
                r_sum = r_sum + r_new(i + 1)
                if (1.0_dp - r_seen < cdf_threshold) exit
            end do

            r_sum = r_sum + r0_mass
            if (r_sum <= 0.0_dp) then
                r_new(1:t) = 0.0_dp
                r_new(1) = 1.0_dp
            else
                r_new(1) = r0_mass
                r_new(1:t) = r_new(1:t) / r_sum
            end if

            call normal_gamma_observe(prior_m, prior_r, prior_s, prior_v, x(t), m_new(1), rr_new(1), s_new(1), v_new(1))
            do i = 1, nold
                call normal_gamma_observe(m_old(i), rr_old(i), s_old(i), v_old(i), x(t), m_new(i + 1), rr_new(i + 1), s_new(i + 1), v_new(i + 1))
            end do

            runlen_argmax(t) = maxloc(r_new(1:t), dim=1) - 1
            r_old(1:t) = r_new(1:t)
            m_old(1:t) = m_new(1:t)
            rr_old(1:t) = rr_new(1:t)
            s_old(1:t) = s_new(1:t)
            v_old(1:t) = v_new(1:t)
        end do

        call map_changepoints_bocpd(runlen_argmax, map_cps)

        deallocate(r_old, r_new, m_old, rr_old, s_old, v_old, m_new, rr_new, s_new, v_new)
    end subroutine solve_bocpd_normal_gamma_1d

    subroutine solve_bocd_student_t_1d(x, lam, map_cps, runlen_argmax, mu0, kappa0, alpha0, beta0)
        !> Exact replication of the local `bocd` package with ConstantHazard and StudentT.
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(in) :: lam
        integer, allocatable, intent(out) :: map_cps(:)
        integer, allocatable, intent(out) :: runlen_argmax(:)
        real(kind=dp), intent(in), optional :: mu0, kappa0, alpha0, beta0

        integer :: n, t, i, nold
        real(kind=dp) :: prior_mu, prior_kappa, prior_alpha, prior_beta, hazard, pp, norm
        real(kind=dp), allocatable :: belief_old(:), belief_new(:)
        real(kind=dp), allocatable :: mu_old(:), kappa_old(:), alpha_old(:), beta_old(:)
        real(kind=dp), allocatable :: mu_new(:), kappa_new(:), alpha_new(:), beta_new(:)

        n = size(x)
        if (n <= 0) then
            allocate(map_cps(0), runlen_argmax(0))
            return
        end if

        prior_mu = 0.0_dp
        prior_kappa = 1.0_dp
        prior_alpha = 1.0_dp
        prior_beta = 1.0_dp
        if (present(mu0)) prior_mu = mu0
        if (present(kappa0)) prior_kappa = kappa0
        if (present(alpha0)) prior_alpha = alpha0
        if (present(beta0)) prior_beta = beta0
        hazard = 1.0_dp / lam

        allocate(runlen_argmax(n))
        allocate(belief_old(n), belief_new(n), mu_old(n), kappa_old(n), alpha_old(n), beta_old(n))
        allocate(mu_new(n), kappa_new(n), alpha_new(n), beta_new(n))

        belief_old = 0.0_dp
        mu_old = 0.0_dp
        kappa_old = 0.0_dp
        alpha_old = 0.0_dp
        beta_old = 0.0_dp
        belief_old(1) = 1.0_dp
        mu_old(1) = prior_mu
        kappa_old(1) = prior_kappa
        alpha_old(1) = prior_alpha
        beta_old(1) = prior_beta

        do t = 1, n
            nold = t
            belief_new(1:t+1) = 0.0_dp
            do i = 1, nold
                pp = student_t_bocd_pdf(x(t), mu_old(i), kappa_old(i), alpha_old(i), beta_old(i))
                belief_new(i + 1) = belief_old(i) * pp * (1.0_dp - hazard)
                belief_new(1) = belief_new(1) + belief_old(i) * pp * hazard
            end do
            norm = sum(belief_new(1:t+1))
            if (norm > 0.0_dp) then
                belief_new(1:t+1) = belief_new(1:t+1) / norm
            else
                belief_new(1:t+1) = 0.0_dp
                belief_new(1) = 1.0_dp
            end if

            mu_new(1) = prior_mu
            kappa_new(1) = prior_kappa
            alpha_new(1) = prior_alpha
            beta_new(1) = prior_beta
            do i = 1, nold
                call student_t_bocd_observe(mu_old(i), kappa_old(i), alpha_old(i), beta_old(i), prior_beta, x(t), &
                                            mu_new(i + 1), kappa_new(i + 1), alpha_new(i + 1), beta_new(i + 1))
            end do

            runlen_argmax(t) = maxloc(belief_new(1:t+1), dim=1) - 1
            belief_old(1:t+1) = belief_new(1:t+1)
            mu_old(1:t+1) = mu_new(1:t+1)
            kappa_old(1:t+1) = kappa_new(1:t+1)
            alpha_old(1:t+1) = alpha_new(1:t+1)
            beta_old(1:t+1) = beta_new(1:t+1)
        end do

        call map_changepoints_bocpd(runlen_argmax, map_cps)
        deallocate(belief_old, belief_new, mu_old, kappa_old, alpha_old, beta_old, mu_new, kappa_new, alpha_new, beta_new)
    end subroutine solve_bocd_student_t_1d

    subroutine solve_bocpd_poisson_gamma_1d(x, lam, map_cps, runlen_argmax, shape0, rate0)
        !> Bayesian online changepoint detection for a univariate count series
        !! with a Gamma prior on the Poisson rate.
        integer, intent(in) :: x(:)
        real(kind=dp), intent(in) :: lam
        integer, allocatable, intent(out) :: map_cps(:)
        integer, allocatable, intent(out) :: runlen_argmax(:)
        real(kind=dp), intent(in), optional :: shape0, rate0

        integer :: n, t, i, nold
        real(kind=dp) :: prior_shape, prior_rate, hazard
        real(kind=dp) :: r0_mass, r_sum, pp, r_seen
        real(kind=dp), parameter :: cdf_threshold = 1.0e-3_dp
        real(kind=dp), allocatable :: r_old(:), r_new(:)
        real(kind=dp), allocatable :: shape_old(:), rate_old(:), shape_new(:), rate_new(:)

        n = size(x)
        if (n <= 0) then
            allocate(map_cps(0), runlen_argmax(0))
            return
        end if

        prior_shape = 1.0_dp
        prior_rate = 1.0_dp
        if (present(shape0)) prior_shape = shape0
        if (present(rate0)) prior_rate = rate0
        hazard = 1.0_dp / lam

        allocate(runlen_argmax(n))
        allocate(r_old(n), r_new(n), shape_old(n), rate_old(n), shape_new(n), rate_new(n))

        r_old = 0.0_dp
        shape_old = 0.0_dp
        rate_old = 0.0_dp

        runlen_argmax(1) = 0
        r_old(1) = 1.0_dp
        call poisson_gamma_observe(prior_shape, prior_rate, x(1), shape_old(1), rate_old(1))

        do t = 2, n
            nold = t - 1
            r_new(1:t) = 0.0_dp
            r0_mass = 0.0_dp
            r_sum = 0.0_dp
            r_seen = 0.0_dp

            do i = nold, 1, -1
                if (r_old(i) <= 0.0_dp) cycle
                pp = exp(poisson_gamma_log_predictive(x(t), shape_old(i), rate_old(i)))
                r_seen = r_seen + r_old(i)
                r_new(i + 1) = r_old(i) * pp * (1.0_dp - hazard)
                r0_mass = r0_mass + r_old(i) * pp * hazard
                r_sum = r_sum + r_new(i + 1)
                if (1.0_dp - r_seen < cdf_threshold) exit
            end do

            r_sum = r_sum + r0_mass
            if (r_sum <= 0.0_dp) then
                r_new(1:t) = 0.0_dp
                r_new(1) = 1.0_dp
            else
                r_new(1) = r0_mass
                r_new(1:t) = r_new(1:t) / r_sum
            end if

            call poisson_gamma_observe(prior_shape, prior_rate, x(t), shape_new(1), rate_new(1))
            do i = 1, nold
                call poisson_gamma_observe(shape_old(i), rate_old(i), x(t), shape_new(i + 1), rate_new(i + 1))
            end do

            runlen_argmax(t) = maxloc(r_new(1:t), dim=1) - 1
            r_old(1:t) = r_new(1:t)
            shape_old(1:t) = shape_new(1:t)
            rate_old(1:t) = rate_new(1:t)
        end do

        call map_changepoints_bocpd(runlen_argmax, map_cps)
        deallocate(r_old, r_new, shape_old, rate_old, shape_new, rate_new)
    end subroutine solve_bocpd_poisson_gamma_1d

    subroutine solve_bocpd_beta_bernoulli_1d(x, lam, map_cps, runlen_argmax, alpha0, beta0)
        !> Bayesian online changepoint detection for a binary series
        !! with a Beta prior on the Bernoulli success probability.
        logical, intent(in) :: x(:)
        real(kind=dp), intent(in) :: lam
        integer, allocatable, intent(out) :: map_cps(:)
        integer, allocatable, intent(out) :: runlen_argmax(:)
        real(kind=dp), intent(in), optional :: alpha0, beta0

        integer :: n, t, i, nold
        real(kind=dp) :: prior_alpha, prior_beta, hazard
        real(kind=dp) :: r0_mass, r_sum, pp, r_seen
        real(kind=dp), parameter :: cdf_threshold = 1.0e-3_dp
        real(kind=dp), allocatable :: r_old(:), r_new(:)
        real(kind=dp), allocatable :: alpha_old(:), beta_old(:), alpha_new(:), beta_new(:)

        n = size(x)
        if (n <= 0) then
            allocate(map_cps(0), runlen_argmax(0))
            return
        end if

        prior_alpha = 0.5_dp
        prior_beta = 0.5_dp
        if (present(alpha0)) prior_alpha = alpha0
        if (present(beta0)) prior_beta = beta0
        hazard = 1.0_dp / lam

        allocate(runlen_argmax(n))
        allocate(r_old(n), r_new(n), alpha_old(n), beta_old(n), alpha_new(n), beta_new(n))

        r_old = 0.0_dp
        alpha_old = 0.0_dp
        beta_old = 0.0_dp

        runlen_argmax(1) = 0
        r_old(1) = 1.0_dp
        call beta_bernoulli_observe(prior_alpha, prior_beta, x(1), alpha_old(1), beta_old(1))

        do t = 2, n
            nold = t - 1
            r_new(1:t) = 0.0_dp
            r0_mass = 0.0_dp
            r_sum = 0.0_dp
            r_seen = 0.0_dp

            do i = nold, 1, -1
                if (r_old(i) <= 0.0_dp) cycle
                pp = exp(beta_bernoulli_log_predictive(x(t), alpha_old(i), beta_old(i)))
                r_seen = r_seen + r_old(i)
                r_new(i + 1) = r_old(i) * pp * (1.0_dp - hazard)
                r0_mass = r0_mass + r_old(i) * pp * hazard
                r_sum = r_sum + r_new(i + 1)
                if (1.0_dp - r_seen < cdf_threshold) exit
            end do

            r_sum = r_sum + r0_mass
            if (r_sum <= 0.0_dp) then
                r_new(1:t) = 0.0_dp
                r_new(1) = 1.0_dp
            else
                r_new(1) = r0_mass
                r_new(1:t) = r_new(1:t) / r_sum
            end if

            call beta_bernoulli_observe(prior_alpha, prior_beta, x(t), alpha_new(1), beta_new(1))
            do i = 1, nold
                call beta_bernoulli_observe(alpha_old(i), beta_old(i), x(t), alpha_new(i + 1), beta_new(i + 1))
            end do

            runlen_argmax(t) = maxloc(r_new(1:t), dim=1) - 1
            r_old(1:t) = r_new(1:t)
            alpha_old(1:t) = alpha_new(1:t)
            beta_old(1:t) = beta_new(1:t)
        end do

        call map_changepoints_bocpd(runlen_argmax, map_cps)
        deallocate(r_old, r_new, alpha_old, beta_old, alpha_new, beta_new)
    end subroutine solve_bocpd_beta_bernoulli_1d

    subroutine solve_bocpd_ar_1d(x, lam_logit_h, map_cps, runlen_argmax, max_lag, hazard_a, hazard_b, epsilon)
        !> Online changepoint detection using an AR predictive model and logistic hazard.
        !! This is a pragmatic AR approximation to ArgpCpd rather than a full GP implementation.
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp), intent(in) :: lam_logit_h
        integer, allocatable, intent(out) :: map_cps(:)
        integer, allocatable, intent(out) :: runlen_argmax(:)
        integer, intent(in), optional :: max_lag
        real(kind=dp), intent(in), optional :: hazard_a, hazard_b, epsilon

        integer :: n, t, i, nold, lag_order
        real(kind=dp) :: h_a, h_b, eps
        real(kind=dp) :: r0_mass, r_sum, pp, r_seen, hazard
        real(kind=dp), allocatable :: r_old(:), r_new(:)

        n = size(x)
        if (n <= 0) then
            allocate(map_cps(0), runlen_argmax(0))
            return
        end if

        lag_order = 3
        if (present(max_lag)) lag_order = max_lag
        h_a = 1.0_dp
        if (present(hazard_a)) h_a = hazard_a
        h_b = 1.0_dp
        if (present(hazard_b)) h_b = hazard_b
        eps = 1.0e-10_dp
        if (present(epsilon)) eps = epsilon

        allocate(runlen_argmax(n), r_old(n), r_new(n))
        r_old = 0.0_dp
        r_old(1) = 1.0_dp
        runlen_argmax(1) = 0

        do t = 2, n
            nold = t - 1
            r_new(1:t) = 0.0_dp
            r0_mass = 0.0_dp
            r_sum = 0.0_dp
            r_seen = 0.0_dp

            do i = nold, 1, -1
                if (r_old(i) <= 0.0_dp) cycle
                pp = exp(ar_segment_log_predictive(x, t - i, t - 1, x(t), lag_order))
                r_seen = r_seen + r_old(i)
                hazard = logistic_hazard_value(real(i, dp), lam_logit_h, h_a, h_b)
                r_new(i + 1) = r_old(i) * pp * (1.0_dp - hazard)
                r0_mass = r0_mass + r_old(i) * pp * hazard
                r_sum = r_sum + r_new(i + 1)
                if (1.0_dp - r_seen < eps) exit
            end do

            r_sum = r_sum + r0_mass
            if (r_sum <= 0.0_dp) then
                r_new(1:t) = 0.0_dp
                r_new(1) = 1.0_dp
            else
                r_new(1) = r0_mass
                r_new(1:t) = r_new(1:t) / r_sum
            end if

            runlen_argmax(t) = maxloc(r_new(1:t), dim=1) - 1
            r_old(1:t) = r_new(1:t)
        end do

        call map_changepoints_bocpd(runlen_argmax, map_cps)
        deallocate(r_old, r_new)
    end subroutine solve_bocpd_ar_1d

    subroutine solve_argpcpd_gp_1d(x, map_cps, runlen_argmax, scale, length_scale, noise_level, max_lag, alpha0, beta0, logistic_hazard_h, logistic_hazard_a, logistic_hazard_b, epsilon, max_train)
        !> Dense GP approximation to ArgpCpd using precomputed lag features and kernel matrix.
        !! This follows the Rust implementation more closely by using one GP factorization
        !! per step for the active suffix and predictive probabilities from NLML differences.
        real(kind=dp), intent(in) :: x(:)
        integer, allocatable, intent(out) :: map_cps(:)
        integer, allocatable, intent(out) :: runlen_argmax(:)
        real(kind=dp), intent(in), optional :: scale, length_scale, noise_level, alpha0, beta0
        real(kind=dp), intent(in), optional :: logistic_hazard_h, logistic_hazard_a, logistic_hazard_b, epsilon
        integer, intent(in), optional :: max_lag, max_train

        integer :: n, lag_order, max_train_, t, mrc, i, m_used
        real(kind=dp) :: scale_, length_scale_, noise_level_, alpha0_, beta0_
        real(kind=dp) :: hazard_h_, hazard_a_, hazard_b_, eps
        real(kind=dp) :: r0_mass, r_sum, pred_prob, hazard
        real(kind=dp), allocatable :: lagmat(:,:), kfull(:,:), nlml_prev(:), nlml_cur(:), pred_probs(:)
        real(kind=dp), allocatable :: r_old(:), r_new(:), r_work(:)

        n = size(x)
        if (n <= 0) then
            allocate(map_cps(0), runlen_argmax(0))
            return
        end if

        scale_ = 3.0_dp
        if (present(scale)) scale_ = scale
        length_scale_ = 10.0_dp
        if (present(length_scale)) length_scale_ = length_scale
        noise_level_ = 0.01_dp
        if (present(noise_level)) noise_level_ = noise_level
        lag_order = 12
        if (present(max_lag)) lag_order = max_lag
        alpha0_ = 2.0_dp
        if (present(alpha0)) alpha0_ = alpha0
        beta0_ = 1.0_dp
        if (present(beta0)) beta0_ = beta0
        hazard_h_ = -5.0_dp
        if (present(logistic_hazard_h)) hazard_h_ = logistic_hazard_h
        hazard_a_ = 1.0_dp
        if (present(logistic_hazard_a)) hazard_a_ = logistic_hazard_a
        hazard_b_ = 1.0_dp
        if (present(logistic_hazard_b)) hazard_b_ = logistic_hazard_b
        eps = 1.0e-10_dp
        if (present(epsilon)) eps = epsilon
        max_train_ = 40
        if (present(max_train)) max_train_ = max_train

        allocate(runlen_argmax(n), lagmat(n, lag_order), kfull(n, n))
        allocate(r_old(n), r_new(n), r_work(n + 1))

        call build_arg_lag_matrix_1d(x, lag_order, lagmat)
        call build_argp_kernel_matrix(lagmat, scale_, length_scale_, noise_level_, kfull)

        r_old = 0.0_dp
        r_old(1) = 1.0_dp
        mrc = 1
        runlen_argmax(1) = 0

        do t = 2, n
            m_used = min(mrc, min(max_train_, t))
            allocate(pred_probs(m_used))
            if (m_used > 1) then
                allocate(nlml_prev(m_used - 1))
                call argpcpd_nlml_vector(x, kfull, t - 1, m_used - 1, alpha0_, beta0_, nlml_prev)
            end if
            allocate(nlml_cur(m_used))
            call argpcpd_nlml_vector(x, kfull, t, m_used, alpha0_, beta0_, nlml_cur)

            pred_probs(1) = exp(-nlml_cur(1))
            do i = 2, m_used
                pred_probs(i) = exp(nlml_prev(i - 1) - nlml_cur(i))
            end do

            r_new(1:mrc + 1) = 0.0_dp
            r0_mass = 0.0_dp
            r_sum = 0.0_dp

            do i = 1, mrc
                if (r_old(i) <= 0.0_dp) cycle
                pred_prob = pred_probs(min(i, m_used))
                hazard = logistic_hazard_value(real(i, dp), hazard_h_, hazard_a_, hazard_b_)
                r_new(i + 1) = r_old(i) * pred_prob * (1.0_dp - hazard)
                r0_mass = r0_mass + r_old(i) * pred_prob * hazard
                r_sum = r_sum + r_new(i + 1)
            end do

            r_new(1) = r0_mass
            r_sum = r_sum + r0_mass
            if (r_sum > 0.0_dp) then
                r_new(1:mrc + 1) = r_new(1:mrc + 1) / r_sum
            else
                r_new(1:mrc + 1) = 0.0_dp
                r_new(1) = 1.0_dp
            end if

            call truncate_run_probs(r_new(1:mrc + 1), eps, r_work, mrc)
            r_old(1:mrc) = r_work(1:mrc)
            runlen_argmax(t) = maxloc(r_old(1:mrc), dim=1) - 1

            if (allocated(nlml_prev)) deallocate(nlml_prev)
            deallocate(nlml_cur, pred_probs)
        end do

        call map_changepoints_bocpd(runlen_argmax, map_cps)
        deallocate(lagmat, kfull, r_old, r_new, r_work)
    end subroutine solve_argpcpd_gp_1d

    subroutine solve_argpcpd_gp_rust_1d(x, map_cps, runlen_argmax, scale, length_scale, noise_level, max_lag, alpha0, beta0, logistic_hazard_h, logistic_hazard_a, logistic_hazard_b, epsilon, max_train)
        !> Closer port of the Rust ArgpCpd step update.
        !! This keeps the incremental GP state across steps, using the same
        !! run-length recursion, Cholesky downdate, and NLML-difference logic.
        !! `max_train`, if provided, caps the active run-length state count; if
        !! omitted, all active states are retained subject only to epsilon truncation.
        real(kind=dp), intent(in) :: x(:)
        integer, allocatable, intent(out) :: map_cps(:)
        integer, allocatable, intent(out) :: runlen_argmax(:)
        real(kind=dp), intent(in), optional :: scale, length_scale, noise_level, alpha0, beta0
        real(kind=dp), intent(in), optional :: logistic_hazard_h, logistic_hazard_a, logistic_hazard_b, epsilon
        integer, intent(in), optional :: max_lag, max_train

        integer :: n, lag_order, t, i, mrc, mrc_keep, u_n, mrc_cap
        real(kind=dp) :: scale_, length_scale_, noise_level_, alpha0_, beta0_
        real(kind=dp) :: hazard_h_, hazard_a_, hazard_b_, eps
        real(kind=dp) :: kss, z, hazard, log_const
        real(kind=dp), allocatable :: lagmat(:,:), kfull(:,:), r_pr(:), next_r(:)
        real(kind=dp), allocatable :: u(:,:), u_new(:,:), alpha(:), alpha_t(:), beta_t(:), last_nlml(:)
        real(kind=dp), allocatable :: kstar(:), kstar_work(:), rev_y(:), csum_alpha2(:), csum_logdiag(:), nlml_cur(:)
        logical :: ok

        n = size(x)
        if (n <= 0) then
            allocate(map_cps(0), runlen_argmax(0))
            return
        end if

        scale_ = 3.0_dp
        if (present(scale)) scale_ = scale
        length_scale_ = 10.0_dp
        if (present(length_scale)) length_scale_ = length_scale
        noise_level_ = 0.01_dp
        if (present(noise_level)) noise_level_ = noise_level
        lag_order = 12
        if (present(max_lag)) lag_order = max_lag
        alpha0_ = 2.0_dp
        if (present(alpha0)) alpha0_ = alpha0
        beta0_ = 1.0_dp
        if (present(beta0)) beta0_ = beta0
        hazard_h_ = -5.0_dp
        if (present(logistic_hazard_h)) hazard_h_ = logistic_hazard_h
        hazard_a_ = 1.0_dp
        if (present(logistic_hazard_a)) hazard_a_ = logistic_hazard_a
        hazard_b_ = 1.0_dp
        if (present(logistic_hazard_b)) hazard_b_ = logistic_hazard_b
        eps = 1.0e-10_dp
        if (present(epsilon)) eps = epsilon
        mrc_cap = n + 1
        if (present(max_train)) mrc_cap = min(max(1, max_train), n + 1)

        log_const = log_gamma(alpha0_) + 0.5_dp * log(2.0_dp * acos(-1.0_dp) * beta0_)

        allocate(runlen_argmax(n), lagmat(n, lag_order), kfull(n, n))
        allocate(r_pr(n + 1), next_r(n + 1))
        allocate(u(n, n), u_new(n, n), alpha(n), alpha_t(n), beta_t(n), last_nlml(n))
        allocate(kstar(n), kstar_work(n), rev_y(n), csum_alpha2(n), csum_logdiag(n), nlml_cur(n))

        call build_arg_lag_matrix_1d(x, lag_order, lagmat)
        call build_argp_kernel_matrix(lagmat, scale_, length_scale_, noise_level_, kfull)

        r_pr = 0.0_dp
        r_pr(1) = 1.0_dp
        next_r = 0.0_dp
        u = 0.0_dp
        if (n >= 1) u(1, 1) = 1.0_dp
        last_nlml = 0.0_dp
        alpha = 0.0_dp
        alpha_t = 0.0_dp
        beta_t = 0.0_dp
        mrc = 1

        do t = 1, n
            if (mrc > mrc_cap) mrc = mrc_cap
            u_n = max(0, mrc - 1)

            kss = sqrt(max(kfull(t, t), 1.0e-12_dp))
            u_new(1:mrc, 1:mrc) = 0.0_dp
            u_new(1, 1) = kss

            if (u_n > 0) then
                do i = 1, u_n
                    kstar(i) = kfull(t - i, t)
                end do
                kstar_work(1:u_n) = kstar(1:u_n) / kss
                u_new(2:mrc, 1) = kstar_work(1:u_n)
                call rank_one_update_lower_real(u(1:u_n, 1:u_n), kstar_work(1:u_n), -1.0_dp, ok)
                if (ok) then
                    u_new(2:mrc, 2:mrc) = u(1:u_n, 1:u_n)
                else
                    call build_argpcpd_reversed_chol(kfull, t, mrc, u_new(1:mrc, 1:mrc), ok)
                    if (.not. ok) then
                        call map_changepoints_bocpd(runlen_argmax(1:max(1, t - 1)), map_cps)
                        return
                    end if
                end if
            end if
            u(1:mrc, 1:mrc) = u_new(1:mrc, 1:mrc)

            do i = 1, mrc
                rev_y(i) = x(t - i + 1)
            end do
            call solve_lower_triangular(u(1:mrc, 1:mrc), rev_y(1:mrc), alpha(1:mrc))

            alpha_t(1) = alpha0_ + 0.5_dp
            beta_t(1) = beta0_ + 0.5_dp * alpha(1) * alpha(1)
            csum_alpha2(1) = alpha(1) * alpha(1)
            csum_logdiag(1) = log(u(1, 1))
            nlml_cur(1) = alpha_t(1) * log(beta_t(1) / beta0_) + csum_logdiag(1) - log_gamma(alpha_t(1)) + log_const
            do i = 2, mrc
                csum_alpha2(i) = csum_alpha2(i - 1) + alpha(i) * alpha(i)
                csum_logdiag(i) = csum_logdiag(i - 1) + log(u(i, i))
                alpha_t(i) = alpha0_ + 0.5_dp * real(i, dp)
                beta_t(i) = beta0_ + 0.5_dp * csum_alpha2(i)
                nlml_cur(i) = alpha_t(i) * log(beta_t(i) / beta0_) + csum_logdiag(i) - log_gamma(alpha_t(i)) + &
                              log_gamma(alpha0_) + 0.5_dp * real(i, dp) * log(2.0_dp * acos(-1.0_dp) * beta0_)
            end do

            next_r(1:mrc + 1) = 0.0_dp
            do i = 1, mrc
                if (i == 1) then
                    z = exp(-nlml_cur(1))
                else
                    z = exp(last_nlml(i - 1) - nlml_cur(i))
                end if
                hazard = logistic_hazard_value(real(i, dp), hazard_h_, hazard_a_, hazard_b_)
                next_r(1) = next_r(1) + r_pr(i) * z * hazard
                next_r(i + 1) = r_pr(i) * z * (1.0_dp - hazard)
            end do

            z = sum(next_r(1:mrc + 1))
            if (z > 0.0_dp) then
                next_r(1:mrc + 1) = next_r(1:mrc + 1) / z
            else
                next_r(1:mrc + 1) = 0.0_dp
                next_r(1) = 1.0_dp
            end if

            call truncate_run_probs(next_r(1:mrc + 1), eps, r_pr, mrc_keep)
            if (mrc_keep > mrc_cap) then
                r_pr(1:mrc_cap) = r_pr(1:mrc_cap) / sum(r_pr(1:mrc_cap))
                mrc_keep = mrc_cap
            end if
            mrc = mrc_keep
            runlen_argmax(t) = maxloc(r_pr(1:mrc), dim=1) - 1

            if (mrc <= 1) then
                last_nlml(1) = 0.0_dp
            else
                last_nlml(1:mrc - 1) = nlml_cur(1:mrc - 1)
            end if
        end do

        call map_changepoints_bocpd(runlen_argmax, map_cps)
        deallocate(lagmat, kfull, r_pr, next_r, u, u_new, alpha, alpha_t, beta_t, last_nlml)
        deallocate(kstar, kstar_work, rev_y, csum_alpha2, csum_logdiag, nlml_cur)
    end subroutine solve_argpcpd_gp_rust_1d

    subroutine map_changepoints_bocpd(runlen_argmax, map_cps)
        !> Reconstruct changepoints from the per-time MAP run lengths.
        integer, intent(in) :: runlen_argmax(:)
        integer, allocatable, intent(out) :: map_cps(:)

        integer :: n, s, rl, ncp, i
        integer, allocatable :: tmp(:)

        n = size(runlen_argmax)
        if (n <= 0) then
            allocate(map_cps(0))
            return
        end if

        allocate(tmp(n))
        ncp = 0
        s = n
        do while (s /= 1)
            rl = runlen_argmax(s)
            if (rl == 0) then
                if (ncp == 0 .or. tmp(ncp) /= s - 1) then
                    ncp = ncp + 1
                    tmp(ncp) = s - 1
                end if
                s = s - 1
            else
                s = max(1, s - rl)
                ncp = ncp + 1
                tmp(ncp) = s - 1
            end if
        end do

        allocate(map_cps(ncp))
        do i = 1, ncp
            map_cps(i) = tmp(ncp - i + 1)
        end do
        deallocate(tmp)
    end subroutine map_changepoints_bocpd

    subroutine reduce_changepoints_normal_1d(x, cps_in, target_ncp, cps_out, min_seg_len, max_iter)
        !> Reduce a changepoint set to target_ncp by repeated greedy merges followed by
        !! local one-breakpoint-at-a-time reoptimization under the univariate normal profile cost.
        real(kind=dp), intent(in) :: x(:)
        integer, intent(in) :: cps_in(:)
        integer, intent(in) :: target_ncp
        integer, allocatable, intent(out) :: cps_out(:)
        integer, intent(in), optional :: min_seg_len, max_iter

        integer :: n, min_len, max_iter_, ncp, i, j, best_j, iter
        integer, allocatable :: work(:)
        real(kind=dp), allocatable :: sz(:), szz(:)
        real(kind=dp) :: best_delta, delta
        logical :: changed

        n = size(x)
        min_len = 2
        if (present(min_seg_len)) min_len = min_seg_len
        max_iter_ = 4
        if (present(max_iter)) max_iter_ = max_iter

        if (size(cps_in) <= target_ncp .or. target_ncp < 1) then
            allocate(cps_out(size(cps_in)))
            cps_out = cps_in
            return
        end if

        allocate(sz(0:n), szz(0:n))
        call build_mean_shift_prefix_1d(x, sz, szz)

        ncp = size(cps_in)
        allocate(work(ncp))
        work = cps_in

        do while (ncp > target_ncp)
            best_delta = huge(1.0_dp)
            best_j = -1
            do j = 2, ncp
                delta = merge_delta_normal_1d(sz, szz, work, j, n, min_len)
                if (delta < best_delta) then
                    best_delta = delta
                    best_j = j
                end if
            end do
            if (best_j < 2) exit

            if (best_j < ncp) work(best_j:ncp-1) = work(best_j+1:ncp)
            ncp = ncp - 1

            do i = max(2, best_j - 1), min(ncp, best_j)
                call refine_one_normal_cp_1d(sz, szz, work, ncp, i, n, min_len)
            end do

            do iter = 1, max_iter_
                changed = .false.
                do i = 2, ncp
                    j = work(i)
                    call refine_one_normal_cp_1d(sz, szz, work, ncp, i, n, min_len)
                    if (work(i) /= j) changed = .true.
                end do
                if (.not. changed) exit
            end do
        end do

        allocate(cps_out(ncp))
        cps_out = work(1:ncp)
        deallocate(sz, szz, work)
    end subroutine reduce_changepoints_normal_1d

    pure subroutine normal_gamma_observe(m, r, s, v, x, m_new, r_new, s_new, v_new)
        real(kind=dp), intent(in) :: m, r, s, v, x
        real(kind=dp), intent(out) :: m_new, r_new, s_new, v_new

        r_new = r + 1.0_dp
        v_new = v + 1.0_dp
        m_new = (r * m + x) / r_new
        s_new = s + x * x + r * m * m - r_new * m_new * m_new
    end subroutine normal_gamma_observe

    pure subroutine student_t_bocd_observe(mu, kappa, alpha, beta, beta0, x, mu_new, kappa_new, alpha_new, beta_new)
        real(kind=dp), intent(in) :: mu, kappa, alpha, beta, beta0, x
        real(kind=dp), intent(out) :: mu_new, kappa_new, alpha_new, beta_new

        mu_new = (kappa * mu + x) / (kappa + 1.0_dp)
        kappa_new = kappa + 1.0_dp
        alpha_new = alpha + 0.5_dp
        beta_new = kappa + (kappa * (x - mu)**2) / (2.0_dp * (kappa + 1.0_dp))
        if (beta_new <= 0.0_dp .or. beta_new >= huge(beta_new)) beta_new = beta0
    end subroutine student_t_bocd_observe

    pure subroutine poisson_gamma_observe(shape, rate, x, shape_new, rate_new)
        real(kind=dp), intent(in) :: shape, rate
        integer, intent(in) :: x
        real(kind=dp), intent(out) :: shape_new, rate_new

        shape_new = shape + real(x, dp)
        rate_new = rate + 1.0_dp
    end subroutine poisson_gamma_observe

    pure subroutine beta_bernoulli_observe(alpha, beta, x, alpha_new, beta_new)
        real(kind=dp), intent(in) :: alpha, beta
        logical, intent(in) :: x
        real(kind=dp), intent(out) :: alpha_new, beta_new

        alpha_new = alpha
        beta_new = beta
        if (x) then
            alpha_new = alpha_new + 1.0_dp
        else
            beta_new = beta_new + 1.0_dp
        end if
    end subroutine beta_bernoulli_observe

    pure function normal_gamma_log_predictive(x, m, r, s, v) result(logp)
        real(kind=dp), intent(in) :: x, m, r, s, v
        real(kind=dp) :: logp
        real(kind=dp) :: rn, mn, sn, half_v, term

        rn = r + 1.0_dp
        mn = (r * m + x) / rn
        sn = rn * mn * mn
        sn = (r * m) * m + x * x + s - sn
        half_v = 0.5_dp * v
        term = 0.5_dp * log(2.0_dp) - 0.5_dp * log(2.0_dp * acos(-1.0_dp)) + &
               0.5_dp * (log(r / rn) + v * log(s)) + &
               (log_gamma(half_v + 0.5_dp) - log_gamma(half_v))
        logp = term - 0.5_dp * (v + 1.0_dp) * log(sn)
    end function normal_gamma_log_predictive

    pure function student_t_bocd_pdf(x, mu, kappa, alpha, beta) result(pdf)
        real(kind=dp), intent(in) :: x, mu, kappa, alpha, beta
        real(kind=dp) :: pdf
        real(kind=dp) :: df, scale, z, logp

        df = 2.0_dp * alpha
        scale = sqrt(beta * (kappa + 1.0_dp) / (alpha * kappa))
        z = (x - mu) / scale
        logp = log_gamma(0.5_dp * (df + 1.0_dp)) - log_gamma(0.5_dp * df) - 0.5_dp * log(df * acos(-1.0_dp)) - &
               log(scale) - 0.5_dp * (df + 1.0_dp) * log(1.0_dp + (z * z) / df)
        pdf = exp(logp)
    end function student_t_bocd_pdf

    pure function poisson_gamma_log_predictive(x, shape, rate) result(logp)
        integer, intent(in) :: x
        real(kind=dp), intent(in) :: shape, rate
        real(kind=dp) :: logp
        real(kind=dp) :: xdp

        xdp = real(x, dp)
        logp = log_gamma(shape + xdp) - log_gamma(shape) - log_gamma(xdp + 1.0_dp) + &
               shape * log(rate / (rate + 1.0_dp)) + xdp * log(1.0_dp / (rate + 1.0_dp))
    end function poisson_gamma_log_predictive

    pure function beta_bernoulli_log_predictive(x, alpha, beta) result(logp)
        logical, intent(in) :: x
        real(kind=dp), intent(in) :: alpha, beta
        real(kind=dp) :: logp

        if (x) then
            logp = log(alpha / (alpha + beta))
        else
            logp = log(beta / (alpha + beta))
        end if
    end function beta_bernoulli_log_predictive

    function ar_segment_log_predictive(x, i1, i2, x_next, max_lag) result(logp)
        !> One-step Gaussian predictive density from an AR fit on x(i1:i2).
        !! Indices i1:i2 are 0-based segment endpoints in BOCPD convention.
        real(kind=dp), intent(in) :: x(:), x_next
        integer, intent(in) :: i1, i2, max_lag
        real(kind=dp) :: logp

        integer :: nseg, order
        real(kind=dp), allocatable :: seg(:), centered(:), ar(:), pred(:)
        real(kind=dp) :: mu, xms, sigma2, resid, xbar

        nseg = i2 - i1 + 1
        if (nseg <= 1) then
            sigma2 = max(1.0_dp, variance_small(x))
            logp = gaussian_logpdf(x_next, x(i2), sigma2)
            return
        end if

        allocate(seg(nseg))
        seg = x(i1 + 1:i2 + 1)
        xbar = sum(seg) / real(nseg, dp)

        if (nseg <= max_lag + 1) then
            sigma2 = max(variance_small(seg), 1.0e-6_dp)
            logp = gaussian_logpdf(x_next, xbar, sigma2)
            deallocate(seg)
            return
        end if

        order = min(max_lag, max(1, nseg / 4))
        allocate(centered(nseg), ar(order), pred(nseg))
        centered = seg - xbar
        call ar_fit(centered, xms, ar)
        pred = ar_pred(centered, ar)
        mu = xbar + pred(nseg)
        resid = centered(nseg) - pred(nseg)
        sigma2 = max(xms, resid * resid, 1.0e-6_dp)
        logp = gaussian_logpdf(x_next, mu, sigma2)

        deallocate(seg, centered, ar, pred)
    end function ar_segment_log_predictive

    pure function logistic_hazard_value(run_length, h, a, b) result(hazard)
        real(kind=dp), intent(in) :: run_length, h, a, b
        real(kind=dp) :: hazard
        hazard = logistic_value(h) * logistic_value(a * run_length + b)
    end function logistic_hazard_value

    pure function logistic_value(x) result(y)
        real(kind=dp), intent(in) :: x
        real(kind=dp) :: y
        y = 1.0_dp / (1.0_dp + exp(-x))
    end function logistic_value

    pure function gaussian_logpdf(x, mu, sigma2) result(logp)
        real(kind=dp), intent(in) :: x, mu, sigma2
        real(kind=dp) :: logp, var_
        var_ = max(sigma2, 1.0e-12_dp)
        logp = -0.5_dp * (log(2.0_dp * acos(-1.0_dp) * var_) + (x - mu) * (x - mu) / var_)
    end function gaussian_logpdf

    pure function variance_small(x) result(v)
        real(kind=dp), intent(in) :: x(:)
        real(kind=dp) :: v, mu
        integer :: n
        n = size(x)
        if (n <= 1) then
            v = 1.0_dp
            return
        end if
        mu = sum(x) / real(n, dp)
        v = sum((x - mu) ** 2) / real(n, dp)
    end function variance_small

    subroutine build_arg_lag_matrix_1d(x, max_lag, lagmat)
        real(kind=dp), intent(in) :: x(:)
        integer, intent(in) :: max_lag
        real(kind=dp), intent(out) :: lagmat(size(x), max_lag)

        integer :: i, j, n
        n = size(x)
        lagmat = 0.0_dp
        do i = 1, n
            do j = 1, max_lag
                if (i > j) lagmat(i, j) = x(i - j)
            end do
        end do
    end subroutine build_arg_lag_matrix_1d

    subroutine build_argp_kernel_matrix(lagmat, scale, length_scale, noise_level, kfull)
        real(kind=dp), intent(in) :: lagmat(:, :), scale, length_scale, noise_level
        real(kind=dp), intent(out) :: kfull(size(lagmat, 1), size(lagmat, 1))

        integer :: n, i, j
        real(kind=dp) :: d2, kval

        n = size(lagmat, 1)
        do i = 1, n
            do j = i, n
                d2 = sum((lagmat(i, :) - lagmat(j, :)) ** 2)
                kval = scale * exp(-0.5_dp * d2 / max(length_scale * length_scale, 1.0e-12_dp))
                if (i == j) kval = kval + noise_level
                kfull(i, j) = kval
                kfull(j, i) = kval
            end do
        end do
    end subroutine build_argp_kernel_matrix

    subroutine argpcpd_nlml_vector(x, kfull, t_end, len_max, alpha0, beta0, nlml)
        !> Compute the vector of negative log marginal likelihoods for suffix lengths
        !! 1:len_max ending at t_end, following the structure in the Rust Argpcpd code.
        real(kind=dp), intent(in) :: x(:), kfull(:, :), alpha0, beta0
        integer, intent(in) :: t_end, len_max
        real(kind=dp), intent(out) :: nlml(len_max)

        real(kind=dp), allocatable :: krev(:, :), lmat(:, :), yrev(:), alpha(:), csum_alpha2(:), csum_logdiag(:)
        logical :: ok
        integer :: m, i, idx_i, idx_j

        if (len_max <= 0) return
        allocate(krev(len_max, len_max), lmat(len_max, len_max), yrev(len_max), alpha(len_max), csum_alpha2(len_max), csum_logdiag(len_max))
        do i = 1, len_max
            idx_i = t_end - i + 1
            yrev(i) = x(idx_i)
            do m = 1, len_max
                idx_j = t_end - m + 1
                krev(i, m) = kfull(idx_i, idx_j)
            end do
        end do

        call cholesky_lower_spd(krev, lmat, ok)
        if (.not. ok) then
            do i = 1, len_max
                nlml(i) = -gaussian_logpdf(yrev(i), 0.0_dp, 1.0_dp)
            end do
            deallocate(krev, lmat, yrev, alpha, csum_alpha2, csum_logdiag)
            return
        end if

        call solve_lower_triangular(lmat, yrev, alpha)
        csum_alpha2(1) = alpha(1) * alpha(1)
        csum_logdiag(1) = log(lmat(1, 1))
        do i = 2, len_max
            csum_alpha2(i) = csum_alpha2(i - 1) + alpha(i) * alpha(i)
            csum_logdiag(i) = csum_logdiag(i - 1) + log(lmat(i, i))
        end do

        do i = 1, len_max
            nlml(i) = (alpha0 + 0.5_dp * real(i, dp)) * log((beta0 + 0.5_dp * csum_alpha2(i)) / beta0) + &
                      log_gamma(alpha0) + csum_logdiag(i) - log_gamma(alpha0 + 0.5_dp * real(i, dp)) + &
                      0.5_dp * real(i, dp) * log(2.0_dp * acos(-1.0_dp) * beta0)
        end do

        deallocate(krev, lmat, yrev, alpha, csum_alpha2, csum_logdiag)
    end subroutine argpcpd_nlml_vector

    subroutine cholesky_lower_spd(a, l, ok)
        real(kind=dp), intent(in) :: a(:, :)
        real(kind=dp), intent(out) :: l(size(a, 1), size(a, 2))
        logical, intent(out) :: ok

        integer :: n, i, j, k
        real(kind=dp) :: s

        n = size(a, 1)
        l = 0.0_dp
        ok = .true.
        do i = 1, n
            do j = 1, i
                s = a(i, j)
                do k = 1, j - 1
                    s = s - l(i, k) * l(j, k)
                end do
                if (i == j) then
                    if (s <= 1.0e-12_dp) then
                        ok = .false.
                        return
                    end if
                    l(i, j) = sqrt(s)
                else
                    l(i, j) = s / l(j, j)
                end if
            end do
        end do
    end subroutine cholesky_lower_spd

    subroutine build_argpcpd_reversed_chol(kfull, t, mrc, lmat, ok)
        real(kind=dp), intent(in) :: kfull(:, :)
        integer, intent(in) :: t, mrc
        real(kind=dp), intent(out) :: lmat(mrc, mrc)
        logical, intent(out) :: ok

        real(kind=dp), allocatable :: krev(:,:)
        integer :: i, j

        allocate(krev(mrc, mrc))
        do i = 1, mrc
            do j = 1, mrc
                krev(i, j) = kfull(t - i + 1, t - j + 1)
            end do
        end do
        call cholesky_lower_spd(krev, lmat, ok)
        deallocate(krev)
    end subroutine build_argpcpd_reversed_chol

    subroutine rank_one_update_lower_real(chol, x, sigma, ok)
        !> Real-valued rank-one update/downdate of a lower-triangular Cholesky factor.
        real(kind=dp), intent(inout) :: chol(:, :)
        real(kind=dp), intent(inout) :: x(:)
        real(kind=dp), intent(in) :: sigma
        logical, intent(out) :: ok

        integer :: n, j
        real(kind=dp) :: beta, diag, diag2, xj, sigma_xj2, gamma, new_diag2, new_diag

        n = size(x)
        ok = .true.
        beta = 1.0_dp

        do j = 1, n
            diag = chol(j, j)
            diag2 = diag * diag
            xj = x(j)
            sigma_xj2 = xj * xj * sigma
            gamma = diag2 * beta + sigma_xj2
            new_diag2 = diag2 + sigma_xj2 / beta
            if (new_diag2 <= 1.0e-12_dp .or. abs(diag) <= 1.0e-12_dp .or. abs(gamma) <= 1.0e-18_dp) then
                ok = .false.
                return
            end if
            new_diag = sqrt(new_diag2)
            chol(j, j) = new_diag

            if (j < n) then
                x(j + 1:n) = x(j + 1:n) - (xj / diag) * chol(j + 1:n, j)
                chol(j + 1:n, j) = (new_diag / diag) * chol(j + 1:n, j) + (new_diag * sigma / gamma) * xj * x(j + 1:n)
            end if
            beta = beta + sigma_xj2 / diag2
        end do
    end subroutine rank_one_update_lower_real

    subroutine solve_lower_triangular(l, b, x)
        real(kind=dp), intent(in) :: l(:, :), b(:)
        real(kind=dp), intent(out) :: x(size(b))
        integer :: n, i

        n = size(b)
        do i = 1, n
            if (i == 1) then
                x(i) = b(i) / l(i, i)
            else
                x(i) = (b(i) - dot_product(l(i, 1:i-1), x(1:i-1))) / l(i, i)
            end if
        end do
    end subroutine solve_lower_triangular

    subroutine truncate_run_probs(r_in, epsilon, r_out, m_out)
        real(kind=dp), intent(in) :: r_in(:), epsilon
        real(kind=dp), intent(out) :: r_out(:)
        integer, intent(out) :: m_out

        integer :: i, last_keep
        real(kind=dp) :: z

        last_keep = size(r_in)
        do i = size(r_in), 1, -1
            if (r_in(i) > epsilon) then
                last_keep = i
                exit
            end if
        end do
        m_out = last_keep
        r_out(1:m_out) = r_in(1:m_out)
        z = sum(r_out(1:m_out))
        if (z > 0.0_dp) r_out(1:m_out) = r_out(1:m_out) / z
    end subroutine truncate_run_probs

    subroutine add_crops_solution(z, sz, szz, penalty, min_len, cps_store, numcps, penalties, seg_costs, n_store, store_cap, idx_out)
        real(kind=dp), intent(in) :: z(:), sz(:), szz(:), penalty
        integer, intent(in) :: min_len
        integer, allocatable, intent(inout) :: cps_store(:,:), numcps(:)
        real(kind=dp), allocatable, intent(inout) :: penalties(:), seg_costs(:)
        integer, intent(inout) :: n_store, store_cap
        integer, intent(out) :: idx_out

        integer :: i, n_found
        integer, allocatable :: bkps(:)

        do i = 1, n_store
            if (abs(penalties(i) - penalty) <= 1.0e-12_dp * max(1.0_dp, abs(penalty))) then
                idx_out = i
                return
            end if
        end do

        call solve_pelt_mean_shift_1d(z, penalty, bkps, n_found, min_seg_len=min_len, jump=1)
        if (n_store + 1 > store_cap) call grow_crops_storage(cps_store, numcps, penalties, seg_costs, store_cap)
        n_store = n_store + 1
        idx_out = n_store
        penalties(idx_out) = penalty
        numcps(idx_out) = n_found
        cps_store(:, idx_out) = 0
        if (n_found > 0) cps_store(1:n_found, idx_out) = bkps(1:n_found)
        seg_costs(idx_out) = segmentation_sse_1d(sz, szz, bkps)
        deallocate(bkps)
    end subroutine add_crops_solution

    subroutine grow_crops_storage(cps_store, numcps, penalties, seg_costs, store_cap)
        integer, allocatable, intent(inout) :: cps_store(:,:), numcps(:)
        real(kind=dp), allocatable, intent(inout) :: penalties(:), seg_costs(:)
        integer, intent(inout) :: store_cap

        integer :: n, new_cap
        integer, allocatable :: cps_tmp(:,:), num_tmp(:)
        real(kind=dp), allocatable :: pen_tmp(:), cost_tmp(:)

        n = size(cps_store, 1)
        new_cap = 2 * store_cap
        allocate(cps_tmp(n, new_cap), num_tmp(new_cap), pen_tmp(new_cap), cost_tmp(new_cap))
        cps_tmp = 0
        num_tmp = 0
        pen_tmp = 0.0_dp
        cost_tmp = 0.0_dp
        cps_tmp(:, 1:store_cap) = cps_store
        num_tmp(1:store_cap) = numcps
        pen_tmp(1:store_cap) = penalties
        cost_tmp(1:store_cap) = seg_costs
        call move_alloc(cps_tmp, cps_store)
        call move_alloc(num_tmp, numcps)
        call move_alloc(pen_tmp, penalties)
        call move_alloc(cost_tmp, seg_costs)
        store_cap = new_cap
    end subroutine grow_crops_storage

    subroutine grow_int_pair_storage(a, b, cap)
        integer, allocatable, intent(inout) :: a(:), b(:)
        integer, intent(inout) :: cap

        integer :: new_cap
        integer, allocatable :: a_tmp(:), b_tmp(:)

        new_cap = 2 * cap
        allocate(a_tmp(new_cap), b_tmp(new_cap))
        a_tmp = 0
        b_tmp = 0
        a_tmp(1:cap) = a
        b_tmp(1:cap) = b
        call move_alloc(a_tmp, a)
        call move_alloc(b_tmp, b)
        cap = new_cap
    end subroutine grow_int_pair_storage

    pure function segmentation_sse_1d(sz, szz, bkps) result(cost)
        real(kind=dp), intent(in) :: sz(:), szz(:)
        integer, intent(in) :: bkps(:)
        real(kind=dp) :: cost

        integer :: start0, i, end0, n

        n = size(sz) - 1
        cost = 0.0_dp
        start0 = 0
        do i = 1, size(bkps)
            end0 = bkps(i)
            cost = cost + mean_shift_segment_sse_half_open_1d(sz, szz, start0, end0)
            start0 = end0
        end do
        cost = cost + mean_shift_segment_sse_half_open_1d(sz, szz, start0, n)
    end function segmentation_sse_1d

    subroutine build_seeded_intervals(n, min_length, max_length, growth_factor, starts, ends)
        integer, intent(in) :: n, min_length, max_length
        real(kind=dp), intent(in) :: growth_factor
        integer, allocatable, intent(out) :: starts(:), ends(:)

        integer :: max_len, n_lengths, i, j, interval_len, step, n_steps, total, pos, last_len
        integer, allocatable :: lens(:), starts_tmp(:), ends_tmp(:)
        real(kind=dp) :: step_factor, frac, log_ratio

        max_len = min(max_length, n)
        if (max_len < min_length) then
            allocate(starts(0), ends(0))
            return
        end if

        if (max_len == min_length) then
            n_lengths = 1
        else
            log_ratio = log(real(max_len, dp) / real(min_length, dp)) / log(growth_factor)
            n_lengths = max(1, ceiling(log_ratio))
        end if

        allocate(lens(n_lengths))
        if (n_lengths == 1) then
            lens(1) = min_length
        else
            do i = 1, n_lengths
                frac = real(i - 1, dp) / real(n_lengths - 1, dp)
                lens(i) = nint(exp(log(real(min_length, dp)) + frac * log(real(max_len, dp) / real(min_length, dp))))
            end do
            lens(1) = min_length
            lens(n_lengths) = max_len
        end if

        total = 0
        last_len = -1
        step_factor = 1.0_dp - 1.0_dp / growth_factor
        do i = 1, n_lengths
            interval_len = lens(i)
            if (interval_len == last_len) cycle
            last_len = interval_len
            step = max(1, nint(step_factor * real(interval_len, dp)))
            n_steps = ceiling(real(n - interval_len, dp) / real(step, dp))
            total = total + n_steps + 1
        end do

        allocate(starts_tmp(total), ends_tmp(total))
        pos = 0
        last_len = -1
        do i = 1, n_lengths
            interval_len = lens(i)
            if (interval_len == last_len) cycle
            last_len = interval_len
            step = max(1, nint(step_factor * real(interval_len, dp)))
            n_steps = ceiling(real(n - interval_len, dp) / real(step, dp))
            do while (n_steps < 0)
                n_steps = 0
            end do
            do j = 0, n_steps
                pos = pos + 1
                starts_tmp(pos) = j * step
                ends_tmp(pos) = min(j * step + interval_len, n)
                if (ends_tmp(pos) - starts_tmp(pos) < min_length) starts_tmp(pos) = n - min_length
            end do
        end do
        allocate(starts(pos), ends(pos))
        starts = starts_tmp(1:pos)
        ends = ends_tmp(1:pos)
        deallocate(lens, starts_tmp, ends_tmp)
    end subroutine build_seeded_intervals

    pure function cusum_weighted_difference_1d(sz, start0, split0, end0) result(val)
        real(kind=dp), intent(in) :: sz(:)
        integer, intent(in) :: start0, split0, end0
        real(kind=dp) :: val
        real(kind=dp) :: n, before_n, after_n, before_sum, after_sum

        n = real(end0 - start0, dp)
        before_n = real(split0 - start0, dp)
        after_n = real(end0 - split0, dp)
        before_sum = sz(split0) - sz(start0)
        after_sum = sz(end0) - sz(split0)
        val = sqrt(after_n / (n * before_n)) * before_sum - sqrt(before_n / (n * after_n)) * after_sum
    end function cusum_weighted_difference_1d

    pure function narrowest_seeded_interval(starts, ends, active) result(best_idx)
        integer, intent(in) :: starts(:), ends(:)
        logical, intent(in) :: active(:)
        integer :: best_idx
        integer :: i, best_width

        best_idx = 1
        best_width = huge(1)
        do i = 1, size(starts)
            if (.not. active(i)) cycle
            if (ends(i) - starts(i) < best_width) then
                best_width = ends(i) - starts(i)
                best_idx = i
            end if
        end do
    end function narrowest_seeded_interval

    pure function merge_delta_normal_1d(sz, szz, cps, j, n, min_len) result(delta)
        real(kind=dp), intent(in) :: sz(:), szz(:)
        integer, intent(in) :: cps(:), j, n, min_len
        real(kind=dp) :: delta
        integer :: left_end, mid_end, right_end
        real(kind=dp), parameter :: huge_cost = 1.0e20_dp

        left_end = cps(j - 1)
        mid_end = cps(j)
        if (j < size(cps)) then
            right_end = cps(j + 1)
        else
            right_end = n
        end if

        delta = normal_segment_cost_1d(sz, szz, left_end + 1, right_end, min_len, huge_cost) - &
                normal_segment_cost_1d(sz, szz, left_end + 1, mid_end, min_len, huge_cost) - &
                normal_segment_cost_1d(sz, szz, mid_end + 1, right_end, min_len, huge_cost)
    end function merge_delta_normal_1d

    subroutine refine_one_normal_cp_1d(sz, szz, cps, ncp, idx, n, min_len)
        real(kind=dp), intent(in) :: sz(:), szz(:)
        integer, intent(inout) :: cps(:)
        integer, intent(in) :: ncp, idx, n, min_len

        integer :: left_end, right_end, k, best_k
        real(kind=dp) :: best_cost, this_cost
        real(kind=dp), parameter :: huge_cost = 1.0e20_dp

        if (idx < 2 .or. idx > ncp) return
        left_end = cps(idx - 1)
        if (idx < ncp) then
            right_end = cps(idx + 1)
        else
            right_end = n
        end if

        best_k = cps(idx)
        best_cost = normal_segment_cost_1d(sz, szz, left_end + 1, best_k, min_len, huge_cost) + &
                    normal_segment_cost_1d(sz, szz, best_k + 1, right_end, min_len, huge_cost)

        do k = left_end + min_len, right_end - min_len
            this_cost = normal_segment_cost_1d(sz, szz, left_end + 1, k, min_len, huge_cost) + &
                        normal_segment_cost_1d(sz, szz, k + 1, right_end, min_len, huge_cost)
            if (this_cost < best_cost) then
                best_cost = this_cost
                best_k = k
            end if
        end do

        cps(idx) = best_k
    end subroutine refine_one_normal_cp_1d

end module changepoint_mod
