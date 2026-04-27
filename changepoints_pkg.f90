module changepoints_pkg_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: solve_wbs_univar_1d, solve_wbs_univar_rob_1d, solve_arc_univar_1d, &
              solve_aarc_univar_1d, solve_online_univar_1d, solve_bs_cov_1d, solve_wbsip_cov_1d, solve_dp_regression_1d, solve_dpdu_regression_1d, solve_cv_dp_regression_1d, solve_local_refine_regression_1d, threshold_wbs_tree_1d, &
              local_refine_wbs_univar_1d, solve_dp_univar_1d, solve_cv_dp_univar_1d, solve_dp_poly_1d, solve_cv_dp_poly_1d, solve_local_refine_poly_1d, solve_dp_var1_1d

contains

    subroutine solve_dp_univar_1d(y, gamma, delta, partition, yhat, cpt_hat)
        !> Dynamic programming for univariate mean change points with l0 penalty, matching changepoints::DP.univar.
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), intent(in) :: gamma
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)
        real(kind=dp), allocatable, intent(out) :: yhat(:)

        integer :: n, r, l, rr, ll, t, ncpt
        real(kind=dp) :: b, dist2
        real(kind=dp), allocatable :: bestvalue(:), cs(:), css(:)
        integer, allocatable :: part_work(:), cps_work(:)
        real(kind=dp) :: seg_sum, seg_ss, seg_mean

        n = size(y)
        allocate(bestvalue(0:n), part_work(0:n), yhat(n), cs(0:n), css(0:n))
        call prefix_sum_1d(y, cs)
        call prefix_sum_sq_1d(y, css)

        bestvalue(0) = -gamma
        part_work(0) = 0
        do r = 1, n
            bestvalue(r) = huge(1.0_dp)
            part_work(r) = 0
            do l = 1, r
                if (r - l > 2 * delta) then
                    seg_sum = cs(r) - cs(l - 1)
                    seg_ss = css(r) - css(l - 1)
                    dist2 = seg_ss - seg_sum * seg_sum / real(r - l + 1, dp)
                    b = bestvalue(l - 1) + gamma + dist2
                else
                    b = huge(1.0_dp)
                end if
                if (b < bestvalue(r)) then
                    bestvalue(r) = b
                    part_work(r) = l - 1
                end if
            end do
        end do

        rr = n
        ll = part_work(rr)
        yhat = 0.0_dp
        do while (rr > 0)
            seg_mean = sum(y(ll + 1:rr)) / real(rr - ll, dp)
            do t = ll + 1, rr
                yhat(t) = seg_mean
            end do
            rr = ll
            ll = part_work(rr)
        end do

        allocate(partition(n))
        partition = part_work(1:n)

        allocate(cps_work(n))
        ncpt = 0
        rr = n
        ll = part_work(rr)
        do while (rr > 0)
            ncpt = ncpt + 1
            cps_work(ncpt) = ll
            rr = ll
            if (rr > 0) ll = part_work(rr)
        end do
        if (ncpt <= 1) then
            allocate(cpt_hat(0))
        else
            allocate(cpt_hat(ncpt - 1))
            do t = 1, ncpt - 1
                cpt_hat(t) = cps_work(ncpt - t)
            end do
        end if

        deallocate(bestvalue, part_work, cs, css, cps_work)
    end subroutine solve_dp_univar_1d

    subroutine solve_dp_poly_1d(y, r, gamma, delta, partition, yhat, cpt_hat)
        !> Dynamic programming for univariate polynomial changepoints, matching changepoints::DP.poly.
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: r, delta
        real(kind=dp), intent(in) :: gamma
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)
        real(kind=dp), allocatable, intent(out) :: yhat(:)

        integer :: n, i, l, rr, ll, t, ncpt
        real(kind=dp) :: b, dist2
        real(kind=dp), allocatable :: bestvalue(:), seg_fit(:)
        integer, allocatable :: part_work(:), cps_work(:)

        n = size(y)
        allocate(bestvalue(0:n), part_work(0:n), yhat(n), cps_work(n))
        bestvalue(0) = -gamma
        part_work(0) = 0

        do i = 1, n
            bestvalue(i) = huge(1.0_dp)
            part_work(i) = 0
            do l = 1, i
                if (i - l > 2 * delta) then
                    call polynomial_segment_fit_1d(y, n, l, i, r, seg_fit, dist2)
                    b = bestvalue(l - 1) + gamma + dist2
                    deallocate(seg_fit)
                else
                    b = huge(1.0_dp)
                end if
                if (b < bestvalue(i)) then
                    bestvalue(i) = b
                    part_work(i) = l - 1
                end if
            end do
        end do

        rr = n
        ll = part_work(rr)
        yhat = 0.0_dp
        do while (rr > 0)
            call polynomial_segment_fit_1d(y, n, ll + 1, rr, r, seg_fit, dist2)
            yhat(ll + 1:rr) = seg_fit
            deallocate(seg_fit)
            rr = ll
            ll = part_work(rr)
        end do

        allocate(partition(n))
        partition = part_work(1:n)

        ncpt = 0
        rr = n
        ll = part_work(rr)
        do while (rr > 0)
            ncpt = ncpt + 1
            cps_work(ncpt) = ll
            rr = ll
            if (rr > 0) ll = part_work(rr)
        end do
        if (ncpt <= 1) then
            allocate(cpt_hat(0))
        else
            allocate(cpt_hat(ncpt - 1))
            do t = 1, ncpt - 1
                cpt_hat(t) = cps_work(ncpt - t)
            end do
        end if

        deallocate(bestvalue, part_work, cps_work)
    end subroutine solve_dp_poly_1d

    subroutine solve_cv_dp_poly_1d(y, r, gamma, delta, cpt_hat, k_hat, test_error, train_error)
        !> Sample-splitting cross-validation wrapper matching changepoints::CV.DP.poly.
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: r, delta
        real(kind=dp), intent(in) :: gamma
        integer, allocatable, intent(out) :: cpt_hat(:)
        integer, intent(out) :: k_hat
        real(kind=dp), intent(out) :: test_error, train_error

        integer :: n, n_train
        real(kind=dp), allocatable :: train_y(:), validation_y(:), yhat_train(:)
        integer, allocatable :: partition(:), cpt_train(:), init_cpt_train_long(:), diff_point(:)

        n = size(y)
        n_train = (n + 1) / 2
        allocate(train_y(n_train), validation_y(n / 2))
        train_y = y(1:n:2)
        validation_y = y(2:n:2)

        call solve_dp_poly_1d(train_y, r, gamma, delta, partition, yhat_train, cpt_train)

        allocate(init_cpt_train_long(size(cpt_train) + 2), diff_point(size(cpt_train) + 1))
        init_cpt_train_long(1) = 0
        if (size(cpt_train) > 0) init_cpt_train_long(2:size(cpt_train) + 1) = cpt_train
        init_cpt_train_long(size(init_cpt_train_long)) = size(train_y)
        diff_point = init_cpt_train_long(2:) - init_cpt_train_long(:size(init_cpt_train_long) - 1)

        if (any(diff_point == 1)) then
            k_hat = size(cpt_train)
            allocate(cpt_hat(k_hat))
            if (k_hat > 0) cpt_hat = 2 * cpt_train - 1
            test_error = huge(1.0_dp)
            train_error = huge(1.0_dp)
        else
            k_hat = size(cpt_train)
            allocate(cpt_hat(k_hat))
            if (k_hat > 0) cpt_hat = 2 * cpt_train - 1
            train_error = sqrt(sum((train_y - yhat_train) ** 2))
            test_error = sqrt(sum((validation_y - yhat_train(1:size(validation_y))) ** 2))
        end if

        deallocate(train_y, validation_y, yhat_train, partition, cpt_train, init_cpt_train_long, diff_point)
    end subroutine solve_cv_dp_poly_1d

    subroutine solve_local_refine_poly_1d(cpt_init, y, r, delta_lr, cpt_refined)
        !> Local refinement for univariate polynomial changepoints, matching changepoints::local.refine.poly.
        integer, intent(in) :: cpt_init(:), r, delta_lr
        real(kind=dp), intent(in) :: y(:)
        integer, allocatable, intent(out) :: cpt_refined(:)

        real(kind=dp), parameter :: w = 0.9_dp
        integer :: n, ncp, k, lower, upper, eta, best_eta
        integer, allocatable :: cpt_ext(:)
        real(kind=dp) :: s, e, best_obj, obj

        n = size(y)
        ncp = size(cpt_init)
        allocate(cpt_refined(ncp))
        if (ncp == 0) return

        allocate(cpt_ext(ncp + 2))
        cpt_ext(1) = 0
        cpt_ext(2:ncp + 1) = cpt_init
        cpt_ext(ncp + 2) = n

        do k = 1, ncp
            s = w * real(cpt_ext(k), dp) + (1.0_dp - w) * real(cpt_ext(k + 1), dp)
            e = (1.0_dp - w) * real(cpt_ext(k + 1), dp) + w * real(cpt_ext(k + 2), dp)
            lower = ceiling(s) + 1
            upper = floor(e) - 1
            best_obj = huge(1.0_dp)
            best_eta = lower
            do eta = lower, upper
                obj = obj_func_lr_poly_1d(y, s, e, eta, r, delta_lr)
                if (obj < best_obj) then
                    best_obj = obj
                    best_eta = eta
                end if
            end do
            cpt_refined(k) = ceiling(s) + (best_eta - lower + 1)
        end do

        deallocate(cpt_ext)
    end subroutine solve_local_refine_poly_1d

    subroutine solve_dp_var1_1d(x_futu, x_curr, gamma, lambda, delta, partition, cpt_hat)
        !> Dynamic programming for multivariate VAR(1) changepoints, matching changepoints::DP.VAR1 for scalar lambda = 0.
        real(kind=dp), intent(in) :: x_futu(:, :), x_curr(:, :)
        real(kind=dp), intent(in) :: gamma, lambda
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)

        integer :: n, i, l, rr, ll, ncpt, t
        real(kind=dp) :: b, mse
        real(kind=dp), allocatable :: bestvalue(:)
        integer, allocatable :: part_work(:), cps_work(:)

        if (abs(lambda) > 1.0e-12_dp) error stop "solve_dp_var1_1d currently supports lambda = 0 only"
        n = size(x_futu, 2)
        allocate(bestvalue(0:n), part_work(0:n), cps_work(n))
        bestvalue(0) = -gamma
        part_work(0) = 0

        do i = 1, n
            bestvalue(i) = huge(1.0_dp)
            part_work(i) = 0
            do l = 1, i
                call error_pred_seg_var1_1d(x_futu, x_curr, l, i, lambda, delta, mse)
                b = bestvalue(l - 1) + gamma + mse
                if (b < bestvalue(i)) then
                    bestvalue(i) = b
                    part_work(i) = l - 1
                end if
            end do
        end do

        allocate(partition(n))
        partition = part_work(1:n)

        ncpt = 0
        rr = n
        ll = part_work(rr)
        do while (rr > 0)
            ncpt = ncpt + 1
            cps_work(ncpt) = ll
            rr = ll
            if (rr > 0) ll = part_work(rr)
        end do
        if (ncpt <= 1) then
            allocate(cpt_hat(0))
        else
            allocate(cpt_hat(ncpt - 1))
            do t = 1, ncpt - 1
                cpt_hat(t) = cps_work(ncpt - t)
            end do
        end if

        deallocate(bestvalue, part_work, cps_work)
    end subroutine solve_dp_var1_1d

    subroutine solve_cv_dp_univar_1d(y, gamma, delta, cpt_hat, k_hat, test_error, train_error)
        !> Sample-splitting cross-validation wrapper matching changepoints::CV.DP.univar for even-length series.
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), intent(in) :: gamma
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: cpt_hat(:)
        integer, intent(out) :: k_hat
        real(kind=dp), intent(out) :: test_error, train_error

        integer :: n, n_half, seg_start, seg_end, j
        real(kind=dp) :: seg_mean
        real(kind=dp), allocatable :: train_y(:), validation_y(:), yhat_train(:)
        integer, allocatable :: partition(:), cpt_train(:)

        n = size(y)
        if (mod(n, 2) /= 0) error stop "solve_cv_dp_univar_1d requires an even-length series"

        n_half = n / 2
        allocate(train_y(n_half), validation_y(n_half))
        train_y = y(1:n:2)
        validation_y = y(2:n:2)

        call solve_dp_univar_1d(train_y, gamma, delta, partition, yhat_train, cpt_train)

        k_hat = size(cpt_train)
        allocate(cpt_hat(k_hat))
        if (k_hat > 0) cpt_hat = 2 * cpt_train - 1

        train_error = 0.0_dp
        test_error = 0.0_dp
        seg_start = 1
        do j = 1, k_hat + 1
            if (j <= k_hat) then
                seg_end = cpt_train(j)
            else
                seg_end = n_half
            end if
            seg_mean = sum(train_y(seg_start:seg_end)) / real(seg_end - seg_start + 1, dp)
            train_error = train_error + sum((train_y(seg_start:seg_end) - seg_mean) ** 2)
            test_error = test_error + sum((validation_y(seg_start:seg_end) - seg_mean) ** 2)
            seg_start = seg_end + 1
        end do

        deallocate(train_y, validation_y, partition, yhat_train, cpt_train)
    end subroutine solve_cv_dp_univar_1d

    subroutine solve_wbs_univar_1d(y, alpha, beta, s, e, delta, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        !> Wild binary segmentation for univariate mean changes, matching changepoints::WBS.univar.
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: alpha(:), beta(:)
        integer, intent(in) :: s, e, delta
        integer, intent(out) :: n_nodes
        integer, allocatable, intent(out) :: s_out(:), level_out(:), parent_out(:, :), parent_idx_out(:)
        real(kind=dp), allocatable, intent(out) :: dval_out(:)

        integer :: nmax
        real(kind=dp), allocatable :: cs(:), work_s(:), work_dval(:)
        integer, allocatable :: work_level(:), work_parent(:, :), work_parent_idx(:)

        nmax = max(size(y), 1)
        allocate(cs(0:size(y)))
        call prefix_sum_1d(y, cs)

        allocate(work_s(nmax), work_dval(nmax), work_level(nmax), work_parent(2, nmax), work_parent_idx(nmax))
        n_nodes = 0
        call wbs_univar_rec(y, cs, alpha, beta, s, e, delta, 0, 0, n_nodes, work_s, work_dval, work_level, work_parent, work_parent_idx)

        allocate(s_out(n_nodes), dval_out(n_nodes), level_out(n_nodes), parent_out(2, n_nodes), parent_idx_out(n_nodes))
        if (n_nodes > 0) then
            s_out = int(work_s(1:n_nodes))
            dval_out = work_dval(1:n_nodes)
            level_out = work_level(1:n_nodes)
            parent_out = work_parent(:, 1:n_nodes)
            parent_idx_out = work_parent_idx(1:n_nodes)
        end if

        deallocate(cs, work_s, work_dval, work_level, work_parent, work_parent_idx)
    end subroutine solve_wbs_univar_1d

    subroutine solve_wbs_univar_rob_1d(y, alpha, beta, s, e, k_huber, delta, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        !> Robust wild binary segmentation for univariate mean changes, matching changepoints::WBS.uni.rob.
        real(kind=dp), intent(in) :: y(:), k_huber
        integer, intent(in) :: alpha(:), beta(:)
        integer, intent(in) :: s, e, delta
        integer, intent(out) :: n_nodes
        integer, allocatable, intent(out) :: s_out(:), level_out(:), parent_out(:, :), parent_idx_out(:)
        real(kind=dp), allocatable, intent(out) :: dval_out(:)

        integer :: nmax
        real(kind=dp), allocatable :: work_s(:), work_dval(:)
        integer, allocatable :: work_level(:), work_parent(:, :), work_parent_idx(:)

        nmax = max(size(y), 1)
        allocate(work_s(nmax), work_dval(nmax), work_level(nmax), work_parent(2, nmax), work_parent_idx(nmax))
        n_nodes = 0
        call wbs_univar_rob_rec(y, alpha, beta, s, e, k_huber, delta, 0, 0, n_nodes, work_s, work_dval, work_level, work_parent, work_parent_idx)

        allocate(s_out(n_nodes), dval_out(n_nodes), level_out(n_nodes), parent_out(2, n_nodes), parent_idx_out(n_nodes))
        if (n_nodes > 0) then
            s_out = int(work_s(1:n_nodes))
            dval_out = work_dval(1:n_nodes)
            level_out = work_level(1:n_nodes)
            parent_out = work_parent(:, 1:n_nodes)
            parent_idx_out = work_parent_idx(1:n_nodes)
        end if

        deallocate(work_s, work_dval, work_level, work_parent, work_parent_idx)
    end subroutine solve_wbs_univar_rob_1d

    subroutine solve_arc_univar_1d(y, h, block_num, epsilon, gaussian, cpt_hat, lambda, cusum)
        !> Adversarially robust mean change detection matching changepoints::ARC.
        real(kind=dp), intent(in) :: y(:), epsilon
        integer, intent(in) :: h, block_num
        logical, intent(in) :: gaussian
        integer, allocatable, intent(out) :: cpt_hat(:)
        real(kind=dp), intent(out) :: lambda
        real(kind=dp), allocatable, intent(out) :: cusum(:)

        integer :: n, s_h, i, j, n_cusum, k_trim, n_keep
        real(kind=dp) :: sp, local_max
        integer, allocatable :: cpt_work(:)

        n = size(y)
        sp = (1.0_dp - 2.0_dp * epsilon - 0.1_dp) / 2.0_dp
        k_trim = floor(sp * real(h, dp))
        s_h = h * block_num
        n_cusum = n - 2 * h
        if (n_cusum <= 0) then
            allocate(cusum(0), cpt_hat(0))
            lambda = 0.0_dp
            return
        end if

        allocate(cusum(n_cusum))
        do i = h + 1, n - h
            cusum(i - h) = abs(rume_1d(y(i - h:i), k_trim) - rume_1d(y(i:i + h), k_trim))
        end do

        if (gaussian) then
            lambda = max(8.0_dp * epsilon, 0.6_dp * sqrt(40.0_dp * log(real(n, dp)) / real(h, dp)))
        else
            lambda = max(8.0_dp * sqrt(epsilon), 0.6_dp * sqrt(40.0_dp * log(real(n, dp)) / real(h, dp)))
        end if

        if (n_cusum - 2 * s_h <= 0) then
            allocate(cpt_hat(0))
            return
        end if

        allocate(cpt_work(n_cusum))
        n_keep = 0
        do j = s_h + 1, n_cusum - s_h
            local_max = maxval(cusum(j - s_h:j + s_h))
            if (cusum(j) == local_max .and. cusum(j) > lambda) then
                n_keep = n_keep + 1
                cpt_work(n_keep) = j + h
            end if
        end do

        allocate(cpt_hat(n_keep))
        if (n_keep > 0) cpt_hat = cpt_work(1:n_keep)
        deallocate(cpt_work)
    end subroutine solve_arc_univar_1d

    subroutine solve_aarc_univar_1d(y, t_dat, guess_true, h, block_num, cpt_hat, eps_hat, arc_epsilon, lambda)
        !> Automatic adversarially robust mean change detection matching changepoints::aARC.
        real(kind=dp), intent(in) :: y(:), t_dat(:), guess_true
        integer, intent(in) :: h, block_num
        integer, allocatable, intent(out) :: cpt_hat(:)
        real(kind=dp), intent(out) :: eps_hat, arc_epsilon, lambda

        integer, parameter :: n_grid = 201
        real(kind=dp), allocatable :: epsilon_grid(:), est_e(:), ho(:), roots(:), cusum(:)
        integer :: e1, e2, ind
        real(kind=dp) :: a, b, cc

        allocate(epsilon_grid(n_grid), est_e(n_grid), ho(n_grid), roots(2))
        do e1 = 1, n_grid
            epsilon_grid(e1) = 0.25_dp + 0.001_dp * real(e1 - 1, dp)
            est_e(e1) = rume_1d(t_dat, floor(epsilon_grid(e1) * real(size(t_dat), dp)))
        end do

        ho = 0.0_dp
        do e1 = 1, n_grid
            do e2 = 1, n_grid
                a = min(est_e(e1), est_e(e2))
                b = max(est_e(e1), est_e(e2))
                if (a == est_e(e1)) then
                    ho(e1) = ho(e1) + real(rtest_1d(t_dat, a, b), dp)
                else if (b == est_e(e1)) then
                    ho(e1) = ho(e1) + abs(1.0_dp - real(rtest_1d(t_dat, a, b), dp))
                end if
            end do
        end do

        ind = minloc(ho, dim=1)
        eps_hat = epsilon_grid(ind)
        cc = 0.96_dp - 2.0_dp * eps_hat
        roots = ((-0.4_dp + [1.0_dp, -1.0_dp] * sqrt(0.4_dp**2 + 8.0_dp * cc)) / 4.0_dp) ** 2
        if (abs(roots(1) - guess_true) <= abs(roots(2) - guess_true)) then
            arc_epsilon = roots(1)
        else
            arc_epsilon = roots(2)
        end if

        call solve_arc_univar_1d(y, h, block_num, arc_epsilon, .true., cpt_hat, lambda, cusum)
        deallocate(epsilon_grid, est_e, ho, roots, cusum)
    end subroutine solve_aarc_univar_1d

    subroutine solve_online_univar_1d(y, b_vec, cpt_hat)
        !> Online univariate mean change detection for a supplied threshold vector, matching changepoints::online.univar.
        real(kind=dp), intent(in) :: y(:), b_vec(:)
        integer, intent(out) :: cpt_hat

        integer :: n, t, s
        real(kind=dp) :: cusum, max_cusum, left_mean, right_mean
        logical :: flag

        n = size(y)
        if (size(b_vec) /= n - 1) error stop "solve_online_univar_1d requires size(b_vec) = size(y) - 1"

        t = 1
        flag = .false.
        do while ((.not. flag) .and. t < n)
            t = t + 1
            max_cusum = -1.0_dp
            do s = 1, t - 1
                left_mean = sum(y(1:s)) / real(s, dp)
                right_mean = sum(y(s + 1:t)) / real(t - s, dp)
                cusum = sqrt(real((t - s) * s, dp) / real(t, dp)) * abs(left_mean - right_mean)
                if (cusum > max_cusum) max_cusum = cusum
            end do
            flag = max_cusum > b_vec(t - 1)
        end do
        cpt_hat = t
    end subroutine solve_online_univar_1d

    subroutine solve_bs_cov_1d(x, s, e, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        !> Binary segmentation for covariance changes through the covariance CUSUM operator norm, matching changepoints::BS.cov.
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: s, e
        integer, intent(out) :: n_nodes
        integer, allocatable, intent(out) :: s_out(:), level_out(:), parent_out(:, :), parent_idx_out(:)
        real(kind=dp), allocatable, intent(out) :: dval_out(:)

        integer :: nmax
        integer, allocatable :: work_s(:), work_level(:), work_parent(:, :), work_parent_idx(:)
        real(kind=dp), allocatable :: work_dval(:)

        nmax = max(size(x, 2), 1)
        allocate(work_s(nmax), work_level(nmax), work_parent(2, nmax), work_parent_idx(nmax), work_dval(nmax))
        n_nodes = 0
        call bs_cov_rec(x, s, e, 0, 0, n_nodes, work_s, work_dval, work_level, work_parent, work_parent_idx)

        allocate(s_out(n_nodes), dval_out(n_nodes), level_out(n_nodes), parent_out(2, n_nodes), parent_idx_out(n_nodes))
        if (n_nodes > 0) then
            s_out = work_s(1:n_nodes)
            dval_out = work_dval(1:n_nodes)
            level_out = work_level(1:n_nodes)
            parent_out = work_parent(:, 1:n_nodes)
            parent_idx_out = work_parent_idx(1:n_nodes)
        end if

        deallocate(work_s, work_dval, work_level, work_parent, work_parent_idx)
    end subroutine solve_bs_cov_1d

    subroutine solve_wbsip_cov_1d(x, x_prime, alpha, beta, s, e, delta, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        !> Wild binary segmentation for covariance changes through independent projection, matching changepoints::WBSIP.cov.
        real(kind=dp), intent(in) :: x(:, :), x_prime(:, :)
        integer, intent(in) :: alpha(:), beta(:)
        integer, intent(in) :: s, e, delta
        integer, intent(out) :: n_nodes
        integer, allocatable, intent(out) :: s_out(:), level_out(:), parent_out(:, :), parent_idx_out(:)
        real(kind=dp), allocatable, intent(out) :: dval_out(:)

        integer :: nmax
        integer, allocatable :: work_s(:), work_level(:), work_parent(:, :), work_parent_idx(:)
        real(kind=dp), allocatable :: work_dval(:)

        nmax = max(size(x, 2), 1)
        allocate(work_s(nmax), work_level(nmax), work_parent(2, nmax), work_parent_idx(nmax), work_dval(nmax))
        n_nodes = 0
        call wbsip_cov_rec(x, x_prime, alpha, beta, s, e, delta, 0, 0, n_nodes, work_s, work_dval, work_level, work_parent, work_parent_idx)

        allocate(s_out(n_nodes), dval_out(n_nodes), level_out(n_nodes), parent_out(2, n_nodes), parent_idx_out(n_nodes))
        if (n_nodes > 0) then
            s_out = work_s(1:n_nodes)
            dval_out = work_dval(1:n_nodes)
            level_out = work_level(1:n_nodes)
            parent_out = work_parent(:, 1:n_nodes)
            parent_idx_out = work_parent_idx(1:n_nodes)
        end if

        deallocate(work_s, work_dval, work_level, work_parent, work_parent_idx)
    end subroutine solve_wbsip_cov_1d

    subroutine solve_dp_regression_1d(y, x, gamma, lambda, delta, eps, partition, cpt_hat)
        !> Dynamic programming for regression changepoints with lasso segment fits, matching changepoints::DP.regression for scalar lambda.
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: gamma, lambda, eps
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)

        integer :: n, p, i, l, rr, ll, ncpt, t
        real(kind=dp) :: b, penalty_scale, mse
        real(kind=dp), allocatable :: bestvalue(:)
        integer, allocatable :: part_work(:), cps_work(:)

        n = size(y)
        p = size(x, 2)
        penalty_scale = gamma * log(real(max(n, p), dp))
        allocate(bestvalue(0:n), part_work(0:n))
        bestvalue(0) = -penalty_scale
        part_work(0) = 0

        do i = 1, n
            bestvalue(i) = huge(1.0_dp)
            part_work(i) = 0
            do l = 1, i
                call error_pred_seg_regression_1d(y, x, l, i, lambda, delta, eps, mse)
                b = bestvalue(l - 1) + penalty_scale + mse
                if (b < bestvalue(i)) then
                    bestvalue(i) = b
                    part_work(i) = l - 1
                end if
            end do
        end do

        allocate(partition(n))
        partition = part_work(1:n)

        allocate(cps_work(n))
        ncpt = 0
        rr = n
        ll = part_work(rr)
        do while (rr > 0)
            ncpt = ncpt + 1
            cps_work(ncpt) = ll
            rr = ll
            if (rr > 0) ll = part_work(rr)
        end do
        if (ncpt <= 1) then
            allocate(cpt_hat(0))
        else
            allocate(cpt_hat(ncpt - 1))
            do t = 1, ncpt - 1
                cpt_hat(t) = cps_work(ncpt - t)
            end do
        end if

        deallocate(bestvalue, part_work, cps_work)
    end subroutine solve_dp_regression_1d

    subroutine solve_dpdu_regression_1d(y, x, lambda, zeta, eps, partition, cpt_hat, beta_mat)
        !> Dynamic programming with dynamic updates for regression changepoints, matching changepoints::DPDU.regression.
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: lambda, eps
        integer, intent(in) :: zeta
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)
        real(kind=dp), allocatable, intent(out) :: beta_mat(:, :)

        integer :: n, p, i, l, n_seg, rr, ll, ncpt, t, n_keep
        real(kind=dp) :: b, ymean_seg, lambda_scaled, loss
        real(kind=dp), allocatable :: bestvalue(:), ymean_new(:), ymean_old(:), weights(:), vtilde(:), beta_start(:), beta_std(:), beta_hat_seg(:)
        real(kind=dp), allocatable :: xmeans_new(:, :), xmeans_old(:, :), v_new(:, :), v_old(:, :)
        real(kind=dp), allocatable :: m_new(:, :, :), m_old(:, :, :), mcentered(:, :), mtilde(:, :)
        real(kind=dp), allocatable :: beta_temp(:, :), beta_best(:, :), cps_work(:)
        integer, allocatable :: part_work(:), cpt_all(:)

        n = size(y)
        p = size(x, 2)
        allocate(bestvalue(0:n), part_work(0:n), ymean_new(n), ymean_old(n), xmeans_new(n, p), xmeans_old(n, p), &
                 v_new(n, p), v_old(n, p), m_new(p, p, n), m_old(p, p, n), mcentered(p, p), mtilde(p, p), weights(p), &
                 vtilde(p), beta_start(p), beta_std(p), beta_hat_seg(p), beta_temp(p + 1, n), beta_best(p + 1, n), cps_work(n))

        bestvalue = 0.0_dp
        part_work = 0
        bestvalue(0) = -real(zeta, dp)
        ymean_new = 0.0_dp
        ymean_old = 0.0_dp
        xmeans_new = 0.0_dp
        xmeans_old = 0.0_dp
        v_new = 0.0_dp
        v_old = 0.0_dp
        m_new = 0.0_dp
        m_old = 0.0_dp
        beta_temp = 0.0_dp
        beta_best = 0.0_dp
        beta_start = 0.0_dp

        do i = 1, n
            bestvalue(i) = huge(1.0_dp)
            do l = 1, i
                n_seg = i - l + 1
                ymean_new(l) = (ymean_old(l) * real(n_seg - 1, dp) + y(i)) / real(n_seg, dp)
                xmeans_new(l, :) = (xmeans_old(l, :) * real(n_seg - 1, dp) + x(i, :)) / real(n_seg, dp)
                m_new(:, :, l) = m_old(:, :, l) + outer_product_row_1d(x(i, :), x(i, :))
                v_new(l, :) = v_old(l, :) + y(i) * x(i, :)
                mcentered = m_new(:, :, l) - real(n_seg, dp) * outer_product_row_1d(xmeans_new(l, :), xmeans_new(l, :))
                do t = 1, p
                    weights(t) = sqrt(max(mcentered(t, t) / real(n_seg, dp), 0.0_dp))
                    if (weights(t) <= 1.0e-12_dp) weights(t) = 1.0_dp
                end do
                mtilde = scale_symmetric_matrix_1d(mcentered, weights)
                vtilde = (v_new(l, :) - real(n_seg, dp) * ymean_new(l) * xmeans_new(l, :)) / weights
                if (n_seg >= zeta) then
                    lambda_scaled = lambda * sqrt(max(log(real(max(n, p), dp)), real(n_seg - 1, dp))) * &
                                    sqrt(log(real(max(n, p), dp))) / real(n_seg, dp)
                    call lasso_dpdu_standardized_1d(mtilde, vtilde, beta_start, n_seg, lambda_scaled, eps, beta_std)
                    beta_start = beta_std
                    loss = dot_product(beta_std, matmul(mtilde, beta_std)) - 2.0_dp * dot_product(vtilde, beta_std)
                    beta_hat_seg = beta_std / weights
                    ymean_seg = ymean_new(l) - dot_product(xmeans_new(l, :), beta_hat_seg)
                    beta_temp(1, l) = ymean_seg
                    beta_temp(2:p + 1, l) = beta_hat_seg
                    b = bestvalue(l - 1) + real(zeta, dp) + loss
                else
                    b = bestvalue(l - 1) + real(zeta, dp)
                end if
                if (b < bestvalue(i)) then
                    bestvalue(i) = b
                    part_work(i) = l - 1
                    beta_best(:, i) = beta_temp(:, l)
                end if
            end do
            m_old = m_new
            v_old = v_new
            ymean_old = ymean_new
            xmeans_old = xmeans_new
        end do

        allocate(partition(n))
        partition = part_work(1:n)
        allocate(beta_mat(p + 1, n))
        beta_mat = 0.0_dp
        rr = n
        ll = part_work(rr)
        do while (rr > 0)
            do t = ll + 1, rr
                beta_mat(:, t) = beta_best(:, rr)
            end do
            rr = ll
            if (rr > 0) ll = part_work(rr)
        end do

        ncpt = 0
        rr = n
        ll = part_work(rr)
        do while (rr > 0)
            ncpt = ncpt + 1
            cps_work(ncpt) = real(ll, dp)
            rr = ll
            if (rr > 0) ll = part_work(rr)
        end do
        if (ncpt <= 1) then
            allocate(cpt_all(0))
        else
            allocate(cpt_all(ncpt - 1))
            do t = 1, ncpt - 1
                cpt_all(t) = int(cps_work(ncpt - t))
            end do
        end if

        if (size(cpt_all) > 0) then
            n_keep = count(cpt_all >= zeta .and. (n - cpt_all) >= zeta)
            allocate(cpt_hat(n_keep))
            if (n_keep > 0) cpt_hat = pack(cpt_all, cpt_all >= zeta .and. (n - cpt_all) >= zeta)
        else
            allocate(cpt_hat(0))
        end if

        deallocate(bestvalue, part_work, ymean_new, ymean_old, xmeans_new, xmeans_old, v_new, v_old, m_new, m_old, &
                   mcentered, mtilde, weights, vtilde, beta_start, beta_std, beta_hat_seg, beta_temp, beta_best, cps_work, cpt_all)
    end subroutine solve_dpdu_regression_1d

    subroutine solve_local_refine_regression_1d(cpt_init, y, x, zeta, cpt_refined)
        !> Local refinement for regression changepoints, matching changepoints::local.refine.regression.
        integer, intent(in) :: cpt_init(:)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: zeta
        integer, allocatable, intent(out) :: cpt_refined(:)

        real(kind=dp), parameter :: w = 0.9_dp
        integer :: n, ncp, k, lower, upper, eta, best_eta
        integer, allocatable :: cpt_ext(:)
        real(kind=dp) :: s, e, best_obj, obj

        n = size(y)
        ncp = size(cpt_init)
        allocate(cpt_refined(ncp))
        if (ncp == 0) return

        allocate(cpt_ext(ncp + 2))
        cpt_ext = [0, cpt_init, n]
        do k = 1, ncp
            s = w * real(cpt_ext(k), dp) + (1.0_dp - w) * real(cpt_ext(k + 1), dp)
            e = (1.0_dp - w) * real(cpt_ext(k + 1), dp) + w * real(cpt_ext(k + 2), dp)
            lower = ceiling(s) + 1
            upper = floor(e) - 1
            best_obj = huge(1.0_dp)
            best_eta = lower
            do eta = lower, upper
                call obj_local_refine_regression_1d(ceiling(s), floor(e), eta, y, x, zeta, obj)
                if (obj < best_obj) then
                    best_obj = obj
                    best_eta = eta
                end if
            end do
            cpt_refined(k) = best_eta
        end do
        deallocate(cpt_ext)
    end subroutine solve_local_refine_regression_1d

    subroutine solve_cv_dp_regression_1d(y, x, gamma, lambda, delta, eps, cpt_hat, k_hat, test_error, train_error)
        !> Cross-validation for DP.regression by odd/even sample splitting, matching changepoints::CV.DP.regression for scalar lambda.
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: gamma, lambda, eps
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: cpt_hat(:)
        integer, intent(out) :: k_hat
        real(kind=dp), intent(out) :: test_error, train_error

        integer :: n, n_train, n_valid, len, j, p, seg_start, seg_end
        logical :: has_consecutive
        integer, allocatable :: odd_idx(:), even_idx(:), cpt_train(:), cpt_train_long(:), init_cpt(:), interval(:, :), partition(:)
        real(kind=dp), allocatable :: train_y(:), validation_y(:), train_x(:, :), validation_x(:, :), beta_hat(:)

        n = size(y)
        if (mod(n, 2) /= 0) error stop "solve_cv_dp_regression_1d requires an even number of observations"
        p = size(x, 2)
        n_train = n / 2
        n_valid = n / 2
        allocate(odd_idx(n_train), even_idx(n_valid))
        do j = 1, n_train
            odd_idx(j) = 2 * j - 1
            even_idx(j) = 2 * j
        end do

        allocate(train_y(n_train), validation_y(n_valid), train_x(n_train, p), validation_x(n_valid, p))
        train_y = y(odd_idx)
        validation_y = y(even_idx)
        do j = 1, p
            train_x(:, j) = x(odd_idx, j)
            validation_x(:, j) = x(even_idx, j)
        end do

        call solve_dp_regression_1d(train_y, train_x, gamma, lambda, delta, eps, partition, cpt_train)
        len = size(cpt_train)
        allocate(init_cpt(len))
        if (len > 0) init_cpt = odd_idx(cpt_train)
        k_hat = len
        allocate(cpt_hat(len))
        if (len > 0) cpt_hat = init_cpt

        if (len > 0) then
            has_consecutive = .false.
            do j = 1, len - 1
                if (cpt_train(j + 1) - cpt_train(j) == 1) then
                    has_consecutive = .true.
                    exit
                end if
            end do
            if (has_consecutive) then
                test_error = huge(1.0_dp)
                train_error = huge(1.0_dp)
                deallocate(odd_idx, even_idx, train_y, validation_y, train_x, validation_x, partition, cpt_train, init_cpt)
                return
            end if
        end if

        allocate(cpt_train_long(len + 1), interval(len + 1, 2))
        if (len > 0) then
            cpt_train_long(1:len) = cpt_train
        end if
        cpt_train_long(len + 1) = n_train
        interval(1, :) = [1, cpt_train_long(1)]
        if (len > 0) then
            do j = 2, len + 1
                interval(j, :) = [cpt_train_long(j - 1) + 1, cpt_train_long(j)]
            end do
        end if

        train_error = 0.0_dp
        test_error = 0.0_dp
        do j = 1, len + 1
            seg_start = interval(j, 1)
            seg_end = interval(j, 2)
            call fit_lasso_segment_regression_1d(train_y, train_x, seg_start, seg_end, lambda, eps, beta_hat)
            train_error = train_error + sum((train_y(seg_start:seg_end) - matmul(train_x(seg_start:seg_end, :), beta_hat)) ** 2)
            test_error = test_error + sum((validation_y(seg_start:seg_end) - matmul(validation_x(seg_start:seg_end, :), beta_hat)) ** 2)
            deallocate(beta_hat)
        end do

        deallocate(odd_idx, even_idx, train_y, validation_y, train_x, validation_x, partition, cpt_train, init_cpt, cpt_train_long, interval)
    end subroutine solve_cv_dp_regression_1d

    recursive subroutine wbs_univar_rec(y, cs, alpha, beta, s, e, delta, level, parent_idx, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        real(kind=dp), intent(in) :: y(:), cs(0:)
        integer, intent(in) :: alpha(:), beta(:)
        integer, intent(in) :: s, e, delta, level, parent_idx
        integer, intent(inout) :: n_nodes
        real(kind=dp), intent(inout) :: s_out(:), dval_out(:)
        integer, intent(inout) :: level_out(:), parent_out(:, :), parent_idx_out(:)

        integer :: m, m_star, n_valid, idx_node, level_here, best_t
        integer, allocatable :: alpha_new(:), beta_new(:)
        real(kind=dp), allocatable :: a(:), b(:)
        real(kind=dp) :: best_value

        call valid_intervals(alpha, beta, s, e, delta, alpha_new, beta_new, n_valid)
        if (n_valid == 0) return

        level_here = level + 1
        allocate(a(n_valid), b(n_valid))
        do m = 1, n_valid
            call best_cusum_on_interval(cs, alpha_new(m), beta_new(m), delta, best_t, best_value)
            a(m) = best_value
            b(m) = real(best_t, dp)
        end do
        m_star = maxloc(a, dim=1)

        n_nodes = n_nodes + 1
        idx_node = n_nodes
        s_out(idx_node) = b(m_star)
        dval_out(idx_node) = a(m_star)
        level_out(idx_node) = level_here
        parent_out(:, idx_node) = [s, e]
        parent_idx_out(idx_node) = parent_idx

        call wbs_univar_rec(y, cs, alpha, beta, s, int(b(m_star)) - 1, delta, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        call wbs_univar_rec(y, cs, alpha, beta, int(b(m_star)), e, delta, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)

        deallocate(alpha_new, beta_new, a, b)
    end subroutine wbs_univar_rec

    recursive subroutine wbs_univar_rob_rec(y, alpha, beta, s, e, k_huber, delta, level, parent_idx, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        real(kind=dp), intent(in) :: y(:), k_huber
        integer, intent(in) :: alpha(:), beta(:)
        integer, intent(in) :: s, e, delta, level, parent_idx
        integer, intent(inout) :: n_nodes
        real(kind=dp), intent(inout) :: s_out(:), dval_out(:)
        integer, intent(inout) :: level_out(:), parent_out(:, :), parent_idx_out(:)

        integer :: m, m_star, n_valid, idx_node, level_here, best_t
        integer, allocatable :: alpha_new(:), beta_new(:)
        real(kind=dp), allocatable :: a(:), b(:)
        real(kind=dp) :: best_value

        call valid_intervals(alpha, beta, s, e, delta, alpha_new, beta_new, n_valid)
        if (n_valid == 0) return

        level_here = level + 1
        allocate(a(n_valid), b(n_valid))
        do m = 1, n_valid
            call best_huber_on_interval(y, alpha_new(m), beta_new(m), k_huber, best_t, best_value)
            a(m) = best_value
            b(m) = real(best_t, dp)
        end do
        m_star = maxloc(a, dim=1)

        n_nodes = n_nodes + 1
        idx_node = n_nodes
        s_out(idx_node) = b(m_star)
        dval_out(idx_node) = a(m_star)
        level_out(idx_node) = level_here
        parent_out(:, idx_node) = [s, e]
        parent_idx_out(idx_node) = parent_idx

        call wbs_univar_rob_rec(y, alpha, beta, s, int(b(m_star)) - 1, k_huber, delta, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        call wbs_univar_rob_rec(y, alpha, beta, int(b(m_star)), e, k_huber, delta, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)

        deallocate(alpha_new, beta_new, a, b)
    end subroutine wbs_univar_rob_rec

    recursive subroutine bs_cov_rec(x, s, e, level, parent_idx, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: s, e, level, parent_idx
        integer, intent(inout) :: n_nodes
        integer, intent(inout) :: s_out(:), level_out(:), parent_out(:, :), parent_idx_out(:)
        real(kind=dp), intent(inout) :: dval_out(:)

        integer :: p, n, delta, s_star, e_star, t, best_t, idx_node, level_here
        real(kind=dp) :: best_value, stat

        p = size(x, 1)
        n = size(x, 2)
        delta = int(2.0_dp * real(p, dp) * log(real(n, dp)) + 1.0_dp)
        if (e - s <= delta) return

        level_here = level + 1
        s_star = ceiling(real(s, dp) + real(p, dp) * log(real(n, dp)))
        e_star = floor(real(e, dp) - real(p, dp) * log(real(n, dp)))
        if (s_star > e_star) return

        best_value = -1.0_dp
        best_t = s_star
        do t = s_star, e_star
            stat = covariance_cusum_norm_1d(x, s, e, t)
            if (stat > best_value) then
                best_value = stat
                best_t = t
            end if
        end do

        n_nodes = n_nodes + 1
        idx_node = n_nodes
        s_out(idx_node) = best_t
        dval_out(idx_node) = best_value
        level_out(idx_node) = level_here
        parent_out(:, idx_node) = [s, e]
        parent_idx_out(idx_node) = parent_idx

        call bs_cov_rec(x, s, best_t - 1, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        call bs_cov_rec(x, best_t, e, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
    end subroutine bs_cov_rec

    recursive subroutine wbsip_cov_rec(x, x_prime, alpha, beta, s, e, delta, level, parent_idx, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        real(kind=dp), intent(in) :: x(:, :), x_prime(:, :)
        integer, intent(in) :: alpha(:), beta(:)
        integer, intent(in) :: s, e, delta, level, parent_idx
        integer, intent(inout) :: n_nodes
        integer, intent(inout) :: s_out(:), level_out(:), parent_out(:, :), parent_idx_out(:)
        real(kind=dp), intent(inout) :: dval_out(:)

        integer :: n, m, m_star, idx_node, level_here, best_t, alpha_m, beta_m, s_star, e_star
        integer :: n_valid, p
        integer, allocatable :: b(:)
        logical, allocatable :: valid(:)
        real(kind=dp), allocatable :: a(:), u_mat(:, :), proj_sq(:, :), proj_cs(:, :)
        real(kind=dp) :: stat, left_mean, right_mean, best_value, logn

        n = size(x, 2)
        p = size(x, 1)
        logn = log(real(n, dp))

        allocate(u_mat(p, size(alpha)))
        call pc_cov_1d(x_prime, alpha, beta, u_mat)

        allocate(proj_sq(size(alpha), n), proj_cs(size(alpha), 0:n))
        do m = 1, size(alpha)
            do idx_node = 1, n
                proj_sq(m, idx_node) = dot_product(u_mat(:, m), x(:, idx_node)) ** 2
            end do
            proj_cs(m, 0) = 0.0_dp
            do idx_node = 1, n
                proj_cs(m, idx_node) = proj_cs(m, idx_node - 1) + proj_sq(m, idx_node)
            end do
        end do

        allocate(a(size(alpha)), b(size(alpha)), valid(size(alpha)))
        a = -1.0_dp
        b = 0
        valid = .false.
        do m = 1, size(alpha)
            alpha_m = ceiling(real(max(alpha(m), s) + delta, dp))
            beta_m = floor(real(min(beta(m), e) - delta, dp))
            if (beta_m - alpha_m >= ceiling(2.0_dp * logn)) then
                s_star = ceiling(real(alpha_m, dp) + logn)
                e_star = floor(real(beta_m, dp) - logn)
                if (s_star <= e_star) then
                    best_value = -1.0_dp
                    best_t = s_star
                    do idx_node = s_star, e_star
                        left_mean = (proj_cs(m, idx_node) - proj_cs(m, alpha_m)) / real(idx_node - alpha_m, dp)
                        right_mean = (proj_cs(m, beta_m) - proj_cs(m, idx_node)) / real(beta_m - idx_node, dp)
                        stat = sqrt(real((idx_node - alpha_m) * (beta_m - idx_node), dp) / real(beta_m - alpha_m, dp)) * abs(left_mean - right_mean)
                        if (stat > best_value) then
                            best_value = stat
                            best_t = idx_node
                        end if
                    end do
                    a(m) = best_value
                    b(m) = best_t
                    valid(m) = .true.
                end if
            end if
        end do

        n_valid = count(valid)
        if (n_valid == 0) then
            deallocate(u_mat, proj_sq, proj_cs, a, b, valid)
            return
        end if

        m_star = maxloc(a, dim=1)
        level_here = level + 1
        n_nodes = n_nodes + 1
        idx_node = n_nodes
        s_out(idx_node) = b(m_star)
        dval_out(idx_node) = a(m_star)
        level_out(idx_node) = level_here
        parent_out(:, idx_node) = [s, e]
        parent_idx_out(idx_node) = parent_idx

        call wbsip_cov_rec(x, x_prime, alpha, beta, s, b(m_star) - 1, delta, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)
        call wbsip_cov_rec(x, x_prime, alpha, beta, b(m_star), e, delta, level_here, idx_node, n_nodes, s_out, dval_out, level_out, parent_out, parent_idx_out)

        deallocate(u_mat, proj_sq, proj_cs, a, b, valid)
    end subroutine wbsip_cov_rec

    subroutine error_pred_seg_regression_1d(y, x, s, e, lambda, delta, eps, mse)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        integer, intent(in) :: s, e, delta
        real(kind=dp), intent(in) :: lambda, eps
        real(kind=dp), intent(out) :: mse
        real(kind=dp), allocatable :: beta_hat(:)

        if (e - s <= 2 * delta) then
            mse = huge(1.0_dp)
            return
        end if

        call fit_lasso_segment_regression_1d(y, x, s, e, lambda, eps, beta_hat)
        mse = sum((y(s:e) - matmul(x(s:e, :), beta_hat)) ** 2)
        deallocate(beta_hat)
    end subroutine error_pred_seg_regression_1d

    subroutine fit_lasso_segment_regression_1d(y, x, s, e, lambda, eps, beta_hat)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        integer, intent(in) :: s, e
        real(kind=dp), intent(in) :: lambda, eps
        real(kind=dp), allocatable, intent(out) :: beta_hat(:)

        integer :: seg_len, n_all, p
        real(kind=dp) :: lambda_scaled

        n_all = size(y)
        p = size(x, 2)
        seg_len = e - s + 1
        lambda_scaled = lambda * sqrt(max(log(real(max(n_all, p), dp)), real(e - s, dp))) * &
                        sqrt(log(real(max(n_all, p), dp))) / real(e - s, dp)
        call lasso_seq_single_1d(x(s:e, :), y(s:e), lambda_scaled, eps, beta_hat)
    end subroutine fit_lasso_segment_regression_1d

    subroutine obj_local_refine_regression_1d(s_extra, e_extra, eta, y, x, zeta, obj)
        integer, intent(in) :: s_extra, e_extra, eta
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: zeta
        real(kind=dp), intent(out) :: obj

        integer :: n_all, p, n_seg
        real(kind=dp) :: lambda_lr
        real(kind=dp), allocatable :: x_convert(:, :), y_convert(:), beta(:), beta1(:), beta2(:), fit(:)

        n_all = size(x, 1)
        p = size(x, 2)
        n_seg = e_extra - s_extra + 1
        allocate(x_convert(n_seg, 2 * p), y_convert(n_seg), beta(2 * p), beta1(p), beta2(p), fit(n_seg))
        call x_glasso_converter_regression_1d(x(s_extra:e_extra, :), eta, s_extra, x_convert)
        y_convert = y(s_extra:e_extra)
        lambda_lr = zeta * sqrt(log(real(max(n_all, p), dp)))
        call group_lasso_two_block_1d(x_convert, y_convert, lambda_lr / real(n_seg, dp), 0.001_dp, beta)
        beta1 = beta(1:p)
        beta2 = beta(p + 1:2 * p)
        fit = matmul(x_convert, beta)
        obj = sum((y_convert - fit) ** 2) + lambda_lr * sum(sqrt(beta1 ** 2 + beta2 ** 2))
        deallocate(x_convert, y_convert, beta, beta1, beta2, fit)
    end subroutine obj_local_refine_regression_1d

    subroutine threshold_wbs_tree_1d(n_nodes, s_in, dval_in, parent_idx_in, tau, n_keep, cpt_hat, dval_hat)
        !> Threshold a WBS tree, keeping nodes above tau and pruning descendants of removed nodes.
        integer, intent(in) :: n_nodes
        integer, intent(in) :: s_in(:), parent_idx_in(:)
        real(kind=dp), intent(in) :: dval_in(:), tau
        integer, intent(out) :: n_keep
        integer, allocatable, intent(out) :: cpt_hat(:)
        real(kind=dp), allocatable, intent(out) :: dval_hat(:)

        integer :: i
        logical, allocatable :: keep(:), parent_kept(:)

        if (tau <= 0.0_dp) error stop "threshold_wbs_tree_1d requires tau > 0"

        allocate(keep(n_nodes), parent_kept(n_nodes))
        keep = .false.
        parent_kept = .false.
        do i = 1, n_nodes
            if (parent_idx_in(i) == 0) then
                parent_kept(i) = .true.
            else
                parent_kept(i) = keep(parent_idx_in(i))
            end if
            keep(i) = parent_kept(i) .and. (dval_in(i) > tau)
        end do

        n_keep = count(keep)
        allocate(cpt_hat(n_keep), dval_hat(n_keep))
        if (n_keep > 0) then
            cpt_hat = pack(s_in, keep)
            dval_hat = pack(dval_in, keep)
        end if

        deallocate(keep, parent_kept)
    end subroutine threshold_wbs_tree_1d

    subroutine local_refine_wbs_univar_1d(cpt_init, y, cpt_refined)
        !> Local refinement matching changepoints::local.refine.univar.
        integer, intent(in) :: cpt_init(:)
        real(kind=dp), intent(in) :: y(:)
        integer, allocatable, intent(out) :: cpt_refined(:)

        real(kind=dp), parameter :: w = 0.9_dp
        integer :: n, k, ncp, lower, upper, eta, best_eta
        integer, allocatable :: cpt_ext(:)
        real(kind=dp) :: s, e, best_cost, cost

        n = size(y)
        ncp = size(cpt_init)
        allocate(cpt_refined(ncp))
        if (ncp == 0) return

        allocate(cpt_ext(ncp + 2))
        cpt_ext = [0, cpt_init, n]
        do k = 1, ncp
            s = w * real(cpt_ext(k), dp) + (1.0_dp - w) * real(cpt_ext(k + 1), dp)
            e = (1.0_dp - w) * real(cpt_ext(k + 1), dp) + w * real(cpt_ext(k + 2), dp)
            lower = ceiling(s) + 1
            upper = floor(e) - 1
            best_cost = huge(1.0_dp)
            best_eta = lower
            do eta = lower, upper
                cost = segment_sse_univar(y, ceiling(s), eta) + segment_sse_univar(y, eta + 1, floor(e))
                if (cost < best_cost) then
                    best_cost = cost
                    best_eta = eta
                end if
            end do
            cpt_refined(k) = best_eta
        end do
        deallocate(cpt_ext)
    end subroutine local_refine_wbs_univar_1d

    subroutine prefix_sum_1d(y, cs)
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), intent(out) :: cs(0:)
        integer :: i
        cs(0) = 0.0_dp
        do i = 1, size(y)
            cs(i) = cs(i - 1) + y(i)
        end do
    end subroutine prefix_sum_1d

    subroutine prefix_sum_sq_1d(y, css)
        real(kind=dp), intent(in) :: y(:)
        real(kind=dp), intent(out) :: css(0:)
        integer :: i
        css(0) = 0.0_dp
        do i = 1, size(y)
            css(i) = css(i - 1) + y(i) * y(i)
        end do
    end subroutine prefix_sum_sq_1d

    subroutine valid_intervals(alpha, beta, s, e, delta, alpha_new, beta_new, n_valid)
        integer, intent(in) :: alpha(:), beta(:), s, e, delta
        integer, allocatable, intent(out) :: alpha_new(:), beta_new(:)
        integer, intent(out) :: n_valid

        integer :: i, a, b
        integer, allocatable :: tmp_a(:), tmp_b(:)

        allocate(tmp_a(size(alpha)), tmp_b(size(alpha)))
        n_valid = 0
        do i = 1, size(alpha)
            a = max(alpha(i), s)
            b = min(beta(i), e)
            if (b - a > 2 * delta) then
                n_valid = n_valid + 1
                tmp_a(n_valid) = a
                tmp_b(n_valid) = b
            end if
        end do
        allocate(alpha_new(n_valid), beta_new(n_valid))
        if (n_valid > 0) then
            alpha_new = tmp_a(1:n_valid)
            beta_new = tmp_b(1:n_valid)
        end if
        deallocate(tmp_a, tmp_b)
    end subroutine valid_intervals

    subroutine best_cusum_on_interval(cs, alpha, beta, delta, best_t, best_value)
        real(kind=dp), intent(in) :: cs(0:)
        integer, intent(in) :: alpha, beta, delta
        integer, intent(out) :: best_t
        real(kind=dp), intent(out) :: best_value

        integer :: t
        real(kind=dp) :: left_mean, right_mean, stat

        best_value = -1.0_dp
        best_t = alpha + delta
        do t = alpha + delta, beta - delta
            left_mean = (cs(t) - cs(alpha)) / real(t - alpha, dp)
            right_mean = (cs(beta) - cs(t)) / real(beta - t, dp)
            stat = sqrt(real((t - alpha) * (beta - t), dp) / real(beta - alpha, dp)) * abs(left_mean - right_mean)
            if (stat > best_value) then
                best_value = stat
                best_t = t
            end if
        end do
    end subroutine best_cusum_on_interval

    subroutine best_huber_on_interval(y, alpha, beta, k_huber, best_t, best_value)
        real(kind=dp), intent(in) :: y(:), k_huber
        integer, intent(in) :: alpha, beta
        integer, intent(out) :: best_t
        real(kind=dp), intent(out) :: best_value

        integer :: n_seg, i
        real(kind=dp) :: mu_huber, partial_sum, score
        real(kind=dp), allocatable :: x(:), r1(:)

        n_seg = beta - alpha
        allocate(x(n_seg), r1(n_seg))
        x = y(alpha + 1:beta)
        mu_huber = huber_mean_1d(x, k_huber)
        r1 = max(min(x - mu_huber, k_huber), -k_huber)

        best_value = -1.0_dp
        best_t = alpha + 1
        partial_sum = 0.0_dp
        do i = 1, n_seg - 1
            partial_sum = partial_sum + r1(i)
            score = partial_sum * partial_sum * real(n_seg, dp) / real((n_seg - i) * i, dp)
            if (score > best_value) then
                best_value = score
                best_t = alpha + i
            end if
        end do

        deallocate(x, r1)
    end subroutine best_huber_on_interval

    real(kind=dp) function covariance_cusum_norm_1d(x, s, e, t) result(val)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: s, e, t

        integer :: p, n_st, n_te, n_se
        real(kind=dp) :: scale
        real(kind=dp), allocatable :: cov_st(:, :), cov_te(:, :), diff(:, :)

        p = size(x, 1)
        n_st = t - s
        n_te = e - t
        n_se = e - s

        allocate(cov_st(p, p), cov_te(p, p), diff(p, p))
        call sample_cov_cols_1d(x, s + 1, t, cov_st)
        call sample_cov_cols_1d(x, t + 1, e, cov_te)
        scale = sqrt(real(n_st * n_te, dp) / real(n_se, dp))
        diff = scale * (cov_st - cov_te)
        val = symmetric_operator_norm_1d(diff)
        deallocate(cov_st, cov_te, diff)
    end function covariance_cusum_norm_1d

    subroutine covariance_cusum_matrix_1d(x, s, e, t, diff)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: s, e, t
        real(kind=dp), intent(out) :: diff(:, :)

        integer :: n_st, n_te, n_se, p
        real(kind=dp) :: scale
        real(kind=dp), allocatable :: cov_st(:, :), cov_te(:, :)

        p = size(x, 1)
        n_st = t - s
        n_te = e - t
        n_se = e - s
        allocate(cov_st(p, p), cov_te(p, p))
        call sample_cov_cols_1d(x, s + 1, t, cov_st)
        call sample_cov_cols_1d(x, t + 1, e, cov_te)
        scale = sqrt(real(n_st * n_te, dp) / real(n_se, dp))
        diff = scale * (cov_st - cov_te)
        deallocate(cov_st, cov_te)
    end subroutine covariance_cusum_matrix_1d

    subroutine pc_cov_1d(x, alpha, beta, u_mat)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: alpha(:), beta(:)
        real(kind=dp), intent(out) :: u_mat(:, :)

        integer :: p, n, m, delta_cov, best_t, s_star, e_star, t
        real(kind=dp) :: best_value, stat, logn
        real(kind=dp), allocatable :: diff(:, :)

        p = size(x, 1)
        n = size(x, 2)
        logn = log(real(n, dp))
        delta_cov = int(2.0_dp * real(p, dp) * logn + 1.0_dp)
        u_mat = 0.0_dp
        allocate(diff(p, p))
        do m = 1, size(alpha)
            if (beta(m) - alpha(m) > delta_cov) then
                s_star = ceiling(real(alpha(m), dp) + real(p, dp) * logn)
                e_star = floor(real(beta(m), dp) - real(p, dp) * logn)
                if (s_star <= e_star) then
                    best_value = -1.0_dp
                    best_t = s_star
                    do t = s_star, e_star
                        stat = covariance_cusum_norm_1d(x, alpha(m), beta(m), t)
                        if (stat > best_value) then
                            best_value = stat
                            best_t = t
                        end if
                    end do
                    call covariance_cusum_matrix_1d(x, alpha(m), beta(m), best_t, diff)
                    call dominant_symmetric_abs_evec_1d(diff, u_mat(:, m))
                end if
            end if
        end do
        deallocate(diff)
    end subroutine pc_cov_1d

    subroutine sample_cov_cols_1d(x, i1, i2, cov_mat)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: i1, i2
        real(kind=dp), intent(out) :: cov_mat(:, :)

        integer :: p, m, i, j, k
        real(kind=dp), allocatable :: mu(:)

        p = size(x, 1)
        m = i2 - i1 + 1
        if (m <= 1) error stop "sample_cov_cols_1d requires at least two columns"

        allocate(mu(p))
        do i = 1, p
            mu(i) = sum(x(i, i1:i2)) / real(m, dp)
        end do

        cov_mat = 0.0_dp
        do k = i1, i2
            do i = 1, p
                do j = 1, p
                    cov_mat(i, j) = cov_mat(i, j) + (x(i, k) - mu(i)) * (x(j, k) - mu(j))
                end do
            end do
        end do
        cov_mat = cov_mat / real(m - 1, dp)
        deallocate(mu)
    end subroutine sample_cov_cols_1d

    real(kind=dp) function symmetric_operator_norm_1d(a) result(val)
        real(kind=dp), intent(in) :: a(:, :)

        integer :: n, p, q, i, j, iter, max_iter
        real(kind=dp) :: app, aqq, apq, phi, c, s, aip, aiq, max_off
        real(kind=dp), allocatable :: work(:, :)

        n = size(a, 1)
        allocate(work(n, n))
        work = a
        max_iter = max(50 * n * n, 1)

        do iter = 1, max_iter
            max_off = 0.0_dp
            p = 1
            q = 1
            do i = 1, n - 1
                do j = i + 1, n
                    if (abs(work(i, j)) > max_off) then
                        max_off = abs(work(i, j))
                        p = i
                        q = j
                    end if
                end do
            end do
            if (max_off <= 1.0e-12_dp) exit

            app = work(p, p)
            aqq = work(q, q)
            apq = work(p, q)
            phi = 0.5_dp * atan2(2.0_dp * apq, aqq - app)
            c = cos(phi)
            s = sin(phi)

            do i = 1, n
                if (i /= p .and. i /= q) then
                    aip = work(i, p)
                    aiq = work(i, q)
                    work(i, p) = c * aip - s * aiq
                    work(p, i) = work(i, p)
                    work(i, q) = s * aip + c * aiq
                    work(q, i) = work(i, q)
                end if
            end do

            work(p, p) = c * c * app - 2.0_dp * s * c * apq + s * s * aqq
            work(q, q) = s * s * app + 2.0_dp * s * c * apq + c * c * aqq
            work(p, q) = 0.0_dp
            work(q, p) = 0.0_dp
        end do

        val = 0.0_dp
        do i = 1, n
            val = max(val, abs(work(i, i)))
        end do
        deallocate(work)
    end function symmetric_operator_norm_1d

    subroutine dominant_symmetric_abs_evec_1d(a, u)
        real(kind=dp), intent(in) :: a(:, :)
        real(kind=dp), intent(out) :: u(:)

        integer :: n, iter, max_iter
        real(kind=dp) :: norm_v
        real(kind=dp), allocatable :: v(:), w(:)

        n = size(a, 1)
        allocate(v(n), w(n))
        v = 1.0_dp / sqrt(real(n, dp))
        max_iter = max(200, 20 * n)
        do iter = 1, max_iter
            w = matmul(a, matmul(a, v))
            norm_v = sqrt(sum(w * w))
            if (norm_v <= 1.0e-14_dp) exit
            v = w / norm_v
        end do
        if (sqrt(sum(v * v)) <= 1.0e-14_dp) then
            v = 0.0_dp
            v(1) = 1.0_dp
        end if
        u = v
        deallocate(v, w)
    end subroutine dominant_symmetric_abs_evec_1d

    subroutine lasso_seq_single_1d(x, y, lambda, eps, beta_hat)
        real(kind=dp), intent(in) :: x(:, :), y(:)
        real(kind=dp), intent(in) :: lambda, eps
        real(kind=dp), allocatable, intent(out) :: beta_hat(:)

        integer :: n, p
        real(kind=dp) :: y_mean
        real(kind=dp), allocatable :: x_means(:), weights(:), xtilde(:, :), ytilde(:), beta_std(:)

        n = size(x, 1)
        p = size(x, 2)
        allocate(x_means(p), weights(p), xtilde(n, p), ytilde(n), beta_std(p), beta_hat(p))

        y_mean = sum(y) / real(n, dp)
        ytilde = y - y_mean
        x_means = 0.0_dp
        weights = 0.0_dp
        call standardize_xy_1d(x, y, x_means, weights, xtilde, ytilde)
        call lasso_standardized_1d(xtilde, ytilde, lambda, eps, beta_std)
        beta_hat = beta_std / weights

        deallocate(x_means, weights, xtilde, ytilde, beta_std)
    end subroutine lasso_seq_single_1d

    subroutine x_glasso_converter_regression_1d(x, eta, s_ceil, xx)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: eta, s_ceil
        real(kind=dp), intent(out) :: xx(:, :)

        integer :: n_seg, p, t

        n_seg = size(x, 1)
        p = size(x, 2)
        t = eta - s_ceil + 1
        xx(:, 1:p) = x
        xx(:, p + 1:2 * p) = x
        if (t < n_seg) xx(t + 1:n_seg, 1:p) = 0.0_dp
        if (t >= 1) xx(1:t, p + 1:2 * p) = 0.0_dp
        xx(:, 1:p) = xx(:, 1:p) / sqrt(real(t - 1, dp))
        xx(:, p + 1:2 * p) = xx(:, p + 1:2 * p) / sqrt(real(n_seg - t, dp))
    end subroutine x_glasso_converter_regression_1d

    subroutine group_lasso_two_block_1d(x, y, lambda, eps, beta)
        real(kind=dp), intent(in) :: x(:, :), y(:)
        real(kind=dp), intent(in) :: lambda, eps
        real(kind=dp), intent(out) :: beta(:)

        integer :: n, p2, g, iter, max_iter, j1, j2
        real(kind=dp) :: step, obj_old, obj_new, diff, normg
        real(kind=dp), allocatable :: grad(:), beta_old(:), beta_next(:), xbeta(:), tempg(:)

        n = size(x, 1)
        p2 = size(x, 2)
        allocate(grad(p2), beta_old(p2), beta_next(p2), xbeta(n), tempg(2))
        beta = 0.0_dp
        step = 1.0_dp / max(2.0_dp * symmetric_operator_norm_1d(matmul(transpose(x), x)), 1.0_dp)
        max_iter = 5000
        obj_old = huge(1.0_dp)
        do iter = 1, max_iter
            beta_old = beta
            xbeta = matmul(x, beta_old)
            grad = 2.0_dp * matmul(transpose(x), xbeta - y)
            beta_next = beta_old - step * grad
            do g = 1, p2 / 2
                j1 = g
                j2 = g + p2 / 2
                tempg = beta_next([j1, j2])
                normg = sqrt(sum(tempg ** 2))
                if (normg > step * lambda) then
                    tempg = (1.0_dp - step * lambda / normg) * tempg
                else
                    tempg = 0.0_dp
                end if
                beta_next(j1) = tempg(1)
                beta_next(j2) = tempg(2)
            end do
            beta = beta_next
            obj_new = sum((y - matmul(x, beta)) ** 2)
            do g = 1, p2 / 2
                obj_new = obj_new + lambda * sqrt(beta(g) ** 2 + beta(g + p2 / 2) ** 2)
            end do
            diff = abs(obj_old - obj_new)
            if (diff < eps) exit
            obj_old = obj_new
        end do
        deallocate(grad, beta_old, beta_next, xbeta, tempg)
    end subroutine group_lasso_two_block_1d

    subroutine standardize_xy_1d(x, y, x_means, weights, xtilde, ytilde)
        real(kind=dp), intent(in) :: x(:, :), y(:)
        real(kind=dp), intent(out) :: x_means(:), weights(:), xtilde(:, :), ytilde(:)

        integer :: n, p, j
        real(kind=dp), allocatable :: x_centered(:, :)

        n = size(x, 1)
        p = size(x, 2)
        allocate(x_centered(n, p))
        ytilde = y - sum(y) / real(n, dp)
        do j = 1, p
            x_means(j) = sum(x(:, j)) / real(n, dp)
            x_centered(:, j) = x(:, j) - x_means(j)
            weights(j) = sqrt(sum(x_centered(:, j) ** 2) / real(n, dp))
            if (weights(j) <= 1.0e-12_dp) weights(j) = 1.0_dp
            xtilde(:, j) = x_centered(:, j) / weights(j)
        end do
        deallocate(x_centered)
    end subroutine standardize_xy_1d

    subroutine lasso_standardized_1d(xtilde, ytilde, lambda, eps, beta_out)
        real(kind=dp), intent(in) :: xtilde(:, :), ytilde(:)
        real(kind=dp), intent(in) :: lambda, eps
        real(kind=dp), intent(out) :: beta_out(:)

        integer :: n, p, i
        real(kind=dp) :: loss_old, loss_new, loss_diff, tmp
        real(kind=dp), allocatable :: beta_last(:), beta_new(:), res_vec(:)

        n = size(xtilde, 1)
        p = size(xtilde, 2)
        allocate(beta_last(p), beta_new(p), res_vec(n))
        beta_new = 0.0_dp
        beta_last = 0.0_dp
        res_vec = ytilde
        loss_diff = 100.0_dp
        do while (loss_diff >= eps)
            beta_last = beta_new
            loss_old = lasso_standardized_obj_1d(xtilde, ytilde, beta_last, lambda)
            do i = 1, p
                tmp = beta_last(i) + dot_product(xtilde(:, i), res_vec) / real(n, dp)
                beta_new(i) = soft_threshold_scalar_1d(tmp, lambda)
                res_vec = res_vec + xtilde(:, i) * (beta_last(i) - beta_new(i))
            end do
            loss_new = lasso_standardized_obj_1d(xtilde, ytilde, beta_new, lambda)
            loss_diff = loss_old - loss_new
            if (loss_diff < 0.0_dp .and. abs(loss_diff) < 1.0e-12_dp) exit
        end do
        beta_out = beta_new
        deallocate(beta_last, beta_new, res_vec)
    end subroutine lasso_standardized_1d

    subroutine lasso_dpdu_standardized_1d(mtilde, vtilde, beta_start, n, lambda, eps, beta_out)
        real(kind=dp), intent(in) :: mtilde(:, :), vtilde(:), beta_start(:)
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: lambda, eps
        real(kind=dp), intent(out) :: beta_out(:)

        integer :: p, i
        real(kind=dp) :: loss_old, loss_new, loss_diff, tmp
        real(kind=dp), allocatable :: beta_last(:), beta_new(:)

        p = size(vtilde)
        allocate(beta_last(p), beta_new(p))
        beta_last = beta_start
        beta_new = beta_start
        loss_diff = 100.0_dp
        do while (loss_diff >= eps)
            beta_last = beta_new
            loss_old = lasso_dpdu_standardized_obj_1d(mtilde, vtilde, beta_last, n, lambda)
            do i = 1, p
                tmp = beta_last(i) + (vtilde(i) - dot_product(mtilde(i, :), beta_new)) / real(n, dp)
                beta_new(i) = soft_threshold_scalar_1d(tmp, lambda)
            end do
            loss_new = lasso_dpdu_standardized_obj_1d(mtilde, vtilde, beta_new, n, lambda)
            loss_diff = loss_old - loss_new
            if (loss_diff < 0.0_dp .and. abs(loss_diff) < 1.0e-12_dp) exit
        end do
        beta_out = beta_new
        deallocate(beta_last, beta_new)
    end subroutine lasso_dpdu_standardized_1d

    real(kind=dp) function lasso_standardized_obj_1d(xtilde, ytilde, beta, lambda) result(obj)
        real(kind=dp), intent(in) :: xtilde(:, :), ytilde(:), beta(:)
        real(kind=dp), intent(in) :: lambda
        integer :: n

        n = size(xtilde, 1)
        obj = sum((ytilde - matmul(xtilde, beta)) ** 2) / (2.0_dp * real(n, dp)) + lambda * sum(abs(beta))
    end function lasso_standardized_obj_1d

    real(kind=dp) function lasso_dpdu_standardized_obj_1d(mtilde, vtilde, beta, n, lambda) result(obj)
        real(kind=dp), intent(in) :: mtilde(:, :), vtilde(:), beta(:)
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: lambda

        obj = (dot_product(beta, matmul(mtilde, beta)) - 2.0_dp * dot_product(vtilde, beta)) / real(n, dp) + lambda * sum(abs(beta))
    end function lasso_dpdu_standardized_obj_1d

    real(kind=dp) function soft_threshold_scalar_1d(x, lambda) result(val)
        real(kind=dp), intent(in) :: x, lambda

        if (x > lambda) then
            val = x - lambda
        else if (x < -lambda) then
            val = x + lambda
        else
            val = 0.0_dp
        end if
    end function soft_threshold_scalar_1d

    function outer_product_row_1d(x, y) result(mat)
        real(kind=dp), intent(in) :: x(:), y(:)
        real(kind=dp) :: mat(size(x), size(y))
        integer :: i, j

        do i = 1, size(x)
            do j = 1, size(y)
                mat(i, j) = x(i) * y(j)
            end do
        end do
    end function outer_product_row_1d

    function scale_symmetric_matrix_1d(a, weights) result(b)
        real(kind=dp), intent(in) :: a(:, :), weights(:)
        real(kind=dp) :: b(size(a, 1), size(a, 2))
        integer :: i, j

        do i = 1, size(a, 1)
            do j = 1, size(a, 2)
                b(i, j) = a(i, j) / (weights(i) * weights(j))
            end do
        end do
    end function scale_symmetric_matrix_1d

    real(kind=dp) function huber_mean_1d(x, tau) result(mu_new)
        real(kind=dp), intent(in) :: x(:), tau
        real(kind=dp), parameter :: eps = 1.0e-8_dp
        real(kind=dp) :: mu_old
        real(kind=dp), allocatable :: r(:), weights(:)

        allocate(r(size(x)), weights(size(x)))
        mu_new = sum(x) / real(size(x), dp)
        mu_old = 0.0_dp
        weights = 1.0_dp
        do while (abs(mu_new - mu_old) > eps)
            mu_old = mu_new
            r = x - mu_new
            where (abs(r) > tau)
                weights = tau / abs(r)
            end where
            mu_new = sum(x * weights) / sum(weights)
        end do
        deallocate(r, weights)
    end function huber_mean_1d

    real(kind=dp) function rume_1d(y, k_trim) result(mu_hat)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: k_trim
        integer :: n, n_even, n_odd, i, best_i, n_ans
        real(kind=dp) :: lower, upper
        real(kind=dp), allocatable :: x(:), e(:), sec(:), ans(:)

        n = size(y)
        n_even = n / 2
        n_odd = n - n_even
        allocate(x(n_even), e(n_odd))
        x = y(2:n:2)
        e = y(1:n:2)
        call sort_real_inplace(x)

        allocate(sec(n_even - k_trim))
        do i = 1, n_even - k_trim
            sec(i) = x(i + k_trim) - x(i)
        end do
        best_i = minloc(sec, dim=1)
        lower = x(best_i)
        upper = x(best_i + k_trim)

        n_ans = count(e >= lower .and. e <= upper)
        if (n_ans == 0) then
            mu_hat = real(count(e < lower), dp) / real(n_odd, dp) * lower + &
                     real(count(e > upper), dp) / real(n_odd, dp) * upper
        else
            allocate(ans(n_ans))
            ans = pack(e, e >= lower .and. e <= upper)
            mu_hat = sum(ans) / real(n_ans, dp)
            deallocate(ans)
        end if

        deallocate(x, e, sec)
    end function rume_1d

    integer function rtest_1d(dat1, th1, th2) result(flag)
        real(kind=dp), intent(in) :: dat1(:), th1, th2
        real(kind=dp) :: a, b, c, midpoint

        a = real(count(exp(-0.5_dp * (dat1 - th1) ** 2) > exp(-0.5_dp * (dat1 - th2) ** 2)), dp) / real(size(dat1), dp)
        midpoint = 0.5_dp * (th1 + th2)
        b = normal_cdf_1d(midpoint, th1)
        c = normal_cdf_1d(midpoint, th2)
        if (abs(a - b) > abs(a - c)) then
            flag = 1
        else
            flag = 0
        end if
    end function rtest_1d

    real(kind=dp) function normal_cdf_1d(x, mu) result(val)
        real(kind=dp), intent(in) :: x, mu
        val = 0.5_dp * (1.0_dp + erf((x - mu) / sqrt(2.0_dp)))
    end function normal_cdf_1d

    subroutine sort_real_inplace(x)
        real(kind=dp), intent(inout) :: x(:)
        integer :: i, j
        real(kind=dp) :: key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1 .and. x(j) > key)
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_real_inplace

    real(kind=dp) function segment_sse_univar(y, s, e) result(sse)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: s, e
        real(kind=dp) :: mu
        integer :: n
        n = e - s + 1
        if (n <= 0) then
            sse = 0.0_dp
            return
        end if
        mu = sum(y(s:e)) / real(n, dp)
        sse = sum((y(s:e) - mu) ** 2)
    end function segment_sse_univar

    subroutine polynomial_segment_fit_1d(y, n_total, s, e, r, yfit, dist2)
        real(kind=dp), intent(in) :: y(:)
        integer, intent(in) :: n_total, s, e, r
        real(kind=dp), allocatable, intent(out) :: yfit(:)
        real(kind=dp), intent(out) :: dist2

        real(kind=dp), allocatable :: basis(:, :), xtx(:, :), xty(:), coef(:), resid(:)

        call basis_poly_1d(n_total, s, e, r, basis)
        allocate(xtx(r + 1, r + 1), xty(r + 1), coef(r + 1), resid(e - s + 1))
        xtx = matmul(transpose(basis), basis)
        xty = matmul(transpose(basis), y(s:e))
        call solve_linear_system_1d(xtx, xty, coef)
        allocate(yfit(e - s + 1))
        yfit = matmul(basis, coef)
        resid = y(s:e) - yfit
        dist2 = dot_product(resid, resid)
        deallocate(basis, xtx, xty, coef, resid)
    end subroutine polynomial_segment_fit_1d

    subroutine error_pred_seg_var1_1d(x_futu, x_curr, s, e, lambda, delta, mse)
        real(kind=dp), intent(in) :: x_futu(:, :), x_curr(:, :)
        integer, intent(in) :: s, e, delta
        real(kind=dp), intent(in) :: lambda
        real(kind=dp), intent(out) :: mse

        integer :: p, m
        real(kind=dp), allocatable :: xseg(:, :), yseg(:, :), gram(:, :), cross(:, :), tran_hat(:, :), tran_t(:, :), resid(:, :)

        if (abs(lambda) > 1.0e-12_dp) error stop "error_pred_seg_var1_1d currently supports lambda = 0 only"
        if (e - s <= 2 * delta) then
            mse = huge(1.0_dp)
            return
        end if

        p = size(x_curr, 1)
        m = e - s + 1
        allocate(xseg(p, m), yseg(p, m), gram(p, p), cross(p, p), tran_hat(p, p), tran_t(p, p), resid(p, m))
        xseg = x_curr(:, s:e)
        yseg = x_futu(:, s:e)
        gram = matmul(xseg, transpose(xseg))
        cross = matmul(yseg, transpose(xseg))
        call solve_linear_matrix_1d(transpose(gram), transpose(cross), tran_t)
        tran_hat = transpose(tran_t)
        resid = matmul(tran_hat, xseg) - yseg
        mse = sum(resid ** 2)
        deallocate(xseg, yseg, gram, cross, tran_hat, tran_t, resid)
    end subroutine error_pred_seg_var1_1d

    real(kind=dp) function obj_func_lr_poly_1d(y, s_inter, e_inter, eta, r, delta_lr) result(val)
        real(kind=dp), intent(in) :: y(:), s_inter, e_inter
        integer, intent(in) :: eta, r, delta_lr

        integer :: n, s_star, e_star
        real(kind=dp) :: dist1, dist2
        real(kind=dp), allocatable :: fit1(:), fit2(:)

        n = size(y)
        s_star = ceiling(s_inter)
        e_star = floor(e_inter)
        if ((eta - s_star < 2 * delta_lr) .or. (e_star - eta + 1 < 2 * delta_lr)) then
            val = huge(1.0_dp)
        else
            call polynomial_segment_fit_1d(y, n, s_star, eta - 1, r, fit1, dist1)
            call polynomial_segment_fit_1d(y, n, eta, e_star, r, fit2, dist2)
            val = dist1 + dist2
            deallocate(fit1, fit2)
        end if
    end function obj_func_lr_poly_1d

    subroutine basis_poly_1d(n_total, s, e, r, basis)
        integer, intent(in) :: n_total, s, e, r
        real(kind=dp), allocatable, intent(out) :: basis(:, :)

        integer :: i, j, nseg
        real(kind=dp) :: xval

        nseg = e - s + 1
        allocate(basis(nseg, r + 1))
        do i = 1, nseg
            xval = real(s + i - 1, dp) / real(n_total, dp)
            basis(i, 1) = 1.0_dp
            do j = 2, r + 1
                basis(i, j) = basis(i, j - 1) * xval
            end do
        end do
    end subroutine basis_poly_1d

    subroutine solve_linear_system_1d(a, b, x)
        real(kind=dp), intent(in) :: a(:, :), b(:)
        real(kind=dp), intent(out) :: x(:)

        integer :: n, i, k, pivot
        real(kind=dp) :: factor, tmp, maxval_abs
        real(kind=dp), allocatable :: aa(:, :), bb(:), rowtmp(:)

        n = size(b)
        allocate(aa(n, n), bb(n), rowtmp(n))
        aa = a
        bb = b

        do k = 1, n - 1
            pivot = k
            maxval_abs = abs(aa(k, k))
            do i = k + 1, n
                if (abs(aa(i, k)) > maxval_abs) then
                    maxval_abs = abs(aa(i, k))
                    pivot = i
                end if
            end do
            if (pivot /= k) then
                rowtmp = aa(k, :)
                aa(k, :) = aa(pivot, :)
                aa(pivot, :) = rowtmp
                tmp = bb(k)
                bb(k) = bb(pivot)
                bb(pivot) = tmp
            end if
            do i = k + 1, n
                factor = aa(i, k) / aa(k, k)
                aa(i, k:n) = aa(i, k:n) - factor * aa(k, k:n)
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do

        x(n) = bb(n) / aa(n, n)
        do i = n - 1, 1, -1
            x(i) = (bb(i) - dot_product(aa(i, i + 1:n), x(i + 1:n))) / aa(i, i)
        end do

        deallocate(aa, bb, rowtmp)
    end subroutine solve_linear_system_1d

    subroutine solve_linear_matrix_1d(a, b, x)
        real(kind=dp), intent(in) :: a(:, :), b(:, :)
        real(kind=dp), intent(out) :: x(:, :)

        integer :: j
        real(kind=dp), allocatable :: rhs(:), sol(:)

        allocate(rhs(size(b, 1)), sol(size(b, 1)))
        do j = 1, size(b, 2)
            rhs = b(:, j)
            call solve_linear_system_1d(a, rhs, sol)
            x(:, j) = sol
        end do
        deallocate(rhs, sol)
    end subroutine solve_linear_matrix_1d

end module changepoints_pkg_mod
