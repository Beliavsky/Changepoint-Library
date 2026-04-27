module changepoints_pkg_dpdu_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: solve_dpdu_regression_1d, solve_dpdu2_regression_1d, solve_cv_dpdu_regression_1d, solve_local_refine_dpdu_regression_1d, &
              solve_ci_dpdu_regression_1d

contains

    subroutine solve_dpdu_regression_1d(y, x, lambda, zeta, eps, partition, cpt_hat, beta_mat)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: lambda, eps
        integer, intent(in) :: zeta
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)
        real(kind=dp), allocatable, intent(out) :: beta_mat(:, :)

        integer :: n, p, i, l, n_seg, rr, ll, ncpt, t, n_keep
        real(kind=dp) :: b, beta0, lambda_scaled, loss
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
                    beta0 = ymean_new(l) - dot_product(xmeans_new(l, :), beta_hat_seg)
                    beta_temp(1, l) = beta0
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

    subroutine solve_dpdu2_regression_1d(y, x, lambda, zeta, eps, partition, cpt_hat, beta_mat)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: lambda, eps
        integer, intent(in) :: zeta
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)
        real(kind=dp), allocatable, intent(out) :: beta_mat(:, :)

        integer :: n, p, i, l, n_seg, rr, ll, ncpt, t, n_keep
        real(kind=dp) :: b, beta0, ymean_new, ymean_old, lambda_scaled, loss
        real(kind=dp), allocatable :: bestvalue(:), weights(:), vtilde(:), beta_start(:), beta_std(:), beta_hat_seg(:)
        real(kind=dp), allocatable :: xmeans_new(:), xmeans_old(:), v_new(:), v_old(:)
        real(kind=dp), allocatable :: m_new(:, :), m_old(:, :), mcentered(:, :), mtilde(:, :)
        real(kind=dp), allocatable :: beta_temp(:, :), beta_best(:, :), cps_work(:)
        integer, allocatable :: part_work(:), cpt_all(:)

        n = size(y)
        p = size(x, 2)
        allocate(bestvalue(0:n), part_work(0:n), weights(p), vtilde(p), beta_start(p), beta_std(p), beta_hat_seg(p), &
                 xmeans_new(p), xmeans_old(p), v_new(p), v_old(p), m_new(p, p), m_old(p, p), mcentered(p, p), mtilde(p, p), &
                 beta_temp(p + 1, n), beta_best(p + 1, n), cps_work(n))

        bestvalue = 0.0_dp
        part_work = 0
        bestvalue(0) = -real(zeta, dp)
        beta_temp = 0.0_dp
        beta_best = 0.0_dp

        do i = 1, n
            bestvalue(i) = huge(1.0_dp)
            m_old = 0.0_dp
            v_old = 0.0_dp
            ymean_old = 0.0_dp
            xmeans_old = 0.0_dp
            beta_start = 0.0_dp
            do l = i, 1, -1
                n_seg = i - l + 1
                ymean_new = (ymean_old * real(n_seg - 1, dp) + y(l)) / real(n_seg, dp)
                xmeans_new = (xmeans_old * real(n_seg - 1, dp) + x(l, :)) / real(n_seg, dp)
                m_new = m_old + outer_product_row_1d(x(l, :), x(l, :))
                v_new = v_old + y(l) * x(l, :)
                mcentered = m_new - real(n_seg, dp) * outer_product_row_1d(xmeans_new, xmeans_new)
                do t = 1, p
                    weights(t) = sqrt(max(mcentered(t, t) / real(n_seg, dp), 0.0_dp))
                    if (weights(t) <= 1.0e-12_dp) weights(t) = 1.0_dp
                end do
                mtilde = scale_symmetric_matrix_1d(mcentered, weights)
                vtilde = (v_new - real(n_seg, dp) * ymean_new * xmeans_new) / weights
                if (n_seg >= zeta) then
                    lambda_scaled = lambda * sqrt(max(log(real(max(n, p), dp)), real(n_seg - 1, dp))) * &
                                    sqrt(log(real(max(n, p), dp))) / real(n_seg, dp)
                    call lasso_dpdu_standardized_1d(mtilde, vtilde, beta_start, n_seg, lambda_scaled, eps, beta_std)
                    beta_start = beta_std
                    loss = dot_product(beta_std, matmul(mtilde, beta_std)) - 2.0_dp * dot_product(vtilde, beta_std)
                    beta_hat_seg = beta_std / weights
                    beta0 = ymean_new - dot_product(xmeans_new, beta_hat_seg)
                    beta_temp(1, l) = beta0
                    beta_temp(2:p + 1, l) = beta_hat_seg
                    b = bestvalue(l - 1) + real(zeta, dp) + loss
                else
                    b = bestvalue(l - 1) + real(zeta, dp)
                end if
                if (b <= bestvalue(i)) then
                    bestvalue(i) = b
                    part_work(i) = l - 1
                    beta_best(:, i) = beta_temp(:, l)
                end if
                m_old = m_new
                v_old = v_new
                ymean_old = ymean_new
                xmeans_old = xmeans_new
            end do
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

        deallocate(bestvalue, part_work, weights, vtilde, beta_start, beta_std, beta_hat_seg, xmeans_new, xmeans_old, &
                   v_new, v_old, m_new, m_old, mcentered, mtilde, beta_temp, beta_best, cps_work, cpt_all)
    end subroutine solve_dpdu2_regression_1d

    subroutine solve_cv_dpdu_regression_1d(y, x, lambda, zeta, eps, cpt_hat, k_hat, test_error, train_error)
        real(kind=dp), intent(in) :: y(:), x(:, :)
        real(kind=dp), intent(in) :: lambda, eps
        integer, intent(in) :: zeta
        integer, allocatable, intent(out) :: cpt_hat(:)
        integer, intent(out) :: k_hat
        real(kind=dp), intent(out) :: test_error, train_error

        integer :: n, p, n_train, n_valid, i, seg_start, seg_stop, len
        integer, allocatable :: odd_indexes(:), even_indexes(:), train_cpt(:), partition_train(:), init_train_cpt_long(:), init_test_cpt_long(:)
        real(kind=dp), allocatable :: train_y(:), valid_y(:), train_x(:, :), valid_x(:, :), beta_mat(:, :), beta_seg(:, :), xseg(:, :)

        n = size(y)
        p = size(x, 2)
        n_train = (n + 1) / 2
        n_valid = n / 2

        allocate(odd_indexes(n_train), even_indexes(n_valid))
        do i = 1, n_train
            odd_indexes(i) = 2 * i - 1
        end do
        do i = 1, n_valid
            even_indexes(i) = 2 * i
        end do

        allocate(train_y(n_train), train_x(n_train, p), valid_y(n_valid), valid_x(n_valid, p))
        train_y = y(odd_indexes)
        train_x = x(odd_indexes, :)
        valid_y = y(even_indexes)
        valid_x = x(even_indexes, :)

        call solve_dpdu_regression_1d(train_y, train_x, lambda, zeta, eps, partition_train, train_cpt, beta_mat)

        if (size(train_cpt) >= 1) then
            allocate(init_train_cpt_long(size(train_cpt) + 2), init_test_cpt_long(size(train_cpt) + 2))
            init_train_cpt_long(1) = 0
            init_train_cpt_long(2:size(train_cpt) + 1) = train_cpt
            init_train_cpt_long(size(init_train_cpt_long)) = n_train
            init_test_cpt_long(1) = 0
            init_test_cpt_long(2:size(train_cpt) + 1) = train_cpt
            init_test_cpt_long(size(init_test_cpt_long)) = n_valid

            len = size(init_train_cpt_long) - 1
            allocate(beta_seg(p + 1, len))
            beta_seg(:, 1:size(train_cpt)) = beta_mat(:, train_cpt)
            beta_seg(:, len) = beta_mat(:, n_train)

            train_error = 0.0_dp
            test_error = 0.0_dp
            do i = 1, len
                seg_start = init_train_cpt_long(i) + 1
                seg_stop = init_train_cpt_long(i + 1)
                allocate(xseg(seg_stop - seg_start + 1, p + 1))
                xseg(:, 1) = 1.0_dp
                xseg(:, 2:p + 1) = train_x(seg_start:seg_stop, :)
                train_error = train_error + lasso_dpdu_error_1d(train_y(seg_start:seg_stop), xseg, beta_seg(:, i))
                deallocate(xseg)

                seg_start = init_test_cpt_long(i) + 1
                seg_stop = init_test_cpt_long(i + 1)
                allocate(xseg(seg_stop - seg_start + 1, p + 1))
                xseg(:, 1) = 1.0_dp
                xseg(:, 2:p + 1) = valid_x(seg_start:seg_stop, :)
                test_error = test_error + lasso_dpdu_error_1d(valid_y(seg_start:seg_stop), xseg, beta_seg(:, i))
                deallocate(xseg)
            end do

            k_hat = len - 1
            allocate(cpt_hat(size(train_cpt)))
            cpt_hat = odd_indexes(train_cpt)
            deallocate(init_train_cpt_long, init_test_cpt_long, beta_seg)
        else
            k_hat = 0
            allocate(cpt_hat(0))
            allocate(xseg(n_train, p + 1))
            xseg(:, 1) = 1.0_dp
            xseg(:, 2:p + 1) = train_x
            train_error = lasso_dpdu_error_1d(train_y, xseg, beta_mat(:, 1))
            deallocate(xseg)

            allocate(xseg(n_valid, p + 1))
            xseg(:, 1) = 1.0_dp
            xseg(:, 2:p + 1) = valid_x
            test_error = lasso_dpdu_error_1d(valid_y, xseg, beta_mat(:, 1))
            deallocate(xseg)
        end if

        deallocate(odd_indexes, even_indexes, train_y, train_x, valid_y, valid_x, train_cpt, partition_train, beta_mat)
    end subroutine solve_cv_dpdu_regression_1d

    subroutine solve_local_refine_dpdu_regression_1d(cpt_init, beta_hat, y, x, w, cpt_refined)
        integer, intent(in) :: cpt_init(:)
        real(kind=dp), intent(in) :: beta_hat(:, :), y(:), x(:, :)
        real(kind=dp), intent(in) :: w
        integer, allocatable, intent(out) :: cpt_refined(:)

        integer :: n, p, k, cpt_init_numb, lower, upper, eta, best_eta
        integer, allocatable :: cpt_init_long(:)
        real(kind=dp) :: s, e, best_value, value
        real(kind=dp), allocatable :: xseg(:, :)

        n = size(x, 1)
        p = size(x, 2)
        cpt_init_numb = size(cpt_init)
        allocate(cpt_refined(cpt_init_numb))
        if (cpt_init_numb == 0) return

        allocate(cpt_init_long(cpt_init_numb + 2))
        cpt_init_long(1) = 0
        cpt_init_long(2:cpt_init_numb + 1) = cpt_init
        cpt_init_long(cpt_init_numb + 2) = n

        do k = 1, cpt_init_numb
            s = w * real(cpt_init_long(k), dp) + (1.0_dp - w) * real(cpt_init_long(k + 1), dp)
            e = (1.0_dp - w) * real(cpt_init_long(k + 1), dp) + w * real(cpt_init_long(k + 2), dp)
            lower = ceiling(s) + 2
            upper = floor(e) - 2
            best_value = huge(1.0_dp)
            best_eta = lower
            do eta = lower, upper
                allocate(xseg(eta - ceiling(s) + 1, p + 1))
                xseg(:, 1) = 1.0_dp
                xseg(:, 2:p + 1) = x(ceiling(s):eta, :)
                value = lasso_dpdu_error_1d(y(ceiling(s):eta), xseg, beta_hat(:, k))
                deallocate(xseg)

                allocate(xseg(floor(e) - eta, p + 1))
                xseg(:, 1) = 1.0_dp
                xseg(:, 2:p + 1) = x(eta + 1:floor(e), :)
                value = value + lasso_dpdu_error_1d(y(eta + 1:floor(e)), xseg, beta_hat(:, k + 1))
                deallocate(xseg)

                if (value < best_value) then
                    best_value = value
                    best_eta = eta
                end if
            end do
            cpt_refined(k) = ceiling(s) + (best_eta - lower + 1)
        end do

        deallocate(cpt_init_long)
    end subroutine solve_local_refine_dpdu_regression_1d

    subroutine solve_ci_dpdu_regression_1d(cpt_init, cpt_lr, beta_hat, y, x, w, alpha_vec, u_draws, rounding, &
                                           ci_array, block_size, lrv_hat, kappa2_hat, drift_hat)
        integer, intent(in) :: cpt_init(:), cpt_lr(:)
        real(kind=dp), intent(in) :: beta_hat(:, :), y(:), x(:, :)
        real(kind=dp), intent(in) :: w, alpha_vec(:), u_draws(:, :)
        logical, intent(in) :: rounding
        real(kind=dp), allocatable, intent(out) :: ci_array(:, :, :)
        integer, intent(out) :: block_size
        real(kind=dp), allocatable, intent(out) :: lrv_hat(:), kappa2_hat(:), drift_hat(:)

        integer :: n, p1, k, j, m_min, seg_len, start_idx, end_idx, pair_numb, idx
        integer, allocatable :: cpt_init_long(:)
        real(kind=dp) :: s, e, qlo, qhi
        real(kind=dp), allocatable :: interval_refine_mat(:, :), x_full(:, :), diffbeta(:, :), xtx(:, :), z_vec(:), z_block(:), zrev_block(:)

        n = size(y)
        p1 = size(beta_hat, 1)
        if (size(cpt_init) /= size(cpt_lr)) stop "cpt_init and cpt_lr must have the same length"

        allocate(ci_array(size(cpt_init), 2, size(alpha_vec)), lrv_hat(size(cpt_init)), kappa2_hat(size(cpt_init)), drift_hat(size(cpt_init)))
        if (size(cpt_init) == 0) return

        allocate(cpt_init_long(size(cpt_init) + 2), interval_refine_mat(size(cpt_init), 2), x_full(n, p1), diffbeta(p1, size(cpt_init)), xtx(p1, p1))
        cpt_init_long(1) = 0
        cpt_init_long(2:size(cpt_init) + 1) = cpt_init
        cpt_init_long(size(cpt_init_long)) = n
        x_full(:, 1) = 1.0_dp
        x_full(:, 2:p1) = x
        xtx = matmul(transpose(x_full), x_full)

        m_min = huge(1)
        do k = 1, size(cpt_init)
            interval_refine_mat(k, 1) = w * real(cpt_init_long(k), dp) + (1.0_dp - w) * real(cpt_init_long(k + 1), dp)
            interval_refine_mat(k, 2) = (1.0_dp - w) * real(cpt_init_long(k + 1), dp) + w * real(cpt_init_long(k + 2), dp)
            seg_len = floor(interval_refine_mat(k, 2)) - ceiling(interval_refine_mat(k, 1))
            m_min = min(m_min, seg_len)
        end do
        block_size = ceiling((real(m_min, dp) ** (2.0_dp / 5.0_dp)) / 2.0_dp)

        diffbeta = beta_hat(:, 1:size(cpt_init)) - beta_hat(:, 2:size(cpt_init) + 1)
        do k = 1, size(cpt_init)
            kappa2_hat(k) = dot_product(diffbeta(:, k), diffbeta(:, k))
            drift_hat(k) = dot_product(diffbeta(:, k), matmul(xtx, diffbeta(:, k))) / (real(n, dp) * kappa2_hat(k))
        end do

        do k = 1, size(cpt_init)
            s = interval_refine_mat(k, 1)
            e = interval_refine_mat(k, 2)
            start_idx = ceiling(s)
            end_idx = floor(e)
            seg_len = end_idx - start_idx + 1
            allocate(z_vec(seg_len))
            do idx = start_idx, end_idx
                z_vec(idx - start_idx + 1) = (2.0_dp * y(idx) - dot_product(x_full(idx, :), beta_hat(:, k) + beta_hat(:, k + 1))) * &
                                             dot_product(x_full(idx, :), beta_hat(:, k + 1) - beta_hat(:, k))
            end do
            pair_numb = seg_len / (2 * block_size)
            allocate(z_block(2 * pair_numb), zrev_block(2 * pair_numb))
            do idx = 1, 2 * pair_numb
                z_block(idx) = sum(z_vec((idx - 1) * block_size + 1:idx * block_size))
                zrev_block(idx) = sum(z_vec(seg_len - idx * block_size + 1:seg_len - (idx - 1) * block_size))
            end do
            lrv_hat(k) = (mean_pairdiff_sq_1d(z_block, block_size) + mean_pairdiff_sq_1d(zrev_block, block_size)) / (2.0_dp * kappa2_hat(k))
            deallocate(z_vec, z_block, zrev_block)
        end do

        do k = 1, size(cpt_init)
            do j = 1, size(alpha_vec)
                qlo = quantile_type7_1d(u_draws(k, :), alpha_vec(j) / 2.0_dp) / kappa2_hat(k) + real(cpt_lr(k), dp)
                qhi = quantile_type7_1d(u_draws(k, :), 1.0_dp - alpha_vec(j) / 2.0_dp) / kappa2_hat(k) + real(cpt_lr(k), dp)
                if (rounding) then
                    ci_array(k, 1, j) = floor(qlo)
                    ci_array(k, 2, j) = ceiling(qhi)
                else
                    ci_array(k, 1, j) = qlo
                    ci_array(k, 2, j) = qhi
                end if
            end do
        end do

        deallocate(cpt_init_long, interval_refine_mat, x_full, diffbeta, xtx)
    end subroutine solve_ci_dpdu_regression_1d

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

    real(kind=dp) function lasso_dpdu_standardized_obj_1d(mtilde, vtilde, beta, n, lambda) result(obj)
        real(kind=dp), intent(in) :: mtilde(:, :), vtilde(:), beta(:)
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: lambda

        obj = (dot_product(beta, matmul(mtilde, beta)) - 2.0_dp * dot_product(vtilde, beta)) / real(n, dp) + lambda * sum(abs(beta))
    end function lasso_dpdu_standardized_obj_1d

    real(kind=dp) function lasso_dpdu_error_1d(y, x, beta_hat) result(error)
        real(kind=dp), intent(in) :: y(:), x(:, :), beta_hat(:)

        real(kind=dp), allocatable :: resid(:)

        allocate(resid(size(y)))
        resid = y - matmul(x, beta_hat)
        error = dot_product(resid, resid)
        deallocate(resid)
    end function lasso_dpdu_error_1d

    real(kind=dp) function mean_pairdiff_sq_1d(zsum, block_size) result(val)
        real(kind=dp), intent(in) :: zsum(:)
        integer, intent(in) :: block_size
        integer :: i, pair_numb

        pair_numb = size(zsum) / 2
        val = 0.0_dp
        do i = 1, pair_numb
            val = val + (zsum(2 * i - 1) - zsum(2 * i)) ** 2 / (2.0_dp * real(block_size, dp))
        end do
        val = val / real(pair_numb, dp)
    end function mean_pairdiff_sq_1d

    real(kind=dp) function quantile_type7_1d(x, prob) result(q)
        real(kind=dp), intent(in) :: x(:), prob
        real(kind=dp), allocatable :: work(:)
        real(kind=dp) :: h, g
        integer :: n, j

        n = size(x)
        allocate(work(n))
        work = x
        call sort_real_1d(work)
        if (prob <= 0.0_dp) then
            q = work(1)
        else if (prob >= 1.0_dp) then
            q = work(n)
        else
            h = (real(n - 1, dp) * prob) + 1.0_dp
            j = floor(h)
            g = h - real(j, dp)
            if (j >= n) then
                q = work(n)
            else
                q = (1.0_dp - g) * work(j) + g * work(j + 1)
            end if
        end if
        deallocate(work)
    end function quantile_type7_1d

    subroutine sort_real_1d(a)
        real(kind=dp), intent(inout) :: a(:)
        integer :: i, j
        real(kind=dp) :: key

        do i = 2, size(a)
            key = a(i)
            j = i - 1
            do while (j >= 1 .and. a(j) > key)
                a(j + 1) = a(j)
                j = j - 1
            end do
            a(j + 1) = key
        end do
    end subroutine sort_real_1d

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

end module changepoints_pkg_dpdu_mod
