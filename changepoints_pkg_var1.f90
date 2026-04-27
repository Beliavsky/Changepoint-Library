module changepoints_pkg_var1_mod
use kind_mod, only: dp
implicit none
private
public :: solve_dp_var1_1d, solve_cv_dp_var1_1d, solve_local_refine_cv_var1_1d, probe_local_refine_cv_var1_case_1d

contains

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

    real(kind=dp) function lasso_standardized_obj_1d(xtilde, ytilde, beta, lambda_val) result(obj)
        real(kind=dp), intent(in) :: xtilde(:, :), ytilde(:), beta(:), lambda_val
        real(kind=dp), allocatable :: resid(:)
        integer :: nseg

        nseg = size(xtilde, 1)
        allocate(resid(nseg))
        resid = ytilde - matmul(xtilde, beta)
        obj = sum(resid ** 2) / (2.0_dp * real(nseg, dp)) + lambda_val * sum(abs(beta))
        deallocate(resid)
    end function lasso_standardized_obj_1d

    subroutine lasso_standardized_1d(xtilde, ytilde, lambda_val, eps, beta)
        real(kind=dp), intent(in) :: xtilde(:, :), ytilde(:), lambda_val, eps
        real(kind=dp), intent(out) :: beta(:)

        integer :: nseg, p, i
        real(kind=dp), allocatable :: beta_last(:), beta_new(:), res_vec(:)
        real(kind=dp) :: loss_diff, loss_old, loss_new

        nseg = size(xtilde, 1)
        p = size(xtilde, 2)
        allocate(beta_last(p), beta_new(p), res_vec(nseg))
        beta_last = 0.0_dp
        beta_new = 0.0_dp
        res_vec = ytilde
        loss_diff = 100.0_dp

        do while (loss_diff >= eps)
            beta_last = beta_new
            loss_old = lasso_standardized_obj_1d(xtilde, ytilde, beta_last, lambda_val)
            do i = 1, p
                beta_new(i) = soft_threshold_scalar_1d(beta_last(i) + dot_product(xtilde(:, i), res_vec) / real(nseg, dp), lambda_val)
                res_vec = res_vec + xtilde(:, i) * (beta_last(i) - beta_new(i))
            end do
            loss_new = lasso_standardized_obj_1d(xtilde, ytilde, beta_new, lambda_val)
            loss_diff = loss_old - loss_new
        end do

        beta = beta_new
        deallocate(beta_last, beta_new, res_vec)
    end subroutine lasso_standardized_1d

    subroutine solve_dp_var1_1d(x_futu, x_curr, gamma, lambda, delta, partition, cpt_hat)
        real(kind=dp), intent(in) :: x_futu(:, :), x_curr(:, :)
        real(kind=dp), intent(in) :: gamma, lambda
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: partition(:), cpt_hat(:)

        integer :: n, i, l, r, k
        real(kind=dp), allocatable :: bestvalue(:)
        integer, allocatable :: part_work(:), cps_work(:)
        real(kind=dp) :: mse, b

        n = size(x_futu, 2)
        allocate(bestvalue(0:n), part_work(0:n), cps_work(n))
        bestvalue(0) = -gamma
        part_work(0) = 0
        do i = 1, n
            bestvalue(i) = huge(1.0_dp)
            part_work(i) = 0
            do l = 1, i
                call fit_var1_segment_1d(x_futu, x_curr, l, i, lambda, delta, mse=mse)
                b = bestvalue(l - 1) + gamma + mse
                if (b < bestvalue(i)) then
                    bestvalue(i) = b
                    part_work(i) = l - 1
                end if
            end do
        end do

        allocate(partition(n))
        partition = part_work(1:n)
        k = 0
        r = n
        l = partition(r)
        do while (r > 0)
            if (l > 0) then
                k = k + 1
                cps_work(k) = l
            end if
            r = l
            if (r <= 0) exit
            l = partition(r)
        end do
        allocate(cpt_hat(k))
        do i = 1, k
            cpt_hat(i) = cps_work(k - i + 1)
        end do

        deallocate(bestvalue, part_work, cps_work)
    end subroutine solve_dp_var1_1d

    subroutine solve_cv_dp_var1_1d(data, gamma, lambda, delta, cpt_hat, k_hat, test_error, train_error)
        real(kind=dp), intent(in) :: data(:, :)
        real(kind=dp), intent(in) :: gamma, lambda
        integer, intent(in) :: delta
        integer, allocatable, intent(out) :: cpt_hat(:)
        integer, intent(out) :: k_hat
        real(kind=dp), intent(out) :: test_error, train_error

        real(kind=dp), allocatable :: data_temp(:, :), x_train(:, :), x_test(:, :), x_curr_train(:, :), x_futu_train(:, :)
        real(kind=dp), allocatable :: x_curr_test(:, :), x_futu_test(:, :), tran_hat(:, :), train_losses(:)
        integer, allocatable :: partition_train(:), cpt_train(:), init_cpt_long(:), lower(:), upper(:)
        integer :: n_total, ntemp, p, len, j, nseg, lower_j, upper_j
        real(kind=dp) :: mse

        n_total = size(data, 2)
        p = size(data, 1)
        if (mod(n_total, 2) /= 0) then
            allocate(data_temp(p, n_total - 1))
            data_temp = data(:, 2:n_total)
        else
            allocate(data_temp(p, n_total))
            data_temp = data
        end if
        ntemp = size(data_temp, 2)

        allocate(x_train(p, ntemp / 2), x_test(p, ntemp / 2))
        x_train = data_temp(:, 1:ntemp:2)
        x_test = data_temp(:, 2:ntemp:2)
        allocate(x_curr_train(p, ntemp / 2 - 1), x_futu_train(p, ntemp / 2 - 1), x_curr_test(p, ntemp / 2 - 1), x_futu_test(p, ntemp / 2 - 1))
        x_curr_train = x_train(:, 1:ntemp / 2 - 1)
        x_futu_train = x_train(:, 2:ntemp / 2)
        x_curr_test = x_test(:, 1:ntemp / 2 - 1)
        x_futu_test = x_test(:, 2:ntemp / 2)

        call solve_dp_var1_1d(x_futu_train, x_curr_train, gamma, lambda, delta, partition_train, cpt_train)
        len = size(cpt_train)
        allocate(cpt_hat(len))
        if (len > 0) cpt_hat = 2 * cpt_train
        k_hat = len

        allocate(init_cpt_long(len + 1))
        if (len > 0) init_cpt_long(1:len) = cpt_train
        init_cpt_long(len + 1) = size(x_curr_train, 2)
        nseg = len + 1
        allocate(lower(nseg), upper(nseg), train_losses(nseg))
        lower(1) = 1
        upper(1) = init_cpt_long(1)
        if (len > 0) then
            do j = 2, nseg
                lower(j) = init_cpt_long(j - 1) + 1
                upper(j) = init_cpt_long(j)
            end do
        end if

        train_error = 0.0_dp
        test_error = 0.0_dp
        do j = 1, nseg
            lower_j = lower(j)
            upper_j = upper(j)
            call fit_var1_segment_1d(x_futu_train, x_curr_train, lower_j, upper_j, lambda, delta, mse, tran_hat)
            train_losses(j) = mse
            train_error = train_error + mse
            test_error = test_error + error_test_var1_1d(x_futu_test, x_curr_test, lower_j, upper_j, tran_hat)
            if (allocated(tran_hat)) deallocate(tran_hat)
        end do

        deallocate(data_temp, x_train, x_test, x_curr_train, x_futu_train, x_curr_test, x_futu_test, partition_train, cpt_train, init_cpt_long, lower, upper, train_losses)
    end subroutine solve_cv_dp_var1_1d

    subroutine solve_local_refine_cv_var1_1d(cpt_init, data, zeta_set, delta_local, cpt_hat, zeta_best)
        integer, intent(in) :: cpt_init(:)
        real(kind=dp), intent(in) :: data(:, :), zeta_set(:)
        integer, intent(in) :: delta_local
        integer, allocatable, intent(out) :: cpt_hat(:)
        real(kind=dp), intent(out) :: zeta_best

        real(kind=dp), parameter :: w = 0.9_dp
        real(kind=dp), allocatable :: data_temp(:, :), x_curr(:, :), x_futu(:, :), x_curr_train(:, :), x_curr_test(:, :), x_futu_train(:, :), x_futu_test(:, :)
        integer, allocatable :: cpt_hat_train_ext(:)
        integer :: n_total, ntemp, ncurr, p, kk, kk_count, temp_estimate
        real(kind=dp) :: s_inter, e_inter

        n_total = size(data, 2)
        p = size(data, 1)
        if (mod(n_total, 2) == 0) then
            allocate(data_temp(p, n_total - 1))
            data_temp = data(:, 2:n_total)
        else
            allocate(data_temp(p, n_total))
            data_temp = data
        end if
        ntemp = size(data_temp, 2)
        ncurr = ntemp - 1

        allocate(x_curr(p, ncurr), x_futu(p, ncurr))
        x_curr = data_temp(:, 1:ncurr)
        x_futu = data_temp(:, 2:ntemp)
        allocate(x_curr_train(p, (ncurr + 1) / 2), x_curr_test(p, ncurr / 2), x_futu_train(p, (ncurr + 1) / 2), x_futu_test(p, ncurr / 2))
        x_curr_train = x_curr(:, 1:ncurr:2)
        x_curr_test = x_curr(:, 2:ncurr:2)
        x_futu_train = x_futu(:, 1:ncurr:2)
        x_futu_test = x_futu(:, 2:ncurr:2)

        kk_count = size(cpt_init)
        allocate(cpt_hat_train_ext(kk_count + 2))
        cpt_hat_train_ext(1) = 1
        do kk = 1, kk_count
            cpt_hat_train_ext(kk + 1) = round_to_int_even_1d(real(cpt_init(kk), dp) / 2.0_dp)
        end do
        cpt_hat_train_ext(kk_count + 2) = floor_div_int_1d(ncurr, 2)

        zeta_best = zeta_set(1)
        if (kk_count > 0) then
            do kk = 1, kk_count
                call glasso_error_test_var1_1d(cpt_hat_train_ext(kk), cpt_hat_train_ext(kk + 2), x_futu_train, x_curr_train, x_futu_test, x_curr_test, delta_local, zeta_set, zeta_best)
                s_inter = w * real(cpt_hat_train_ext(kk), dp) + (1.0_dp - w) * real(cpt_hat_train_ext(kk + 1), dp)
                e_inter = (1.0_dp - w) * real(cpt_hat_train_ext(kk + 1), dp) + w * real(cpt_hat_train_ext(kk + 2), dp)
                call find_one_change_grouplasso_var1_1d(ceiling_int_1d(s_inter * 2.0_dp), floor_int_1d(e_inter * 2.0_dp), x_futu, x_curr, delta_local, zeta_best, temp_estimate)
                cpt_hat_train_ext(kk + 1) = round_to_int_even_1d(real(temp_estimate, dp) / 2.0_dp)
            end do
        end if

        allocate(cpt_hat(kk_count))
        do kk = 1, kk_count
            cpt_hat(kk) = 2 * cpt_hat_train_ext(kk + 1)
        end do

        deallocate(data_temp, x_curr, x_futu, x_curr_train, x_curr_test, x_futu_train, x_futu_test, cpt_hat_train_ext)
    end subroutine solve_local_refine_cv_var1_1d

    subroutine probe_local_refine_cv_var1_case_1d(data)
        real(kind=dp), intent(in) :: data(:, :)

        integer, parameter :: s = 20, e = 60
        integer, parameter :: eta_list(4) = [25, 38, 45, 48]
        real(kind=dp), parameter :: zeta_list(2) = [0.01_dp, 0.05_dp]
        real(kind=dp), allocatable :: data_temp(:, :), x_curr(:, :), x_futu(:, :), x_curr_train(:, :), x_curr_test(:, :), x_futu_train(:, :), x_futu_test(:, :)
        real(kind=dp), allocatable :: x_convert(:, :), x_test_convert(:, :), y_vec(:), y_test(:), beta(:), resid(:)
        integer :: n_total, ntemp, ncurr, p, ii, jj, m, eta_abs, eta_rel, len_seg
        real(kind=dp) :: zeta_val, obj_term, test_term

        n_total = size(data, 2)
        p = size(data, 1)
        if (mod(n_total, 2) == 0) then
            allocate(data_temp(p, n_total - 1))
            data_temp = data(:, 2:n_total)
        else
            allocate(data_temp(p, n_total))
            data_temp = data
        end if
        ntemp = size(data_temp, 2)
        ncurr = ntemp - 1

        allocate(x_curr(p, ncurr), x_futu(p, ncurr))
        x_curr = data_temp(:, 1:ncurr)
        x_futu = data_temp(:, 2:ntemp)
        allocate(x_curr_train(p, (ncurr + 1) / 2), x_curr_test(p, ncurr / 2), x_futu_train(p, (ncurr + 1) / 2), x_futu_test(p, ncurr / 2))
        x_curr_train = x_curr(:, 1:ncurr:2)
        x_curr_test = x_curr(:, 2:ncurr:2)
        x_futu_train = x_futu(:, 1:ncurr:2)
        x_futu_test = x_futu(:, 2:ncurr:2)

        len_seg = e - s + 1
        allocate(y_vec(len_seg), y_test(len_seg), beta(2 * p), resid(len_seg))

        write(*, '(A,I0)') 'probe s            = ', s
        write(*, '(A,I0)') 'probe e            = ', e
        do ii = 1, size(eta_list)
            eta_abs = eta_list(ii)
            eta_rel = eta_abs - s + 2
            do jj = 1, size(zeta_list)
                zeta_val = zeta_list(jj)
                write(*, '(A)') 'probe mode         = objective'
                write(*, '(A,I0)') 'probe eta          = ', eta_abs
                write(*, '(A,F0.8)') 'probe zeta         = ', zeta_val
                allocate(x_convert(len_seg, 2 * p))
                call x_glasso_converter_var1_1d(x_curr_train(:, s:e), eta_abs, s, x_convert)
                do m = 1, p
                    y_vec = x_futu_train(m, s:e)
                    call solve_group_lasso_var1_1d(x_convert, y_vec, zeta_val / real(len_seg, dp), build_group_labels_rep_1_to_p_twice_1d(p), beta)
                    resid = y_vec - matmul(x_convert, beta)
                    obj_term = sum(resid ** 2) + zeta_val * sqrt(sum(beta ** 2))
                    write(*, '(A,I0)') 'probe response     = ', m
                    write(*, '(A,*(1X,ES24.16E3))') 'probe beta         = ', beta
                    write(*, '(A,ES24.16E3)') 'probe obj term     = ', obj_term
                end do
                deallocate(x_convert)

                write(*, '(A)') 'probe mode         = test'
                write(*, '(A,I0)') 'probe eta          = ', eta_abs
                write(*, '(A,I0)') 'probe eta rel      = ', eta_rel
                write(*, '(A,F0.8)') 'probe zeta         = ', zeta_val
                allocate(x_convert(len_seg, 2 * p), x_test_convert(len_seg, 2 * p))
                call x_glasso_converter_var1_1d(x_curr_train(:, s:e), eta_rel, 1, x_convert)
                call x_glasso_converter_var1_1d(x_curr_test(:, s:e), eta_rel, 1, x_test_convert)
                do m = 1, p
                    y_vec = x_futu_train(m, s:e)
                    y_test = x_futu_test(m, s:e)
                    call solve_group_lasso_var1_1d(x_convert, y_vec, zeta_val / real(len_seg, dp), build_group_labels_each2_1d(p), beta)
                    resid = y_test - matmul(x_test_convert, beta)
                    test_term = sum(resid ** 2) + zeta_val * sqrt(sum(beta ** 2))
                    write(*, '(A,I0)') 'probe response     = ', m
                    write(*, '(A,*(1X,ES24.16E3))') 'probe beta         = ', beta
                    write(*, '(A,ES24.16E3)') 'probe test term    = ', test_term
                end do
                deallocate(x_convert, x_test_convert)
            end do
        end do

        deallocate(data_temp, x_curr, x_futu, x_curr_train, x_curr_test, x_futu_train, x_futu_test, y_vec, y_test, beta, resid)
    end subroutine probe_local_refine_cv_var1_case_1d

    subroutine fit_var1_segment_1d(x_futu, x_curr, s, e, lambda, delta, mse, tran_hat_out)
        real(kind=dp), intent(in) :: x_futu(:, :), x_curr(:, :)
        integer, intent(in) :: s, e, delta
        real(kind=dp), intent(in) :: lambda
        real(kind=dp), intent(out) :: mse
        real(kind=dp), allocatable, intent(out), optional :: tran_hat_out(:, :)

        integer :: p, nseg, j
        real(kind=dp), allocatable :: xseg(:, :), yseg(:, :), tran_hat(:, :), resid(:, :)
        real(kind=dp), allocatable :: xobs(:, :), xobs_center(:, :), xtilde(:, :), yvec(:), ytilde(:), beta_std(:), beta_raw(:), weights(:), xmeans(:)
        real(kind=dp) :: lambda_val, eps

        if (e - s <= 2 * delta) then
            mse = huge(1.0_dp)
            if (present(tran_hat_out)) allocate(tran_hat_out(0, 0))
            return
        end if

        p = size(x_curr, 1)
        nseg = e - s + 1
        eps = 0.001_dp
        lambda_val = lambda / sqrt(real(e - s, dp))
        allocate(xseg(p, nseg), yseg(p, nseg), tran_hat(p, p), resid(p, nseg))
        allocate(xobs(nseg, p), xobs_center(nseg, p), xtilde(nseg, p), yvec(nseg), ytilde(nseg), beta_std(p), beta_raw(p), weights(p), xmeans(p))
        xseg = x_curr(:, s:e)
        yseg = x_futu(:, s:e)
        xobs = transpose(xseg)
        xmeans = sum(xobs, dim=1) / real(nseg, dp)
        xobs_center = xobs
        do j = 1, p
            xobs_center(:, j) = xobs_center(:, j) - xmeans(j)
        end do
        weights = sqrt(sum(xobs_center ** 2, dim=1) / real(nseg, dp))
        do j = 1, p
            if (weights(j) < 1.0e-12_dp) weights(j) = 1.0_dp
            xtilde(:, j) = xobs_center(:, j) / weights(j)
        end do
        tran_hat = 0.0_dp
        do j = 1, p
            yvec = yseg(j, :)
            ytilde = yvec - sum(yvec) / real(nseg, dp)
            call lasso_standardized_1d(xtilde, ytilde, lambda_val, eps, beta_std)
            beta_raw = beta_std / weights
            tran_hat(j, :) = beta_raw
        end do
        resid = matmul(tran_hat, xseg) - yseg
        mse = sum(resid ** 2)
        if (present(tran_hat_out)) then
            allocate(tran_hat_out(p, p))
            tran_hat_out = tran_hat
        end if
        deallocate(xseg, yseg, tran_hat, resid, xobs, xobs_center, xtilde, yvec, ytilde, beta_std, beta_raw, weights, xmeans)
    end subroutine fit_var1_segment_1d

    real(kind=dp) function error_test_var1_1d(x_futu, x_curr, lower, upper, tran_hat) result(res)
        real(kind=dp), intent(in) :: x_futu(:, :), x_curr(:, :)
        integer, intent(in) :: lower, upper
        real(kind=dp), intent(in) :: tran_hat(:, :)
        real(kind=dp), allocatable :: diff(:, :)
        if (size(tran_hat, 1) == 0 .or. size(tran_hat, 2) == 0) then
            res = huge(1.0_dp)
            return
        end if
        allocate(diff(size(x_futu, 1), upper - lower + 1))
        diff = x_futu(:, lower:upper) - matmul(tran_hat, x_curr(:, lower:upper))
        res = sum(diff ** 2)
        deallocate(diff)
    end function error_test_var1_1d

    subroutine glasso_error_test_var1_1d(s, e, y_train, x_train, y_test, x_test, delta_local, zeta_set, zeta_best)
        integer, intent(in) :: s, e, delta_local
        real(kind=dp), intent(in) :: y_train(:, :), x_train(:, :), y_test(:, :), x_test(:, :), zeta_set(:)
        real(kind=dp), intent(out) :: zeta_best

        integer :: p, lenz, ll, m, estimate
        real(kind=dp), allocatable :: test_error_temp(:)
        real(kind=dp) :: res, best_val

        p = size(y_train, 1)
        lenz = size(zeta_set)
        allocate(test_error_temp(lenz))
        do ll = 1, lenz
            call find_one_change_grouplasso_var1_1d(s, e, y_train, x_train, delta_local, zeta_set(ll), estimate)
            test_error_temp(ll) = 0.0_dp
            do m = 1, p
                res = test_res_glasso_var1_1d(p, estimate - s + 2, y_train(m, s:e), x_train(:, s:e), y_test(m, s:e), x_test(:, s:e), zeta_set(ll))
                test_error_temp(ll) = test_error_temp(ll) + res
            end do
        end do
        best_val = test_error_temp(1)
        zeta_best = zeta_set(1)
        do ll = 2, lenz
            if (test_error_temp(ll) < best_val) then
                best_val = test_error_temp(ll)
                zeta_best = zeta_set(ll)
            end if
        end do
        deallocate(test_error_temp)
    end subroutine glasso_error_test_var1_1d

    subroutine find_one_change_grouplasso_var1_1d(s, e, y_train, x_train, delta_local, zeta_group, estimate)
        integer, intent(in) :: s, e, delta_local
        real(kind=dp), intent(in) :: y_train(:, :), x_train(:, :)
        real(kind=dp), intent(in) :: zeta_group
        integer, intent(out) :: estimate

        integer :: t
        real(kind=dp) :: best_obj, obj

        estimate = (s + e) / 2
        if (e - s > 2 * delta_local) then
            best_obj = huge(1.0_dp)
            do t = s + delta_local, e - delta_local
                obj = obj_lr_var1_1d(s, e, t, y_train, x_train, zeta_group)
                if (obj < best_obj) then
                    best_obj = obj
                    estimate = t
                end if
            end do
        end if
    end subroutine find_one_change_grouplasso_var1_1d

    real(kind=dp) function obj_lr_var1_1d(s_inter_ceil, e_inter_floor, eta, x_futu, x_curr, zeta_group) result(val)
        integer, intent(in) :: s_inter_ceil, e_inter_floor, eta
        real(kind=dp), intent(in) :: x_futu(:, :), x_curr(:, :)
        real(kind=dp), intent(in) :: zeta_group

        integer :: p, len_inter, m
        real(kind=dp), allocatable :: x_convert(:, :), y_convert(:, :), y_vec(:), beta(:), resid(:)

        p = size(x_futu, 1)
        len_inter = e_inter_floor - s_inter_ceil + 1
        allocate(x_convert(len_inter, 2 * p), y_convert(p, len_inter), y_vec(len_inter), beta(2 * p), resid(len_inter))
        call x_glasso_converter_var1_1d(x_curr(:, s_inter_ceil:e_inter_floor), eta, s_inter_ceil, x_convert)
        y_convert = x_futu(:, s_inter_ceil:e_inter_floor)
        val = 0.0_dp
        do m = 1, p
            y_vec = y_convert(m, :)
            call solve_group_lasso_var1_1d(x_convert, y_vec, zeta_group / real(len_inter, dp), build_group_labels_each2_1d(p), beta)
            resid = y_vec - matmul(x_convert, beta)
            val = val + sum(resid ** 2) + zeta_group * sqrt(sum(beta ** 2))
        end do
        deallocate(x_convert, y_convert, y_vec, beta, resid)
    end function obj_lr_var1_1d

    real(kind=dp) function test_res_glasso_var1_1d(p, eta, y_train, x_train, y_test, x_test, zeta_group) result(res)
        integer, intent(in) :: p, eta
        real(kind=dp), intent(in) :: y_train(:), x_train(:, :), y_test(:), x_test(:, :)
        real(kind=dp), intent(in) :: zeta_group

        real(kind=dp), allocatable :: x_convert(:, :), x_test_convert(:, :), beta(:), resid(:)
        integer :: nseg

        nseg = size(x_train, 2)
        allocate(x_convert(nseg, 2 * p), x_test_convert(size(x_test, 2), 2 * p), beta(2 * p), resid(size(y_test)))
        call x_glasso_converter_var1_1d(x_train, eta, 1, x_convert)
        call x_glasso_converter_var1_1d(x_test, eta, 1, x_test_convert)
        call solve_group_lasso_var1_1d(x_convert, y_train, zeta_group / real(nseg, dp), build_group_labels_each2_1d(p), beta)
        resid = y_test - matmul(x_test_convert, beta)
        res = sum(resid ** 2) + zeta_group * sqrt(sum(beta ** 2))
        deallocate(x_convert, x_test_convert, beta, resid)
    end function test_res_glasso_var1_1d

    subroutine x_glasso_converter_var1_1d(x, eta, s_ceil, xxx)
        real(kind=dp), intent(in) :: x(:, :)
        integer, intent(in) :: eta, s_ceil
        real(kind=dp), intent(out) :: xxx(:, :)

        integer :: n, p, t, pp
        real(kind=dp), allocatable :: xx1(:, :), xx2(:, :), xx(:, :)

        n = size(x, 2)
        p = size(x, 1)
        t = eta - s_ceil + 1
        allocate(xx1(n, p), xx2(n, p), xx(n, 2 * p))
        xx1 = transpose(x)
        xx2 = xx1
        if (t + 1 <= n) xx1(t + 1:n, :) = 0.0_dp
        if (t >= 1) xx2(1:t, :) = 0.0_dp
        xx(:, 1:p) = xx1 / sqrt(real(t - 1, dp))
        xx(:, p + 1:2 * p) = xx2 / sqrt(real(n - t, dp))
        do pp = 1, p
            xxx(:, 2 * pp - 1) = xx(:, pp)
            xxx(:, 2 * pp) = xx(:, pp + p)
        end do
        deallocate(xx1, xx2, xx)
    end subroutine x_glasso_converter_var1_1d

    function build_group_labels_rep_1_to_p_twice_1d(p) result(groups)
        integer, intent(in) :: p
        integer :: groups(2 * p)
        integer :: j
        do j = 1, 2 * p
            groups(j) = mod(j - 1, p) + 1
        end do
    end function build_group_labels_rep_1_to_p_twice_1d

    function build_group_labels_each2_1d(p) result(groups)
        integer, intent(in) :: p
        integer :: groups(2 * p)
        integer :: j
        do j = 1, p
            groups(2 * j - 1) = j
            groups(2 * j) = j
        end do
    end function build_group_labels_each2_1d

    subroutine solve_group_lasso_var1_1d(x, y, lambda_val, groups, beta)
        real(kind=dp), intent(in) :: x(:, :), y(:), lambda_val
        integer, intent(in) :: groups(:)
        real(kind=dp), intent(out) :: beta(:)

        integer :: n, d, g, iter, ng, idx_count, j
        real(kind=dp) :: diff, gamma_g, tau_g
        real(kind=dp), allocatable :: beta_old(:), resid(:), old_group(:), new_group(:), score(:)
        real(kind=dp), allocatable :: xg(:, :)
        integer, allocatable :: idx(:)

        n = size(x, 1)
        d = size(x, 2)
        ng = maxval(groups)
        allocate(beta_old(d), resid(n))
        beta = 0.0_dp
        resid = y
        do iter = 1, 1000
            beta_old = beta
            do g = 1, ng
                idx_count = count(groups == g)
                allocate(idx(idx_count), xg(n, idx_count), old_group(idx_count), new_group(idx_count), score(idx_count))
                idx_count = 0
                do j = 1, d
                    if (groups(j) == g) then
                        idx_count = idx_count + 1
                        idx(idx_count) = j
                    end if
                end do
                xg = x(:, idx)
                old_group = beta(idx)
                gamma_g = max(group_lipschitz_var1_1d(xg), 1.0e-8_dp)
                score = old_group + matmul(transpose(xg), resid) / (real(n, dp) * gamma_g)
                tau_g = lambda_val * sqrt(real(size(idx), dp)) / gamma_g
                call apply_group_shrink_vector_1d(score, tau_g, new_group)
                beta(idx) = new_group
                resid = resid - matmul(xg, new_group - old_group)
                deallocate(idx, xg, old_group, new_group, score)
            end do
            diff = maxval(abs(beta - beta_old))
            if (diff < 1.0e-8_dp) exit
        end do
        deallocate(beta_old, resid)
    end subroutine solve_group_lasso_var1_1d

    real(kind=dp) function group_lipschitz_var1_1d(xg) result(val)
        real(kind=dp), intent(in) :: xg(:, :)

        integer :: n, k, iter
        real(kind=dp), allocatable :: gram(:, :), v(:), w(:)
        real(kind=dp) :: nv

        n = size(xg, 1)
        k = size(xg, 2)
        allocate(gram(k, k), v(k), w(k))
        gram = matmul(transpose(xg), xg) / real(n, dp)
        v = 1.0_dp / sqrt(real(k, dp))
        do iter = 1, 50
            w = matmul(gram, v)
            nv = sqrt(sum(w ** 2))
            if (nv <= 1.0e-16_dp) exit
            v = w / nv
        end do
        val = dot_product(v, matmul(gram, v))
        deallocate(gram, v, w)
    end function group_lipschitz_var1_1d

    subroutine apply_group_shrink_vector_1d(v, tau, out)
        real(kind=dp), intent(in) :: v(:), tau
        real(kind=dp), intent(out) :: out(:)
        real(kind=dp) :: nv, scale
        nv = sqrt(sum(v ** 2))
        if (nv <= tau .or. nv <= 1.0e-16_dp) then
            out = 0.0_dp
        else
            scale = 1.0_dp - tau / nv
            out = scale * v
        end if
    end subroutine apply_group_shrink_vector_1d

    integer function ceiling_int_1d(x) result(v)
        real(kind=dp), intent(in) :: x
        integer :: i
        i = int(x)
        if (real(i, dp) < x) then
            v = i + 1
        else
            v = i
        end if
    end function ceiling_int_1d

    integer function floor_int_1d(x) result(v)
        real(kind=dp), intent(in) :: x
        v = floor_div_int_1d(int(floor(x)), 1)
    end function floor_int_1d

    integer function floor_div_int_1d(a, b) result(v)
        integer, intent(in) :: a, b
        v = a / b
    end function floor_div_int_1d

    integer function round_to_int_even_1d(x) result(v)
        real(kind=dp), intent(in) :: x
        real(kind=dp) :: frac
        integer :: base
        base = int(floor(x))
        frac = x - real(base, dp)
        if (frac < 0.5_dp) then
            v = base
        else if (frac > 0.5_dp) then
            v = base + 1
        else
            if (mod(base, 2) == 0) then
                v = base
            else
                v = base + 1
            end if
        end if
    end function round_to_int_even_1d

end module changepoints_pkg_var1_mod
