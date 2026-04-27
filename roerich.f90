module roerich_mod
use kind_mod, only: dp
implicit none
private
public :: solve_sliding_windows_energy_1d, solve_cpdc_qda_klsym_1d, solve_cpdc_qda_klsym_cv_1d, solve_rulsif_linear_pesym_1d

contains

subroutine solve_sliding_windows_energy_1d(z, window_size, cp, best_score, periods, step)
    !> Raw sliding-window energy-distance detector matching the core of roerich.SlidingWindows.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_size
    integer, intent(out) :: cp
    real(kind=dp), intent(out) :: best_score
    integer, intent(in), optional :: periods, step

    integer :: p, jump, n, n_eval, i, t_idx, t_best
    real(kind=dp), allocatable :: x_auto(:, :), x_ref(:, :), x_test(:, :)
    real(kind=dp) :: score

    p = 1
    if (present(periods)) p = periods
    jump = 1
    if (present(step)) jump = step

    n = size(z)
    cp = -1
    best_score = -huge(1.0_dp)
    if (n < 2 * window_size) return

    allocate(x_auto(n, p))
    call autoregression_matrix_1d(z, p, x_auto)

    n_eval = 0
    do t_idx = 2 * window_size, n, jump
        n_eval = n_eval + 1
    end do

    t_best = -1
    do t_idx = 2 * window_size, n, jump
        allocate(x_ref(window_size, p), x_test(window_size, p))
        do i = 1, window_size
            x_ref(i, :) = x_auto(t_idx - 2 * window_size + i, :)
            x_test(i, :) = x_auto(t_idx - window_size + i, :)
        end do
        score = energy_distance_mv(x_ref, x_test)
        if (score > best_score) then
            best_score = score
            t_best = t_idx - 1
        end if
        deallocate(x_ref, x_test)
    end do

    if (t_best >= 0) cp = t_best - window_size
    deallocate(x_auto)
end subroutine solve_sliding_windows_energy_1d

subroutine solve_cpdc_qda_klsym_1d(z, window_size, cp, best_score, periods, step, seed)
    !> Deterministic roerich ChangePointDetectionClassifier-style detector with QDA and KL_sym.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_size
    integer, intent(out) :: cp
    real(kind=dp), intent(out) :: best_score
    integer, intent(in), optional :: periods, step, seed

    integer :: p, jump, seed_use, n, t_idx, t_best, i
    real(kind=dp), allocatable :: x_auto(:, :), x_ref(:, :), x_test(:, :)
    real(kind=dp) :: score

    p = 1
    if (present(periods)) p = periods
    jump = 1
    if (present(step)) jump = step
    seed_use = 2357
    if (present(seed)) seed_use = seed

    n = size(z)
    cp = -1
    best_score = -huge(1.0_dp)
    if (n < 2 * window_size) return

    allocate(x_auto(n, p))
    call autoregression_matrix_1d(z, p, x_auto)

    t_best = -1
    do t_idx = 2 * window_size, n, jump
        allocate(x_ref(window_size, p), x_test(window_size, p))
        do i = 1, window_size
            x_ref(i, :) = x_auto(t_idx - 2 * window_size + i, :)
            x_test(i, :) = x_auto(t_idx - window_size + i, :)
        end do
        call cpdc_qda_klsym_score(x_ref, x_test, seed_use + t_idx, score)
        if (score > best_score) then
            best_score = score
            t_best = t_idx - 1
        end if
        deallocate(x_ref, x_test)
    end do

    if (t_best >= 0) cp = t_best - window_size
    deallocate(x_auto)
end subroutine solve_cpdc_qda_klsym_1d

subroutine solve_cpdc_qda_klsym_cv_1d(z, window_size, cp, best_score, periods, step, seed, n_splits)
    !> Deterministic roerich ChangePointDetectionClassifierCV-style detector with QDA and KL_sym.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_size
    integer, intent(out) :: cp
    real(kind=dp), intent(out) :: best_score
    integer, intent(in), optional :: periods, step, seed, n_splits

    integer :: p, jump, seed_use, nfold, n, t_idx, t_best, i
    real(kind=dp), allocatable :: x_auto(:, :), x_ref(:, :), x_test(:, :)
    real(kind=dp) :: score

    p = 1
    if (present(periods)) p = periods
    jump = 1
    if (present(step)) jump = step
    seed_use = 2357
    if (present(seed)) seed_use = seed
    nfold = 5
    if (present(n_splits)) nfold = n_splits

    n = size(z)
    cp = -1
    best_score = -huge(1.0_dp)
    if (n < 2 * window_size) return

    allocate(x_auto(n, p))
    call autoregression_matrix_1d(z, p, x_auto)

    t_best = -1
    do t_idx = 2 * window_size, n, jump
        allocate(x_ref(window_size, p), x_test(window_size, p))
        do i = 1, window_size
            x_ref(i, :) = x_auto(t_idx - 2 * window_size + i, :)
            x_test(i, :) = x_auto(t_idx - window_size + i, :)
        end do
        call cpdc_qda_klsym_cv_score(x_ref, x_test, seed_use + t_idx, nfold, score)
        if (score > best_score) then
            best_score = score
            t_best = t_idx - 1
        end if
        deallocate(x_ref, x_test)
    end do

    if (t_best >= 0) cp = t_best - window_size
    deallocate(x_auto)
end subroutine solve_cpdc_qda_klsym_cv_1d

subroutine solve_rulsif_linear_pesym_1d(z, window_size, cp, best_score, periods, step, seed, alpha, l2)
    !> Deterministic linear RuLSIF-style detector with symmetric PE score.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_size
    integer, intent(out) :: cp
    real(kind=dp), intent(out) :: best_score
    integer, intent(in), optional :: periods, step, seed
    real(kind=dp), intent(in), optional :: alpha, l2

    integer :: p, jump, seed_use, n, t_idx, t_best, i
    real(kind=dp) :: alpha_use, l2_use, score
    real(kind=dp), allocatable :: x_auto(:, :), x_ref(:, :), x_test(:, :)

    p = 1
    if (present(periods)) p = periods
    jump = 1
    if (present(step)) jump = step
    seed_use = 2357
    if (present(seed)) seed_use = seed
    alpha_use = 0.05_dp
    if (present(alpha)) alpha_use = alpha
    l2_use = 1.0e-3_dp
    if (present(l2)) l2_use = l2

    n = size(z)
    cp = -1
    best_score = -huge(1.0_dp)
    if (n < 2 * window_size) return

    allocate(x_auto(n, p))
    call autoregression_matrix_1d(z, p, x_auto)

    t_best = -1
    do t_idx = 2 * window_size, n, jump
        allocate(x_ref(window_size, p), x_test(window_size, p))
        do i = 1, window_size
            x_ref(i, :) = x_auto(t_idx - 2 * window_size + i, :)
            x_test(i, :) = x_auto(t_idx - window_size + i, :)
        end do
        call rulsif_linear_pesym_score(x_ref, x_test, seed_use + t_idx, alpha_use, l2_use, score)
        if (score > best_score) then
            best_score = score
            t_best = t_idx - 1
        end if
        deallocate(x_ref, x_test)
    end do

    if (t_best >= 0) cp = t_best - window_size
    deallocate(x_auto)
end subroutine solve_rulsif_linear_pesym_1d

subroutine autoregression_matrix_1d(z, periods, x_auto)
    !> Build the autoregression feature matrix used by roerich for a univariate series.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: periods
    real(kind=dp), intent(out) :: x_auto(size(z), periods)

    integer :: i, lag

    do lag = 1, periods
        do i = 1, size(z)
            if (i - lag + 1 >= 1) then
                x_auto(i, lag) = z(i - lag + 1)
            else
                x_auto(i, lag) = 0.0_dp
            end if
        end do
    end do
end subroutine autoregression_matrix_1d

real(kind=dp) function energy_distance_mv(x_ref, x_test) result(score)
    !> Energy distance between two window matrices.
    real(kind=dp), intent(in) :: x_ref(:, :), x_test(:, :)
    real(kind=dp) :: a, b, c

    a = mean_pairwise_distance_between(x_ref, x_test)
    b = mean_pairwise_distance_within(x_test)
    c = mean_pairwise_distance_within(x_ref)
    score = 2.0_dp * a - b - c
end function energy_distance_mv

real(kind=dp) function mean_pairwise_distance_between(x, y) result(mean_dist)
    real(kind=dp), intent(in) :: x(:, :), y(:, :)
    integer :: i, j
    real(kind=dp) :: total

    total = 0.0_dp
    do i = 1, size(x, 1)
        do j = 1, size(y, 1)
            total = total + euclidean_distance(x(i, :), y(j, :))
        end do
    end do
    mean_dist = total / real(size(x, 1) * size(y, 1), dp)
end function mean_pairwise_distance_between

real(kind=dp) function mean_pairwise_distance_within(x) result(mean_dist)
    real(kind=dp), intent(in) :: x(:, :)
    integer :: i, j
    real(kind=dp) :: total

    total = 0.0_dp
    do i = 1, size(x, 1)
        do j = 1, size(x, 1)
            total = total + euclidean_distance(x(i, :), x(j, :))
        end do
    end do
    mean_dist = total / real(size(x, 1) * size(x, 1), dp)
end function mean_pairwise_distance_within

real(kind=dp) function euclidean_distance(x, y) result(dist)
    real(kind=dp), intent(in) :: x(:), y(:)
    dist = sqrt(sum((x - y) ** 2))
end function euclidean_distance

subroutine cpdc_qda_klsym_score(x_ref, x_test, seed, score)
    !> Deterministic QDA/KL_sym score for one reference/test window pair.
    real(kind=dp), intent(in) :: x_ref(:, :), x_test(:, :)
    integer, intent(in) :: seed
    real(kind=dp), intent(out) :: score

    integer :: n_ref, n_test, p, n_total, i, j, ntrain0, ntest0, ntrain1, ntest1
    real(kind=dp), allocatable :: x(:,:), x_scaled(:,:), mu(:), sig(:)
    integer, allocatable :: y(:), idx0(:), idx1(:), perm0(:), perm1(:), train_idx(:), test_idx(:)
    real(kind=dp), allocatable :: prob1(:), ref_preds(:), test_preds(:)

    n_ref = size(x_ref, 1)
    n_test = size(x_test, 1)
    p = size(x_ref, 2)
    n_total = n_ref + n_test

    allocate(x(n_total, p), y(n_total), mu(p), sig(p), x_scaled(n_total, p))
    x(1:n_ref, :) = x_ref
    x(n_ref + 1:n_total, :) = x_test
    y(1:n_ref) = 0
    y(n_ref + 1:n_total) = 1

    do j = 1, p
        mu(j) = sum(x(:, j)) / real(n_total, dp)
        sig(j) = sqrt(max(sum((x(:, j) - mu(j)) ** 2) / real(n_total, dp), 0.0_dp))
        if (sig(j) < 1.0e-12_dp) sig(j) = 1.0_dp
        x_scaled(:, j) = (x(:, j) - mu(j)) / sig(j)
    end do

    allocate(idx0(n_ref), idx1(n_test))
    do i = 1, n_ref
        idx0(i) = i
    end do
    do i = 1, n_test
        idx1(i) = n_ref + i
    end do
    call random_permutation_indices_pm(n_ref, seed + 37, perm0)
    call random_permutation_indices_pm(n_test, seed + 74, perm1)

    ntest0 = n_ref / 2
    ntrain0 = n_ref - ntest0
    ntest1 = n_test / 2
    ntrain1 = n_test - ntest1
    allocate(train_idx(ntrain0 + ntrain1), test_idx(ntest0 + ntest1))
    do i = 1, ntest0
        test_idx(i) = idx0(perm0(i))
    end do
    do i = 1, ntrain0
        train_idx(i) = idx0(perm0(ntest0 + i))
    end do
    do i = 1, ntest1
        test_idx(ntest0 + i) = idx1(perm1(i))
    end do
    do i = 1, ntrain1
        train_idx(ntrain0 + i) = idx1(perm1(ntest1 + i))
    end do

    call qda_predict_proba_binary(x_scaled, y, train_idx, test_idx, 0.01_dp, prob1)

    allocate(ref_preds(ntest0), test_preds(ntest1))
    j = 0
    do i = 1, size(test_idx)
        if (y(test_idx(i)) == 0) then
            j = j + 1
            ref_preds(j) = prob1(i)
        end if
    end do
    j = 0
    do i = 1, size(test_idx)
        if (y(test_idx(i)) == 1) then
            j = j + 1
            test_preds(j) = prob1(i)
        end if
    end do
    score = klsym_metric(ref_preds, test_preds)

    deallocate(x, y, mu, sig, x_scaled, idx0, idx1, perm0, perm1, train_idx, test_idx, prob1, ref_preds, test_preds)
end subroutine cpdc_qda_klsym_score

subroutine cpdc_qda_klsym_cv_score(x_ref, x_test, seed, n_splits, score)
    !> Deterministic stratified-KFold QDA/KL_sym score for one reference/test pair.
    real(kind=dp), intent(in) :: x_ref(:, :), x_test(:, :)
    integer, intent(in) :: seed, n_splits
    real(kind=dp), intent(out) :: score

    integer :: n_ref, n_test, p, n_total, i, nfold
    real(kind=dp), allocatable :: x(:,:), x_scaled(:,:), mu(:), sig(:)
    integer, allocatable :: y(:), idx0(:), idx1(:), perm0(:), perm1(:), fold0(:), fold1(:), train_idx(:), test_idx(:), y_test_all(:)
    real(kind=dp), allocatable :: prob_all(:), prob1(:)

    n_ref = size(x_ref, 1)
    n_test = size(x_test, 1)
    p = size(x_ref, 2)
    n_total = n_ref + n_test
    nfold = min(max(2, n_splits), min(n_ref, n_test))

    allocate(x(n_total, p), y(n_total), mu(p), sig(p), x_scaled(n_total, p))
    x(1:n_ref, :) = x_ref
    x(n_ref + 1:n_total, :) = x_test
    y(1:n_ref) = 0
    y(n_ref + 1:n_total) = 1

    do i = 1, p
        mu(i) = sum(x(:, i)) / real(n_total, dp)
        sig(i) = sqrt(max(sum((x(:, i) - mu(i)) ** 2) / real(n_total, dp), 0.0_dp))
        if (sig(i) < 1.0e-12_dp) sig(i) = 1.0_dp
        x_scaled(:, i) = (x(:, i) - mu(i)) / sig(i)
    end do

    allocate(idx0(n_ref), idx1(n_test))
    do i = 1, n_ref
        idx0(i) = i
    end do
    do i = 1, n_test
        idx1(i) = n_ref + i
    end do
    call random_permutation_indices_pm(n_ref, seed + 37, perm0)
    call random_permutation_indices_pm(n_test, seed + 74, perm1)
    allocate(fold0(n_ref), fold1(n_test))
    do i = 1, n_ref
        fold0(i) = modulo(i - 1, nfold) + 1
    end do
    do i = 1, n_test
        fold1(i) = modulo(i - 1, nfold) + 1
    end do

    allocate(prob_all(n_total), y_test_all(n_total))
    prob_all = 0.0_dp
    y_test_all = -1

    do i = 1, nfold
        call build_cv_indices(idx0, perm0, fold0, idx1, perm1, fold1, i, train_idx, test_idx)
        call qda_predict_proba_binary(x_scaled, y, train_idx, test_idx, 0.01_dp, prob1)
        prob_all(test_idx) = prob1
        y_test_all(test_idx) = y(test_idx)
        deallocate(train_idx, test_idx, prob1)
    end do

    call score_klsym_from_predictions(prob_all, y_test_all, score)
    deallocate(x, y, mu, sig, x_scaled, idx0, idx1, perm0, perm1, fold0, fold1, prob_all, y_test_all)
end subroutine cpdc_qda_klsym_cv_score

subroutine qda_predict_proba_binary(x, y, train_idx, test_idx, reg_param, prob1)
    !> Binary QDA posterior probabilities for class 1 on test_idx.
    real(kind=dp), intent(in) :: x(:, :), reg_param
    integer, intent(in) :: y(:), train_idx(:), test_idx(:)
    real(kind=dp), allocatable, intent(out) :: prob1(:)

    integer :: p, i, n0, n1
    real(kind=dp), allocatable :: mu0(:), mu1(:), cov0(:,:), cov1(:,:), inv0(:,:), inv1(:,:)
    real(kind=dp), allocatable :: diff(:)
    real(kind=dp) :: prior0, prior1, logp0, logp1, det0, det1

    p = size(x, 2)
    n0 = count(y(train_idx) == 0)
    n1 = count(y(train_idx) == 1)
    prior0 = real(n0, dp) / real(size(train_idx), dp)
    prior1 = real(n1, dp) / real(size(train_idx), dp)

    allocate(mu0(p), mu1(p), cov0(p, p), cov1(p, p), inv0(p, p), inv1(p, p), diff(p), prob1(size(test_idx)))
    call class_mean_cov(x, y, train_idx, 0, reg_param, mu0, cov0)
    call class_mean_cov(x, y, train_idx, 1, reg_param, mu1, cov1)
    call invert_spd(cov0, inv0, det0)
    call invert_spd(cov1, inv1, det1)

    do i = 1, size(test_idx)
        diff = x(test_idx(i), :) - mu0
        logp0 = log(prior0) - 0.5_dp * log(max(det0, 1.0e-12_dp)) - 0.5_dp * quad_form(inv0, diff)
        diff = x(test_idx(i), :) - mu1
        logp1 = log(prior1) - 0.5_dp * log(max(det1, 1.0e-12_dp)) - 0.5_dp * quad_form(inv1, diff)
        prob1(i) = 1.0_dp / (1.0_dp + exp(logp0 - logp1))
    end do

    deallocate(mu0, mu1, cov0, cov1, inv0, inv1, diff)
end subroutine qda_predict_proba_binary

subroutine build_cv_indices(idx0, perm0, fold0, idx1, perm1, fold1, fold_id, train_idx, test_idx)
    integer, intent(in) :: idx0(:), perm0(:), fold0(:), idx1(:), perm1(:), fold1(:), fold_id
    integer, allocatable, intent(out) :: train_idx(:), test_idx(:)
    integer :: i, ntrain, ntest, p

    ntrain = 0
    ntest = 0
    do i = 1, size(idx0)
        if (fold0(i) == fold_id) then
            ntest = ntest + 1
        else
            ntrain = ntrain + 1
        end if
    end do
    do i = 1, size(idx1)
        if (fold1(i) == fold_id) then
            ntest = ntest + 1
        else
            ntrain = ntrain + 1
        end if
    end do

    allocate(train_idx(ntrain), test_idx(ntest))
    ntrain = 0
    ntest = 0
    do i = 1, size(idx0)
        if (fold0(i) == fold_id) then
            ntest = ntest + 1
            test_idx(ntest) = idx0(perm0(i))
        else
            ntrain = ntrain + 1
            train_idx(ntrain) = idx0(perm0(i))
        end if
    end do
    do i = 1, size(idx1)
        if (fold1(i) == fold_id) then
            ntest = ntest + 1
            test_idx(ntest) = idx1(perm1(i))
        else
            ntrain = ntrain + 1
            train_idx(ntrain) = idx1(perm1(i))
        end if
    end do
end subroutine build_cv_indices

subroutine score_klsym_from_predictions(prob_all, y_test_all, score)
    real(kind=dp), intent(in) :: prob_all(:)
    integer, intent(in) :: y_test_all(:)
    real(kind=dp), intent(out) :: score
    integer :: n0, n1, i, j0, j1
    real(kind=dp), allocatable :: ref_preds(:), test_preds(:)

    n0 = count(y_test_all == 0)
    n1 = count(y_test_all == 1)
    allocate(ref_preds(n0), test_preds(n1))
    j0 = 0
    j1 = 0
    do i = 1, size(prob_all)
        if (y_test_all(i) == 0) then
            j0 = j0 + 1
            ref_preds(j0) = prob_all(i)
        else if (y_test_all(i) == 1) then
            j1 = j1 + 1
            test_preds(j1) = prob_all(i)
        end if
    end do
    score = klsym_metric(ref_preds, test_preds)
    deallocate(ref_preds, test_preds)
end subroutine score_klsym_from_predictions

subroutine rulsif_linear_pesym_score(x_ref, x_test, seed, alpha, l2, score)
    !> Deterministic two-sided linear RuLSIF score for one reference/test pair.
    real(kind=dp), intent(in) :: x_ref(:, :), x_test(:, :), alpha, l2
    integer, intent(in) :: seed
    real(kind=dp), intent(out) :: score

    integer :: n_ref, n_test, p, n_total, i, ntest0, ntest1, ntrain0, ntrain1
    real(kind=dp), allocatable :: x(:,:), x_scaled(:,:), mu(:), sig(:)
    integer, allocatable :: y(:), idx0(:), idx1(:), perm0(:), perm1(:), train_idx(:), test_idx(:)
    real(kind=dp) :: score_right, score_left

    n_ref = size(x_ref, 1)
    n_test = size(x_test, 1)
    p = size(x_ref, 2)
    n_total = n_ref + n_test

    allocate(x(n_total, p), y(n_total), mu(p), sig(p), x_scaled(n_total, p))
    x(1:n_ref, :) = x_ref
    x(n_ref + 1:n_total, :) = x_test
    y(1:n_ref) = 0
    y(n_ref + 1:n_total) = 1

    do i = 1, p
        mu(i) = sum(x(:, i)) / real(n_total, dp)
        sig(i) = sqrt(max(sum((x(:, i) - mu(i)) ** 2) / real(n_total, dp), 0.0_dp))
        if (sig(i) < 1.0e-12_dp) sig(i) = 1.0_dp
        x_scaled(:, i) = (x(:, i) - mu(i)) / sig(i)
    end do

    allocate(idx0(n_ref), idx1(n_test))
    do i = 1, n_ref
        idx0(i) = i
    end do
    do i = 1, n_test
        idx1(i) = n_ref + i
    end do
    call random_permutation_indices_pm(n_ref, seed + 37, perm0)
    call random_permutation_indices_pm(n_test, seed + 74, perm1)

    ntest0 = n_ref / 2
    ntrain0 = n_ref - ntest0
    ntest1 = n_test / 2
    ntrain1 = n_test - ntest1
    allocate(train_idx(ntrain0 + ntrain1), test_idx(ntest0 + ntest1))
    do i = 1, ntest0
        test_idx(i) = idx0(perm0(i))
    end do
    do i = 1, ntrain0
        train_idx(i) = idx0(perm0(ntest0 + i))
    end do
    do i = 1, ntest1
        test_idx(ntest0 + i) = idx1(perm1(i))
    end do
    do i = 1, ntrain1
        train_idx(ntrain0 + i) = idx1(perm1(ntest1 + i))
    end do

    call one_side_rulsif_score(x_scaled, y, train_idx, test_idx, alpha, l2, score_right)
    call one_side_rulsif_score(x_scaled, 1 - y, test_idx, train_idx, alpha, l2, score_left)
    score = score_right + score_left

    deallocate(x, y, mu, sig, x_scaled, idx0, idx1, perm0, perm1, train_idx, test_idx)
end subroutine rulsif_linear_pesym_score

subroutine one_side_rulsif_score(x, y, train_idx, eval_idx, alpha, l2, score)
    !> One directional PE score from a linear RuLSIF fit.
    real(kind=dp), intent(in) :: x(:, :), alpha, l2
    integer, intent(in) :: y(:), train_idx(:), eval_idx(:)
    real(kind=dp), intent(out) :: score

    integer :: p, ntrain, neval, i, j, k, n0, n1
    real(kind=dp), allocatable :: phi_train(:, :), phi_eval(:, :), phi0(:, :), phi1(:, :), h(:), hmat(:, :), theta(:), ratios(:)

    p = size(x, 2)
    ntrain = size(train_idx)
    neval = size(eval_idx)
    allocate(phi_train(ntrain, p + 1), phi_eval(neval, p + 1))
    phi_train(:, 1) = 1.0_dp
    phi_eval(:, 1) = 1.0_dp
    do i = 1, ntrain
        phi_train(i, 2:p + 1) = x(train_idx(i), :)
    end do
    do i = 1, neval
        phi_eval(i, 2:p + 1) = x(eval_idx(i), :)
    end do

    n0 = count(y(train_idx) == 0)
    n1 = count(y(train_idx) == 1)
    allocate(phi0(n0, p + 1), phi1(n1, p + 1))
    j = 0
    k = 0
    do i = 1, ntrain
        if (y(train_idx(i)) == 0) then
            j = j + 1
            phi0(j, :) = phi_train(i, :)
        else
            k = k + 1
            phi1(k, :) = phi_train(i, :)
        end if
    end do

    allocate(hmat(p + 1, p + 1), h(p + 1), theta(p + 1))
    hmat = (1.0_dp - alpha) / real(n0, dp) * matmul(transpose(phi0), phi0) + &
           alpha / real(n1, dp) * matmul(transpose(phi1), phi1)
    do i = 1, p + 1
        hmat(i, i) = hmat(i, i) + l2
    end do
    h = 0.0_dp
    do i = 1, n1
        h = h + phi1(i, :)
    end do
    h = h / real(n1, dp)
    call solve_linear_system(hmat, h, theta)

    allocate(ratios(neval))
    ratios = matmul(phi_eval, theta)
    where (ratios < 0.0_dp) ratios = 0.0_dp
    score = 0.5_dp * sum(ratios, mask=y(eval_idx) == 1) / real(count(y(eval_idx) == 1), dp) - 0.5_dp

    deallocate(phi_train, phi_eval, phi0, phi1, hmat, h, theta, ratios)
end subroutine one_side_rulsif_score

subroutine solve_linear_system(a, b, x)
    !> Solve A x = b by Gauss-Jordan elimination.
    real(kind=dp), intent(in) :: a(:, :), b(:)
    real(kind=dp), intent(out) :: x(:)

    integer :: n, i, j
    real(kind=dp), allocatable :: aug(:, :)
    real(kind=dp) :: pivot, factor

    n = size(b)
    allocate(aug(n, n + 1))
    aug(:, 1:n) = a
    aug(:, n + 1) = b
    do i = 1, n
        pivot = aug(i, i)
        if (abs(pivot) < 1.0e-12_dp) pivot = sign(1.0e-12_dp, pivot + 1.0e-12_dp)
        aug(i, :) = aug(i, :) / pivot
        do j = 1, n
            if (j == i) cycle
            factor = aug(j, i)
            aug(j, :) = aug(j, :) - factor * aug(i, :)
        end do
    end do
    x = aug(:, n + 1)
    deallocate(aug)
end subroutine solve_linear_system

subroutine class_mean_cov(x, y, train_idx, class_label, reg_param, mu, cov)
    !> Class-specific mean and covariance with sklearn-style QDA regularization.
    real(kind=dp), intent(in) :: x(:, :), reg_param
    integer, intent(in) :: y(:), train_idx(:), class_label
    real(kind=dp), intent(out) :: mu(:), cov(:, :)

    integer :: p, i, j, k, ncls
    real(kind=dp), allocatable :: centered(:,:)

    p = size(x, 2)
    ncls = count(y(train_idx) == class_label)
    allocate(centered(ncls, p))
    mu = 0.0_dp
    k = 0
    do i = 1, size(train_idx)
        if (y(train_idx(i)) /= class_label) cycle
        k = k + 1
        mu = mu + x(train_idx(i), :)
        centered(k, :) = x(train_idx(i), :)
    end do
    mu = mu / real(ncls, dp)
    do i = 1, ncls
        centered(i, :) = centered(i, :) - mu
    end do
    cov = 0.0_dp
    do i = 1, ncls
        do j = 1, p
            do k = 1, p
                cov(j, k) = cov(j, k) + centered(i, j) * centered(i, k)
            end do
        end do
    end do
    if (ncls > 1) then
        cov = cov / real(ncls - 1, dp)
    else
        cov = 0.0_dp
    end if
    cov = (1.0_dp - reg_param) * cov
    do i = 1, p
        cov(i, i) = cov(i, i) + reg_param
    end do
    deallocate(centered)
end subroutine class_mean_cov

subroutine invert_spd(a, inv_a, det_a)
    !> Invert a small SPD matrix with Gauss-Jordan elimination and return its determinant.
    real(kind=dp), intent(in) :: a(:, :)
    real(kind=dp), intent(out) :: inv_a(:, :)
    real(kind=dp), intent(out) :: det_a

    integer :: n, i, j
    real(kind=dp), allocatable :: aug(:, :)
    real(kind=dp) :: pivot, factor

    n = size(a, 1)
    allocate(aug(n, 2 * n))
    aug(:, 1:n) = a
    aug(:, n + 1:2 * n) = 0.0_dp
    do i = 1, n
        aug(i, n + i) = 1.0_dp
    end do

    det_a = 1.0_dp
    do i = 1, n
        pivot = aug(i, i)
        if (abs(pivot) < 1.0e-12_dp) pivot = sign(1.0e-12_dp, pivot + 1.0e-12_dp)
        det_a = det_a * pivot
        aug(i, :) = aug(i, :) / pivot
        do j = 1, n
            if (j == i) cycle
            factor = aug(j, i)
            aug(j, :) = aug(j, :) - factor * aug(i, :)
        end do
    end do

    inv_a = aug(:, n + 1:2 * n)
    deallocate(aug)
end subroutine invert_spd

real(kind=dp) function quad_form(a, x) result(val)
    real(kind=dp), intent(in) :: a(:, :), x(:)
    val = dot_product(x, matmul(a, x))
end function quad_form

real(kind=dp) function klsym_metric(ref_preds, test_preds) result(score)
    real(kind=dp), intent(in) :: ref_preds(:), test_preds(:)
    score = mean_log(test_preds + 1.0e-3_dp) - mean_log(1.0_dp - test_preds + 1.0e-3_dp) + &
            mean_log(1.0_dp - ref_preds + 1.0e-3_dp) - mean_log(ref_preds + 1.0e-3_dp)
end function klsym_metric

real(kind=dp) function mean_log(x) result(val)
    real(kind=dp), intent(in) :: x(:)
    val = sum(log(x)) / real(size(x), dp)
end function mean_log

subroutine random_permutation_indices_pm(n, seed, perm)
    integer, intent(in) :: n, seed
    integer, allocatable, intent(out) :: perm(:)
    integer :: i, j, state

    allocate(perm(n))
    do i = 1, n
        perm(i) = i
    end do
    state = modulo(abs(seed), 2147483646) + 1
    do i = n, 2, -1
        state = park_miller_next(state)
        j = 1 + modulo(state - 1, i)
        call swap_int_local(perm(i), perm(j))
    end do
end subroutine random_permutation_indices_pm

integer function park_miller_next(state) result(next_state)
    integer, intent(in) :: state
    integer, parameter :: a = 16807
    integer, parameter :: m = 2147483647
    integer, parameter :: q = 127773
    integer, parameter :: r = 2836
    integer :: hi, lo, test

    hi = state / q
    lo = modulo(state, q)
    test = a * lo - r * hi
    if (test > 0) then
        next_state = test
    else
        next_state = test + m
    end if
end function park_miller_next

subroutine swap_int_local(a, b)
    integer, intent(inout) :: a, b
    integer :: tmp
    tmp = a
    a = b
    b = tmp
end subroutine swap_int_local

end module roerich_mod
