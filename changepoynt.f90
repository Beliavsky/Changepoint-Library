module changepoynt_mod
use kind_mod, only: dp
use pca_jacobi_mod, only: jacobi_eigen_sym
implicit none
private
public :: solve_sst_naive_1d, solve_esst_exact_1d, solve_rulsif_gaussian_1d, solve_bocpd_gaussian_mean_1d

contains

subroutine solve_sst_naive_1d(z, window_length, n_windows, lag, rank, score, cp, best_score, scoring_step, scale_signal)
    !> Score a univariate series with changepoynt SST using the exact slow "naive" subspace-overlap path.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_length
    integer, intent(in) :: n_windows
    integer, intent(in) :: lag
    integer, intent(in) :: rank
    real(kind=dp), intent(out) :: score(size(z))
    integer, intent(out) :: cp
    real(kind=dp), intent(out) :: best_score
    integer, intent(in), optional :: scoring_step
    logical, intent(in), optional :: scale_signal

    integer :: n
    integer :: step
    integer :: start_idx
    integer :: offset
    integer :: idx
    integer :: left_idx
    integer :: right_idx
    integer :: cp_pos(1)
    logical :: do_scale
    real(kind=dp), allocatable :: x(:)
    real(kind=dp), allocatable :: hankel_past(:, :)
    real(kind=dp), allocatable :: hankel_future(:, :)
    real(kind=dp) :: local_score

    n = size(z)
    score = 0.0_dp
    cp = -1
    best_score = -huge(1.0_dp)

    if (window_length <= 0 .or. n_windows <= 0 .or. lag <= 0 .or. rank <= 0) return
    if (rank > min(window_length, n_windows)) return

    step = 1
    if (present(scoring_step)) step = scoring_step
    do_scale = .true.
    if (present(scale_signal)) do_scale = scale_signal

    start_idx = window_length + n_windows + lag
    if (start_idx >= n) return
    offset = n_windows / 2 + lag

    allocate(x(n))
    x = z
    if (do_scale) call min_max_scale_1d(x, 1.0_dp, 2.0_dp)

    allocate(hankel_past(window_length, n_windows), hankel_future(window_length, n_windows))
    do idx = start_idx, n - 1, step
        call compile_hankel_1d(x, idx - lag, window_length, n_windows, hankel_past)
        call compile_hankel_1d(x, idx, window_length, n_windows, hankel_future)
        call sst_naive_score(hankel_past, hankel_future, rank, local_score)

        left_idx = idx - offset - step / 2 + 1
        right_idx = idx - offset + (step + 1) / 2
        left_idx = max(left_idx, 1)
        right_idx = min(right_idx, n)
        if (left_idx <= right_idx) score(left_idx:right_idx) = local_score
    end do

    cp_pos = maxloc(score)
    cp = cp_pos(1) - 1
    best_score = score(cp_pos(1))

    deallocate(hankel_past, hankel_future, x)
end subroutine solve_sst_naive_1d

subroutine solve_esst_exact_1d(z, window_length, n_windows, lag, rank, score, cp, best_score, scoring_step, scale_signal)
    !> Score a univariate series with a deterministic exact-SVD ESST analog matching the package score formula.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_length
    integer, intent(in) :: n_windows
    integer, intent(in) :: lag
    integer, intent(in) :: rank
    real(kind=dp), intent(out) :: score(size(z))
    integer, intent(out) :: cp
    real(kind=dp), intent(out) :: best_score
    integer, intent(in), optional :: scoring_step
    logical, intent(in), optional :: scale_signal

    integer :: n
    integer :: step
    integer :: start_idx
    integer :: offset
    integer :: idx
    integer :: left_idx
    integer :: right_idx
    integer :: cp_pos(1)
    logical :: do_scale
    real(kind=dp), allocatable :: x(:)
    real(kind=dp), allocatable :: hankel_past(:, :)
    real(kind=dp), allocatable :: hankel_future(:, :)
    real(kind=dp), allocatable :: hankel(:, :)
    real(kind=dp) :: local_score

    n = size(z)
    score = 0.0_dp
    cp = -1
    best_score = -huge(1.0_dp)

    if (window_length <= 0 .or. n_windows <= 0 .or. lag <= 0 .or. rank <= 0) return

    step = 1
    if (present(scoring_step)) step = scoring_step
    do_scale = .true.
    if (present(scale_signal)) do_scale = scale_signal

    start_idx = window_length + n_windows + lag
    if (start_idx >= n) return
    offset = n_windows + lag

    allocate(x(n))
    x = z
    if (do_scale) call min_max_scale_1d(x, 1.0_dp, 2.0_dp)

    allocate(hankel_past(window_length, n_windows), hankel_future(window_length, n_windows), &
             hankel(window_length, 2 * n_windows))
    do idx = start_idx, n - 1, step
        call compile_hankel_1d(x, idx - lag, window_length, n_windows, hankel_past)
        call compile_hankel_1d(x, idx, window_length, n_windows, hankel_future)
        hankel(:, 1:n_windows) = hankel_past
        hankel(:, n_windows + 1:2 * n_windows) = hankel_future
        call esst_exact_score(hankel, rank, local_score)

        left_idx = idx - offset - step / 2 + 1
        right_idx = idx - offset + (step + 1) / 2
        left_idx = max(left_idx, 1)
        right_idx = min(right_idx, n)
        if (left_idx <= right_idx) score(left_idx:right_idx) = local_score
    end do

    cp_pos = maxloc(score)
    cp = cp_pos(1) - 1
    best_score = score(cp_pos(1))

    deallocate(hankel_past, hankel_future, hankel, x)
end subroutine solve_esst_exact_1d

subroutine solve_rulsif_gaussian_1d(z, window_length, n_windows, lag, alpha, sigma, lambda_reg, score, cp, best_score, &
                                    scoring_step, symmetric)
    !> Score a univariate series with a deterministic Gaussian-kernel RuLSIF PE-divergence detector.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_length
    integer, intent(in) :: n_windows
    integer, intent(in) :: lag
    real(kind=dp), intent(in) :: alpha
    real(kind=dp), intent(in) :: sigma
    real(kind=dp), intent(in) :: lambda_reg
    real(kind=dp), intent(out) :: score(size(z))
    integer, intent(out) :: cp
    real(kind=dp), intent(out) :: best_score
    integer, intent(in), optional :: scoring_step
    logical, intent(in), optional :: symmetric

    integer :: step
    logical :: do_symmetric
    real(kind=dp), allocatable :: score_forward(:), score_backward(:), z_rev(:)
    integer :: cp_pos(1)

    step = 1
    if (present(scoring_step)) step = scoring_step
    do_symmetric = .true.
    if (present(symmetric)) do_symmetric = symmetric

    allocate(score_forward(size(z)))
    call rulsif_gaussian_one_direction(z, window_length, n_windows, lag, alpha, sigma, lambda_reg, score_forward, step)

    if (do_symmetric) then
        allocate(score_backward(size(z)), z_rev(size(z)))
        z_rev = z(size(z):1:-1)
        call rulsif_gaussian_one_direction(z_rev, window_length, n_windows, lag, alpha, sigma, lambda_reg, score_backward, step)
        score = score_forward + score_backward(size(z):1:-1)
        deallocate(score_backward, z_rev)
    else
        score = score_forward
    end if

    cp_pos = maxloc(score)
    cp = cp_pos(1) - 1
    best_score = score(cp_pos(1))
    deallocate(score_forward)
end subroutine solve_rulsif_gaussian_1d

subroutine solve_bocpd_gaussian_mean_1d(x, run_length, score, cp_raw, cp_post, best_score, prior_mean, prior_var, signal_var, &
                                        change_length_threshold)
    !> Exact replication of changepoynt.BOCPD with GaussianUnknownMean and constant hazard 1/run_length.
    real(kind=dp), intent(in) :: x(:)
    integer, intent(in) :: run_length
    real(kind=dp), intent(out) :: score(size(x))
    integer, intent(out) :: cp_raw
    integer, intent(out) :: cp_post
    real(kind=dp), intent(out) :: best_score
    real(kind=dp), intent(out) :: prior_mean
    real(kind=dp), intent(out) :: prior_var
    real(kind=dp), intent(out) :: signal_var
    integer, intent(in), optional :: change_length_threshold

    integer :: n
    integer :: n_windows
    integer :: t
    integer :: thr
    integer :: imax(1)
    integer :: i
    integer :: ipost(1)
    real(kind=dp) :: hazard
    real(kind=dp) :: log_hazard
    real(kind=dp) :: log_one_minus_hazard
    real(kind=dp) :: prior_prec
    real(kind=dp) :: log_cp_prob
    real(kind=dp) :: norm_const
    real(kind=dp), parameter :: tiny_var = 1.0e-12_dp
    real(kind=dp), parameter :: neg_inf = -huge(1.0_dp) / 8.0_dp
    real(kind=dp), allocatable :: prefix(:)
    real(kind=dp), allocatable :: prefix2(:)
    real(kind=dp), allocatable :: window_means(:)
    real(kind=dp), allocatable :: window_vars(:)
    real(kind=dp), allocatable :: mean_params(:)
    real(kind=dp), allocatable :: prec_params(:)
    real(kind=dp), allocatable :: post_var(:)
    real(kind=dp), allocatable :: log_message(:)
    real(kind=dp), allocatable :: log_pis(:)
    real(kind=dp), allocatable :: log_growth(:)
    real(kind=dp), allocatable :: new_log_joint(:)
    real(kind=dp), allocatable :: new_prec(:)
    real(kind=dp), allocatable :: new_mean(:)
    real(kind=dp), allocatable :: row_prob(:)

    n = size(x)
    score = 0.0_dp
    cp_raw = -1
    cp_post = -1
    best_score = 0.0_dp
    prior_mean = 0.0_dp
    prior_var = 0.0_dp
    signal_var = 0.0_dp
    if (run_length <= 0 .or. n < run_length) return

    n_windows = n - run_length + 1
    allocate(prefix(0:n), prefix2(0:n), window_means(n_windows), window_vars(n_windows))
    prefix(0) = 0.0_dp
    prefix2(0) = 0.0_dp
    do t = 1, n
        prefix(t) = prefix(t - 1) + x(t)
        prefix2(t) = prefix2(t - 1) + x(t) * x(t)
    end do
    do t = 1, n_windows
        call sliding_mean_var(prefix, prefix2, t, run_length, window_means(t), window_vars(t))
    end do

    prior_mean = median_copy(window_means)
    prior_var = population_variance(window_means)
    signal_var = median_copy(window_vars)
    prior_var = max(prior_var, tiny_var)
    signal_var = max(signal_var, tiny_var)
    prior_prec = 1.0_dp / prior_var

    hazard = 1.0_dp / real(run_length, dp)
    log_hazard = log(hazard)
    log_one_minus_hazard = log(1.0_dp - hazard)
    thr = max(0, int(real(run_length, dp) * 0.1_dp))
    if (present(change_length_threshold)) thr = max(0, change_length_threshold)

    allocate(mean_params(n + 1), prec_params(n + 1), post_var(n + 1))
    allocate(log_message(1))
    mean_params = 0.0_dp
    prec_params = 0.0_dp
    mean_params(1) = prior_mean
    prec_params(1) = prior_prec
    log_message(1) = 0.0_dp

    do t = 1, n - 1
        post_var(1:t) = 1.0_dp / prec_params(1:t) + signal_var

        allocate(log_pis(t), log_growth(t), new_log_joint(t + 1), row_prob(t + 1))
        do i = 1, t
            log_pis(i) = normal_logpdf(x(t), mean_params(i), sqrt(post_var(i)))
        end do
        log_growth = log_pis + log_message + log_one_minus_hazard
        log_cp_prob = logsumexp_vec(log_pis + log_message + log_hazard)
        new_log_joint(1) = log_cp_prob
        new_log_joint(2:t + 1) = log_growth
        norm_const = logsumexp_vec(new_log_joint)
        row_prob = exp(new_log_joint - norm_const)
        score(t) = sum(row_prob(1:min(thr + 1, size(row_prob))))

        allocate(new_prec(t), new_mean(t))
        new_prec = prec_params(1:t) + 1.0_dp / signal_var
        new_mean = (mean_params(1:t) * prec_params(1:t) + x(t) / signal_var) / new_prec
        prec_params(2:t + 1) = new_prec
        prec_params(1) = prior_prec
        mean_params(2:t + 1) = new_mean
        mean_params(1) = prior_mean

        deallocate(log_message)
        allocate(log_message(t + 1))
        log_message = new_log_joint

        deallocate(log_pis, log_growth, new_log_joint, new_prec, new_mean, row_prob)
    end do
    score(n) = 0.0_dp

    imax = maxloc(score)
    cp_raw = imax(1) - 1
    best_score = score(imax(1))
    if (run_length <= n) then
        ipost = maxloc(score(run_length:n))
        cp_post = run_length + ipost(1) - 2
    end if

    deallocate(prefix, prefix2, window_means, window_vars, mean_params, prec_params, post_var, log_message)
end subroutine solve_bocpd_gaussian_mean_1d

subroutine min_max_scale_1d(x, min_val, max_val)
    !> Apply the package's min-max scaling convention to a 1D series.
    real(kind=dp), intent(in out) :: x(:)
    real(kind=dp), intent(in) :: min_val
    real(kind=dp), intent(in) :: max_val

    real(kind=dp) :: xmin, xmax

    xmin = minval(x)
    xmax = maxval(x)
    if (xmax == xmin) then
        x = x - xmin
    else
        x = (x - xmin) / (xmax - xmin)
    end if
    x = x * (max_val - min_val) + min_val
end subroutine min_max_scale_1d

subroutine compile_hankel_1d(x, end_index, window_length, n_windows, hankel)
    !> Build the slow explicit Hankel matrix with the same column ordering as changepoynt.utils.linalg.compile_hankel.
    real(kind=dp), intent(in) :: x(:)
    integer, intent(in) :: end_index
    integer, intent(in) :: window_length
    integer, intent(in) :: n_windows
    real(kind=dp), intent(out) :: hankel(window_length, n_windows)

    integer :: cx
    integer :: col
    integer :: i1
    integer :: i2

    do cx = 0, n_windows - 1
        col = n_windows - cx
        i1 = end_index - window_length - cx + 1
        i2 = end_index - cx
        hankel(:, col) = x(i1:i2)
    end do
end subroutine compile_hankel_1d

subroutine rulsif_gaussian_one_direction(z, window_length, n_windows, lag, alpha, sigma, lambda_reg, score, scoring_step)
    !> Compute one directional RuLSIF score path.
    real(kind=dp), intent(in) :: z(:)
    integer, intent(in) :: window_length
    integer, intent(in) :: n_windows
    integer, intent(in) :: lag
    real(kind=dp), intent(in) :: alpha
    real(kind=dp), intent(in) :: sigma
    real(kind=dp), intent(in) :: lambda_reg
    real(kind=dp), intent(out) :: score(size(z))
    integer, intent(in) :: scoring_step

    integer :: idx
    integer :: start_idx
    integer :: offset
    integer :: left_idx
    integer :: right_idx
    real(kind=dp), allocatable :: hankel(:, :)
    real(kind=dp) :: local_score

    score = 0.0_dp
    start_idx = window_length + n_windows + lag
    if (start_idx >= size(z)) return
    offset = n_windows

    allocate(hankel(window_length, 2 * n_windows))
    do idx = start_idx, size(z) - 1, scoring_step
        call compile_hankel_1d(z, idx, window_length, 2 * n_windows, hankel)
        call rulsif_gaussian_pair_score(hankel(:, 1:n_windows), hankel(:, n_windows + 1:2 * n_windows), alpha, sigma, lambda_reg, local_score)
        left_idx = idx - offset - scoring_step / 2 + 1
        right_idx = idx - offset + (scoring_step + 1) / 2
        left_idx = max(left_idx, 1)
        right_idx = min(right_idx, size(z))
        if (left_idx <= right_idx) score(left_idx:right_idx) = local_score
    end do
    deallocate(hankel)
end subroutine rulsif_gaussian_one_direction

subroutine rulsif_gaussian_pair_score(reference_samples, test_samples, alpha, sigma, lambda_reg, score)
    !> Deterministic one-pair RuLSIF score using reference columns as Gaussian centers.
    real(kind=dp), intent(in) :: reference_samples(:, :)
    real(kind=dp), intent(in) :: test_samples(:, :)
    real(kind=dp), intent(in) :: alpha
    real(kind=dp), intent(in) :: sigma
    real(kind=dp), intent(in) :: lambda_reg
    real(kind=dp), intent(out) :: score

    integer :: d, n_ref, n_test, n_centers, i
    real(kind=dp), allocatable :: ref_scaled(:, :), test_scaled(:, :), all_samples(:, :)
    real(kind=dp), allocatable :: std_row(:), centers(:, :)
    real(kind=dp), allocatable :: k_ref(:, :), k_test(:, :)
    real(kind=dp), allocatable :: h_hat(:), hmat(:, :), theta(:), ratios_test(:)

    d = size(reference_samples, 1)
    n_ref = size(reference_samples, 2)
    n_test = size(test_samples, 2)
    n_centers = n_ref

    allocate(ref_scaled(d, n_ref), test_scaled(d, n_test), all_samples(d, n_ref + n_test), std_row(d))
    all_samples(:, 1:n_ref) = reference_samples
    all_samples(:, n_ref + 1:n_ref + n_test) = test_samples
    do i = 1, d
        std_row(i) = sqrt(sum((all_samples(i, :) - sum(all_samples(i, :)) / real(n_ref + n_test, dp)) ** 2) / real(n_ref + n_test, dp))
        std_row(i) = std_row(i) + epsilon(1.0_dp)
        if (std_row(i) < 1.0e-12_dp) std_row(i) = 1.0_dp
        ref_scaled(i, :) = reference_samples(i, :) / std_row(i)
        test_scaled(i, :) = test_samples(i, :) / std_row(i)
    end do

    allocate(centers(d, n_centers))
    centers = ref_scaled
    allocate(k_ref(n_centers, n_ref), k_test(n_centers, n_test))
    call gaussian_kernel_matrix(ref_scaled, centers, sigma, k_ref)
    call gaussian_kernel_matrix(test_scaled, centers, sigma, k_test)

    allocate(hmat(n_centers, n_centers), h_hat(n_centers), theta(n_centers))
    hmat = alpha / real(n_ref, dp) * matmul(k_ref, transpose(k_ref)) + &
           (1.0_dp - alpha) / real(n_test, dp) * matmul(k_test, transpose(k_test))
    do i = 1, n_centers
        hmat(i, i) = hmat(i, i) + lambda_reg
    end do
    h_hat = sum(k_ref, dim=2) / real(n_ref, dp)
    call solve_linear_system_local(hmat, h_hat, theta)

    allocate(ratios_test(n_test))
    ratios_test = matmul(transpose(k_test), theta)
    do i = 1, n_test
        if (ratios_test(i) < 0.0_dp) ratios_test(i) = 0.0_dp
    end do
    score = 0.5_dp * sum(ratios_test) / real(n_test, dp) - 0.5_dp

    deallocate(ref_scaled, test_scaled, all_samples, std_row, centers, k_ref, k_test, hmat, h_hat, theta, ratios_test)
end subroutine rulsif_gaussian_pair_score

subroutine gaussian_kernel_matrix(samples, centers, sigma, kernel)
    !> Build the Gaussian kernel matrix K(samples, centers) with center-major row layout.
    real(kind=dp), intent(in) :: samples(:, :)
    real(kind=dp), intent(in) :: centers(:, :)
    real(kind=dp), intent(in) :: sigma
    real(kind=dp), intent(out) :: kernel(size(centers, 2), size(samples, 2))

    integer :: i, j
    real(kind=dp) :: sqdist

    do i = 1, size(centers, 2)
        do j = 1, size(samples, 2)
            sqdist = sum((samples(:, j) - centers(:, i)) ** 2)
            kernel(i, j) = exp(-sqdist / (2.0_dp * sigma * sigma))
        end do
    end do
end subroutine gaussian_kernel_matrix

subroutine sliding_mean_var(prefix, prefix2, start_idx, window_length, mean_val, var_val)
    !> Population mean and variance for one sliding window via prefix sums.
    real(kind=dp), intent(in) :: prefix(0:)
    real(kind=dp), intent(in) :: prefix2(0:)
    integer, intent(in) :: start_idx
    integer, intent(in) :: window_length
    real(kind=dp), intent(out) :: mean_val
    real(kind=dp), intent(out) :: var_val

    integer :: j1, j2
    real(kind=dp) :: s1, s2

    j1 = start_idx - 1
    j2 = start_idx + window_length - 1
    s1 = prefix(j2) - prefix(j1)
    s2 = prefix2(j2) - prefix2(j1)
    mean_val = s1 / real(window_length, dp)
    var_val = s2 / real(window_length, dp) - mean_val * mean_val
    if (var_val < 0.0_dp .and. abs(var_val) < 1.0e-12_dp) var_val = 0.0_dp
end subroutine sliding_mean_var

real(kind=dp) function median_copy(x) result(med)
    !> Median of a vector using a sorted copy.
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp), allocatable :: work(:)
    integer :: n

    n = size(x)
    allocate(work(n))
    work = x
    call sort_real_ascending(work)
    if (mod(n, 2) == 1) then
        med = work((n + 1) / 2)
    else
        med = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
    end if
    deallocate(work)
end function median_copy

real(kind=dp) function population_variance(x) result(var_val)
    !> Population variance matching numpy.var default ddof=0.
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: mu

    mu = sum(x) / real(size(x), dp)
    var_val = sum((x - mu) ** 2) / real(size(x), dp)
end function population_variance

subroutine sort_real_ascending(x)
    !> In-place ascending sort for small vectors.
    real(kind=dp), intent(in out) :: x(:)
    integer :: i, j, k, n
    real(kind=dp) :: tmp

    n = size(x)
    do i = 1, n - 1
        k = i
        do j = i + 1, n
            if (x(j) < x(k)) k = j
        end do
        if (k /= i) then
            tmp = x(i)
            x(i) = x(k)
            x(k) = tmp
        end if
    end do
end subroutine sort_real_ascending

real(kind=dp) function normal_logpdf(x, mu, sigma) result(lp)
    !> Univariate Gaussian log-density.
    real(kind=dp), intent(in) :: x
    real(kind=dp), intent(in) :: mu
    real(kind=dp), intent(in) :: sigma
    real(kind=dp) :: z

    z = (x - mu) / sigma
    lp = -0.5_dp * log(2.0_dp * acos(-1.0_dp)) - log(sigma) - 0.5_dp * z * z
end function normal_logpdf

real(kind=dp) function logsumexp_vec(x) result(val)
    !> Stable log-sum-exp for a vector.
    real(kind=dp), intent(in) :: x(:)
    real(kind=dp) :: xmax

    xmax = maxval(x)
    if (xmax < -huge(1.0_dp) / 16.0_dp) then
        val = xmax
    else
        val = xmax + log(sum(exp(x - xmax)))
    end if
end function logsumexp_vec

subroutine solve_linear_system_local(a, b, x)
    !> Solve A x = b by Gauss-Jordan elimination for small dense systems.
    real(kind=dp), intent(in) :: a(:, :)
    real(kind=dp), intent(in) :: b(:)
    real(kind=dp), intent(out) :: x(:)

    integer :: n, i, j, pivot_row
    real(kind=dp), allocatable :: aug(:, :)
    real(kind=dp) :: pivot, factor, best
    real(kind=dp), allocatable :: tmp(:)

    n = size(b)
    allocate(aug(n, n + 1), tmp(n + 1))
    aug(:, 1:n) = a
    aug(:, n + 1) = b
    do i = 1, n
        pivot_row = i
        best = abs(aug(i, i))
        do j = i + 1, n
            if (abs(aug(j, i)) > best) then
                best = abs(aug(j, i))
                pivot_row = j
            end if
        end do
        if (pivot_row /= i) then
            tmp = aug(i, :)
            aug(i, :) = aug(pivot_row, :)
            aug(pivot_row, :) = tmp
        end if
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
    deallocate(aug, tmp)
end subroutine solve_linear_system_local

subroutine sst_naive_score(hankel_past, hankel_future, rank, score)
    !> Compute the naive SST score 1 - s_max(U_p^T U_f).
    real(kind=dp), intent(in) :: hankel_past(:, :)
    real(kind=dp), intent(in) :: hankel_future(:, :)
    integer, intent(in) :: rank
    real(kind=dp), intent(out) :: score

    real(kind=dp), allocatable :: cov_past(:, :)
    real(kind=dp), allocatable :: cov_future(:, :)
    real(kind=dp), allocatable :: evals_past(:), evals_future(:)
    real(kind=dp), allocatable :: evecs_past(:, :), evecs_future(:, :)
    real(kind=dp), allocatable :: cross(:, :)
    real(kind=dp), allocatable :: gram(:, :)
    real(kind=dp), allocatable :: gram_evals(:), gram_evecs(:, :)
    integer :: r
    real(kind=dp) :: sigma_max

    r = min(rank, min(size(hankel_past, 1), size(hankel_past, 2)))

    allocate(cov_past(size(hankel_past, 1), size(hankel_past, 1)))
    allocate(cov_future(size(hankel_future, 1), size(hankel_future, 1)))
    cov_past = matmul(hankel_past, transpose(hankel_past))
    cov_future = matmul(hankel_future, transpose(hankel_future))

    call jacobi_eigen_sym(cov_past, evals_past, evecs_past)
    call jacobi_eigen_sym(cov_future, evals_future, evecs_future)

    allocate(cross(r, r), gram(r, r))
    cross = matmul(transpose(evecs_past(:, 1:r)), evecs_future(:, 1:r))
    gram = matmul(transpose(cross), cross)
    call jacobi_eigen_sym(gram, gram_evals, gram_evecs)

    sigma_max = sqrt(max(gram_evals(1), 0.0_dp))
    score = 1.0_dp - sigma_max
    if (score < 0.0_dp .and. abs(score) < 1.0e-12_dp) score = 0.0_dp

    deallocate(cov_past, cov_future, evals_past, evals_future, evecs_past, evecs_future, cross, gram, gram_evals, gram_evecs)
end subroutine sst_naive_score

subroutine esst_exact_score(hankel, rank, score)
    !> Compute the deterministic ESST score from the top right singular vectors of the combined Hankel matrix.
    real(kind=dp), intent(in) :: hankel(:, :)
    integer, intent(in) :: rank
    real(kind=dp), intent(out) :: score

    real(kind=dp), allocatable :: cov_right(:, :)
    real(kind=dp), allocatable :: evals(:), evecs(:, :)
    real(kind=dp), allocatable :: right_vectors(:, :)
    real(kind=dp), allocatable :: singular_values(:)
    real(kind=dp), allocatable :: skew(:)
    real(kind=dp) :: wsum
    integer :: r
    integer :: i
    integer :: nhalf

    r = min(rank, min(size(hankel, 1), size(hankel, 2)))
    nhalf = size(hankel, 2) / 2

    allocate(cov_right(size(hankel, 2), size(hankel, 2)))
    cov_right = matmul(transpose(hankel), hankel)
    call jacobi_eigen_sym(cov_right, evals, evecs)

    allocate(right_vectors(r, size(hankel, 2)), singular_values(r), skew(r))
    right_vectors = transpose(evecs(:, 1:r))
    singular_values = sqrt(max(evals(1:r), 0.0_dp))

    do i = 1, r
        right_vectors(i, :) = right_vectors(i, :) - minval(right_vectors(i, :)) + 1.0_dp
        right_vectors(i, :) = right_vectors(i, :) / sum(right_vectors(i, :))
        skew(i) = abs(sum(right_vectors(i, 1:nhalf) - right_vectors(i, nhalf + 1:2 * nhalf)) / real(nhalf, dp))
    end do

    wsum = sum(singular_values)
    if (wsum > 0.0_dp) then
        score = sum(singular_values * skew) / wsum
    else
        score = 0.0_dp
    end if

    deallocate(cov_right, evals, evecs, right_vectors, singular_values, skew)
end subroutine esst_exact_score

end module changepoynt_mod
