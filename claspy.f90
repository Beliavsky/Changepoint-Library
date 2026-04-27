module claspy_mod
    use kind_mod, only: dp
    use util_mod, only: sort_int
    implicit none
    private
    public :: solve_binary_clasp_1d, solve_binary_clasp_changepoynt_1d, solve_clasp_ensemble_1d, &
              solve_agglomerative_clap_1d, collapse_segment_process_1d, solve_class_1d, &
              solve_streaming_class_1d, solve_clap_centroid_1d

contains

    subroutine solve_binary_clasp_1d(z, window_size, n_segments, cps, k_neighbours, excl_radius)
        !> Fixed-K binary segmentation using a univariate ClaSP profile with z-normalized distances.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: window_size, n_segments
        integer, allocatable, intent(out) :: cps(:)
        integer, intent(in), optional :: k_neighbours, excl_radius

        integer :: n, k, excl, min_seg_size, queue_cap, n_queue, ncp, idx
        integer, allocatable :: queue_l(:), queue_u(:), queue_cp(:), tmp(:)
        real(kind=dp), allocatable :: queue_score(:)
        integer :: cp_best
        real(kind=dp) :: score_best, priority
        logical :: valid

        n = size(z)
        k = 3
        if (present(k_neighbours)) k = k_neighbours
        excl = 5
        if (present(excl_radius)) excl = excl_radius
        min_seg_size = window_size * excl

        if (n_segments <= 1 .or. n < 2 * min_seg_size) then
            allocate(cps(0))
            return
        end if

        queue_cap = max(8, 4 * n_segments)
        allocate(queue_l(queue_cap), queue_u(queue_cap), queue_cp(queue_cap), queue_score(queue_cap), tmp(n_segments - 1))
        n_queue = 0
        ncp = 0

        call local_clasp_candidate_1d(z, 0, n, window_size, k, min_seg_size, cp_best, score_best, valid)
        if (valid) then
            priority = clasp_priority(score_best, 0, n, n)
            call push_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, queue_cap, 0, n, cp_best, priority)
        end if

        do while (n_queue > 0 .and. ncp < n_segments - 1)
            call pop_best_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, idx)
            ncp = ncp + 1
            tmp(ncp) = queue_cp(idx)

            call local_clasp_candidate_1d(z, queue_l(idx), queue_cp(idx), window_size, k, min_seg_size, cp_best, score_best, valid)
            if (valid) then
                priority = clasp_priority(score_best, queue_l(idx), queue_cp(idx), n)
                call push_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, queue_cap, queue_l(idx), queue_cp(idx), cp_best, priority)
            end if

            call local_clasp_candidate_1d(z, queue_cp(idx), queue_u(idx), window_size, k, min_seg_size, cp_best, score_best, valid)
            if (valid) then
                priority = clasp_priority(score_best, queue_cp(idx), queue_u(idx), n)
                call push_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, queue_cap, queue_cp(idx), queue_u(idx), cp_best, priority)
            end if
        end do

        allocate(cps(ncp))
        if (ncp > 0) then
            cps = tmp(1:ncp)
            if (ncp > 1) call sort_int(cps)
        end if

        deallocate(queue_l, queue_u, queue_cp, queue_score, tmp)
    end subroutine solve_binary_clasp_1d

    subroutine solve_binary_clasp_changepoynt_1d(z, cps, window_size_used, n_segments_used, &
            k_neighbours, excl_radius, n_estimators, threshold, seed)
        !> BinaryClaSPSegmentation close to changepoynt/claspy defaults:
        !> SuSS window selection, learned segment cap, ClaSP ensemble candidates,
        !> and significance-test validation.
        real(kind=dp), intent(in) :: z(:)
        integer, allocatable, intent(out) :: cps(:)
        integer, intent(out) :: window_size_used, n_segments_used
        integer, intent(in), optional :: k_neighbours, excl_radius, n_estimators, seed
        real(kind=dp), intent(in), optional :: threshold

        integer :: n, k, excl, nest, seed_use, min_seg_size, queue_cap, n_queue, ncp, idx
        integer :: cp_best
        integer, allocatable :: queue_l(:), queue_u(:), queue_cp(:), tmp(:)
        real(kind=dp), allocatable :: queue_score(:)
        real(kind=dp) :: score_best, threshold_use
        logical :: valid

        n = size(z)
        k = 3
        if (present(k_neighbours)) k = k_neighbours
        excl = 5
        if (present(excl_radius)) excl = excl_radius
        nest = 10
        if (present(n_estimators)) nest = n_estimators
        seed_use = 2357
        if (present(seed)) seed_use = seed
        threshold_use = 1.0e-15_dp
        if (present(threshold)) threshold_use = threshold

        window_size_used = max(3, suss_window_size_1d(z) / 2)
        min_seg_size = window_size_used * excl

        if (n < 2 * min_seg_size) then
            window_size_used = min(window_size_used, max(1, n / 2))
            n_segments_used = 1
            allocate(cps(0))
            return
        end if

        n_segments_used = n / min_seg_size
        if (n_segments_used <= 1) then
            allocate(cps(0))
            return
        end if

        queue_cap = max(8, 4 * n_segments_used)
        allocate(queue_l(queue_cap), queue_u(queue_cap), queue_cp(queue_cap), queue_score(queue_cap), tmp(n_segments_used - 1))
        n_queue = 0
        ncp = 0

        call clasp_segment_candidate_1d(z, 0, n, window_size_used, k, excl, nest, threshold_use, seed_use, &
            tmp, ncp, cp_best, score_best, valid)
        if (valid) call push_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, queue_cap, 0, n, cp_best, score_best)

        do while (n_queue > 0 .and. ncp < n_segments_used - 1)
            call pop_best_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, idx)
            ncp = ncp + 1
            tmp(ncp) = queue_cp(idx)

            call clasp_segment_candidate_1d(z, queue_l(idx), queue_cp(idx), window_size_used, k, excl, nest, threshold_use, seed_use, &
                tmp, ncp, cp_best, score_best, valid)
            if (valid) then
                call push_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, queue_cap, &
                    queue_l(idx), queue_cp(idx), cp_best, score_best)
            end if

            call clasp_segment_candidate_1d(z, queue_cp(idx), queue_u(idx), window_size_used, k, excl, nest, threshold_use, seed_use, &
                tmp, ncp, cp_best, score_best, valid)
            if (valid) then
                call push_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, queue_cap, &
                    queue_cp(idx), queue_u(idx), cp_best, score_best)
            end if
        end do

        allocate(cps(ncp))
        if (ncp > 0) then
            cps = tmp(1:ncp)
            if (ncp > 1) call sort_int(cps)
        end if

        deallocate(queue_l, queue_u, queue_cp, queue_score, tmp)
    end subroutine solve_binary_clasp_changepoynt_1d

    subroutine solve_clasp_ensemble_1d(z, window_size, cp, score, found, n_estimators, k_neighbours, excl_radius, early_stopping, seed)
        !> Deterministic ClaSPEnsemble-style single-split selection with shared temporal constraints.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: window_size, n_estimators
        integer, intent(out) :: cp
        real(kind=dp), intent(out) :: score
        logical, intent(out) :: found
        integer, intent(in), optional :: k_neighbours, excl_radius, seed
        logical, intent(in), optional :: early_stopping

        integer :: k, excl, seed_use
        logical :: early_stop

        k = 3
        if (present(k_neighbours)) k = k_neighbours
        excl = 5
        if (present(excl_radius)) excl = excl_radius
        seed_use = 2357
        if (present(seed)) seed_use = seed
        early_stop = .true.
        if (present(early_stopping)) early_stop = early_stopping

        cp = -1
        score = -huge(1.0_dp)
        found = .false.
        call ensemble_best_candidate_1d(z, 0, size(z), window_size, n_estimators, k, excl, early_stop, seed_use, cp, score, found)
        if (found) score = clasp_priority(score, 0, size(z), size(z))
    end subroutine solve_clasp_ensemble_1d

    subroutine solve_agglomerative_clap_1d(z, true_bkps, window_size, n_splits, segment_labels, gain, sample_cap, seed)
        !> Agglomerative recurring-state detection using the deterministic CLaP-style classifier.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: true_bkps(:), window_size, n_splits
        integer, allocatable, intent(out) :: segment_labels(:)
        real(kind=dp), intent(out) :: gain
        integer, intent(in), optional :: sample_cap, seed

        integer :: nseg, i, j, cap_use, seed_use, label1, label2, best_i, n_active
        integer, allocatable :: labels(:), state_labels(:), y_true(:), y_pred(:), active_labels(:), &
                                conf(:,:), conf_idx(:), row_sum(:), order(:), y_true_m(:), y_pred_m(:)
        real(kind=dp) :: score, current_gain, new_gain, best_loss
        real(kind=dp), allocatable :: conf_loss(:)
        logical :: merged
        logical, allocatable :: ignore_pair(:,:)

        cap_use = 1000
        if (present(sample_cap)) cap_use = sample_cap
        seed_use = 2357
        if (present(seed)) seed_use = seed

        nseg = size(true_bkps) + 1
        allocate(labels(nseg), state_labels(size(z)), ignore_pair(nseg, nseg))
        do i = 1, nseg
            labels(i) = i
        end do
        ignore_pair = .false.

        merged = .true.
        do while (merged)
            call create_state_labels_from_bkps(true_bkps, labels, size(z), state_labels)
            call solve_clap_centroid_1d(z, state_labels, window_size, n_splits, y_true, y_pred, score, cap_use, seed_use)
            current_gain = classification_gain_labels(y_true, y_pred)
            call unique_labels(labels, active_labels)
            n_active = size(active_labels)
            if (n_active <= 1) exit

            allocate(conf(n_active, n_active), conf_idx(n_active), row_sum(n_active), order(n_active), conf_loss(n_active))
            conf = 0
            do i = 1, size(y_true)
                do j = 1, n_active
                    if (y_true(i) == active_labels(j)) exit
                end do
                if (j > n_active) cycle
                best_i = 1
                do best_i = 1, n_active
                    if (y_pred(i) == active_labels(best_i)) exit
                end do
                if (best_i <= n_active) conf(j, best_i) = conf(j, best_i) + 1
            end do

            do i = 1, n_active
                row_sum(i) = sum(conf(i, :))
                best_loss = -1.0_dp
                conf_idx(i) = i
                do j = 1, n_active
                    if (j == i) cycle
                    if (real(conf(i, j), dp) > best_loss) then
                        best_loss = real(conf(i, j), dp)
                        conf_idx(i) = j
                    end if
                end do
                if (row_sum(i) > 0) then
                    conf_loss(i) = real(conf(i, conf_idx(i)), dp) / real(row_sum(i), dp)
                else
                    conf_loss(i) = 0.0_dp
                end if
                order(i) = i
            end do

            do i = 1, n_active - 1
                do j = i + 1, n_active
                    if (conf_loss(order(j)) > conf_loss(order(i))) call swap_int(order(i), order(j))
                end do
            end do

            merged = .false.
            do i = 1, n_active
                label1 = active_labels(order(i))
                label2 = active_labels(conf_idx(order(i)))
                if (label1 == label2) cycle
                if (ignore_pair(min(label1, label2), max(label1, label2))) cycle

                allocate(y_true_m(size(y_true)), y_pred_m(size(y_pred)))
                y_true_m = y_true
                y_pred_m = y_pred
                where (y_true_m == label2) y_true_m = label1
                where (y_pred_m == label2) y_pred_m = label1
                new_gain = classification_gain_labels(y_true_m, y_pred_m)
                deallocate(y_true_m, y_pred_m)

                if (current_gain > new_gain) then
                    ignore_pair(min(label1, label2), max(label1, label2)) = .true.
                    cycle
                end if

                if (label2 > label1) call swap_int(label1, label2)
                where (labels == label2) labels = label1
                merged = .true.
                exit
            end do

            deallocate(y_true, y_pred, active_labels, conf, conf_idx, row_sum, order, conf_loss)
            if (.not. merged) exit
        end do

        call unique_labels(labels, active_labels)
        do i = 1, size(active_labels)
            where (labels == active_labels(i)) labels = i
        end do

        call create_state_labels_from_bkps(true_bkps, labels, size(z), state_labels)
        call solve_clap_centroid_1d(z, state_labels, window_size, n_splits, y_true, y_pred, score, cap_use, seed_use)
        gain = classification_gain_labels(y_true, y_pred)

        allocate(segment_labels(nseg))
        segment_labels = labels

        deallocate(labels, state_labels, ignore_pair, active_labels, y_true, y_pred)
    end subroutine solve_agglomerative_clap_1d

    subroutine collapse_segment_process_1d(true_bkps, segment_labels, change_points_out, labels_out)
        !> Collapse consecutive equal segment labels into the sparse recurring-state process.
        integer, intent(in) :: true_bkps(:), segment_labels(:)
        integer, allocatable, intent(out) :: change_points_out(:), labels_out(:)

        integer :: i, nchange, pos

        nchange = 0
        do i = 2, size(segment_labels)
            if (segment_labels(i) /= segment_labels(i - 1)) nchange = nchange + 1
        end do

        allocate(change_points_out(nchange), labels_out(nchange + 1))
        labels_out(1) = segment_labels(1)
        pos = 1
        nchange = 0
        do i = 2, size(segment_labels)
            if (segment_labels(i) /= segment_labels(i - 1)) then
                nchange = nchange + 1
                change_points_out(nchange) = true_bkps(i - 1)
                pos = pos + 1
                labels_out(pos) = segment_labels(i)
            end if
        end do
    end subroutine collapse_segment_process_1d

    subroutine local_clasp_candidate_1d(z, lbound, ubound, window_size, k_neighbours, min_seg_size, cp_best, score_best, valid)
        !> Best local ClaSP split on [lbound, ubound) in 0-based raw time indices.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: lbound, ubound, window_size, k_neighbours, min_seg_size
        integer, intent(out) :: cp_best
        real(kind=dp), intent(out) :: score_best
        logical, intent(out) :: valid

        integer :: n_local, m, split_idx, start_cp, end_cp
        real(kind=dp), allocatable :: x(:), means(:), stds(:)
        integer, allocatable :: offsets(:, :), y_true(:), y_pred(:)
        real(kind=dp) :: score

        n_local = ubound - lbound
        if (n_local < 2 * min_seg_size .or. n_local < window_size + 1) then
            valid = .false.
            cp_best = -1
            score_best = -huge(1.0_dp)
            return
        end if

        allocate(x(n_local))
        x = z(lbound + 1:ubound)
        m = n_local - window_size + 1
        start_cp = min_seg_size
        end_cp = m - min_seg_size + window_size - 1
        if (m <= 0 .or. end_cp < start_cp) then
            valid = .false.
            cp_best = -1
            score_best = -huge(1.0_dp)
            deallocate(x)
            return
        end if

        allocate(means(m), stds(m), offsets(m, k_neighbours), y_true(m), y_pred(m))
        call compute_subseq_mean_std_1d(x, window_size, means, stds)
        call knn_offsets_znormed_1d(x, window_size, k_neighbours, means, stds, offsets)

        score_best = -huge(1.0_dp)
        cp_best = start_cp
        do split_idx = start_cp, end_cp
            call cross_val_labels_1d(offsets, split_idx, window_size, y_true, y_pred)
            score = roc_auc_binary_score(y_true, y_pred)
            if (score > score_best) then
                score_best = score
                cp_best = split_idx
            end if
        end do
        cp_best = lbound + cp_best
        valid = .true.

        deallocate(x, means, stds, offsets, y_true, y_pred)
    end subroutine local_clasp_candidate_1d

    subroutine clasp_segment_candidate_1d(z, lbound, ubound, window_size, k_neighbours, excl_radius, n_estimators, threshold, &
            seed, change_points, n_change_points, cp_best, score_best, valid)
        !> One validated BinaryClaSPSegmentation candidate on [lbound, ubound).
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: lbound, ubound, window_size, k_neighbours, excl_radius, n_estimators, seed
        real(kind=dp), intent(in) :: threshold
        integer, intent(in) :: change_points(:), n_change_points
        integer, intent(out) :: cp_best
        real(kind=dp), intent(out) :: score_best
        logical, intent(out) :: valid

        integer :: min_seg_size
        logical :: significant

        min_seg_size = window_size * excl_radius
        cp_best = -1
        score_best = -huge(1.0_dp)
        valid = .false.
        if (ubound - lbound < 2 * min_seg_size) return

        call ensemble_best_candidate_1d(z, lbound, ubound, window_size, n_estimators, k_neighbours, excl_radius, &
            .true., seed, cp_best, score_best, valid)
        if (.not. valid) return
        if (.not. cp_is_valid_1d(cp_best, change_points, n_change_points, min_seg_size, size(z))) then
            cp_best = -1
            score_best = -huge(1.0_dp)
            valid = .false.
            return
        end if

        significant = clasp_significance_valid_1d(z, lbound, ubound, window_size, k_neighbours, cp_best, threshold)
        if (.not. significant) then
            cp_best = -1
            score_best = -huge(1.0_dp)
            valid = .false.
        end if
    end subroutine clasp_segment_candidate_1d

    subroutine ensemble_best_candidate_1d(z, lbound_in, ubound_in, window_size, n_estimators, k_neighbours, excl_radius, &
            early_stopping, seed, cp_best, score_best, found)
        !> Deterministic ClaSPEnsemble candidate with the raw ClaSP profile score.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: lbound_in, ubound_in, window_size, n_estimators, k_neighbours, excl_radius, seed
        logical, intent(in) :: early_stopping
        integer, intent(out) :: cp_best
        real(kind=dp), intent(out) :: score_best
        logical, intent(out) :: found

        integer :: n, min_seg_size, seed_use, ntc, i
        integer :: lbound, area, ubound, state, cp_local
        integer, allocatable :: tc_l(:), tc_u(:)
        logical :: valid
        real(kind=dp) :: local_score

        n = ubound_in - lbound_in
        min_seg_size = window_size * excl_radius
        cp_best = -1
        score_best = -huge(1.0_dp)
        found = .false.
        if (n < 2 * min_seg_size) return

        allocate(tc_l(max(1, n_estimators)), tc_u(max(1, n_estimators)))
        ntc = 1
        tc_l(1) = lbound_in
        tc_u(1) = ubound_in

        seed_use = modulo(abs(seed), 2147483646) + 1
        state = seed_use
        do while (ntc < n_estimators .and. n > 3 * min_seg_size)
            state = park_miller_next(state)
            lbound = lbound_in + modulo(state - 1, n)
            state = park_miller_next(state)
            area = modulo(state - 1, n)
            if (ubound_in - lbound < area) area = ubound_in - lbound
            ubound = lbound + area
            if (ubound - lbound < 2 * min_seg_size) cycle
            ntc = ntc + 1
            tc_l(ntc) = lbound
            tc_u(ntc) = ubound
        end do

        call sort_temporal_constraints_desc(tc_l, tc_u, ntc)

        do i = 1, ntc
            call local_clasp_candidate_1d(z, tc_l(i), tc_u(i), window_size, k_neighbours, min_seg_size, cp_local, local_score, valid)
            if (.not. valid) cycle
            if (local_score > score_best .or. (cp_best < 0 .and. i == ntc)) then
                score_best = local_score
                cp_best = cp_local
                found = .true.
            else
                if (early_stopping) exit
            end if
        end do

        deallocate(tc_l, tc_u)
    end subroutine ensemble_best_candidate_1d

    logical function cp_is_valid_1d(candidate, change_points, n_change_points, min_seg_size, n_total) result(valid)
        !> Mirror BinaryClaSPSegmentation._cp_is_valid for a global raw change point.
        integer, intent(in) :: candidate, change_points(:), n_change_points, min_seg_size, n_total

        integer :: i

        valid = .false.
        if (candidate < min_seg_size) return
        if (candidate >= n_total - min_seg_size) return

        do i = 1, n_change_points
            if (candidate >= change_points(i) - min_seg_size .and. candidate < change_points(i) + min_seg_size) return
        end do

        valid = .true.
    end function cp_is_valid_1d

    logical function clasp_significance_valid_1d(z, lbound, ubound, window_size, k_neighbours, cp_global, threshold) result(valid)
        !> Port of claspy.validation.significance_test for one univariate local segment.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: lbound, ubound, window_size, k_neighbours, cp_global
        real(kind=dp), intent(in) :: threshold

        integer :: n_local, m, cp_local
        integer, allocatable :: offsets(:, :), y_true(:), y_pred(:)
        real(kind=dp), allocatable :: x(:), means(:), stds(:)
        real(kind=dp) :: p_value

        valid = .false.
        n_local = ubound - lbound
        cp_local = cp_global - lbound
        m = n_local - window_size + 1
        if (m <= 0) return
        if (cp_local <= 0 .or. cp_local >= m) return

        allocate(x(n_local), means(m), stds(m), offsets(m, k_neighbours), y_true(m), y_pred(m))
        x = z(lbound + 1:ubound)
        call compute_subseq_mean_std_1d(x, window_size, means, stds)
        call knn_offsets_znormed_1d(x, window_size, k_neighbours, means, stds, offsets)
        call cross_val_labels_1d(offsets, cp_local, window_size, y_true, y_pred)
        call rank_sums_binary_pvalue(y_pred(1:cp_local), y_pred(cp_local + 1:m), p_value)
        valid = (p_value <= threshold)

        deallocate(x, means, stds, offsets, y_true, y_pred)
    end function clasp_significance_valid_1d

    subroutine rank_sums_binary_pvalue(x, y, p_value)
        !> Two-sided rank-sum p-value for binary vectors, matching claspy.validation._rank_sums_test.
        integer, intent(in) :: x(:), y(:)
        real(kind=dp), intent(out) :: p_value

        integer :: n1, n2, n_total, n_zero
        real(kind=dp) :: zero_mean, one_mean, s, expected, denom, z_score

        n1 = size(x)
        n2 = size(y)
        n_total = n1 + n2
        if (n1 <= 0 .or. n2 <= 0) then
            p_value = 1.0_dp
            return
        end if

        n_zero = count(x == 0) + count(y == 0)
        zero_mean = 0.0_dp
        if (n_zero > 0) zero_mean = 0.5_dp * real(n_zero - 1, dp) + 1.0_dp
        one_mean = 0.0_dp
        if (n_zero < n_total) one_mean = 0.5_dp * real(n_zero + n_total - 1, dp) + 1.0_dp

        s = zero_mean * real(count(x == 0), dp) + one_mean * real(count(x == 1), dp)
        expected = real(n1, dp) * real(n_total + 1, dp) / 2.0_dp
        denom = sqrt(real(n1 * n2 * (n_total + 1), dp) / 12.0_dp)
        if (denom <= 0.0_dp) then
            p_value = 1.0_dp
            return
        end if

        z_score = (s - expected) / denom
        p_value = erfc(abs(z_score) / sqrt(2.0_dp))
    end subroutine rank_sums_binary_pvalue

    integer function suss_window_size_1d(z) result(window_size)
        !> Port of claspy.window_size.suss for a univariate series.
        real(kind=dp), intent(in) :: z(:)

        integer, parameter :: lbound_default = 10
        real(kind=dp), parameter :: threshold = 0.89_dp

        integer :: n, lbound, ubound, expn, w
        real(kind=dp), allocatable :: x(:)
        real(kind=dp) :: zmin, zmax, ts_mean, ts_std, ts_range, max_score, min_score, score

        n = size(z)
        if (n < lbound_default) then
            window_size = n
            return
        end if

        zmin = minval(z)
        zmax = maxval(z)
        if (zmin == zmax) then
            window_size = lbound_default
            return
        end if

        allocate(x(n))
        x = (z - zmin) / (zmax - zmin)
        ts_mean = sum(x) / real(n, dp)
        ts_std = sqrt(max(sum((x - ts_mean) ** 2) / real(n, dp), 0.0_dp))
        ts_range = maxval(x) - minval(x)

        max_score = suss_score_1d(x, 1, ts_mean, ts_std, ts_range)
        min_score = suss_score_1d(x, n - 1, ts_mean, ts_std, ts_range)
        if (min_score == max_score) then
            window_size = lbound_default
            deallocate(x)
            return
        end if

        expn = 0
        do
            w = 2 ** expn
            if (w < lbound_default) then
                expn = expn + 1
                cycle
            end if

            score = 1.0_dp - (suss_score_1d(x, w, ts_mean, ts_std, ts_range) - min_score) / (max_score - min_score)
            if (score > threshold) exit
            expn = expn + 1
        end do

        lbound = max(lbound_default, 2 ** max(expn - 1, 0))
        ubound = 2 ** expn + 1

        do while (lbound <= ubound)
            w = (lbound + ubound) / 2
            score = 1.0_dp - (suss_score_1d(x, w, ts_mean, ts_std, ts_range) - min_score) / (max_score - min_score)
            if (score < threshold) then
                lbound = w + 1
            else if (score > threshold) then
                ubound = w - 1
            else
                exit
            end if
        end do

        window_size = 2 * lbound
        deallocate(x)
    end function suss_window_size_1d

    real(kind=dp) function suss_score_1d(x, window_size, ts_mean, ts_std, ts_range) result(score)
        !> Summary-statistics similarity score used by SuSS.
        real(kind=dp), intent(in) :: x(:), ts_mean, ts_std, ts_range
        integer, intent(in) :: window_size

        integer :: n, start_idx, count_seg
        real(kind=dp) :: seg_mean, seg_std, seg_min, seg_max, dist

        n = size(x)
        score = 0.0_dp
        count_seg = 0
        do start_idx = 2, n - window_size + 1
            seg_mean = sum(x(start_idx:start_idx + window_size - 1)) / real(window_size, dp)
            seg_std = sqrt(max(sum((x(start_idx:start_idx + window_size - 1) - seg_mean) ** 2) / real(window_size, dp), 0.0_dp))
            seg_min = minval(x(start_idx:start_idx + window_size - 1))
            seg_max = maxval(x(start_idx:start_idx + window_size - 1))

            dist = (seg_mean - ts_mean) ** 2 + (seg_std - ts_std) ** 2 + ((seg_max - seg_min) - ts_range) ** 2
            score = score + sqrt(dist) / sqrt(real(window_size, dp))
            count_seg = count_seg + 1
        end do

        if (count_seg > 0) score = score / real(count_seg, dp)
    end function suss_score_1d

    pure real(kind=dp) function clasp_priority(score, lbound, ubound, n_total) result(priority)
        !> Priority score used by BinaryClaSPSegmentation with a single ClaSP ensemble member.
        real(kind=dp), intent(in) :: score
        integer, intent(in) :: lbound, ubound, n_total

        priority = (score + real(ubound - lbound, dp) / real(n_total, dp)) / 2.0_dp
    end function clasp_priority

    subroutine compute_subseq_mean_std_1d(x, window_size, means, stds)
        !> Sliding means and standard deviations for all subsequences of length window_size.
        real(kind=dp), intent(in) :: x(:)
        integer, intent(in) :: window_size
        real(kind=dp), intent(out) :: means(:), stds(:)

        real(kind=dp), allocatable :: s(:), ss(:)
        real(kind=dp) :: seg_sum, seg_sumsq, var
        integer :: n, m, i

        n = size(x)
        m = size(means)
        allocate(s(0:n), ss(0:n))
        s(0) = 0.0_dp
        ss(0) = 0.0_dp
        do i = 1, n
            s(i) = s(i - 1) + x(i)
            ss(i) = ss(i - 1) + x(i) * x(i)
        end do

        do i = 1, m
            seg_sum = s(i + window_size - 1) - s(i - 1)
            seg_sumsq = ss(i + window_size - 1) - ss(i - 1)
            means(i) = seg_sum / real(window_size, dp)
            var = seg_sumsq / real(window_size, dp) - means(i) * means(i)
            stds(i) = sqrt(max(var, 0.0_dp))
            if (abs(stds(i)) < 1.0e-3_dp) stds(i) = 1.0_dp
        end do

        deallocate(s, ss)
    end subroutine compute_subseq_mean_std_1d

    subroutine knn_offsets_znormed_1d(x, window_size, k_neighbours, means, stds, offsets)
        !> K nearest subsequence neighbours using z-normalized Euclidean distance and a half-window exclusion zone.
        real(kind=dp), intent(in) :: x(:), means(:), stds(:)
        integer, intent(in) :: window_size, k_neighbours
        integer, intent(out) :: offsets(:, :)

        integer :: m, i, j, t, excl_zone
        real(kind=dp) :: dot, dist
        real(kind=dp), allocatable :: best_dist(:)
        integer, allocatable :: best_idx(:)

        m = size(means)
        excl_zone = nint(real(window_size, dp) / 2.0_dp)
        allocate(best_dist(k_neighbours), best_idx(k_neighbours))

        do i = 1, m
            best_dist = huge(1.0_dp)
            best_idx = 0
            do j = 1, m
                if (abs(j - i) <= excl_zone) cycle
                dot = 0.0_dp
                do t = 1, window_size
                    dot = dot + x(i + t - 1) * x(j + t - 1)
                end do
                dist = 2.0_dp * real(window_size, dp) * (1.0_dp - &
                       (dot - real(window_size, dp) * means(i) * means(j)) / &
                       (real(window_size, dp) * stds(i) * stds(j)))
                dist = max(dist, 0.0_dp)
                call insert_smallest(dist, j - 1, best_dist, best_idx)
            end do
            offsets(i, :) = best_idx
        end do

        deallocate(best_dist, best_idx)
    end subroutine knn_offsets_znormed_1d

    subroutine insert_smallest(dist, idx, best_dist, best_idx)
        !> Insert one candidate into an ascending top-k distance list.
        real(kind=dp), intent(in) :: dist
        integer, intent(in) :: idx
        real(kind=dp), intent(inout) :: best_dist(:)
        integer, intent(inout) :: best_idx(:)

        integer :: k, pos

        k = size(best_dist)
        do pos = 1, k
            if (dist < best_dist(pos)) then
                if (pos < k) then
                    best_dist(pos + 1:k) = best_dist(pos:k - 1)
                    best_idx(pos + 1:k) = best_idx(pos:k - 1)
                end if
                best_dist(pos) = dist
                best_idx(pos) = idx
                exit
            end if
        end do
    end subroutine insert_smallest

    subroutine cross_val_labels_1d(offsets, split_idx, window_size, y_true, y_pred)
        !> Generate ClaSP true and predicted binary labels for one candidate split.
        integer, intent(in) :: offsets(:, :), split_idx, window_size
        integer, intent(out) :: y_true(:), y_pred(:)

        integer :: n_timepoints, k_neighbours, i, j, ones, zeros, start_excl, end_excl

        n_timepoints = size(offsets, 1)
        k_neighbours = size(offsets, 2)

        do i = 1, n_timepoints
            if (i - 1 < split_idx) then
                y_true(i) = 0
            else
                y_true(i) = 1
            end if
        end do

        do i = 1, n_timepoints
            ones = 0
            do j = 1, k_neighbours
                if (y_true(offsets(i, j) + 1) == 1) ones = ones + 1
            end do
            zeros = k_neighbours - ones
            if (ones > zeros) then
                y_pred(i) = 1
            else
                y_pred(i) = 0
            end if
        end do

        start_excl = max(0, split_idx - window_size)
        end_excl = min(n_timepoints - 1, split_idx - 1)
        if (end_excl >= start_excl) y_pred(start_excl + 1:end_excl + 1) = 1
    end subroutine cross_val_labels_1d

    real(kind=dp) function roc_auc_binary_score(y_score, y_true) result(area)
        !> Port of claspy.scoring.roc_auc_score for binary integer arrays.
        integer, intent(in) :: y_score(:), y_true(:)

        integer :: n, i, n_thresh, idx
        integer, allocatable :: score_rev(:), thresh(:), tps(:), fps(:)
        logical, allocatable :: truth_rev(:)
        integer :: cum_tp
        real(kind=dp), allocatable :: fpr(:), tpr(:)
        real(kind=dp) :: dx

        n = size(y_score)
        allocate(score_rev(n), truth_rev(n))
        do i = 1, n
            score_rev(i) = y_score(n - i + 1)
            truth_rev(i) = (y_true(n - i + 1) == 1)
        end do

        n_thresh = 1
        do i = 1, n - 1
            if (score_rev(i + 1) /= score_rev(i)) n_thresh = n_thresh + 1
        end do
        allocate(thresh(n_thresh))
        idx = 0
        do i = 1, n - 1
            if (score_rev(i + 1) /= score_rev(i)) then
                idx = idx + 1
                thresh(idx) = i
            end if
        end do
        thresh(n_thresh) = n

        allocate(tps(0:n_thresh), fps(0:n_thresh))
        tps(0) = 0
        fps(0) = 0
        cum_tp = 0
        idx = 1
        do i = 1, n
            if (truth_rev(i)) cum_tp = cum_tp + 1
            if (i == thresh(idx)) then
                tps(idx) = cum_tp
                fps(idx) = i - cum_tp
                if (idx < n_thresh) idx = idx + 1
            end if
        end do

        if (fps(n_thresh) <= 0 .or. tps(n_thresh) <= 0) then
            area = -huge(1.0_dp)
            deallocate(score_rev, truth_rev, thresh, tps, fps)
            return
        end if

        allocate(fpr(0:n_thresh), tpr(0:n_thresh))
        do i = 0, n_thresh
            fpr(i) = real(fps(i), dp) / real(fps(n_thresh), dp)
            tpr(i) = real(tps(i), dp) / real(tps(n_thresh), dp)
        end do

        area = 0.0_dp
        do i = 1, n_thresh
            dx = fpr(i) - fpr(i - 1)
            area = area + 0.5_dp * (tpr(i) + tpr(i - 1)) * dx
        end do

        deallocate(score_rev, truth_rev, thresh, tps, fps, fpr, tpr)
    end function roc_auc_binary_score

    subroutine push_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, queue_cap, lbound, ubound, cp, score)
        integer, intent(inout) :: queue_l(:), queue_u(:), queue_cp(:), n_queue, queue_cap
        real(kind=dp), intent(inout) :: queue_score(:)
        integer, intent(in) :: lbound, ubound, cp
        real(kind=dp), intent(in) :: score

        if (n_queue >= queue_cap) error stop "candidate queue capacity exceeded in solve_binary_clasp_1d"

        n_queue = n_queue + 1
        queue_l(n_queue) = lbound
        queue_u(n_queue) = ubound
        queue_cp(n_queue) = cp
        queue_score(n_queue) = score
    end subroutine push_candidate

    subroutine pop_best_candidate(queue_l, queue_u, queue_cp, queue_score, n_queue, idx)
        integer, intent(inout) :: queue_l(:), queue_u(:), queue_cp(:), n_queue
        real(kind=dp), intent(inout) :: queue_score(:)
        integer, intent(out) :: idx

        integer :: i
        real(kind=dp) :: best_score

        idx = 1
        best_score = queue_score(1)
        do i = 2, n_queue
            if (queue_score(i) > best_score) then
                best_score = queue_score(i)
                idx = i
            end if
        end do

        queue_l(idx) = queue_l(n_queue)
        queue_u(idx) = queue_u(n_queue)
        queue_cp(idx) = queue_cp(n_queue)
        queue_score(idx) = queue_score(n_queue)
        n_queue = n_queue - 1
    end subroutine pop_best_candidate

    subroutine solve_class_1d(z, window_size, cp, score, found, k_neighbours, excl_radius, threshold)
        !> Best split from the streaming ClaSS profile on a fixed window.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: window_size
        integer, intent(out) :: cp
        real(kind=dp), intent(out) :: score
        logical, intent(out) :: found
        integer, intent(in), optional :: k_neighbours, excl_radius
        real(kind=dp), intent(in), optional :: threshold

        integer :: k, excl, min_seg_size, m, split_idx
        integer, allocatable :: offsets(:, :)
        real(kind=dp), allocatable :: means(:), stds(:)
        real(kind=dp) :: threshold_use, score_split

        k = 3
        if (present(k_neighbours)) k = k_neighbours
        excl = 5
        if (present(excl_radius)) excl = excl_radius
        threshold_use = 0.5_dp
        if (present(threshold)) threshold_use = threshold
        min_seg_size = window_size * excl

        cp = -1
        score = -huge(1.0_dp)
        found = .false.
        m = size(z) - window_size + 1

        if (size(z) < 2 * min_seg_size .or. m <= 0) return

        allocate(means(m), stds(m), offsets(m, k))
        call compute_subseq_mean_std_1d(z, window_size, means, stds)
        call knn_offsets_znormed_1d(z, window_size, k, means, stds, offsets)

        do split_idx = min_seg_size, m - min_seg_size - 1
            score_split = class_f1_score_1d(offsets, split_idx, window_size)
            if (score_split > score) then
                score = score_split
                cp = split_idx
            end if
        end do

        found = (cp >= 0 .and. score >= threshold_use)
        if (.not. found) cp = -1

        deallocate(means, stds, offsets)
    end subroutine solve_class_1d

    real(kind=dp) function class_f1_score_1d(offsets, split_idx, window_size) result(score)
        !> Direct macro-F1 score for one streaming ClaSS split candidate.
        integer, intent(in) :: offsets(:, :), split_idx, window_size

        integer :: n_timepoints, k_neighbours, i, j, ones, zeros, truth_i, pred_i
        integer :: tp0, fp0, fn0, tn0
        integer :: excl_start, excl_end
        real(kind=dp) :: f0, f1

        n_timepoints = size(offsets, 1)
        k_neighbours = size(offsets, 2)
        excl_start = split_idx
        excl_end = min(n_timepoints - 1, split_idx + window_size - 1)

        tp0 = 0
        fp0 = 0
        fn0 = 0
        tn0 = 0

        do i = 1, n_timepoints
            if (i - 1 >= excl_start .and. i - 1 <= excl_end) cycle

            if (i - 1 < split_idx) then
                truth_i = 0
            else
                truth_i = 1
            end if

            ones = 0
            do j = 1, k_neighbours
                if (offsets(i, j) >= split_idx) ones = ones + 1
            end do
            zeros = k_neighbours - ones
            if (ones > zeros) then
                pred_i = 1
            else
                pred_i = 0
            end if

            if (truth_i == 0 .and. pred_i == 0) then
                tp0 = tp0 + 1
            else if (truth_i == 1 .and. pred_i == 0) then
                fp0 = fp0 + 1
            else if (truth_i == 0 .and. pred_i == 1) then
                fn0 = fn0 + 1
            else
                tn0 = tn0 + 1
            end if
        end do

        if (2 * tp0 + fp0 + fn0 <= 0) then
            score = -huge(1.0_dp)
            return
        end if
        f0 = 2.0_dp * real(tp0, dp) / real(2 * tp0 + fp0 + fn0, dp)

        if (2 * tn0 + fn0 + fp0 <= 0) then
            score = -huge(1.0_dp)
            return
        end if
        f1 = 2.0_dp * real(tn0, dp) / real(2 * tn0 + fn0 + fp0, dp)

        score = 0.5_dp * (f0 + f1)
    end function class_f1_score_1d

    subroutine solve_streaming_class_1d(z, n_timepoints, n_warmup, window_size, jump, cps, last_cp, &
            k_neighbours, excl_radius, threshold)
        !> Streaming ClaSP analog with an incremental streaming k-NN state.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: n_timepoints, n_warmup, window_size, jump
        integer, allocatable, intent(out) :: cps(:)
        integer, intent(out) :: last_cp
        integer, intent(in), optional :: k_neighbours, excl_radius
        real(kind=dp), intent(in), optional :: threshold

        integer :: n, k, excl, min_seg_size, n_windows, insert, lag, excl_zone
        integer :: ncp, i, ingested, lbound, n_filled, knn_filled, cp_local, m_local
        integer, allocatable :: tmp(:), knns(:, :)
        real(kind=dp), allocatable :: buffer(:), means(:), stds(:), dists(:, :), dot_rolled(:), prof(:)
        real(kind=dp) :: threshold_use, score_local
        logical :: found, has_dot

        n = size(z)
        k = 3
        if (present(k_neighbours)) k = k_neighbours
        excl = 5
        if (present(excl_radius)) excl = excl_radius
        threshold_use = 0.5_dp
        if (present(threshold)) threshold_use = threshold

        min_seg_size = 5 * window_size
        n_windows = n_timepoints - window_size + 1
        insert = n_windows - window_size / 2 - k - 1
        lag = window_size + window_size / 2 + k
        excl_zone = window_size / 2

        if (n <= 0 .or. n_timepoints < 2 * min_seg_size .or. insert <= 0) then
            allocate(cps(0))
            last_cp = 0
            return
        end if

        allocate(buffer(n_timepoints), tmp(max(0, n)), means(n_windows), stds(n_windows), &
            dists(n_windows, k), knns(n_windows, k), dot_rolled(n_windows), prof(n_windows))
        buffer = 0.0_dp
        means = 0.0_dp
        stds = 1.0_dp
        dists = huge(1.0_dp)
        knns = -1
        dot_rolled = huge(1.0_dp)
        prof = -huge(1.0_dp)
        ncp = 0
        last_cp = 0
        lbound = 0
        n_filled = 0
        knn_filled = 0
        has_dot = .false.

        do i = 1, n
            call streaming_knn_update_1d(z(i), last_cp, n_windows, window_size, k, insert, excl_zone, &
                buffer, means, stds, dists, knns, dot_rolled, has_dot, lbound, n_filled, knn_filled)

            ingested = i
            if (ingested < max(n_warmup, 2 * min_seg_size)) cycle
            if (insert - knn_filled == 0) last_cp = max(0, last_cp - 1)
            if (insert - lbound < 2 * min_seg_size) cycle
            if (mod(ingested, jump) /= 0) cycle

            call class_profile_from_offsets(knns(lbound + 1:insert, :), lbound, window_size, min_seg_size, prof(1:insert - lbound), &
                cp_local, score_local, found, threshold_use)
            if (found) then
                m_local = insert - lbound
                if (cp_local >= min_seg_size .and. m_local - cp_local >= min_seg_size) then
                    last_cp = last_cp + cp_local
                    ncp = ncp + 1
                    tmp(ncp) = ingested - lag - (insert - lbound) + cp_local + window_size
                end if
            end if
        end do

        allocate(cps(ncp))
        if (ncp > 0) cps = tmp(1:ncp)
        deallocate(buffer, tmp, means, stds, dists, knns, dot_rolled, prof)
    end subroutine solve_streaming_class_1d

    subroutine streaming_knn_update_1d(x_new, change_point, n_windows, window_size, k_neighbours, &
            insert, excl_zone, time_series, means, stds, dists, knns, dot_rolled, has_dot, lbound, n_filled, knn_filled)
        real(kind=dp), intent(in) :: x_new
        integer, intent(in) :: change_point, n_windows, window_size, k_neighbours, insert, excl_zone
        real(kind=dp), intent(inout) :: time_series(:), means(:), stds(:), dists(:, :), dot_rolled(:)
        integer, intent(inout) :: knns(:, :), lbound, n_filled, knn_filled
        logical, intent(inout) :: has_dot

        integer :: lbound_pre
        integer, allocatable :: knn_new(:)
        real(kind=dp), allocatable :: dist(:)
        logical :: first_flag

        lbound_pre = insert - knn_filled + 1 + change_point
        n_filled = min(n_filled + 1, n_windows)

        call roll_left_real(time_series, x_new)
        call compute_streaming_mean_std_1d(time_series, n_filled, window_size, means, stds)

        lbound = lbound_pre
        if (n_filled < window_size + excl_zone + k_neighbours) return

        if (knn_filled > 0) then
            call roll_left_matrix_real(dists, huge(1.0_dp))
            call roll_left_matrix_int(knns, -1)
            if (insert - knn_filled + 1 <= insert) knns(insert - knn_filled + 1:insert, :) = &
                knns(insert - knn_filled + 1:insert, :) - 1
        end if

        allocate(dist(n_windows), knn_new(k_neighbours))
        first_flag = .not. has_dot
        call streaming_knn_core_1d(time_series, means, stds, n_filled, n_windows, window_size, k_neighbours, &
            insert, excl_zone, lbound_pre, dot_rolled, first_flag, dist, knn_new)
        has_dot = .true.
        call streaming_roll_knns_1d(dist, knn_new, insert, k_neighbours, n_windows, lbound_pre, dists, knns, lbound, knn_filled)
        deallocate(dist, knn_new)
    end subroutine streaming_knn_update_1d

    subroutine compute_streaming_mean_std_1d(x, n_filled, window_size, means, stds)
        real(kind=dp), intent(in) :: x(:)
        integer, intent(in) :: n_filled, window_size
        real(kind=dp), intent(inout) :: means(:), stds(:)

        integer :: n_windows, n_valid, start_idx, i
        real(kind=dp), allocatable :: mean_tail(:), std_tail(:)

        n_windows = size(means)
        means = 0.0_dp
        stds = 1.0_dp
        if (n_filled < window_size) return

        n_valid = n_filled - window_size + 1
        start_idx = n_windows - n_valid + 1
        allocate(mean_tail(n_valid), std_tail(n_valid))
        call compute_subseq_mean_std_1d(x(size(x) - n_filled + 1:size(x)), window_size, mean_tail, std_tail)
        means(start_idx:n_windows) = mean_tail
        stds(start_idx:n_windows) = std_tail
        do i = 1, start_idx - 1
            means(i) = 0.0_dp
            stds(i) = 1.0_dp
        end do
        deallocate(mean_tail, std_tail)
    end subroutine compute_streaming_mean_std_1d

    subroutine class_profile_from_offsets(offsets_full, lbound, window_size, min_seg_size, profile, cp, score, found, threshold)
        integer, intent(in) :: offsets_full(:, :), lbound, window_size, min_seg_size
        real(kind=dp), intent(out) :: profile(:), score
        integer, intent(out) :: cp
        logical, intent(out) :: found
        real(kind=dp), intent(in) :: threshold

        integer :: nrows, k, i, j, split_idx
        integer, allocatable :: offsets(:, :)
        real(kind=dp) :: score_split

        nrows = size(offsets_full, 1)
        k = size(offsets_full, 2)
        profile = -huge(1.0_dp)
        cp = -1
        score = -huge(1.0_dp)
        found = .false.
        if (nrows < 2 * min_seg_size) return

        allocate(offsets(nrows, k))
        offsets = offsets_full - lbound
        do i = 1, nrows
            do j = 1, k
                offsets(i, j) = min(max(offsets(i, j), 0), nrows - 1)
            end do
        end do

        do split_idx = min_seg_size, nrows - min_seg_size - 1
            score_split = class_f1_score_1d(offsets, split_idx, window_size)
            profile(split_idx + 1) = score_split
            if (score_split > score) then
                score = score_split
                cp = split_idx
            end if
        end do

        found = (cp >= 0 .and. score >= threshold)
        if (.not. found) cp = -1
        deallocate(offsets)
    end subroutine class_profile_from_offsets

    subroutine solve_clap_centroid_1d(z, state_labels, window_size, n_splits, y_true_out, y_pred_out, score, sample_cap, seed)
        !> CLaP-style window classification with deterministic nearest-centroid cross-validation.
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: state_labels(:), window_size, n_splits
        integer, allocatable, intent(out) :: y_true_out(:), y_pred_out(:)
        real(kind=dp), intent(out) :: score
        integer, intent(in), optional :: sample_cap, seed

        integer :: stride, cap, seed_use
        real(kind=dp), allocatable :: xmat(:, :)
        integer, allocatable :: y(:)

        stride = max(1, window_size / 2)
        cap = 1000
        if (present(sample_cap)) cap = sample_cap
        seed_use = 2357
        if (present(seed)) seed_use = seed

        call create_clap_dataset_1d(z, state_labels, window_size, stride, xmat, y)
        call subselect_clap_dataset_1d(xmat, y, cap, seed_use)
        call crossval_centroid_1d(xmat, y, n_splits, seed_use, y_true_out, y_pred_out)
        score = macro_f1_labels(y_true_out, y_pred_out)

        deallocate(xmat, y)
    end subroutine solve_clap_centroid_1d

    subroutine create_clap_dataset_1d(z, state_labels, window_size, stride, xmat, y)
        real(kind=dp), intent(in) :: z(:)
        integer, intent(in) :: state_labels(:), window_size, stride
        real(kind=dp), allocatable, intent(out) :: xmat(:, :)
        integer, allocatable, intent(out) :: y(:)

        integer :: n, i, cp_count, idx, nwin, start_cp, end_cp
        integer, allocatable :: change_points(:)
        logical, allocatable :: excl(:)

        n = size(z)
        allocate(excl(n))
        excl = .false.
        cp_count = count(state_labels(1:n - 1) /= state_labels(2:n))
        allocate(change_points(cp_count))
        idx = 0
        do i = 1, n - 1
            if (state_labels(i) /= state_labels(i + 1)) then
                idx = idx + 1
                change_points(idx) = i
            end if
        end do

        do i = 1, cp_count
            start_cp = max(1, change_points(i) - window_size / 2 + 1)
            end_cp = max(start_cp - 1, change_points(i) - 1)
            if (end_cp >= start_cp) excl(start_cp:end_cp) = .true.
        end do

        nwin = 0
        do i = 1, n - window_size + 1, stride
            if (.not. excl(i)) nwin = nwin + 1
        end do

        allocate(xmat(window_size, nwin), y(nwin))
        idx = 0
        do i = 1, n - window_size + 1, stride
            if (excl(i)) cycle
            idx = idx + 1
            xmat(:, idx) = z(i:i + window_size - 1)
            y(idx) = state_labels(i)
        end do

        deallocate(excl, change_points)
    end subroutine create_clap_dataset_1d

    subroutine subselect_clap_dataset_1d(xmat, y, sample_cap, seed)
        real(kind=dp), allocatable, intent(inout) :: xmat(:, :)
        integer, allocatable, intent(inout) :: y(:)
        integer, intent(in) :: sample_cap, seed

        integer :: n, i, label_count, keep_n, pos, idx
        integer, allocatable :: labels(:), counts(:), label_idx(:), perm(:), keep(:)
        real(kind=dp), allocatable :: xnew(:, :)
        integer, allocatable :: ynew(:)

        n = size(y)
        call unique_labels(y, labels)
        label_count = size(labels)
        allocate(counts(label_count))
        counts = 0
        do i = 1, label_count
            counts(i) = count(y == labels(i))
        end do

        keep_n = 0
        do i = 1, label_count
            keep_n = keep_n + min(counts(i), sample_cap)
        end do

        allocate(keep(keep_n))
        pos = 0
        do i = 1, label_count
            allocate(label_idx(counts(i)))
            idx = 0
            do n = 1, size(y)
                if (y(n) == labels(i)) then
                    idx = idx + 1
                    label_idx(idx) = n
                end if
            end do
            call random_permutation_indices(counts(i), seed + 37 * i, perm)
            do n = 1, min(counts(i), sample_cap)
                pos = pos + 1
                keep(pos) = label_idx(perm(n))
            end do
            deallocate(label_idx, perm)
        end do

        call random_permutation_indices(keep_n, seed + 911, perm)
        allocate(xnew(size(xmat, 1), keep_n), ynew(keep_n))
        do i = 1, keep_n
            xnew(:, i) = xmat(:, keep(perm(i)))
            ynew(i) = y(keep(perm(i)))
        end do

        call move_alloc(xnew, xmat)
        call move_alloc(ynew, y)
        deallocate(labels, counts, keep, perm)
    end subroutine subselect_clap_dataset_1d

    subroutine crossval_centroid_1d(xmat, y, n_splits, seed, y_true_out, y_pred_out)
        real(kind=dp), intent(in) :: xmat(:, :)
        integer, intent(in) :: y(:), n_splits, seed
        integer, allocatable, intent(out) :: y_true_out(:), y_pred_out(:)

        integer :: n, nfold, i, start_idx, end_idx, fold, ntrain, ntest, pos_train, pos_test
        integer, allocatable :: perm(:), test_mask(:), ytrain(:), test_idx(:), train_idx(:), labels(:), counts(:)
        real(kind=dp), allocatable :: xtrain(:, :), centroids(:, :)
        real(kind=dp) :: best_dist, dist
        integer :: best_label, j, k

        n = size(y)
        nfold = min(max(2, n_splits), n)
        allocate(y_true_out(n), y_pred_out(n))
        y_true_out = y
        y_pred_out = y
        if (n < 2) return

        call random_permutation_indices(n, seed, perm)
        allocate(test_mask(n))
        do fold = 1, nfold
            test_mask = 0
            start_idx = ((fold - 1) * n) / nfold + 1
            end_idx = (fold * n) / nfold
            do i = start_idx, end_idx
                test_mask(perm(i)) = 1
            end do

            ntest = count(test_mask == 1)
            ntrain = n - ntest
            allocate(test_idx(ntest), train_idx(ntrain), xtrain(size(xmat, 1), ntrain), ytrain(ntrain))
            pos_train = 0
            pos_test = 0
            do i = 1, n
                if (test_mask(i) == 1) then
                    pos_test = pos_test + 1
                    test_idx(pos_test) = i
                else
                    pos_train = pos_train + 1
                    train_idx(pos_train) = i
                    xtrain(:, pos_train) = xmat(:, i)
                    ytrain(pos_train) = y(i)
                end if
            end do

            call unique_labels(ytrain, labels)
            allocate(counts(size(labels)), centroids(size(xmat, 1), size(labels)))
            counts = 0
            centroids = 0.0_dp
            do i = 1, ntrain
                do j = 1, size(labels)
                    if (ytrain(i) == labels(j)) then
                        counts(j) = counts(j) + 1
                        centroids(:, j) = centroids(:, j) + xtrain(:, i)
                        exit
                    end if
                end do
            end do
            do j = 1, size(labels)
                if (counts(j) > 0) centroids(:, j) = centroids(:, j) / real(counts(j), dp)
            end do

            do i = 1, ntest
                best_dist = huge(1.0_dp)
                best_label = labels(1)
                do j = 1, size(labels)
                    dist = 0.0_dp
                    do k = 1, size(xmat, 1)
                        dist = dist + (xmat(k, test_idx(i)) - centroids(k, j)) ** 2
                    end do
                    if (dist < best_dist) then
                        best_dist = dist
                        best_label = labels(j)
                    end if
                end do
                y_pred_out(test_idx(i)) = best_label
            end do

            deallocate(test_idx, train_idx, xtrain, ytrain, labels, counts, centroids)
        end do
        deallocate(perm, test_mask)
    end subroutine crossval_centroid_1d

    real(kind=dp) function macro_f1_labels(y_true, y_pred) result(score)
        integer, intent(in) :: y_true(:), y_pred(:)
        integer, allocatable :: labels(:)
        integer :: i, j, tp, fp, fn
        real(kind=dp) :: f1_sum

        call unique_labels(y_true, labels)
        f1_sum = 0.0_dp
        do i = 1, size(labels)
            tp = 0
            fp = 0
            fn = 0
            do j = 1, size(y_true)
                if (y_true(j) == labels(i) .and. y_pred(j) == labels(i)) tp = tp + 1
                if (y_true(j) /= labels(i) .and. y_pred(j) == labels(i)) fp = fp + 1
                if (y_true(j) == labels(i) .and. y_pred(j) /= labels(i)) fn = fn + 1
            end do
            if (2 * tp + fp + fn > 0) then
                f1_sum = f1_sum + 2.0_dp * real(tp, dp) / real(2 * tp + fp + fn, dp)
            end if
        end do
        score = f1_sum / real(size(labels), dp)
        deallocate(labels)
    end function macro_f1_labels

    real(kind=dp) function random_f1_labels(y_true) result(score)
        integer, intent(in) :: y_true(:)
        integer, allocatable :: labels(:)
        integer :: i, n, pos_instances, neg_instances
        real(kind=dp) :: tp, fp, fn, pre, rec

        call unique_labels(y_true, labels)
        n = size(y_true)
        score = 0.0_dp
        do i = 1, size(labels)
            pos_instances = count(y_true == labels(i))
            neg_instances = n - pos_instances
            tp = real(pos_instances * pos_instances, dp) / real(n, dp)
            fn = real(pos_instances * neg_instances, dp) / real(n, dp)
            fp = real(neg_instances * pos_instances, dp) / real(n, dp)
            pre = tp / (tp + fp)
            rec = tp / (tp + fn)
            if (pre + rec > 0.0_dp) score = score + 2.0_dp * pre * rec / (pre + rec)
        end do
        score = score / real(size(labels), dp)
        deallocate(labels)
    end function random_f1_labels

    real(kind=dp) function classification_gain_labels(y_true, y_pred) result(gain)
        integer, intent(in) :: y_true(:), y_pred(:)
        gain = macro_f1_labels(y_true, y_pred) - random_f1_labels(y_true)
    end function classification_gain_labels

    subroutine unique_labels(y, labels)
        integer, intent(in) :: y(:)
        integer, allocatable, intent(out) :: labels(:)
        integer :: i, nuniq
        integer, allocatable :: tmp(:)

        allocate(tmp(size(y)))
        nuniq = 0
        do i = 1, size(y)
            if (.not. any(tmp(1:nuniq) == y(i))) then
                nuniq = nuniq + 1
                tmp(nuniq) = y(i)
            end if
        end do
        allocate(labels(nuniq))
        labels = tmp(1:nuniq)
        deallocate(tmp)
    end subroutine unique_labels

    subroutine create_state_labels_from_bkps(true_bkps, segment_labels, n, state_labels)
        integer, intent(in) :: true_bkps(:), segment_labels(:), n
        integer, intent(out) :: state_labels(n)
        integer :: i, start_idx, end_idx

        start_idx = 1
        do i = 1, size(true_bkps)
            end_idx = true_bkps(i)
            state_labels(start_idx:end_idx) = segment_labels(i)
            start_idx = end_idx + 1
        end do
        state_labels(start_idx:n) = segment_labels(size(segment_labels))
    end subroutine create_state_labels_from_bkps

    subroutine random_permutation_indices(n, seed, perm)
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
            call swap_int(perm(i), perm(j))
        end do
    end subroutine random_permutation_indices

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

    subroutine swap_int(a, b)
        integer, intent(inout) :: a, b
        integer :: tmp
        tmp = a
        a = b
        b = tmp
    end subroutine swap_int

    subroutine sort_temporal_constraints_desc(tc_l, tc_u, ntc)
        integer, intent(inout) :: tc_l(:), tc_u(:)
        integer, intent(in) :: ntc
        integer :: i, j
        integer :: tmp_l, tmp_u

        do i = 1, ntc - 1
            do j = i + 1, ntc
                if (tc_u(j) - tc_l(j) > tc_u(i) - tc_l(i)) then
                    tmp_l = tc_l(i)
                    tmp_u = tc_u(i)
                    tc_l(i) = tc_l(j)
                    tc_u(i) = tc_u(j)
                    tc_l(j) = tmp_l
                    tc_u(j) = tmp_u
                end if
            end do
        end do
    end subroutine sort_temporal_constraints_desc

    subroutine roll_left_real(x, fill_value)
        real(kind=dp), intent(inout) :: x(:)
        real(kind=dp), intent(in) :: fill_value
        if (size(x) > 1) x(1:size(x) - 1) = x(2:size(x))
        x(size(x)) = fill_value
    end subroutine roll_left_real

    subroutine roll_left_matrix_real(x, fill_value)
        real(kind=dp), intent(inout) :: x(:, :)
        real(kind=dp), intent(in) :: fill_value
        if (size(x, 1) > 1) x(1:size(x, 1) - 1, :) = x(2:size(x, 1), :)
        x(size(x, 1), :) = fill_value
    end subroutine roll_left_matrix_real

    subroutine roll_left_matrix_int(x, fill_value)
        integer, intent(inout) :: x(:, :)
        integer, intent(in) :: fill_value
        if (size(x, 1) > 1) x(1:size(x, 1) - 1, :) = x(2:size(x, 1), :)
        x(size(x, 1), :) = fill_value
    end subroutine roll_left_matrix_int

    subroutine sliding_dot_naive(query, series, out)
        real(kind=dp), intent(in) :: query(:), series(:)
        real(kind=dp), intent(out) :: out(:)
        integer :: m, n, i

        m = size(query)
        n = size(series)
        do i = 1, n - m + 1
            out(i) = sum(query * series(i:i + m - 1))
        end do
    end subroutine sliding_dot_naive

    subroutine argkmin_restore(dist, k, lbound, args)
        real(kind=dp), intent(inout) :: dist(:)
        integer, intent(in) :: k, lbound
        integer, intent(out) :: args(:)

        integer :: i, j, best_idx
        real(kind=dp) :: best_val
        real(kind=dp), allocatable :: vals(:)

        allocate(vals(k))
        do i = 1, k
            best_idx = max(1, lbound + 1)
            best_val = huge(1.0_dp)
            do j = max(1, lbound + 1), size(dist)
                if (dist(j) < best_val) then
                    best_val = dist(j)
                    best_idx = j
                end if
            end do
            args(i) = best_idx - 1
            vals(i) = dist(best_idx)
            dist(best_idx) = huge(1.0_dp)
        end do
        do i = 1, k
            dist(args(i) + 1) = vals(i)
        end do
        deallocate(vals)
    end subroutine argkmin_restore

    subroutine streaming_knn_core_1d(time_series, means, stds, n_filled, n_windows, window_size, k_neighbours, &
            insert, excl_zone, lbound_pre, dot_rolled, first_flag, dist, knn_new)
        real(kind=dp), intent(in) :: time_series(:), means(:), stds(:)
        integer, intent(in) :: n_filled, n_windows, window_size, k_neighbours, insert, excl_zone, lbound_pre
        real(kind=dp), intent(inout) :: dot_rolled(:)
        logical, intent(in) :: first_flag
        real(kind=dp), intent(out) :: dist(:)
        integer, intent(out) :: knn_new(:)

        integer :: start_idx, valid_start, valid_len, excl_lo, excl_hi, j
        real(kind=dp) :: dot, maxdist
        real(kind=dp), allocatable :: dots(:)

        start_idx = max(0, lbound_pre - 1)
        valid_start = start_idx + 1
        dist = huge(1.0_dp)

        if (first_flag) then
            valid_len = n_filled - window_size + 1
            allocate(dots(valid_len))
            call sliding_dot_naive(time_series(insert + 1:insert + window_size), &
                time_series(size(time_series) - n_filled + 1:size(time_series)), dots)
            dot_rolled(valid_start:n_windows) = dots
            deallocate(dots)
        else
            dot_rolled = dot_rolled + time_series(insert + window_size) * time_series(window_size:size(time_series))
            dot = sum(time_series(valid_start:valid_start + window_size - 1) * time_series(insert + 1:insert + window_size))
            dot_rolled(valid_start) = dot
        end if

        do j = valid_start, n_windows
            dist(j) = 2.0_dp * real(window_size, dp) * (1.0_dp - &
                (dot_rolled(j) - real(window_size, dp) * means(j) * means(insert + 1)) / &
                (real(window_size, dp) * stds(j) * stds(insert + 1)))
        end do

        maxdist = maxval(dist(valid_start:n_windows))
        excl_lo = max(0, insert - excl_zone) + 1
        excl_hi = min(insert + excl_zone + 1, n_windows)
        dist(excl_lo:excl_hi) = maxdist

        call argkmin_restore(dist, k_neighbours, max(0, lbound_pre), knn_new)
        dot_rolled = dot_rolled - time_series(insert + 1) * time_series(1:n_windows)
    end subroutine streaming_knn_core_1d

    subroutine streaming_roll_knns_1d(dist, knn_new, insert, k_neighbours, n_windows, lbound_pre, dists, knns, lbound, knn_filled)
        real(kind=dp), intent(in) :: dist(:)
        integer, intent(in) :: knn_new(:), insert, k_neighbours, n_windows, lbound_pre
        real(kind=dp), intent(inout) :: dists(:, :)
        integer, intent(inout) :: knns(:, :)
        integer, intent(inout) :: lbound, knn_filled

        integer :: start_idx, pos, col, kk
        logical, allocatable :: change_mask(:)

        dists(insert + 1, :) = dist(knn_new + 1)
        knns(insert + 1, :) = knn_new

        start_idx = max(0, lbound_pre)
        allocate(change_mask(n_windows - start_idx))
        change_mask = .true.

        do col = 1, k_neighbours - 1
            do pos = start_idx, n_windows - 1
                if (.not. change_mask(pos - start_idx + 1)) cycle
                if (dist(pos + 1) < dists(pos + 1, col)) then
                    change_mask(pos - start_idx + 1) = .false.
                    do kk = k_neighbours, col + 1, -1
                        knns(pos + 1, kk) = knns(pos + 1, kk - 1)
                        dists(pos + 1, kk) = dists(pos + 1, kk - 1)
                    end do
                    knns(pos + 1, col) = insert
                    dists(pos + 1, col) = dist(pos + 1)
                end if
            end do
        end do

        lbound = max(0, lbound_pre - 1)
        knn_filled = min(knn_filled + 1, insert)
        deallocate(change_mask)
    end subroutine streaming_roll_knns_1d

end module claspy_mod
