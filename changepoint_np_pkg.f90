module changepoint_np_pkg_mod
use kind_mod, only: dp
implicit none
private
public :: solve_cpt_np_pelt_1d, changepoint_np_penalty_value

contains

    subroutine solve_cpt_np_pelt_1d(x, penalty, minseglen, nquantiles, cpts, lastchangecpts, lastchangelike, numchangecpts, pen_value)
        real(kind=dp), intent(in) :: x(:)
        character(len=*), intent(in) :: penalty
        integer, intent(in) :: minseglen, nquantiles
        integer, allocatable, intent(out) :: cpts(:), lastchangecpts(:), numchangecpts(:)
        real(kind=dp), allocatable, intent(out) :: lastchangelike(:)
        real(kind=dp), intent(in), optional :: pen_value

        real(kind=dp), allocatable :: q(:, :)
        real(kind=dp) :: pen, best_val, cand_val
        real(kind=dp), parameter :: huge_cost = 1.0e100_dp
        integer, allocatable :: parent(:), counts(:), tmp_cpts(:)
        integer :: n, t, tau, best_tau, n_found, i, b

        n = size(x)
        if (n <= 0) then
            allocate(cpts(0), lastchangecpts(0), numchangecpts(0), lastchangelike(0))
            return
        end if

        pen = changepoint_np_penalty_value(trim(penalty), n, 1, pen_value)
        call build_nonparametric_ed_sumstat_1d(x, nquantiles, q)

        allocate(lastchangelike(0:n), parent(0:n), counts(0:n))
        lastchangelike = huge_cost
        parent = -1
        counts = 0
        lastchangelike(0) = -pen
        parent(0) = 0
        counts(0) = 0
        do t = 1, min(minseglen - 1, n)
            lastchangelike(t) = 0.0_dp
            parent(t) = 0
            counts(t) = 0
        end do

        do t = minseglen, n
            best_val = huge_cost
            best_tau = -1
            do tau = 0, t - minseglen
                cand_val = lastchangelike(tau) + nonparametric_ed_segment_cost(q, tau, t, minseglen) + pen
                if (cand_val < best_val) then
                    best_val = cand_val
                    best_tau = tau
                end if
            end do

            lastchangelike(t) = best_val
            parent(t) = best_tau
            if (best_tau >= 0) then
                counts(t) = counts(best_tau) + 1
            else
                counts(t) = 0
            end if
        end do

        allocate(lastchangecpts(0:n), numchangecpts(0:n))
        lastchangecpts = 0
        numchangecpts = counts
        do t = 1, n
            if (parent(t) > 0) lastchangecpts(t) = parent(t)
        end do

        allocate(tmp_cpts(n))
        n_found = 1
        tmp_cpts(1) = n
        b = n
        do while (parent(b) > 0)
            b = parent(b)
            n_found = n_found + 1
            tmp_cpts(n_found) = b
        end do

        allocate(cpts(n_found))
        do i = 1, n_found
            cpts(i) = tmp_cpts(n_found - i + 1)
        end do

        deallocate(q, parent, counts, tmp_cpts)
    end subroutine solve_cpt_np_pelt_1d

    real(kind=dp) function changepoint_np_penalty_value(penalty, n, diffparam, pen_value) result(penalty_out)
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
            if (.not. present(pen_value)) error stop "Manual changepoint.np penalty requires pen_value"
            penalty_out = pen_value
        case default
            error stop "Unsupported changepoint.np penalty"
        end select

        if (penalty_out < 0.0_dp) error stop "Negative changepoint.np penalty is invalid"
    end function changepoint_np_penalty_value

    subroutine build_nonparametric_ed_sumstat_1d(x, nquantiles, q)
        real(kind=dp), intent(in) :: x(:)
        integer, intent(in) :: nquantiles
        real(kind=dp), allocatable, intent(out) :: q(:, :)

        real(kind=dp), allocatable :: xs(:)
        real(kind=dp) :: cval, yk, pk, thresh
        integer :: n, k, i, j, k_use

        n = size(x)
        k_use = min(nquantiles, n)
        allocate(q(k_use, 0:n), xs(n))
        xs = x
        call sort_real_in_place(xs)

        cval = -log(real(2 * n - 1, dp))
        q = 0.0_dp

        do k = 1, k_use
            yk = -1.0_dp + (2.0_dp * real(k, dp) / real(k_use, dp) - 1.0_dp / real(k_use, dp))
            pk = 1.0_dp / (1.0_dp + exp(cval * yk))
            j = int(real(n - 1, dp) * pk + 1.0_dp)
            j = max(1, min(n, j))
            thresh = xs(j)
            do i = 1, n
                q(k, i) = q(k, i - 1)
                if (x(i) < thresh) then
                    q(k, i) = q(k, i) + 1.0_dp
                else if (x(i) == thresh) then
                    q(k, i) = q(k, i) + 0.5_dp
                end if
            end do
        end do

        deallocate(xs)
    end subroutine build_nonparametric_ed_sumstat_1d

    real(kind=dp) function nonparametric_ed_segment_cost(q, s, t, minseglen) result(cost)
        real(kind=dp), intent(in) :: q(:, 0:)
        integer, intent(in) :: s, t, minseglen

        real(kind=dp) :: seg_len, seg_count, factor, cval
        integer :: n, k, k_use

        n = ubound(q, 2)
        seg_len = real(t - s, dp)
        if (t <= s .or. seg_len < real(minseglen, dp)) then
            cost = 1.0e100_dp
            return
        end if

        k_use = size(q, 1)
        cval = -log(real(2 * n - 1, dp))
        factor = -2.0_dp * cval / real(k_use, dp)
        cost = 0.0_dp
        do k = 1, k_use
            seg_count = q(k, t) - q(k, s)
            cost = cost + factor * entropy_term(seg_len, seg_count)
        end do
    end function nonparametric_ed_segment_cost

    pure real(kind=dp) function entropy_term(seg_len, seg_count) result(term)
        real(kind=dp), intent(in) :: seg_len, seg_count
        term = xlogx(seg_len) - xlogx(seg_count) - xlogx(seg_len - seg_count)
    end function entropy_term

    pure real(kind=dp) function xlogx(x) result(val)
        real(kind=dp), intent(in) :: x
        if (x <= 0.0_dp) then
            val = 0.0_dp
        else
            val = x * log(x)
        end if
    end function xlogx

    subroutine sort_real_in_place(x)
        real(kind=dp), intent(in out) :: x(:)
        integer :: i, j, n
        real(kind=dp) :: tmp
        n = size(x)
        do i = 1, n - 1
            do j = i + 1, n
                if (x(i) > x(j)) then
                    tmp = x(i)
                    x(i) = x(j)
                    x(j) = tmp
                end if
            end do
        end do
    end subroutine sort_real_in_place

end module changepoint_np_pkg_mod
