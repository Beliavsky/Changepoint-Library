module compare_sim_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: seed_rng_fixed, simulate_piecewise_normal_1d, simulate_piecewise_normal_mv, &
              simulate_piecewise_poisson_1d, simulate_piecewise_exponential_1d, simulate_piecewise_gamma_1d, simulate_piecewise_bernoulli_1d, simulate_piecewise_ar1_1d, &
              simulate_piecewise_linear_1d, simulate_piecewise_linear_mv_1d, simulate_piecewise_linear_harmonic_1d, &
              simulate_piecewise_sine_1d, &
              simulate_piecewise_linear_irregular_1d, &
              simulate_irregular_time_grid

    logical, save :: has_spare_normal = .false.
    real(kind=dp), save :: spare_normal = 0.0_dp

contains

    subroutine seed_rng_fixed(seed)
        !> Seed the intrinsic RNG reproducibly from a single integer seed.
        integer, intent(in) :: seed

        integer :: n, i
        integer, allocatable :: seed_vec(:)

        call random_seed(size=n)
        allocate(seed_vec(n))
        do i = 1, n
            seed_vec(i) = modulo(seed + 104729 * (i - 1), huge(1) - 1) + 1
        end do
        call random_seed(put=seed_vec)
        has_spare_normal = .false.
        deallocate(seed_vec)
    end subroutine seed_rng_fixed

    subroutine simulate_piecewise_normal_1d(n, regime_starts, means, sds, x)
        !> Simulate a piecewise-Gaussian univariate series.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: means(:), sds(:)
        real(kind=dp), intent(out) :: x(n)

        integer :: i, i1, i2

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            call fill_normal_segment(x(i1:i2), means(i), sds(i))
        end do
    end subroutine simulate_piecewise_normal_1d

    subroutine simulate_piecewise_normal_mv(n, p, regime_starts, means, sds, x)
        !> Simulate a piecewise-Gaussian multivariate series with diagonal covariance.
        integer, intent(in) :: n, p
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: means(:, :), sds(:, :)
        real(kind=dp), intent(out) :: x(n, p)

        integer :: i, i1, i2, j

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            do j = 1, p
                call fill_normal_segment(x(i1:i2, j), means(i, j), sds(i, j))
            end do
        end do
    end subroutine simulate_piecewise_normal_mv

    subroutine simulate_piecewise_poisson_1d(n, regime_starts, rates, x)
        !> Simulate a piecewise-Poisson univariate count series.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: rates(:)
        integer, intent(out) :: x(n)

        integer :: i, j, i1, i2

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            do j = i1, i2
                x(j) = sample_poisson_knuth(rates(i))
            end do
        end do
    end subroutine simulate_piecewise_poisson_1d

    subroutine simulate_piecewise_exponential_1d(n, regime_starts, means, x)
        !> Simulate a piecewise-Exponential univariate series with segment means `means`.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: means(:)
        real(kind=dp), intent(out) :: x(n)

        integer :: i, j, i1, i2
        real(kind=dp) :: u

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            do j = i1, i2
                call random_number(u)
                u = max(u, tiny(1.0_dp))
                x(j) = -means(i) * log(u)
            end do
        end do
    end subroutine simulate_piecewise_exponential_1d

    subroutine simulate_piecewise_gamma_1d(n, regime_starts, shape, means, x)
        !> Simulate a piecewise-Gamma univariate series with fixed shape and segment means `means`.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: shape
        real(kind=dp), intent(in) :: means(:)
        real(kind=dp), intent(out) :: x(n)

        integer :: i, j, i1, i2, k, shape_int
        real(kind=dp) :: u, scale

        shape_int = nint(shape)
        if (abs(shape - real(shape_int, dp)) > 1.0e-12_dp .or. shape_int <= 0) then
            error stop "simulate_piecewise_gamma_1d currently requires a positive integer shape"
        end if

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            scale = means(i) / shape
            do j = i1, i2
                x(j) = 0.0_dp
                do k = 1, shape_int
                    call random_number(u)
                    u = max(u, tiny(1.0_dp))
                    x(j) = x(j) - scale * log(u)
                end do
            end do
        end do
    end subroutine simulate_piecewise_gamma_1d

    subroutine simulate_piecewise_bernoulli_1d(n, regime_starts, probs, x)
        !> Simulate a piecewise-Bernoulli univariate binary series.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: probs(:)
        logical, intent(out) :: x(n)

        integer :: i, j, i1, i2
        real(kind=dp) :: u

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            do j = i1, i2
                call random_number(u)
                x(j) = (u < probs(i))
            end do
        end do
    end subroutine simulate_piecewise_bernoulli_1d

    subroutine simulate_piecewise_ar1_1d(n, regime_starts, phis, sigmas, x)
        !> Simulate a piecewise AR(1) univariate series with regime-local resets.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: phis(:), sigmas(:)
        real(kind=dp), intent(out) :: x(n)

        integer :: i, j, i1, i2
        real(kind=dp) :: prev

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            prev = 0.0_dp
            do j = i1, i2
                prev = phis(i) * prev + sigmas(i) * sample_standard_normal()
                x(j) = prev
            end do
        end do
    end subroutine simulate_piecewise_ar1_1d

    subroutine simulate_piecewise_linear_1d(n, regime_starts, intercepts, slopes, sigma, x)
        !> Simulate a piecewise linear univariate series with Gaussian noise.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: intercepts(:), slopes(:), sigma
        real(kind=dp), intent(out) :: x(n)

        integer :: i, i1, i2, j

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            do j = i1, i2
                x(j) = intercepts(i) + slopes(i) * real(j, dp) + sigma * sample_standard_normal()
            end do
        end do
    end subroutine simulate_piecewise_linear_1d

    subroutine simulate_piecewise_linear_mv_1d(n, p, regime_starts, intercepts, slopes, sigma, x)
        !> Simulate multiple piecewise linear series with Gaussian noise.
        integer, intent(in) :: n, p
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: intercepts(:, :), slopes(:, :), sigma
        real(kind=dp), intent(out) :: x(n, p)

        integer :: j

        do j = 1, p
            call simulate_piecewise_linear_1d(n, regime_starts, intercepts(:, j), slopes(:, j), sigma, x(:, j))
        end do
    end subroutine simulate_piecewise_linear_mv_1d

    subroutine simulate_piecewise_linear_harmonic_1d(n, regime_starts, intercepts, slopes, sigma, period, amp_sin, amp_cos, x)
        !> Simulate a piecewise linear series with a global harmonic seasonal component.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: intercepts(:), slopes(:), sigma, period, amp_sin, amp_cos
        real(kind=dp), intent(out) :: x(n)

        integer :: i, i1, i2, j
        real(kind=dp) :: omega, t

        omega = 2.0_dp * acos(-1.0_dp) / period
        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            do j = i1, i2
                t = real(j, dp)
                x(j) = intercepts(i) + slopes(i) * t + amp_sin * sin(omega * t) + amp_cos * cos(omega * t) + &
                       sigma * sample_standard_normal()
            end do
        end do
    end subroutine simulate_piecewise_linear_harmonic_1d

    subroutine simulate_piecewise_sine_1d(n, regime_starts, periods, sigma, x)
        !> Simulate a piecewise sinusoidal series with regime-specific periods.
        integer, intent(in) :: n
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: periods(:), sigma
        real(kind=dp), intent(out) :: x(n)

        integer :: i, i1, i2, j
        real(kind=dp) :: omega, t

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = n
            end if
            omega = 2.0_dp * acos(-1.0_dp) / periods(i)
            do j = i1, i2
                t = real(j - i1, dp)
                x(j) = sin(omega * t) + sigma * sample_standard_normal()
            end do
        end do
    end subroutine simulate_piecewise_sine_1d

    subroutine simulate_irregular_time_grid(n, start_time, dt_mean, dt_jitter, t)
        !> Simulate a strictly increasing irregular time grid.
        integer, intent(in) :: n
        real(kind=dp), intent(in) :: start_time, dt_mean, dt_jitter
        real(kind=dp), intent(out) :: t(n)

        integer :: i
        real(kind=dp) :: u

        t(1) = start_time
        do i = 2, n
            call random_number(u)
            t(i) = t(i - 1) + dt_mean + dt_jitter * (2.0_dp * u - 1.0_dp)
            if (t(i) <= t(i - 1)) t(i) = t(i - 1) + 0.1_dp * dt_mean
        end do
    end subroutine simulate_irregular_time_grid

    subroutine simulate_piecewise_linear_irregular_1d(t, regime_starts, intercepts, slopes, sigma, x)
        !> Simulate a piecewise linear series observed on an irregular time grid.
        real(kind=dp), intent(in) :: t(:)
        integer, intent(in) :: regime_starts(:)
        real(kind=dp), intent(in) :: intercepts(:), slopes(:), sigma
        real(kind=dp), intent(out) :: x(size(t))

        integer :: i, i1, i2, j

        do i = 1, size(regime_starts)
            i1 = regime_starts(i) + 1
            if (i < size(regime_starts)) then
                i2 = regime_starts(i + 1)
            else
                i2 = size(t)
            end if
            do j = i1, i2
                x(j) = intercepts(i) + slopes(i) * t(j) + sigma * sample_standard_normal()
            end do
        end do
    end subroutine simulate_piecewise_linear_irregular_1d

    subroutine fill_normal_segment(x, mu, sigma)
        real(kind=dp), intent(out) :: x(:)
        real(kind=dp), intent(in) :: mu, sigma

        integer :: i

        do i = 1, size(x)
            x(i) = mu + sigma * sample_standard_normal()
        end do
    end subroutine fill_normal_segment

    function sample_standard_normal() result(z)
        real(kind=dp) :: z
        real(kind=dp) :: u1, u2, r, theta

        if (has_spare_normal) then
            z = spare_normal
            has_spare_normal = .false.
            return
        end if

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, 1.0e-12_dp)
        r = sqrt(-2.0_dp * log(u1))
        theta = 2.0_dp * acos(-1.0_dp) * u2
        z = r * cos(theta)
        spare_normal = r * sin(theta)
        has_spare_normal = .true.
    end function sample_standard_normal

    function sample_poisson_knuth(lambda) result(k)
        real(kind=dp), intent(in) :: lambda
        integer :: k
        real(kind=dp) :: l, p, u

        if (lambda <= 0.0_dp) then
            k = 0
            return
        end if

        l = exp(-lambda)
        k = 0
        p = 1.0_dp
        do
            k = k + 1
            call random_number(u)
            p = p * u
            if (p <= l) exit
        end do
        k = k - 1
    end function sample_poisson_knuth

end module compare_sim_mod
