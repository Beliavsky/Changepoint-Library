program xsim_online_univar
! simulate a univariate mean-change series and detect the changepoint by online.univar with a supplied threshold vector
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_mod, only: solve_online_univar_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 150 ! number of observations
real(kind=dp), parameter :: b_value = 3.5_dp ! constant supplied threshold value

real(kind=dp), allocatable :: y(:), b_vec(:)
integer(kind=long_int) :: t_start
integer :: cpt_hat

call system_clock(t_start)
allocate(y(n), b_vec(n - 1))
call seed_rng_fixed(0)
call fill_normal(y, 0.0_dp, 1.0_dp)
y(101:n) = y(101:n) + 1.5_dp
b_vec = b_value
call solve_online_univar_1d(y, b_vec, cpt_hat)

print *, "n                  =", n
print *, "b_value            =", b_value
print *, "estimated          =", cpt_hat

deallocate(y, b_vec)
call print_wall_time(t_start)

contains

    subroutine fill_normal(x, mu, sigma)
        real(kind=dp), intent(out) :: x(:)
        real(kind=dp), intent(in) :: mu, sigma
        real(kind=dp) :: u1, u2, r, theta
        integer :: i
        i = 1
        do while (i <= size(x))
            call random_number(u1)
            call random_number(u2)
            if (u1 <= 0.0_dp) cycle
            r = sqrt(-2.0_dp * log(u1))
            theta = 2.0_dp * acos(-1.0_dp) * u2
            x(i) = mu + sigma * r * cos(theta)
            if (i + 1 <= size(x)) x(i + 1) = mu + sigma * r * sin(theta)
            i = i + 2
        end do
    end subroutine fill_normal

end program xsim_online_univar
