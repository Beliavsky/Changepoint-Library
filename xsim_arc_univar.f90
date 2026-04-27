program xsim_arc_univar
! simulate a contaminated univariate mean-change series and estimate changepoints by ARC
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_mod, only: solve_arc_univar_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 1000 ! number of observations
integer, parameter :: h = 120 ! block size used by ARC
integer, parameter :: block_num = 1 ! number of blocks used in local maximiser search
real(kind=dp), parameter :: epsilon = 0.1_dp ! contamination proportion
logical, parameter :: gaussian = .true. ! use Gaussian threshold calibration

real(kind=dp), allocatable :: y(:), cusum(:)
integer, allocatable :: cpt_hat(:)
integer, dimension(10) :: outlier_idx
real(kind=dp), dimension(10) :: outlier_shift
real(kind=dp) :: lambda
integer(kind=long_int) :: t_start
integer :: i

call system_clock(t_start)
allocate(y(n))
call seed_rng_fixed(0)
call fill_normal(y(1:500), 0.0_dp, 1.0_dp)
call fill_normal(y(501:n), 1.0_dp, 1.0_dp)
outlier_idx = [80, 150, 220, 310, 420, 560, 640, 730, 820, 910]
outlier_shift = [12.0_dp, -11.0_dp, 10.0_dp, -13.0_dp, 11.0_dp, 12.0_dp, -10.0_dp, 13.0_dp, -12.0_dp, 11.0_dp]
do i = 1, size(outlier_idx)
    y(outlier_idx(i)) = y(outlier_idx(i)) + outlier_shift(i)
end do

call solve_arc_univar_1d(y, h, block_num, epsilon, gaussian, cpt_hat, lambda, cusum)

print *, "n                  =", n
print *, "h                  =", h
print *, "block_num          =", block_num
print *, "epsilon            =", epsilon
print *, "gaussian           =", gaussian
print *, "lambda             =", lambda
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if

deallocate(y, cusum, cpt_hat)
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

    subroutine print_int_list(x)
        integer, intent(in) :: x(:)
        integer :: i
        do i = 1, size(x)
            write (*,'(1X,I0)', advance='no') x(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_arc_univar
