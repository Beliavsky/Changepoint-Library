program xsim_dp_poly
! simulate univariate polynomial changepoints and estimate them by DP.poly
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_mod, only: solve_dp_poly_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 150 ! number of observations
integer, parameter :: r = 2 ! polynomial order
real(kind=dp), parameter :: gamma = 1.0_dp ! l0 penalty parameter
integer, parameter :: delta = 5 ! minimum spacing

real(kind=dp) :: y(n)
real(kind=dp), parameter :: coef1(3) = [0.0_dp, 1.5_dp, -2.5_dp]
real(kind=dp), parameter :: coef2(3) = [2.5_dp, -8.0_dp, 10.0_dp]
real(kind=dp), parameter :: coef3(3) = [-2.0_dp, 5.0_dp, -2.5_dp]
real(kind=dp), allocatable :: yhat(:)
integer, allocatable :: partition(:), cpt_hat(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_response(y)
call solve_dp_poly_1d(y, r, gamma, delta, partition, yhat, cpt_hat)

print *, "n                  = ", n
print *, "r                  = ", r
print *, "gamma              = ", gamma
print *, "delta              = ", delta
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          = "
end if
print *, "yhat checksum      = ", sum(yhat)

deallocate(partition, cpt_hat, yhat)
call print_wall_time(t_start)

contains

    subroutine fill_response(y_out)
        real(kind=dp), intent(out) :: y_out(:)
        integer :: i
        real(kind=dp) :: xval

        do i = 1, 50
            xval = real(i, dp) / real(n, dp)
            y_out(i) = coef1(1) + coef1(2) * xval + coef1(3) * xval * xval + 0.05_dp * sample_standard_normal()
        end do
        do i = 51, 100
            xval = real(i, dp) / real(n, dp)
            y_out(i) = coef2(1) + coef2(2) * xval + coef2(3) * xval * xval + 0.05_dp * sample_standard_normal()
        end do
        do i = 101, 150
            xval = real(i, dp) / real(n, dp)
            y_out(i) = coef3(1) + coef3(2) * xval + coef3(3) * xval * xval + 0.05_dp * sample_standard_normal()
        end do
    end subroutine fill_response

    real(kind=dp) function sample_standard_normal() result(z)
        real(kind=dp) :: u1, u2
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, 1.0e-12_dp)
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
    end function sample_standard_normal

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_dp_poly
