program xsim_reg_amoc
! simulate a regression series and estimate one changepoint by changepoint::cpt.reg(method="AMOC")
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoint_mod, only: solve_amoc_reg_norm_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: seed = 1 ! RNG seed used for the simulated series
integer, parameter :: n = 200 ! number of observations in the simulated series
integer, parameter :: true_cp = 100 ! true changepoint in the regression coefficients
integer, parameter :: minseglen = 5 ! minimum segment length used by cpt.reg AMOC

real(kind=dp), allocatable :: data(:, :)
real(kind=dp), allocatable :: noise(:)
real(kind=dp) :: penalty_value, null_rss, alt_rss
integer :: i, cpt
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(data(n, 3), noise(n))
call seed_rng_fixed(seed)
call random_number(noise)
do i = 1, n
    data(i, 2) = 1.0_dp
    data(i, 3) = real(i, dp)
end do
call fill_normal_noise(noise)
do i = 1, n
    if (i <= true_cp) then
        data(i, 1) = 0.0_dp + 1.0_dp * data(i, 3) + noise(i)
    else
        data(i, 1) = 50.0_dp + 0.25_dp * data(i, 3) + noise(i)
    end if
end do

call solve_amoc_reg_norm_1d(data, cpt, penalty_value, null_rss, alt_rss, penalty="SIC", minseglen=minseglen)

print *, "n                  =", n
print *, "penalty            =", "SIC"
print *, "minseglen          =", minseglen
print *, "true changepoint   =", true_cp
print *, "estimated          =", cpt
print *, "pen.value          =", penalty_value

deallocate(data, noise)
call print_wall_time(t_start)

contains

    subroutine fill_normal_noise(z)
        real(kind=dp), intent(inout) :: z(:)
        integer :: j
        real(kind=dp) :: u1, u2, r, theta

        j = 1
        do while (j <= size(z))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, 1.0e-12_dp)
            r = sqrt(-2.0_dp * log(u1))
            theta = 2.0_dp * acos(-1.0_dp) * u2
            z(j) = r * cos(theta)
            if (j + 1 <= size(z)) z(j + 1) = r * sin(theta)
            j = j + 2
        end do
    end subroutine fill_normal_noise

end program xsim_reg_amoc
