program xsim_dp_univar
! simulate a univariate mean-change series and estimate changepoints by l0-penalized dynamic programming
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed, simulate_piecewise_normal_1d
use changepoints_pkg_mod, only: solve_dp_univar_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 300 ! number of observations
integer, parameter :: delta = 5 ! minimum spacing parameter in DP.univar
real(kind=dp), parameter :: gamma = 5.0_dp ! l0 penalty used by DP.univar

integer, dimension(4) :: regime_starts
real(kind=dp), dimension(4) :: means, sds
real(kind=dp), allocatable :: y(:), yhat(:)
integer, allocatable :: partition(:), cpt_hat(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
allocate(y(n))

call seed_rng_fixed(0)
regime_starts = [0, 20, 50, 170]
means = [0.0_dp, 2.0_dp, 0.0_dp, -2.0_dp]
sds = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
call simulate_piecewise_normal_1d(n, regime_starts, means, sds, y)
call solve_dp_univar_1d(y, gamma, delta, partition, yhat, cpt_hat)

print *, "n                  =", n
print *, "gamma              =", gamma
print *, "delta              =", delta
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
print *, "n changepoints     =", size(cpt_hat)
print *, "partition last     =", partition(size(partition))
print *, "yhat checksum      =", sum(yhat)

deallocate(y, partition, yhat, cpt_hat)
call print_wall_time(t_start)

contains

    subroutine print_int_list(x)
        integer, intent(in) :: x(:)
        integer :: i
        do i = 1, size(x)
            write (*,'(1X,I0)', advance='no') x(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_dp_univar
