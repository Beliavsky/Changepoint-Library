program xcosts_mean_variance
! compare several dynp costs on a univariate series with changes in mean and variance
use kind_mod, only: dp
use changepoint_mod, only: solve_dynp_cost_1d
implicit none

integer, parameter :: n = 600  ! total number of observations
integer, parameter :: n_bkps = 3  ! number of changepoints to estimate
integer, parameter :: max_m = n_bkps + 1  ! number of segments in the fitted partition
integer, parameter :: min_seg_len = 30  ! minimum allowed segment length
integer, parameter :: seed = 11  ! random seed used to simulate the signal
integer, dimension(max_m), parameter :: true_bkps = [150, 300, 450, n]  ! true segment endpoints
real(kind=dp), dimension(max_m), parameter :: means = [0.0_dp, 2.5_dp, -1.0_dp, 1.5_dp]  ! segment means
real(kind=dp), dimension(max_m), parameter :: sds = [1.0_dp, 2.2_dp, 0.6_dp, 1.8_dp]  ! segment standard deviations

real(kind=dp) :: z(n), u, t0, t1, wall_start, wall_end
integer :: i, seed_size
integer, allocatable :: seed_put(:)
integer :: seg_ends_l1(max_m), seg_ends_l2(max_m), seg_ends_normal(max_m), seg_ends_rbf(max_m)

call cpu_time(wall_start)
call random_seed(size=seed_size)
allocate(seed_put(seed_size))
seed_put = seed
call random_seed(put=seed_put)
deallocate(seed_put)

do i = 1, n
    call random_number(u)
    if (u <= 1.0e-12_dp) u = 1.0e-12_dp
    if (i <= true_bkps(1)) then
        z(i) = means(1) + sds(1) * normal_from_uniform(u)
    else if (i <= true_bkps(2)) then
        z(i) = means(2) + sds(2) * normal_from_uniform(u)
    else if (i <= true_bkps(3)) then
        z(i) = means(3) + sds(3) * normal_from_uniform(u)
    else
        z(i) = means(4) + sds(4) * normal_from_uniform(u)
    end if
end do

print *, "n           =", n
print *, "true bkps   =", true_bkps
print *, "segment means =", means
print *, "segment sds   =", sds
print *, "n_bkps      =", n_bkps
print *, "min_seg_len =", min_seg_len
print *

call cpu_time(t0)
call solve_dynp_cost_1d(z, max_m, "l1", seg_ends_l1, min_seg_len)
call cpu_time(t1)
print *, "model       = l1"
print *, "note        = robust location shifts"
print *, "estimated   =", seg_ends_l1
write(*,'(a,f8.3)') "elapsed seconds = ", t1 - t0
print *

call cpu_time(t0)
call solve_dynp_cost_1d(z, max_m, "l2", seg_ends_l2, min_seg_len)
call cpu_time(t1)
print *, "model       = l2"
print *, "note        = mean shifts"
print *, "estimated   =", seg_ends_l2
write(*,'(a,f8.3)') "elapsed seconds = ", t1 - t0
print *

call cpu_time(t0)
call solve_dynp_cost_1d(z, max_m, "normal", seg_ends_normal, min_seg_len)
call cpu_time(t1)
print *, "model       = normal"
print *, "note        = mean and variance shifts"
print *, "estimated   =", seg_ends_normal
write(*,'(a,f8.3)') "elapsed seconds = ", t1 - t0
print *

call cpu_time(t0)
call solve_dynp_cost_1d(z, max_m, "rbf", seg_ends_rbf, min_seg_len)
call cpu_time(t1)
print *, "model       = rbf"
print *, "note        = general distribution shifts"
print *, "estimated   =", seg_ends_rbf
write(*,'(a,f8.3)') "elapsed seconds = ", t1 - t0
print *

call cpu_time(wall_end)
write(*,'(a,f8.3)') "wall time elapsed (s) = ", wall_end - wall_start

contains

real(kind=dp) function normal_from_uniform(u1) result(z1)
real(kind=dp), intent(in) :: u1
real(kind=dp) :: u2

call random_number(u2)
if (u2 <= 1.0e-12_dp) u2 = 1.0e-12_dp
z1 = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)

end function normal_from_uniform

end program xcosts_mean_variance
