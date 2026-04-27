program xsim_mean_bottomup
! simulate a univariate mean-shift series and estimate changepoints by bottom-up segmentation
use kind_mod, only: dp
use changepoint_mod, only: solve_bottomup_mean_shift_1d, refine_bkps_mean_shift_1d
implicit none

integer :: n, n_bkps, min_seg_len, jump
integer, allocatable :: bkps_initial(:), bkps(:)
integer :: i, b1, b2, b3, b4
real(kind=dp), allocatable :: z(:)
real(kind=dp) :: u, sigma
real(kind=dp), dimension(5) :: mu

n = 100000
n_bkps = 4
min_seg_len = max(5, n / 25)
jump = 5
sigma = 1.0_dp
mu = [0.0_dp, 3.0_dp, -1.0_dp, 2.0_dp, 0.5_dp]

b1 = max(2, min(n - 4, nint(0.2_dp * n)))
b2 = max(b1 + 1, min(n - 3, nint(0.4_dp * n)))
b3 = max(b2 + 1, min(n - 2, nint(0.6_dp * n)))
b4 = max(b3 + 1, min(n - 1, nint(0.8_dp * n)))

allocate(z(n))
call random_seed()

do i = 1, n
    call random_number(u)
    if (u <= 1.0e-12_dp) u = 1.0e-12_dp
    if (i <= b1) then
        z(i) = mu(1) + sigma * normal_from_uniform(u)
    else if (i <= b2) then
        z(i) = mu(2) + sigma * normal_from_uniform(u)
    else if (i <= b3) then
        z(i) = mu(3) + sigma * normal_from_uniform(u)
    else if (i <= b4) then
        z(i) = mu(4) + sigma * normal_from_uniform(u)
    else
        z(i) = mu(5) + sigma * normal_from_uniform(u)
    end if
end do

call solve_bottomup_mean_shift_1d(z, n_bkps, bkps_initial, min_seg_len, jump)
call refine_bkps_mean_shift_1d(z, bkps_initial, bkps, min_seg_len)

print *, "n           =", n
print *, "true cps    =", b1, b2, b3, b4
if (size(bkps_initial) > 0) print *, "initial     =", bkps_initial
if (size(bkps) > 0) print *, "estimated   =", bkps

deallocate(z, bkps_initial, bkps)

contains

real(kind=dp) function normal_from_uniform(u1) result(z1)
real(kind=dp), intent(in) :: u1
real(kind=dp) :: u2

call random_number(u2)
if (u2 <= 1.0e-12_dp) u2 = 1.0e-12_dp
z1 = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)

end function normal_from_uniform

end program xsim_mean_bottomup
