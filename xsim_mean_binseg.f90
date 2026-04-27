program xsim_mean_binseg
! simulate a univariate mean-shift series and estimate changepoints by binary segmentation
use kind_mod, only: dp
use changepoint_binseg_mod, only: solve_binseg_mean_shift
implicit none

integer :: n, max_cp, ncp_found, min_seg_len
integer, allocatable :: cp(:)
integer :: i, ios, unit_out
real(kind=dp), allocatable :: z(:)
real(kind=dp) :: u, frac1, frac2, mu1, mu2, mu3, sigma, min_gain
integer :: b1, b2
logical :: write_data
character(len=*), parameter :: outfile = "xsim_mean_binseg_data.txt"

n = 100000
max_cp = 5
min_seg_len = 20
frac1 = 1.0_dp / 3.0_dp
frac2 = 7.0_dp / 12.0_dp
mu1 = 0.0_dp
mu2 = 3.0_dp
mu3 = -1.0_dp
sigma = 1.0_dp
min_gain = 0.0_dp
write_data = .false.

b1 = nint(frac1 * n)
b2 = nint(frac2 * n)
b1 = max(2, min(n - 2, b1))
b2 = max(b1 + 1, min(n - 1, b2))

allocate(z(n), cp(max_cp))

call random_seed()

do i = 1, n
    call random_number(u)
    if (u <= 1.0e-12_dp) u = 1.0e-12_dp
    if (i <= b1) then
        z(i) = mu1 + sigma * normal_from_uniform(u)
    else if (i <= b2) then
        z(i) = mu2 + sigma * normal_from_uniform(u)
    else
        z(i) = mu3 + sigma * normal_from_uniform(u)
    end if
end do

call solve_binseg_mean_shift(z, max_cp, cp, ncp_found, min_seg_len, min_gain)

print *, "n           =", n
print *, "true cps    =", b1, b2
print *, "ncp_found   =", ncp_found
if (ncp_found > 0) print *, "estimated   =", cp(1:ncp_found)

if (write_data) then
    open(newunit=unit_out, file=outfile, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "could not open output file"
    write(unit_out,"(a)") "# i z"
    do i = 1, n
        write(unit_out,"(i0,1x,f12.6)") i, z(i)
    end do
    close(unit_out)
    print *, "wrote data to ", trim(outfile)
end if

deallocate(z, cp)

contains

real(kind=dp) function normal_from_uniform(u1) result(z1)
! transform one uniform random number into one standard normal draw
real(kind=dp), intent(in) :: u1
real(kind=dp) :: u2
call random_number(u2)
if (u2 <= 1.0e-12_dp) u2 = 1.0e-12_dp
z1 = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
end function normal_from_uniform

end program xsim_mean_binseg
