program xsim_bcp_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use bcp_pkg_mod, only: solve_bcp_series_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: burnin = 50, mcmc = 200
integer, parameter :: wh_seed(3) = [2439, 10153, 8035]
real(kind=dp), parameter :: p0 = 0.2_dp, w0 = 0.2_dp
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:), posterior_mean(:), posterior_prob(:)
integer, allocatable :: blocks(:), top_idx(:)
real(kind=dp), allocatable :: top_val(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xbcp_data.txt", data_file)
call read_series_file(data_file, x)

call solve_bcp_series_1d(x, burnin, mcmc, p0, w0, wh_seed, posterior_mean, posterior_prob, blocks)
call top_k_prob(posterior_prob(:size(posterior_prob) - 1), 5, top_idx, top_val)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(x)
print *, "burnin             = ", burnin
print *, "mcmc               = ", mcmc
print *, "p0                 = ", p0
print *, "w0                 = ", w0
print *, "rng                = Wichmann-Hill"
print *, "seed               = ", join_ints(wh_seed)
print *, "pm checksum        = ", sum(posterior_mean)
print *, "pp checksum        = ", sum(posterior_prob(:size(posterior_prob) - 1))
print *, "top idx            = ", join_ints(top_idx)
print *, "top val            = ", join_reals(top_val)

call print_wall_time(t_start)

contains

subroutine top_k_prob(pp, k, idx, vals)
real(kind=dp), intent(in) :: pp(:)
integer, intent(in) :: k
integer, allocatable, intent(out) :: idx(:)
real(kind=dp), allocatable, intent(out) :: vals(:)
logical, allocatable :: used(:)
integer :: i, j, best
real(kind=dp) :: best_val

allocate(idx(k), vals(k), used(size(pp)))
used = .false.
do i = 1, k
    best = 1
    best_val = -huge(1.0_dp)
    do j = 1, size(pp)
        if (used(j)) cycle
        if (pp(j) > best_val) then
            best = j
            best_val = pp(j)
        end if
    end do
    idx(i) = best
    vals(i) = pp(best)
    used(best) = .true.
end do
end subroutine top_k_prob

function join_ints(vals) result(out)
integer, intent(in) :: vals(:)
character(len=256) :: out
integer :: i, pos

out = ""
pos = 1
do i = 1, size(vals)
    write(out(pos:), '(I0)') vals(i)
    pos = len_trim(out) + 1
    if (i < size(vals)) then
        out(pos:pos) = " "
        pos = pos + 1
    end if
end do
end function join_ints

function join_reals(vals) result(out)
real(kind=dp), intent(in) :: vals(:)
character(len=512) :: out
character(len=32) :: buf
integer :: i, pos

out = ""
pos = 1
do i = 1, size(vals)
    write(buf, '(F0.15)') vals(i)
    out(pos:) = trim(adjustl(buf))
    pos = len_trim(out) + 1
    if (i < size(vals)) then
        out(pos:pos) = " "
        pos = pos + 1
    end if
end do
end function join_reals

end program xsim_bcp_file
