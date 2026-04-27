program xsim_probe_bcp_reg_file
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use bcp_pkg_mod, only: solve_bcp_regression_1x_probe
implicit none

integer, parameter :: burnin = 50, mcmc = 200, nreg = 2
integer, parameter :: wh_seed(3) = [2439, 10153, 8035]
real(kind=dp), parameter :: p0 = 0.2_dp, w0(2) = [0.2_dp, 0.2_dp], d = 10.0_dp
character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:), posterior_prob(:)
integer, allocatable :: blocks(:), rho_counts(:), blocks_head(:), top_idx(:)

call get_data_file_arg("xbcp_reg_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1)))
y = dat(:, 1)
x = dat(:, 2)

call solve_bcp_regression_1x_probe(y, x, burnin, mcmc, p0, w0, d, nreg, wh_seed, posterior_prob, blocks, rho_counts, blocks_head)
call top_k_prob(posterior_prob(:size(posterior_prob) - 1), 5, top_idx)

print *, "probe pp checksum =", sum(posterior_prob(:size(posterior_prob) - 1))
print *, "probe rho counts checksum =", real(sum(rho_counts(:size(rho_counts) - 1)), dp)
print *, "probe top idx =", join_ints(top_idx)
print *, "probe blocks head =", join_ints(blocks_head)
print *, "probe rho counts head =", join_ints(rho_counts(:20))
print *, "probe rho counts around 40 =", join_ints(rho_counts(35:45))
print *, "probe rho counts around 80 =", join_ints(rho_counts(75:85))
print *, "probe rho nz =", join_indexed_counts(rho_counts(:size(rho_counts) - 1))

contains

subroutine top_k_prob(pp, k, idx)
real(kind=dp), intent(in) :: pp(:)
integer, intent(in) :: k
integer, allocatable, intent(out) :: idx(:)
logical, allocatable :: used(:)
integer :: i, j, best
real(kind=dp) :: best_val

allocate(idx(k), used(size(pp)))
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
    used(best) = .true.
end do
end subroutine top_k_prob

function join_ints(vals) result(out)
integer, intent(in) :: vals(:)
character(len=512) :: out
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

function join_indexed_counts(vals) result(out)
integer, intent(in) :: vals(:)
character(len=4096) :: out
character(len=32) :: token
integer :: i, pos

out = ""
pos = 1
do i = 1, size(vals)
    if (vals(i) == 0) cycle
    write(token, '(I0, A, I0)') i, ":", vals(i)
    if (pos > 1) then
        out(pos:pos) = " "
        pos = pos + 1
    end if
    out(pos:) = trim(token)
    pos = len_trim(out) + 1
end do
end function join_indexed_counts

end program xsim_probe_bcp_reg_file
