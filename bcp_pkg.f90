module bcp_pkg_mod
use kind_mod, only: dp
implicit none
private
public :: solve_bcp_series_1d, solve_bcp_series_2d, solve_bcp_regression_1x, solve_bcp_regression_1x_probe

integer, parameter :: wh_m1 = 30269, wh_m2 = 30307, wh_m3 = 30323
real(kind=dp), parameter :: tiny_dp = 1.0e-12_dp

contains

subroutine solve_bcp_series_1d(y, burnin, mcmc, p0, w0, wh_seed, posterior_mean, posterior_prob, blocks)
real(kind=dp), intent(in) :: y(:), p0, w0
integer, intent(in) :: burnin, mcmc
integer, intent(in) :: wh_seed(3)
real(kind=dp), allocatable, intent(out) :: posterior_mean(:), posterior_prob(:)
integer, allocatable, intent(out) :: blocks(:)

integer :: n, mm, m, j, k
real(kind=dp), allocatable :: cumy(:), priors(:), bmean(:)
real(kind=dp), allocatable :: rho(:)
integer, allocatable :: bend(:), bsize(:)
real(kind=dp), allocatable :: bz(:)
real(kind=dp) :: ybar, ysqall, b_scalar, w_scalar, wstar, tmp_mean
integer :: b_count
integer :: seed(3)

n = size(y)
mm = burnin + mcmc
seed = wh_seed

allocate(posterior_mean(n), posterior_prob(n), blocks(mm))
allocate(cumy(n), priors(max(1, n - 3)))
allocate(rho(n), bend(n), bsize(n), bz(n), bmean(n))

posterior_mean = 0.0_dp
posterior_prob = 0.0_dp
blocks = 0

cumy = cumulative_sum(y)
ybar = sum(y) / real(n, dp)
ysqall = sum(y * y)

call build_priors(n, p0, priors)
call init_step_1d(cumy, ybar, ysqall, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)

do m = 1, mm
    call bcp_pass_1d(cumy, priors, w0, seed, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)
    blocks(m) = b_count
    if (m <= burnin) cycle
    wstar = compute_wstar_1d(b_scalar, w_scalar, b_count, n, w0)
    k = 1
    do j = 1, n
        posterior_prob(j) = posterior_prob(j) + rho(j)
        tmp_mean = bmean(k) * (1.0_dp - wstar) + ybar * wstar
        posterior_mean(j) = posterior_mean(j) + tmp_mean
        if (nint(rho(j)) == 1) k = k + 1
    end do
end do

posterior_mean = posterior_mean / real(mcmc, dp)
posterior_prob = posterior_prob / real(mcmc, dp)
end subroutine solve_bcp_series_1d

subroutine solve_bcp_series_2d(y, burnin, mcmc, p0, w0, wh_seed, posterior_mean, posterior_prob, blocks)
real(kind=dp), intent(in) :: y(:, :), p0, w0
integer, intent(in) :: burnin, mcmc
integer, intent(in) :: wh_seed(3)
real(kind=dp), allocatable, intent(out) :: posterior_mean(:, :), posterior_prob(:)
integer, allocatable, intent(out) :: blocks(:)

integer :: n, p, mm, m, j, k
real(kind=dp), allocatable :: cumy(:, :), priors(:), bmean(:, :)
real(kind=dp), allocatable :: rho(:)
integer, allocatable :: bend(:), bsize(:)
real(kind=dp), allocatable :: bz(:)
real(kind=dp) :: ybar, ysqall, b_scalar, w_scalar, wstar
integer :: b_count
integer :: seed(3)

n = size(y, 1)
p = size(y, 2)
mm = burnin + mcmc
seed = wh_seed

allocate(posterior_mean(n, p), posterior_prob(n), blocks(mm))
allocate(cumy(n, p), priors(max(1, n - 3)))
allocate(rho(n), bend(n), bsize(n), bz(n), bmean(n, p))

posterior_mean = 0.0_dp
posterior_prob = 0.0_dp
blocks = 0

cumy = cumulative_sum_mat(y)
ybar = sum(y) / real(n * p, dp)
ysqall = sum(y * y)

call build_priors(n, p0, priors)
call init_step_k(cumy, ybar, ysqall, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)

do m = 1, mm
    call bcp_pass_k(cumy, priors, w0, seed, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)
    blocks(m) = b_count
    if (m <= burnin) cycle
    wstar = compute_wstar_k(b_scalar, w_scalar, b_count, n, p, w0)
    k = 1
    do j = 1, n
        posterior_prob(j) = posterior_prob(j) + rho(j)
        posterior_mean(j, :) = posterior_mean(j, :) + bmean(k, :) * (1.0_dp - wstar) + ybar * wstar
        if (nint(rho(j)) == 1) k = k + 1
    end do
end do

posterior_mean = posterior_mean / real(mcmc, dp)
posterior_prob = posterior_prob / real(mcmc, dp)
end subroutine solve_bcp_series_2d

subroutine solve_bcp_regression_1x(y, x, burnin, mcmc, p0, w0_vec, d, nreg, wh_seed, posterior_mean, posterior_prob, blocks)
real(kind=dp), intent(in) :: y(:), x(:), p0, w0_vec(2), d
integer, intent(in) :: burnin, mcmc, nreg
integer, intent(in) :: wh_seed(3)
real(kind=dp), allocatable, intent(out) :: posterior_mean(:), posterior_prob(:)
integer, allocatable, intent(out) :: blocks(:)

integer :: n, mm, m, b_count, seed(3)
real(kind=dp) :: ybar, cumy_end, cumysq_end, bigb, wsum, logc, qsum, ksum, lik, wstar
real(kind=dp), allocatable :: cumy(:), cumysq(:), cumx(:), cumxx(:), cumxy(:), prior_logs(:)
real(kind=dp), allocatable :: rho(:), bz(:), blogc(:), bk(:), bq(:), w(:)
integer, allocatable :: bend(:), bsize(:), btau(:)

n = size(y)
mm = burnin + mcmc
seed = wh_seed

allocate(posterior_mean(n), posterior_prob(n), blocks(mm))
allocate(cumy(n), cumysq(n), cumx(n), cumxx(n), cumxy(n), prior_logs(max(1, n - 3)))
allocate(rho(n), bz(n), blogc(n), bk(n), bq(n), w(2))
allocate(bend(n), bsize(n), btau(n))

posterior_mean = 0.0_dp
posterior_prob = 0.0_dp
blocks = 0

call build_prefix_reg_1x(y, x, cumy, cumysq, cumx, cumxx, cumxy)
call build_priors_reg(n, p0, prior_logs)
ybar = sum(y) / real(n, dp)
cumy_end = cumy(n)
cumysq_end = cumysq(n)

rho = 0.0_dp
rho(n) = 1.0_dp
bend = 0
bsize = 0
btau = 0
bz = 0.0_dp
blogc = 0.0_dp
bk = 0.0_dp
bq = 0.0_dp
w(1) = w0_vec(1)
w(2) = 0.5_dp * w0_vec(2)

bend(1) = n
bsize(1) = n
btau(1) = 0
bz(1) = cumy_end * cumy_end / real(n, dp)
bk(1) = logkcalc_reg(n, 0, d, nreg)
blogc(1) = 0.0_dp
bq(1) = 0.0_dp
b_count = 1
bigb = 0.0_dp
wsum = cumysq_end - bz(1)
logc = 0.0_dp
qsum = 0.0_dp
ksum = bk(1)
lik = likelihood_reg(bigb, wsum, b_count, n, w(1), logc, qsum, ksum, prior_logs)

do m = 1, mm
    call bcp_pass_reg_1x(cumy, cumx, cumxx, cumxy, prior_logs, w0_vec, d, nreg, seed, &
        rho, bend, bsize, btau, bz, blogc, bk, bq, w, bigb, wsum, logc, qsum, ksum, lik, b_count)
    blocks(m) = b_count
    if (m <= burnin) cycle
    wstar = compute_wstar_reg(bigb, wsum - qsum, b_count, n, w(1))
    call accumulate_posterior_reg_1x(x, cumy, cumx, cumxx, cumxy, ybar, w, wstar, rho, bend, bsize, btau, posterior_mean, posterior_prob)
end do

posterior_mean = posterior_mean / real(mcmc, dp)
posterior_prob = posterior_prob / real(mcmc, dp)
end subroutine solve_bcp_regression_1x

subroutine solve_bcp_regression_1x_probe(y, x, burnin, mcmc, p0, w0_vec, d, nreg, wh_seed, &
    posterior_prob, blocks, rho_counts, blocks_head)
real(kind=dp), intent(in) :: y(:), x(:), p0, w0_vec(2), d
integer, intent(in) :: burnin, mcmc, nreg
integer, intent(in) :: wh_seed(3)
real(kind=dp), allocatable, intent(out) :: posterior_prob(:)
integer, allocatable, intent(out) :: blocks(:), rho_counts(:), blocks_head(:)

integer :: n, mm, m, b_count, seed(3), nhead
real(kind=dp) :: ybar, cumy_end, cumysq_end, bigb, wsum, logc, qsum, ksum, lik
real(kind=dp), allocatable :: cumy(:), cumysq(:), cumx(:), cumxx(:), cumxy(:), prior_logs(:)
real(kind=dp), allocatable :: rho(:), bz(:), blogc(:), bk(:), bq(:), w(:)
integer, allocatable :: bend(:), bsize(:), btau(:)

n = size(y)
mm = burnin + mcmc
seed = wh_seed

allocate(posterior_prob(n), blocks(mm), rho_counts(n))
allocate(cumy(n), cumysq(n), cumx(n), cumxx(n), cumxy(n), prior_logs(max(1, n - 3)))
allocate(rho(n), bz(n), blogc(n), bk(n), bq(n), w(2))
allocate(bend(n), bsize(n), btau(n))

posterior_prob = 0.0_dp
blocks = 0
rho_counts = 0

call build_prefix_reg_1x(y, x, cumy, cumysq, cumx, cumxx, cumxy)
call build_priors_reg(n, p0, prior_logs)
ybar = sum(y) / real(n, dp)
cumy_end = cumy(n)
cumysq_end = cumysq(n)

rho = 0.0_dp
rho(n) = 1.0_dp
bend = 0
bsize = 0
btau = 0
bz = 0.0_dp
blogc = 0.0_dp
bk = 0.0_dp
bq = 0.0_dp
w(1) = w0_vec(1)
w(2) = 0.5_dp * w0_vec(2)

bend(1) = n
bsize(1) = n
btau(1) = 0
bz(1) = cumy_end * cumy_end / real(n, dp)
bk(1) = logkcalc_reg(n, 0, d, nreg)
blogc(1) = 0.0_dp
bq(1) = 0.0_dp
b_count = 1
bigb = 0.0_dp
wsum = cumysq_end - bz(1)
logc = 0.0_dp
qsum = 0.0_dp
ksum = bk(1)
lik = likelihood_reg(bigb, wsum, b_count, n, w(1), logc, qsum, ksum, prior_logs)

do m = 1, mm
    call bcp_pass_reg_1x(cumy, cumx, cumxx, cumxy, prior_logs, w0_vec, d, nreg, seed, &
        rho, bend, bsize, btau, bz, blogc, bk, bq, w, bigb, wsum, logc, qsum, ksum, lik, b_count)
    blocks(m) = b_count
    if (m <= burnin) cycle
    posterior_prob = posterior_prob + rho
    rho_counts = rho_counts + nint(rho)
end do

posterior_prob = posterior_prob / real(mcmc, dp)
nhead = min(20, size(blocks))
allocate(blocks_head(nhead))
blocks_head = blocks(:nhead)
end subroutine solve_bcp_regression_1x_probe

subroutine init_step_1d(cumy, ybar, ysqall, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)
real(kind=dp), intent(in) :: cumy(:), ybar, ysqall
real(kind=dp), intent(out) :: rho(:), bz(:), bmean(:), b_scalar, w_scalar
integer, intent(out) :: bend(:), bsize(:), b_count
integer :: n
real(kind=dp) :: z0

n = size(cumy)
rho = 0.0_dp
rho(n) = 1.0_dp
bend = 0
bsize = 0
bz = 0.0_dp
bmean = 0.0_dp

z0 = cumy(n) * cumy(n) / real(n, dp)
bend(1) = n
bsize(1) = n
bz(1) = z0
bmean(1) = cumy(n) / real(n, dp)
b_scalar = z0 - real(n, dp) * ybar * ybar
w_scalar = ysqall - z0
b_count = 1
end subroutine init_step_1d

subroutine init_step_k(cumy, ybar, ysqall, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)
real(kind=dp), intent(in) :: cumy(:, :), ybar, ysqall
real(kind=dp), intent(out) :: rho(:), bz(:), bmean(:, :), b_scalar, w_scalar
integer, intent(out) :: bend(:), bsize(:), b_count
integer :: n, p, j
real(kind=dp) :: z0

n = size(cumy, 1)
p = size(cumy, 2)
rho = 0.0_dp
rho(n) = 1.0_dp
bend = 0
bsize = 0
bz = 0.0_dp
bmean = 0.0_dp

z0 = 0.0_dp
bmean(1, :) = cumy(n, :) / real(n, dp)
do j = 1, p
    z0 = z0 + bmean(1, j) * bmean(1, j) * real(n, dp)
end do
bend(1) = n
bsize(1) = n
bz(1) = z0
b_scalar = z0 - real(n * p, dp) * ybar * ybar
w_scalar = ysqall - z0
b_count = 1
end subroutine init_step_k

subroutine build_prefix_reg_1x(y, x, cumy, cumysq, cumx, cumxx, cumxy)
real(kind=dp), intent(in) :: y(:), x(:)
real(kind=dp), intent(out) :: cumy(:), cumysq(:), cumx(:), cumxx(:), cumxy(:)
integer :: i, n

n = size(y)
cumy(1) = y(1)
cumysq(1) = y(1) * y(1)
cumx(1) = x(1)
cumxx(1) = x(1) * x(1)
cumxy(1) = x(1) * y(1)
do i = 2, n
    cumy(i) = cumy(i - 1) + y(i)
    cumysq(i) = cumysq(i - 1) + y(i) * y(i)
    cumx(i) = cumx(i - 1) + x(i)
    cumxx(i) = cumxx(i - 1) + x(i) * x(i)
    cumxy(i) = cumxy(i - 1) + x(i) * y(i)
end do
end subroutine build_prefix_reg_1x

subroutine build_priors_reg(n, p0, prior_logs)
integer, intent(in) :: n
real(kind=dp), intent(in) :: p0
real(kind=dp), intent(out) :: prior_logs(:)
integer :: j

prior_logs = 0.0_dp
do j = 1, n - 3
    prior_logs(j) = log(reg_incomplete_beta(p0, real(j, dp), real(n - j + 1, dp))) + &
        log_beta(real(j, dp), real(n - j + 1, dp))
end do
end subroutine build_priors_reg

subroutine build_priors(n, p0, priors)
integer, intent(in) :: n
real(kind=dp), intent(in) :: p0
real(kind=dp), intent(out) :: priors(:)
integer :: j

priors = 0.0_dp
do j = 1, n - 3
    priors(j) = exp(log_beta(real(j + 1, dp), real(n - j, dp))) * &
        reg_incomplete_beta(p0, real(j + 1, dp), real(n - j, dp)) / &
        (exp(log_beta(real(j, dp), real(n - j + 1, dp))) * &
         reg_incomplete_beta(p0, real(j, dp), real(n - j + 1, dp)))
end do
end subroutine build_priors

subroutine bcp_pass_1d(cumy, priors, w0, seed, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)
real(kind=dp), intent(in) :: cumy(:), priors(:), w0
integer, intent(inout) :: seed(3), b_count
real(kind=dp), intent(inout) :: rho(:), bz(:), bmean(:), b_scalar, w_scalar
integer, intent(inout) :: bend(:), bsize(:)

integer :: n, i, prevblock, currblock, cp, old_blocks
integer :: thisbend, lastbend, bsize1, bsize2, bsize3
integer :: new_count
real(kind=dp) :: thisblockz, lastblockz, bmean1, bmean2, bmean3, bmeanlast, bz1, bz2, bz3
real(kind=dp) :: tmp, p, u
real(kind=dp) :: bvals(2), wvals(2), bigb(2)
real(kind=dp), allocatable :: rho_old(:), bz_old(:), bmean_old(:)
integer, allocatable :: bend_old(:), bsize_old(:)
real(kind=dp), allocatable :: rho_new(:), bz_new(:), bmean_new(:)
integer, allocatable :: bend_new(:), bsize_new(:)
real(kind=dp) :: new_b, new_w
integer :: new_blocks

n = size(cumy)
allocate(rho_old(n), bz_old(n), bmean_old(n), bend_old(n), bsize_old(n))
allocate(rho_new(n), bz_new(n), bmean_new(n), bend_new(n), bsize_new(n))

rho_old = rho
bz_old = bz
bmean_old = bmean
bend_old = bend
bsize_old = bsize
old_blocks = b_count

rho_new = 0.0_dp
bz_new = 0.0_dp
bmean_new = 0.0_dp
bend_new = 0
bsize_new = 0

prevblock = 1
currblock = 0
thisblockz = bz_old(1)
thisbend = bend_old(1)
lastblockz = 0.0_dp
lastbend = 0
new_b = b_scalar
new_w = w_scalar
new_blocks = b_count

do i = 1, n - 1
    if (i == bend_old(prevblock)) then
        lastblockz = thisblockz
        prevblock = prevblock + 1
        thisbend = bend_old(prevblock)
        thisblockz = bz_old(prevblock)
    end if

    bvals(1) = real(new_blocks, dp)
    if (nint(rho_old(i)) == 0) then
        bigb(1) = new_b
        wvals(1) = new_w
    else
        bvals(1) = real(new_blocks - 1, dp)
        tmp = thisblockz + lastblockz
        if (lastbend > 0) then
            bsize3 = thisbend - lastbend
            bmean3 = (cumy(thisbend) - cumy(lastbend)) / real(bsize3, dp)
        else
            bsize3 = thisbend
            bmean3 = cumy(thisbend) / real(bsize3, dp)
        end if
        bz3 = bmean3 * bmean3 * real(bsize3, dp)
        if (abs(bvals(1) - 1.0_dp) < tiny_dp) then
            bigb(1) = 0.0_dp
        else
            bigb(1) = new_b - tmp + bz3
        end if
        wvals(1) = new_w + tmp - bz3
    end if

    bvals(2) = real(new_blocks, dp)
    if (nint(rho_old(i)) == 1) then
        bigb(2) = new_b
        wvals(2) = new_w
    else
        bvals(2) = real(new_blocks + 1, dp)
        bsize2 = thisbend - i
        if (lastbend > 0) then
            bsize1 = i - lastbend
            bmean1 = (cumy(i) - cumy(lastbend)) / real(bsize1, dp)
        else
            bsize1 = i
            bmean1 = cumy(i) / real(bsize1, dp)
        end if
        bmean2 = (cumy(thisbend) - cumy(i)) / real(bsize2, dp)
        bz1 = bmean1 * bmean1 * real(bsize1, dp)
        bz2 = bmean2 * bmean2 * real(bsize2, dp)
        tmp = thisblockz
        bigb(2) = new_b - tmp + bz1 + bz2
        wvals(2) = new_w + tmp - bz1 - bz2
    end if

    p = getprob_1d(bigb(1), bigb(2), wvals(1), wvals(2), int(nint(bvals(1))), n, w0, priors)
    u = wh_runif(seed)
    if (u < p) then
        cp = 2
    else
        cp = 1
    end if

    new_b = bigb(cp)
    new_w = wvals(cp)
    new_blocks = int(nint(bvals(cp)))

    if (cp /= nint(rho_old(i)) + 1) then
        if (cp == 1) then
            thisblockz = bz3
            if (currblock > 0) then
                lastbend = bend_new(currblock)
                lastblockz = bz_new(currblock)
            else
                lastbend = 0
                lastblockz = 0.0_dp
            end if
        else
            thisblockz = bz2
            lastblockz = bz1
        end if
    end if

    rho_new(i) = real(cp - 1, dp)

    if (nint(rho_new(i)) == 1) then
        if (nint(rho_old(i)) == 1) then
            if (lastbend > 0) then
                bsize1 = i - lastbend
                bmean1 = (cumy(i) - cumy(lastbend)) / real(bsize1, dp)
            else
                bsize1 = i
                bmean1 = cumy(i) / real(bsize1, dp)
            end if
            lastblockz = bmean1 * bmean1 * real(bsize1, dp)
        end if
        currblock = currblock + 1
        bsize_new(currblock) = bsize1
        bend_new(currblock) = i
        bmean_new(currblock) = bmean1
        bz_new(currblock) = lastblockz
        lastbend = i
    end if
end do

new_count = currblock + 1
if (lastbend > 0) then
    bsize_new(new_count) = n - lastbend
    bmeanlast = (cumy(n) - cumy(lastbend)) / real(bsize_new(new_count), dp)
else
    bsize_new(new_count) = n
    bmeanlast = cumy(n) / real(n, dp)
end if
rho_new(n) = 1.0_dp
bend_new(new_count) = n
bmean_new(new_count) = bmeanlast
bz_new(new_count) = bmeanlast * bmeanlast * real(bsize_new(new_count), dp)

rho = rho_new
bz = bz_new
bmean = bmean_new
bend = bend_new
bsize = bsize_new
b_scalar = new_b
w_scalar = new_w
b_count = new_blocks
end subroutine bcp_pass_1d

subroutine bcp_pass_k(cumy, priors, w0, seed, rho, bend, bsize, bz, bmean, b_scalar, w_scalar, b_count)
real(kind=dp), intent(in) :: cumy(:, :), priors(:), w0
integer, intent(inout) :: seed(3), b_count
real(kind=dp), intent(inout) :: rho(:), bz(:), bmean(:, :), b_scalar, w_scalar
integer, intent(inout) :: bend(:), bsize(:)

integer :: n, p, i, prevblock, currblock, cp
integer :: thisbend, lastbend, bsize1, bsize2, bsize3, new_count, new_blocks
real(kind=dp) :: thisblockz, lastblockz, bz1, bz2, bz3, bmeanlast_norm
real(kind=dp) :: tmp, prob, u, new_b, new_w
real(kind=dp) :: bvals(2), wvals(2), bigb(2)
real(kind=dp), allocatable :: rho_old(:), bz_old(:), bmean_old(:, :)
integer, allocatable :: bend_old(:), bsize_old(:)
real(kind=dp), allocatable :: rho_new(:), bz_new(:), bmean_new(:, :)
integer, allocatable :: bend_new(:), bsize_new(:)
real(kind=dp), allocatable :: bmean1(:), bmean2(:), bmean3(:), bmeanlast(:)

n = size(cumy, 1)
p = size(cumy, 2)
allocate(rho_old(n), bz_old(n), bmean_old(n, p), bend_old(n), bsize_old(n))
allocate(rho_new(n), bz_new(n), bmean_new(n, p), bend_new(n), bsize_new(n))
allocate(bmean1(p), bmean2(p), bmean3(p), bmeanlast(p))

rho_old = rho
bz_old = bz
bmean_old = bmean
bend_old = bend
bsize_old = bsize

rho_new = 0.0_dp
bz_new = 0.0_dp
bmean_new = 0.0_dp
bend_new = 0
bsize_new = 0

prevblock = 1
currblock = 0
thisblockz = bz_old(1)
thisbend = bend_old(1)
lastblockz = 0.0_dp
lastbend = 0
new_b = b_scalar
new_w = w_scalar
new_blocks = b_count

do i = 1, n - 1
    if (i == bend_old(prevblock)) then
        lastblockz = thisblockz
        prevblock = prevblock + 1
        thisbend = bend_old(prevblock)
        thisblockz = bz_old(prevblock)
    end if

    bvals(1) = real(new_blocks, dp)
    if (nint(rho_old(i)) == 0) then
        bigb(1) = new_b
        wvals(1) = new_w
    else
        bvals(1) = real(new_blocks - 1, dp)
        tmp = thisblockz + lastblockz
        if (lastbend > 0) then
            bsize3 = thisbend - lastbend
            bmean3 = (cumy(thisbend, :) - cumy(lastbend, :)) / real(bsize3, dp)
        else
            bsize3 = thisbend
            bmean3 = cumy(thisbend, :) / real(bsize3, dp)
        end if
        bz3 = sum(bmean3 * bmean3) * real(bsize3, dp)
        if (abs(bvals(1) - 1.0_dp) < tiny_dp .and. p == 1) then
            bigb(1) = 0.0_dp
        else
            bigb(1) = new_b - tmp + bz3
        end if
        wvals(1) = new_w + tmp - bz3
    end if

    bvals(2) = real(new_blocks, dp)
    if (nint(rho_old(i)) == 1) then
        bigb(2) = new_b
        wvals(2) = new_w
    else
        bvals(2) = real(new_blocks + 1, dp)
        bsize2 = thisbend - i
        if (lastbend > 0) then
            bsize1 = i - lastbend
            bmean1 = (cumy(i, :) - cumy(lastbend, :)) / real(bsize1, dp)
        else
            bsize1 = i
            bmean1 = cumy(i, :) / real(bsize1, dp)
        end if
        bmean2 = (cumy(thisbend, :) - cumy(i, :)) / real(bsize2, dp)
        bz1 = sum(bmean1 * bmean1) * real(bsize1, dp)
        bz2 = sum(bmean2 * bmean2) * real(bsize2, dp)
        tmp = thisblockz
        bigb(2) = new_b - tmp + bz1 + bz2
        wvals(2) = new_w + tmp - bz1 - bz2
    end if

    prob = getprob_k(bigb(1), bigb(2), wvals(1), wvals(2), int(nint(bvals(1))), n, p, w0, priors)
    u = wh_runif(seed)
    if (u < prob) then
        cp = 2
    else
        cp = 1
    end if

    new_b = bigb(cp)
    new_w = wvals(cp)
    new_blocks = int(nint(bvals(cp)))

    if (cp /= nint(rho_old(i)) + 1) then
        if (cp == 1) then
            thisblockz = bz3
            if (currblock > 0) then
                lastbend = bend_new(currblock)
                lastblockz = bz_new(currblock)
            else
                lastbend = 0
                lastblockz = 0.0_dp
            end if
        else
            thisblockz = bz2
            lastblockz = bz1
        end if
    end if

    rho_new(i) = real(cp - 1, dp)

    if (nint(rho_new(i)) == 1) then
        if (nint(rho_old(i)) == 1) then
            if (lastbend > 0) then
                bsize1 = i - lastbend
                bmean1 = (cumy(i, :) - cumy(lastbend, :)) / real(bsize1, dp)
            else
                bsize1 = i
                bmean1 = cumy(i, :) / real(bsize1, dp)
            end if
            lastblockz = sum(bmean1 * bmean1) * real(bsize1, dp)
        end if
        currblock = currblock + 1
        bsize_new(currblock) = bsize1
        bend_new(currblock) = i
        bmean_new(currblock, :) = bmean1
        bz_new(currblock) = lastblockz
        lastbend = i
    end if
end do

new_count = currblock + 1
if (lastbend > 0) then
    bsize_new(new_count) = n - lastbend
    bmeanlast = (cumy(n, :) - cumy(lastbend, :)) / real(bsize_new(new_count), dp)
else
    bsize_new(new_count) = n
    bmeanlast = cumy(n, :) / real(n, dp)
end if
bmeanlast_norm = sum(bmeanlast * bmeanlast) * real(bsize_new(new_count), dp)
rho_new(n) = 1.0_dp
bend_new(new_count) = n
bmean_new(new_count, :) = bmeanlast
bz_new(new_count) = bmeanlast_norm

rho = rho_new
bz = bz_new
bmean = bmean_new
bend = bend_new
bsize = bsize_new
b_scalar = new_b
w_scalar = new_w
b_count = new_blocks
end subroutine bcp_pass_k

subroutine bcp_pass_reg_1x(cumy, cumx, cumxx, cumxy, prior_logs, w0_vec, d, nreg, seed, &
    rho, bend, bsize, btau, bz, blogc, bk, bq, w, bigb, wsum, logc, qsum, ksum, lik, b_count)
real(kind=dp), intent(in) :: cumy(:), cumx(:), cumxx(:), cumxy(:), prior_logs(:), w0_vec(2), d
integer, intent(in) :: nreg
integer, intent(inout) :: seed(3), b_count
real(kind=dp), intent(inout) :: rho(:), bz(:), blogc(:), bk(:), bq(:), w(:), bigb, wsum, logc, qsum, ksum, lik
integer, intent(inout) :: bend(:), bsize(:), btau(:)

integer :: n, i, s, cp, tau, prevblock, currblock, thisbend, lastbend
integer :: bsize1, bsize2, bsize3, thisbtau, lastbtau, new_blocks
real(kind=dp) :: thisblockz, lastblockz, thisblogc, lastblogc, thisbk, lastbk, thisbq, lastbq
real(kind=dp) :: bmean1, bmean2, bmean3, tmp, qval0, logcval0, kval0, maxlik
real(kind=dp) :: bvals(2), wvals(2), bbig(2), likvals(6), qvals(6), logcvals(6), kvals(6)
real(kind=dp) :: bK1, bK2, bK3, bK4, matQ, matLogC, matQ1, matLogC1, matQ2, matLogC2
real(kind=dp), allocatable :: rho_old(:), bz_old(:), blogc_old(:), bk_old(:), bq_old(:), w_old(:)
real(kind=dp), allocatable :: rho_new(:), bz_new(:), blogc_new(:), bk_new(:), bq_new(:)
integer, allocatable :: bend_old(:), bsize_old(:), btau_old(:), bend_new(:), bsize_new(:), btau_new(:)

n = size(cumy)
allocate(rho_old(n), bz_old(n), blogc_old(n), bk_old(n), bq_old(n), w_old(2))
allocate(rho_new(n), bz_new(n), blogc_new(n), bk_new(n), bq_new(n))
allocate(bend_old(n), bsize_old(n), btau_old(n), bend_new(n), bsize_new(n), btau_new(n))

rho_old = rho
bz_old = bz
blogc_old = blogc
bk_old = bk
bq_old = bq
w_old = w
bend_old = bend
bsize_old = bsize
btau_old = btau

rho_new = 0.0_dp
bz_new = 0.0_dp
blogc_new = 0.0_dp
bk_new = 0.0_dp
bq_new = 0.0_dp
bend_new = 0
bsize_new = 0
btau_new = 0

prevblock = 1
currblock = 0
thisblockz = bz_old(1)
thisbend = bend_old(1)
thisbtau = btau_old(1)
thisblogc = blogc_old(1)
thisbk = bk_old(1)
thisbq = bq_old(1)
lastblockz = 0.0_dp
lastbend = 0
lastbtau = 0
lastblogc = 0.0_dp
lastbk = 0.0_dp
lastbq = 0.0_dp
new_blocks = b_count

do i = 1, n - 1
    maxlik = -huge(1.0_dp)
    if (i == bend_old(prevblock)) then
        lastblockz = thisblockz
        lastbtau = thisbtau
        lastblogc = thisblogc
        lastbk = thisbk
        lastbq = thisbq
        prevblock = prevblock + 1
        thisbend = bend_old(prevblock)
        thisblockz = bz_old(prevblock)
        thisbtau = btau_old(prevblock)
        thisblogc = blogc_old(prevblock)
        thisbk = bk_old(prevblock)
        thisbq = bq_old(prevblock)
    end if

    if (nint(rho_old(i)) == 1) then
        qval0 = qsum - thisbq - lastbq
        logcval0 = logc - thisblogc - lastblogc
        kval0 = ksum - thisbk - lastbk
    else
        qval0 = qsum - thisbq
        logcval0 = logc - thisblogc
        kval0 = ksum - thisbk
    end if

    bvals(1) = real(new_blocks - merge(1, 0, nint(rho_old(i)) == 1), dp)
    if (lastbend > 0) then
        bsize3 = thisbend - lastbend
        bmean3 = (cumy(thisbend) - cumy(lastbend)) / real(bsize3, dp)
    else
        bsize3 = thisbend
        bmean3 = cumy(thisbend) / real(bsize3, dp)
    end if
    tmp = 0.0_dp
    if (nint(rho_old(i)) == 1) tmp = thisblockz + lastblockz - bmean3 * bmean3 * real(bsize3, dp)
    if (int(nint(bvals(1))) == 1) then
        bbig(1) = 0.0_dp
    else
        bbig(1) = bigb - tmp
    end if
    wvals(1) = wsum + tmp

    do tau = 0, 1
        if (tau == 1 .and. bsize3 < nreg) then
            likvals(tau + 1) = -huge(1.0_dp)
            exit
        end if
        qvals(tau + 1) = qval0
        logcvals(tau + 1) = logcval0
        if (tau == 1) then
            call matrixcalcs_reg_1x(cumy, cumx, cumxx, cumxy, lastbend + 1, thisbend, w(2), matQ, matLogC)
            qvals(tau + 1) = qvals(tau + 1) + matQ
            logcvals(tau + 1) = logcvals(tau + 1) + matLogC
        end if
        kvals(tau + 1) = kval0 + logkcalc_reg(bsize3, tau, d, nreg)
        likvals(tau + 1) = likelihood_reg(bbig(1), wvals(1), int(nint(bvals(1))), n, w(1), logcvals(tau + 1), qvals(tau + 1), kvals(tau + 1), prior_logs)
        if (likvals(tau + 1) > maxlik) maxlik = likvals(tau + 1)
    end do

    bvals(2) = real(new_blocks + merge(0, 1, nint(rho_old(i)) == 1), dp)
    bsize2 = thisbend - i
    if (lastbend > 0) then
        bsize1 = i - lastbend
        bmean1 = (cumy(i) - cumy(lastbend)) / real(bsize1, dp)
    else
        bsize1 = i
        bmean1 = cumy(i) / real(bsize1, dp)
    end if
    bmean2 = (cumy(thisbend) - cumy(i)) / real(bsize2, dp)
    tmp = 0.0_dp
    if (nint(rho_old(i)) == 0) tmp = thisblockz - bmean1 * bmean1 * real(bsize1, dp) - bmean2 * bmean2 * real(bsize2, dp)
    bbig(2) = bigb - tmp
    wvals(2) = wsum + tmp

    bK1 = logkcalc_reg(bsize1, 0, d, nreg)
    if (bsize1 >= nreg) then
        bK2 = logkcalc_reg(bsize1, 1, d, nreg)
        call matrixcalcs_reg_1x(cumy, cumx, cumxx, cumxy, lastbend + 1, i, w(2), matQ1, matLogC1)
    else
        bK2 = 0.0_dp
        matQ1 = 0.0_dp
        matLogC1 = 0.0_dp
    end if
    bK3 = logkcalc_reg(bsize2, 0, d, nreg)
    if (bsize2 >= nreg) then
        bK4 = logkcalc_reg(bsize2, 1, d, nreg)
        call matrixcalcs_reg_1x(cumy, cumx, cumxx, cumxy, i + 1, thisbend, w(2), matQ2, matLogC2)
    else
        bK4 = 0.0_dp
        matQ2 = 0.0_dp
        matLogC2 = 0.0_dp
    end if

    do tau = 2, 5
        if ((bsize1 < nreg .and. tau > 3) .or. (bsize2 < nreg .and. (tau == 3 .or. tau == 5))) then
            likvals(tau + 1) = -huge(1.0_dp)
            cycle
        end if
        kvals(tau + 1) = kval0 + merge(bK1, bK2, tau < 4) + merge(bK3, bK4, tau == 2 .or. tau == 4)
        qvals(tau + 1) = qval0
        logcvals(tau + 1) = logcval0
        if (tau >= 4) then
            qvals(tau + 1) = qvals(tau + 1) + matQ1
            logcvals(tau + 1) = logcvals(tau + 1) + matLogC1
        end if
        if (tau == 3 .or. tau == 5) then
            qvals(tau + 1) = qvals(tau + 1) + matQ2
            logcvals(tau + 1) = logcvals(tau + 1) + matLogC2
        end if
        likvals(tau + 1) = likelihood_reg(bbig(2), wvals(2), int(nint(bvals(2))), n, w(1), logcvals(tau + 1), qvals(tau + 1), kvals(tau + 1), prior_logs)
        if (likvals(tau + 1) > maxlik) maxlik = likvals(tau + 1)
    end do

    s = sample_from_likelihoods(likvals, maxlik, seed)
    if (s <= 2) then
        cp = 0
    else
        cp = 1
    end if

    if (cp /= nint(rho_old(i))) then
        if (cp == 0) then
            thisblockz = bmean3 * bmean3 * real(bsize3, dp)
            thisbk = thisbk + kvals(s) - ksum + lastbk
            if (currblock > 0) then
                lastbend = bend_new(currblock)
                lastblockz = bz_new(currblock)
                lastblogc = blogc_new(currblock)
                lastbq = bq_new(currblock)
                lastbk = bk_new(currblock)
            else
                lastblockz = 0.0_dp
                lastbend = 0
                lastblogc = 0.0_dp
                lastbq = 0.0_dp
                lastbk = 0.0_dp
            end if
        else
            thisblockz = bmean2 * bmean2 * real(bsize2, dp)
            lastblockz = bmean1 * bmean1 * real(bsize1, dp)
        end if
    end if

    rho_new(i) = real(cp, dp)
    if (nint(rho_new(i)) == 0) then
        if (nint(rho_old(i)) == 0) thisbk = thisbk + kvals(s) - ksum
        if (s == 2) then
            thisbq = matQ
            thisblogc = matLogC
        else
            thisbq = 0.0_dp
            thisblogc = 0.0_dp
        end if
        thisbtau = s - 1
    else
        if (s < 5) then
            lastbtau = 0
            lastbq = 0.0_dp
            lastblogc = 0.0_dp
        else
            lastbtau = 1
            lastbq = matQ1
            lastblogc = matLogC1
        end if
        if (s == 3 .or. s == 5) then
            thisbtau = 0
            thisbq = 0.0_dp
            thisblogc = 0.0_dp
        else
            thisbtau = 1
            thisbq = matQ2
            thisblogc = matLogC2
        end if
        lastbk = merge(bK1, bK2, s == 3 .or. s == 4)
        thisbk = merge(bK3, bK4, s == 3 .or. s == 5)
        currblock = currblock + 1
        bsize_new(currblock) = bsize1
        bend_new(currblock) = i
        bz_new(currblock) = lastblockz
        btau_new(currblock) = lastbtau
        blogc_new(currblock) = lastblogc
        bq_new(currblock) = lastbq
        bk_new(currblock) = lastbk
        lastbend = i
    end if

    lik = likvals(s)
    bigb = bbig(cp + 1)
    wsum = wvals(cp + 1)
    new_blocks = int(nint(bvals(cp + 1)))
    logc = logcvals(s)
    qsum = qvals(s)
    ksum = kvals(s)
end do

if (lastbend == 0) then
    bsize_new(new_blocks) = n
else
    bsize_new(new_blocks) = n - lastbend
end if
 bq_new(new_blocks) = thisbq
 bend_new(new_blocks) = n
 btau_new(new_blocks) = thisbtau
 blogc_new(new_blocks) = thisblogc
 bz_new(new_blocks) = thisblockz
 bk_new(new_blocks) = thisbk

call sample_w_reg_1x(cumy, cumx, cumxx, cumxy, prior_logs, w0_vec, seed, bend_new, btau_new, &
    blogc_new, bq_new, new_blocks, bigb, wsum, ksum, w, logc, qsum, lik)

rho_new(n) = 1.0_dp
rho = rho_new
bz = bz_new
blogc = blogc_new
bk = bk_new
bq = bq_new
bend = bend_new
bsize = bsize_new
btau = btau_new
b_count = new_blocks
end subroutine bcp_pass_reg_1x

subroutine sample_w_reg_1x(cumy, cumx, cumxx, cumxy, prior_logs, w0_vec, seed, bend, btau, blogc, bq, &
    b_count, bigb, wsum, ksum, w, logc, qsum, lik)
real(kind=dp), intent(in) :: cumy(:), cumx(:), cumxx(:), cumxy(:), prior_logs(:), w0_vec(2)
integer, intent(in) :: b_count, bend(:), btau(:)
integer, intent(inout) :: seed(3)
real(kind=dp), intent(inout) :: blogc(:), bq(:), bigb, wsum, ksum, w(:), logc, qsum, lik
integer :: s, lastbend, thisbend
real(kind=dp) :: w2cand, matQ, matLogC, q0cand, logc0cand, likcand, pacc, u
real(kind=dp), allocatable :: logcvec(:), qvec(:)

w2cand = w(2) + (2.0_dp * wh_runif(seed) - 1.0_dp) * 0.05_dp * w0_vec(2)
if (w2cand < 0.0_dp .or. w2cand > w0_vec(2)) return

allocate(logcvec(b_count), qvec(b_count))
q0cand = 0.0_dp
logc0cand = 0.0_dp
lastbend = 0
do s = 1, b_count
    thisbend = bend(s)
    if (btau(s) == 1) then
        call matrixcalcs_reg_1x(cumy, cumx, cumxx, cumxy, lastbend + 1, thisbend, w2cand, matQ, matLogC)
        qvec(s) = matQ
        logcvec(s) = matLogC
        q0cand = q0cand + qvec(s)
        logc0cand = logc0cand + logcvec(s)
    else
        qvec(s) = 0.0_dp
        logcvec(s) = 0.0_dp
    end if
    lastbend = thisbend
end do
likcand = likelihood_reg(bigb, wsum, b_count, size(cumy), w(1), logc0cand, q0cand, ksum, prior_logs)
pacc = exp(likcand - lik)
pacc = pacc / (1.0_dp + pacc)
u = wh_runif(seed)
if (u < pacc) then
    w(2) = w2cand
    blogc(:b_count) = logcvec
    bq(:b_count) = qvec
    logc = logc0cand
    qsum = q0cand
    lik = likcand
end if
end subroutine sample_w_reg_1x

subroutine accumulate_posterior_reg_1x(x, cumy, cumx, cumxx, cumxy, ybar, w, wstar, rho, bend, bsize, btau, posterior_mean, posterior_prob)
real(kind=dp), intent(in) :: x(:), cumy(:), cumx(:), cumxx(:), cumxy(:), ybar, w(:), wstar, rho(:)
integer, intent(in) :: bend(:), bsize(:), btau(:)
real(kind=dp), intent(inout) :: posterior_mean(:), posterior_prob(:)
integer :: start_idx, end_idx, i
real(kind=dp) :: bmean, tmpalpha, bhat, xbar

start_idx = 1
do i = 1, count(nint(rho) == 1)
    end_idx = bend(i)
    if (start_idx > 1) then
        bmean = (cumy(end_idx) - cumy(start_idx - 1)) / real(bsize(i), dp)
    else
        bmean = cumy(end_idx) / real(bsize(i), dp)
    end if
    tmpalpha = bmean * (1.0_dp - wstar) + ybar * wstar
    posterior_mean(start_idx:end_idx) = posterior_mean(start_idx:end_idx) + tmpalpha
    if (btau(i) == 1) then
        bhat = compute_bhat_reg_1x(cumy, cumx, cumxx, cumxy, start_idx, end_idx, w(2))
        xbar = block_mean_x(cumx, start_idx, end_idx)
        posterior_mean(start_idx:end_idx) = posterior_mean(start_idx:end_idx) + (x(start_idx:end_idx) - xbar) * bhat
    end if
    start_idx = end_idx + 1
end do
posterior_prob = posterior_prob + rho
end subroutine accumulate_posterior_reg_1x

real(kind=dp) function logkcalc_reg(bsize, tau, d, nreg) result(val)
integer, intent(in) :: bsize, tau, nreg
real(kind=dp), intent(in) :: d
real(kind=dp) :: kratio, tmp
integer :: canfit

kratio = d / (real(bsize, dp) + d)
canfit = merge(1, 0, bsize >= nreg)
if (tau == 0) then
    tmp = kratio * real(canfit, dp) + real(1 - canfit, dp)
else
    tmp = (1.0_dp - kratio) * real(canfit, dp)
end if
val = log(max(tmp, tiny_dp))
end function logkcalc_reg

subroutine matrixcalcs_reg_1x(cumy, cumx, cumxx, cumxy, start_idx, end_idx, w1, qval, logcval)
real(kind=dp), intent(in) :: cumy(:), cumx(:), cumxx(:), cumxy(:), w1
integer, intent(in) :: start_idx, end_idx
real(kind=dp), intent(out) :: qval, logcval
integer :: n
real(kind=dp) :: sumx, sumxx, sumy, sumxy, sxx, sxy

n = end_idx - start_idx + 1
sumx = block_sum_1d(cumx, start_idx, end_idx)
sumxx = block_sum_1d(cumxx, start_idx, end_idx)
sumy = block_sum_1d(cumy, start_idx, end_idx)
sumxy = block_sum_1d(cumxy, start_idx, end_idx)
sxx = sumxx - sumx * sumx / real(n, dp)
sxy = sumxy - sumx * sumy / real(n, dp)
if (sxx < tiny_dp) then
    qval = 0.0_dp
else
    qval = (1.0_dp - w1) * sxy * sxy / sxx
end if
logcval = 0.5_dp * log(max(w1, tiny_dp))
end subroutine matrixcalcs_reg_1x

real(kind=dp) function compute_bhat_reg_1x(cumy, cumx, cumxx, cumxy, start_idx, end_idx, w1) result(bhat)
real(kind=dp), intent(in) :: cumy(:), cumx(:), cumxx(:), cumxy(:), w1
integer, intent(in) :: start_idx, end_idx
integer :: n
real(kind=dp) :: sumx, sumxx, sumy, sumxy, sxx, sxy

n = end_idx - start_idx + 1
sumx = block_sum_1d(cumx, start_idx, end_idx)
sumxx = block_sum_1d(cumxx, start_idx, end_idx)
sumy = block_sum_1d(cumy, start_idx, end_idx)
sumxy = block_sum_1d(cumxy, start_idx, end_idx)
sxx = sumxx - sumx * sumx / real(n, dp)
sxy = sumxy - sumx * sumy / real(n, dp)
if (sxx < tiny_dp) then
    bhat = 0.0_dp
else
    bhat = (1.0_dp - w1) * sxy / sxx
end if
end function compute_bhat_reg_1x

real(kind=dp) function block_mean_x(cumx, start_idx, end_idx) result(xbar)
real(kind=dp), intent(in) :: cumx(:)
integer, intent(in) :: start_idx, end_idx
xbar = block_sum_1d(cumx, start_idx, end_idx) / real(end_idx - start_idx + 1, dp)
end function block_mean_x

real(kind=dp) function block_sum_1d(csum, start_idx, end_idx) result(s)
real(kind=dp), intent(in) :: csum(:)
integer, intent(in) :: start_idx, end_idx
if (start_idx > 1) then
    s = csum(end_idx) - csum(start_idx - 1)
else
    s = csum(end_idx)
end if
end function block_sum_1d

integer function sample_from_likelihoods(likvals, maxlik, seed) result(idx)
real(kind=dp), intent(in) :: likvals(:), maxlik
integer, intent(inout) :: seed(3)
real(kind=dp) :: weights(size(likvals)), csum, u
integer :: i

csum = 0.0_dp
do i = 1, size(likvals)
    if (likvals(i) <= -huge(1.0_dp)/2) then
        weights(i) = csum
    else
        csum = csum + exp(likvals(i) - maxlik)
        weights(i) = csum
    end if
end do
u = wh_runif(seed) * csum
do i = 1, size(weights)
    if (u < weights(i)) then
        idx = i
        return
    end if
end do
idx = size(weights)
end function sample_from_likelihoods

real(kind=dp) function likelihood_reg(bigb, wsum, b, n, w0, logc, qsum, ksum, prior_logs) result(lik)
real(kind=dp), intent(in) :: bigb, wsum, w0, logc, qsum, ksum, prior_logs(:)
integer, intent(in) :: b, n
real(kind=dp) :: wtilde, xmax

wtilde = max(wsum - qsum, tiny_dp)
if (b == 1) then
    lik = logc + log(max(w0, tiny_dp)) - 0.5_dp * real(n - 1, dp) * log(wtilde)
else
    if (bigb <= tiny_dp) then
        lik = -huge(1.0_dp)
        return
    end if
    xmax = (bigb * w0 / wtilde) / (1.0_dp + bigb * w0 / wtilde)
    xmax = min(1.0_dp - tiny_dp, max(tiny_dp, xmax))
    lik = logc - 0.5_dp * real(b + 1, dp) * log(bigb) - 0.5_dp * real(n - b - 2, dp) * log(wtilde) + &
        log_beta(0.5_dp * real(b + 1, dp), 0.5_dp * real(n - b - 2, dp)) + &
        log(reg_incomplete_beta(xmax, 0.5_dp * real(b + 1, dp), 0.5_dp * real(n - b - 2, dp)))
end if
lik = lik + ksum + prior_logs(b)
end function likelihood_reg

real(kind=dp) function compute_wstar_reg(bigb, wtilde, b, n, w0) result(wstar)
real(kind=dp), intent(in) :: bigb, wtilde, w0
integer, intent(in) :: b, n
real(kind=dp) :: xmax

if (b > 1) then
    if (bigb <= tiny_dp) then
        wstar = 0.5_dp * w0
        return
    end if
    xmax = (bigb * w0 / wtilde) / (1.0_dp + bigb * w0 / wtilde)
    xmax = min(1.0_dp - tiny_dp, max(tiny_dp, xmax))
    wstar = exp(log(wtilde) - log(bigb) + &
        log_beta(0.5_dp * real(b + 3, dp), 0.5_dp * real(n - b - 4, dp)) + &
        log(reg_incomplete_beta(xmax, 0.5_dp * real(b + 3, dp), 0.5_dp * real(n - b - 4, dp))) - &
        log_beta(0.5_dp * real(b + 1, dp), 0.5_dp * real(n - b - 2, dp)) - &
        log(reg_incomplete_beta(xmax, 0.5_dp * real(b + 1, dp), 0.5_dp * real(n - b - 2, dp))))
else
    wstar = 0.5_dp * w0
end if
end function compute_wstar_reg

real(kind=dp) function getprob_1d(b0, b1, w0sum, w1sum, b, n, w0, priors) result(p)
real(kind=dp), intent(in) :: b0, b1, w0sum, w1sum, w0
integer, intent(in) :: b, n
real(kind=dp), intent(in) :: priors(:)
real(kind=dp) :: xmax1, xmax0, ratio, log_ratio, a1, c1, a0, c0

if (b >= n - 4) then
    p = 0.0_dp
    return
end if

a1 = 0.5_dp * real(b + 2, dp)
c1 = 0.5_dp * real(n - b - 3, dp)
xmax1 = safe_xmax(b1, w1sum, w0)
ratio = priors(b) * exp(log_beta(a1, c1)) * reg_incomplete_beta(xmax1, a1, c1)

if (abs(b0) < tiny_dp) then
    log_ratio = 0.5_dp * real(n - 1, dp) * log(max(w0sum, tiny_dp)) + &
        log(0.5_dp * real(b + 1, dp)) - &
        0.5_dp * (real(b + 1, dp) * log(max(w0, tiny_dp)) + &
                  real(b + 2, dp) * log(max(b1, tiny_dp)) + &
                  real(n - b - 3, dp) * log(max(w1sum, tiny_dp)))
    ratio = ratio * exp(log_ratio)
    p = ratio / (1.0_dp + ratio)
    return
end if

xmax0 = safe_xmax(b0, w0sum, w0)
a0 = 0.5_dp * real(b + 1, dp)
c0 = 0.5_dp * real(n - b - 2, dp)
log_ratio = 0.5_dp * (real(n - b - 2, dp) * log(max(w0sum / max(w1sum, tiny_dp), tiny_dp)) + &
                      real(b + 1, dp) * log(max(b0 / max(b1, tiny_dp), tiny_dp)) + &
                      log(max(w1sum / max(b1, tiny_dp), tiny_dp)))
ratio = ratio * exp(log_ratio) / &
    (exp(log_beta(a0, c0)) * reg_incomplete_beta(xmax0, a0, c0))
p = ratio / (1.0_dp + ratio)
if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) p = min(1.0_dp, max(0.0_dp, p))
end function getprob_1d

real(kind=dp) function getprob_k(b0, b1, w0sum, w1sum, b, n, p, w0, priors) result(prob)
real(kind=dp), intent(in) :: b0, b1, w0sum, w1sum, w0
integer, intent(in) :: b, n, p
real(kind=dp), intent(in) :: priors(:)
real(kind=dp) :: xmax1, xmax0, ratio, log_ratio, a1, c1, a0, c0

if (b >= n - 4 / p) then
    prob = 0.0_dp
    return
end if

a1 = 0.5_dp * real(p * (b + 1) + 1, dp)
c1 = 0.5_dp * real((n - b - 1) * p - 2, dp)
xmax1 = safe_xmax(b1, w1sum, w0)
ratio = priors(b) * exp(log_beta(a1, c1)) * reg_incomplete_beta(xmax1, a1, c1)

if (abs(b0) < tiny_dp) then
    log_ratio = 0.5_dp * real(n * p - 1, dp) * log(max(w0sum, tiny_dp)) + &
        log(0.5_dp * real(p * b + 1, dp)) - &
        0.5_dp * (real(p * b + 1, dp) * log(max(w0, tiny_dp)) + &
                  real(p * (b + 1) + 1, dp) * log(max(b1, tiny_dp)) + &
                  real((n - b - 1) * p - 2, dp) * log(max(w1sum, tiny_dp)))
        ratio = ratio * exp(log_ratio)
        prob = ratio / (1.0_dp + ratio)
        return
end if

xmax0 = safe_xmax(b0, w0sum, w0)
a0 = 0.5_dp * real(p * b + 1, dp)
c0 = 0.5_dp * real((n - b) * p - 2, dp)
log_ratio = 0.5_dp * (real((n - b) * p - 2, dp) * log(max(w0sum / max(w1sum, tiny_dp), tiny_dp)) + &
                      real(p * b + 1, dp) * log(max(b0 / max(b1, tiny_dp), tiny_dp)) + &
                      real(p, dp) * log(max(w1sum / max(b1, tiny_dp), tiny_dp)))
ratio = ratio * exp(log_ratio) / (exp(log_beta(a0, c0)) * reg_incomplete_beta(xmax0, a0, c0))
prob = ratio / (1.0_dp + ratio)
if (.not. (prob >= 0.0_dp .and. prob <= 1.0_dp)) prob = min(1.0_dp, max(0.0_dp, prob))
end function getprob_k

real(kind=dp) function compute_wstar_1d(bigb, wsum, b, n, w0) result(wstar)
real(kind=dp), intent(in) :: bigb, wsum, w0
integer, intent(in) :: b, n
real(kind=dp) :: xmax, num_a, num_c, den_a, den_c

if (abs(bigb) < tiny_dp) then
    wstar = w0 * real(b + 1, dp) / real(b + 3, dp)
    return
end if

xmax = safe_xmax(bigb, wsum, w0)
num_a = 0.5_dp * real(b + 3, dp)
num_c = 0.5_dp * real(n - b - 4, dp)
den_a = 0.5_dp * real(b + 1, dp)
den_c = 0.5_dp * real(n - b - 2, dp)

wstar = (wsum / bigb) * exp(log_beta(num_a, num_c) - log_beta(den_a, den_c)) * &
    reg_incomplete_beta(xmax, num_a, num_c) / reg_incomplete_beta(xmax, den_a, den_c)
end function compute_wstar_1d

real(kind=dp) function compute_wstar_k(bigb, wsum, b, n, p, w0) result(wstar)
real(kind=dp), intent(in) :: bigb, wsum, w0
integer, intent(in) :: b, n, p
real(kind=dp) :: xmax, num_a, num_c, den_a, den_c

if (abs(bigb) < tiny_dp) then
    wstar = w0 * real(p * b + 1, dp) / real(p * b + 3, dp)
    return
end if

xmax = safe_xmax(bigb, wsum, w0)
num_a = 0.5_dp * real(p * b + 3, dp)
num_c = 0.5_dp * real((n - b) * p - 4, dp)
den_a = 0.5_dp * real(p * b + 1, dp)
den_c = 0.5_dp * real((n - b) * p - 2, dp)

wstar = (wsum / bigb) * exp(log_beta(num_a, num_c) - log_beta(den_a, den_c)) * &
    reg_incomplete_beta(xmax, num_a, num_c) / reg_incomplete_beta(xmax, den_a, den_c)
end function compute_wstar_k

pure real(kind=dp) function safe_xmax(bigb, wsum, w0) result(x)
real(kind=dp), intent(in) :: bigb, wsum, w0
real(kind=dp) :: denom

denom = 1.0_dp + (bigb * w0) / max(wsum, tiny_dp)
x = (bigb * w0 / max(wsum, tiny_dp)) / max(denom, tiny_dp)
x = min(1.0_dp - tiny_dp, max(tiny_dp, x))
end function safe_xmax

pure function cumulative_sum(x) result(cs)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: cs(size(x))
integer :: i

cs(1) = x(1)
do i = 2, size(x)
    cs(i) = cs(i - 1) + x(i)
end do
end function cumulative_sum

pure function cumulative_sum_mat(x) result(cs)
real(kind=dp), intent(in) :: x(:, :)
real(kind=dp) :: cs(size(x, 1), size(x, 2))
integer :: i

cs(1, :) = x(1, :)
do i = 2, size(x, 1)
    cs(i, :) = cs(i - 1, :) + x(i, :)
end do
end function cumulative_sum_mat

real(kind=dp) function wh_runif(seed) result(u)
integer, intent(inout) :: seed(3)

seed(1) = mod(171 * seed(1), wh_m1)
seed(2) = mod(172 * seed(2), wh_m2)
seed(3) = mod(170 * seed(3), wh_m3)
u = modulo(real(seed(1), dp) / real(wh_m1, dp) + &
           real(seed(2), dp) / real(wh_m2, dp) + &
           real(seed(3), dp) / real(wh_m3, dp), 1.0_dp)
if (u <= 0.0_dp) u = tiny_dp
if (u >= 1.0_dp) u = 1.0_dp - tiny_dp
end function wh_runif

pure real(kind=dp) function log_beta(a, b) result(val)
real(kind=dp), intent(in) :: a, b
val = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
end function log_beta

real(kind=dp) function reg_incomplete_beta(x, a, b) result(val)
real(kind=dp), intent(in) :: x, a, b
real(kind=dp) :: bt

if (x <= 0.0_dp) then
    val = 0.0_dp
    return
end if
if (x >= 1.0_dp) then
    val = 1.0_dp
    return
end if

bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1.0_dp - x))
if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
    val = bt * betacf(a, b, x) / a
else
    val = 1.0_dp - bt * betacf(b, a, 1.0_dp - x) / b
end if
val = min(1.0_dp, max(0.0_dp, val))
end function reg_incomplete_beta

real(kind=dp) function betacf(a, b, x) result(cf)
real(kind=dp), intent(in) :: a, b, x
integer, parameter :: max_iter = 200
real(kind=dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
integer :: m, m2
real(kind=dp) :: aa, c, d, del, qab, qam, qap

qab = a + b
qap = a + 1.0_dp
qam = a - 1.0_dp
c = 1.0_dp
d = 1.0_dp - qab * x / qap
if (abs(d) < fpmin) d = fpmin
d = 1.0_dp / d
cf = d

do m = 1, max_iter
    m2 = 2 * m
    aa = real(m, dp) * (b - real(m, dp)) * x / ((qam + real(m2, dp)) * (a + real(m2, dp)))
    d = 1.0_dp + aa * d
    if (abs(d) < fpmin) d = fpmin
    c = 1.0_dp + aa / c
    if (abs(c) < fpmin) c = fpmin
    d = 1.0_dp / d
    cf = cf * d * c

    aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / ((a + real(m2, dp)) * (qap + real(m2, dp)))
    d = 1.0_dp + aa * d
    if (abs(d) < fpmin) d = fpmin
    c = 1.0_dp + aa / c
    if (abs(c) < fpmin) c = fpmin
    d = 1.0_dp / d
    del = d * c
    cf = cf * del
    if (abs(del - 1.0_dp) <= eps) exit
end do
end function betacf

end module bcp_pkg_mod
