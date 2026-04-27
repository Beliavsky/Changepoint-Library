program xsim_dp_var1
! simulate multivariate VAR(1) changepoints and estimate them by DP.VAR1 with lambda = 0
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_var1_mod, only: solve_dp_var1_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: p = 2 ! number of variables
integer, parameter :: n_trans = 120 ! number of transitions
real(kind=dp), parameter :: gamma = 0.2_dp ! l0 penalty parameter
real(kind=dp), parameter :: lambda = 0.0_dp ! lasso penalty parameter
integer, parameter :: delta = 5 ! minimum spacing

real(kind=dp) :: data(p, n_trans + 1), x_curr(p, n_trans), x_futu(p, n_trans)
real(kind=dp), parameter :: a1(p, p) = reshape([0.6_dp, 0.0_dp, 0.0_dp, 0.4_dp], [p, p])
real(kind=dp), parameter :: a2(p, p) = reshape([-0.4_dp, 0.2_dp, 0.3_dp, 0.5_dp], [p, p])
real(kind=dp), parameter :: a3(p, p) = reshape([0.7_dp, 0.1_dp, -0.2_dp, -0.5_dp], [p, p])
integer, allocatable :: partition(:), cpt_hat(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_data(data)
x_curr = data(:, 1:n_trans)
x_futu = data(:, 2:n_trans + 1)
call solve_dp_var1_1d(x_futu, x_curr, gamma, lambda, delta, partition, cpt_hat)

print *, "p                  = ", p
print *, "n transitions      = ", n_trans
print *, "gamma              = ", gamma
print *, "lambda             = ", lambda
print *, "delta              = ", delta
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          = "
end if
print *, "partition checksum = ", sum(partition)

deallocate(partition, cpt_hat)
call print_wall_time(t_start)

contains

    subroutine fill_data(dat)
        real(kind=dp), intent(out) :: dat(:, :)
        integer :: t
        dat = 0.0_dp
        dat(:, 1) = [0.2_dp, -0.1_dp]
        do t = 1, n_trans
            if (t <= 40) then
                dat(:, t + 1) = matmul(a1, dat(:, t)) + 0.1_dp * sample_normal_vec()
            else if (t <= 80) then
                dat(:, t + 1) = matmul(a2, dat(:, t)) + 0.1_dp * sample_normal_vec()
            else
                dat(:, t + 1) = matmul(a3, dat(:, t)) + 0.1_dp * sample_normal_vec()
            end if
        end do
    end subroutine fill_data

    function sample_normal_vec() result(z)
        real(kind=dp) :: z(p)
        real(kind=dp) :: u1, u2, r, theta
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, 1.0e-12_dp)
        r = sqrt(-2.0_dp * log(u1))
        theta = 2.0_dp * acos(-1.0_dp) * u2
        z(1) = r * cos(theta)
        z(2) = r * sin(theta)
    end function sample_normal_vec

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_dp_var1
