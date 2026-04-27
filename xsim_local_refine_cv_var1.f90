program xsim_local_refine_cv_var1
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_var1_mod, only: solve_cv_dp_var1_1d, solve_local_refine_cv_var1_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: p = 2, n_trans = 120, delta = 5, delta_local = 5
real(kind=dp), parameter :: gamma_set(3) = [0.01_dp, 0.05_dp, 0.1_dp]
real(kind=dp), parameter :: lambda_set(3) = [0.0_dp, 0.05_dp, 0.1_dp]
real(kind=dp), parameter :: zeta_set(4) = [0.01_dp, 0.05_dp, 0.1_dp, 0.2_dp]
real(kind=dp), parameter :: a1(p, p) = reshape([0.6_dp, 0.0_dp, 0.0_dp, 0.4_dp], [p, p])
real(kind=dp), parameter :: a2(p, p) = reshape([-0.4_dp, 0.2_dp, 0.3_dp, 0.5_dp], [p, p])
real(kind=dp), parameter :: a3(p, p) = reshape([0.7_dp, 0.1_dp, -0.2_dp, -0.5_dp], [p, p])
real(kind=dp) :: data(p, n_trans + 1), test_error, train_error, best_test, zeta_best
integer, allocatable :: cpt_hat(:), cpt_init(:), refined(:)
integer :: i, j, k_hat, best_i, best_j
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_data(data)
best_test = huge(1.0_dp)
best_i = 1
best_j = 1

do j = 1, size(lambda_set)
    do i = 1, size(gamma_set)
        call solve_cv_dp_var1_1d(data, gamma_set(i), lambda_set(j), delta, cpt_hat, k_hat, test_error, train_error)
        if (test_error < best_test) then
            best_test = test_error
            best_i = i
            best_j = j
            if (allocated(cpt_init)) deallocate(cpt_init)
            allocate(cpt_init(size(cpt_hat)))
            if (size(cpt_hat) > 0) cpt_init = cpt_hat
        end if
        if (allocated(cpt_hat)) deallocate(cpt_hat)
    end do
end do

call solve_local_refine_cv_var1_1d(cpt_init, data, zeta_set, delta_local, refined, zeta_best)

print *, "p                  = ", p
print *, "n                  = ", n_trans + 1
print *, "gamma_set          = ", gamma_set
print *, "lambda_set         = ", lambda_set
print *, "zeta_set           = ", zeta_set
print *, "delta              = ", delta
print *, "delta_local        = ", delta_local
print *, "best gamma         = ", gamma_set(best_i)
print *, "best lambda        = ", lambda_set(best_j)
if (size(cpt_init) > 0) then
    write (*,'(A)', advance='no') " initial            ="
    call print_int_list(cpt_init)
else
    print *, "initial            = "
end if
if (size(refined) > 0) then
    write (*,'(A)', advance='no') " refined            ="
    call print_int_list(refined)
else
    print *, "refined            = "
end if
print *, "zeta               = ", zeta_best

if (allocated(cpt_init)) deallocate(cpt_init)
if (allocated(refined)) deallocate(refined)
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
        integer :: ii
        do ii = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(ii)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_local_refine_cv_var1
