program xsim_cv_dp_var1
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_var1_mod, only: solve_cv_dp_var1_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: p = 2, n_trans = 120, delta = 5
real(kind=dp), parameter :: gamma_set(3) = [0.01_dp, 0.05_dp, 0.1_dp]
real(kind=dp), parameter :: lambda_set(3) = [0.0_dp, 0.05_dp, 0.1_dp]
real(kind=dp) :: data(p, n_trans + 1), test_error, train_error, best_test
real(kind=dp), parameter :: a1(p, p) = reshape([0.6_dp, 0.0_dp, 0.0_dp, 0.4_dp], [p, p])
real(kind=dp), parameter :: a2(p, p) = reshape([-0.4_dp, 0.2_dp, 0.3_dp, 0.5_dp], [p, p])
real(kind=dp), parameter :: a3(p, p) = reshape([0.7_dp, 0.1_dp, -0.2_dp, -0.5_dp], [p, p])
integer, allocatable :: cpt_hat(:), best_cpt(:)
integer :: i, j, k_hat, best_i, best_j
integer(kind=long_int) :: t_start

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_data(data)
best_test = huge(1.0_dp)
best_i = 1
best_j = 1

print *, "p                  = ", p
print *, "n                  = ", n_trans + 1
print *, "gamma_set          = ", gamma_set
print *, "lambda_set         = ", lambda_set
print *, "delta              = ", delta

do j = 1, size(lambda_set)
    do i = 1, size(gamma_set)
        call solve_cv_dp_var1_1d(data, gamma_set(i), lambda_set(j), delta, cpt_hat, k_hat, test_error, train_error)
        print *, "gamma              = ", gamma_set(i)
        print *, "lambda             = ", lambda_set(j)
        if (size(cpt_hat) > 0) then
            write (*,'(A)', advance='no') " estimated          ="
            call print_int_list(cpt_hat)
        else
            print *, "estimated          = "
        end if
        print *, "n changepoints     = ", k_hat
        print *, "test error         = ", test_error
        print *, "train error        = ", train_error
        if (test_error < best_test) then
            best_test = test_error
            best_i = i
            best_j = j
            if (allocated(best_cpt)) deallocate(best_cpt)
            allocate(best_cpt(size(cpt_hat)))
            if (size(cpt_hat) > 0) best_cpt = cpt_hat
        end if
        if (allocated(cpt_hat)) deallocate(cpt_hat)
    end do
end do

print *, "best gamma         = ", gamma_set(best_i)
print *, "best lambda        = ", lambda_set(best_j)
print *, "best index         = ", best_i, best_j
if (allocated(best_cpt)) then
    if (size(best_cpt) > 0) then
        write (*,'(A)', advance='no') " best changepoints  ="
        call print_int_list(best_cpt)
    else
        print *, "best changepoints  = "
    end if
end if
print *, "best test error    = ", best_test

if (allocated(best_cpt)) deallocate(best_cpt)
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

end program xsim_cv_dp_var1
