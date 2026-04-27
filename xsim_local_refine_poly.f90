program xsim_local_refine_poly
! simulate univariate polynomial changepoints, cross-validate DP.poly, and refine locally
use kind_mod, only: dp, long_int
use compare_sim_mod, only: seed_rng_fixed
use changepoints_pkg_mod, only: solve_cv_dp_poly_1d, solve_local_refine_poly_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: n = 150 ! number of observations
integer, parameter :: r = 2 ! polynomial order
integer, parameter :: n_gamma = 4 ! number of candidate penalties
real(kind=dp), parameter :: gamma_set(n_gamma) = [0.5_dp, 1.0_dp, 2.0_dp, 3.0_dp] ! candidate l0 penalties
integer, parameter :: delta = 5 ! minimum spacing
integer, parameter :: delta_lr = 5 ! local refinement spacing

real(kind=dp) :: y(n), test_error(n_gamma), train_error(n_gamma)
real(kind=dp), parameter :: coef1(3) = [0.0_dp, 1.5_dp, -2.5_dp]
real(kind=dp), parameter :: coef2(3) = [2.5_dp, -8.0_dp, 10.0_dp]
real(kind=dp), parameter :: coef3(3) = [-2.0_dp, 5.0_dp, -2.5_dp]
integer, allocatable :: cpt_hat(:, :), k_hat(:), cpt_one(:), cpt_refined(:)
integer(kind=long_int) :: t_start
integer :: j, max_k_local, min_idx(1)

call system_clock(t_start)
call seed_rng_fixed(123)
call fill_response(y)

allocate(k_hat(n_gamma))
max_k_local = 0
do j = 1, n_gamma
    call solve_cv_dp_poly_1d(y, r, gamma_set(j), delta, cpt_one, k_hat(j), test_error(j), train_error(j))
    max_k_local = max(max_k_local, k_hat(j))
    deallocate(cpt_one)
end do
allocate(cpt_hat(max_k_local, n_gamma))
cpt_hat = 0
do j = 1, n_gamma
    call solve_cv_dp_poly_1d(y, r, gamma_set(j), delta, cpt_one, k_hat(j), test_error(j), train_error(j))
    if (k_hat(j) > 0) cpt_hat(1:k_hat(j), j) = cpt_one
    deallocate(cpt_one)
end do
min_idx = minloc(test_error)
call solve_local_refine_poly_1d(cpt_hat(1:k_hat(min_idx(1)), min_idx(1)), y, r, delta_lr, cpt_refined)

print *, "n                  = ", n
print *, "r                  = ", r
print *, "gamma_set          = ", gamma_set
print *, "delta              = ", delta
print *, "delta_lr           = ", delta_lr
print *, "best gamma         = ", gamma_set(min_idx(1))
if (k_hat(min_idx(1)) > 0) then
    write (*,'(A)', advance='no') " initial            ="
    call print_int_list(cpt_hat(1:k_hat(min_idx(1)), min_idx(1)))
else
    print *, "initial            = "
end if
if (size(cpt_refined) > 0) then
    write (*,'(A)', advance='no') " refined            ="
    call print_int_list(cpt_refined)
else
    print *, "refined            = "
end if

deallocate(cpt_hat, k_hat, cpt_refined)
call print_wall_time(t_start)

contains

    subroutine fill_response(y_out)
        real(kind=dp), intent(out) :: y_out(:)
        integer :: i
        real(kind=dp) :: xval

        do i = 1, 50
            xval = real(i, dp) / real(n, dp)
            y_out(i) = coef1(1) + coef1(2) * xval + coef1(3) * xval * xval + 0.05_dp * sample_standard_normal()
        end do
        do i = 51, 100
            xval = real(i, dp) / real(n, dp)
            y_out(i) = coef2(1) + coef2(2) * xval + coef2(3) * xval * xval + 0.05_dp * sample_standard_normal()
        end do
        do i = 101, 150
            xval = real(i, dp) / real(n, dp)
            y_out(i) = coef3(1) + coef3(2) * xval + coef3(3) * xval * xval + 0.05_dp * sample_standard_normal()
        end do
    end subroutine fill_response

    real(kind=dp) function sample_standard_normal() result(z)
        real(kind=dp) :: u1, u2
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, 1.0e-12_dp)
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
    end function sample_standard_normal

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_local_refine_poly
