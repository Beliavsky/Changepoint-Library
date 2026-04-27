program xsim_local_refine_cv_var1_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_var1_mod, only: solve_cv_dp_var1_1d, solve_local_refine_cv_var1_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: gamma_set(3) = [0.01_dp, 0.05_dp, 0.1_dp]
real(kind=dp), parameter :: lambda_set(3) = [0.0_dp, 0.05_dp, 0.1_dp]
real(kind=dp), parameter :: zeta_set(4) = [0.01_dp, 0.05_dp, 0.1_dp, 0.2_dp]
integer, parameter :: delta = 5, delta_local = 5

character(len=256) :: data_file
real(kind=dp), allocatable :: data(:, :), data_t(:, :)
real(kind=dp) :: test_error, train_error, best_test, zeta_best
integer, allocatable :: cpt_hat(:), cpt_init(:), refined(:)
integer :: i, j, k_hat, best_i, best_j
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xcv_dp_var1_data.txt", data_file)
call read_matrix_file(data_file, data)
allocate(data_t(size(data, 2), size(data, 1)))
data_t = transpose(data)
best_test = huge(1.0_dp)
best_i = 1
best_j = 1

do j = 1, size(lambda_set)
    do i = 1, size(gamma_set)
        call solve_cv_dp_var1_1d(data_t, gamma_set(i), lambda_set(j), delta, cpt_hat, k_hat, test_error, train_error)
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

call solve_local_refine_cv_var1_1d(cpt_init, data_t, zeta_set, delta_local, refined, zeta_best)

print *, "file               = ", trim(data_file)
print *, "p                  = ", size(data_t, 1)
print *, "n                  = ", size(data_t, 2)
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
deallocate(data, data_t)
call print_wall_time(t_start)

contains

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: ii
        do ii = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(ii)
        end do
        write (*,*)
    end subroutine print_int_list

end program xsim_local_refine_cv_var1_file
