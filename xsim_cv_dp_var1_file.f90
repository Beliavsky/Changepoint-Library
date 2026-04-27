program xsim_cv_dp_var1_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_var1_mod, only: solve_cv_dp_var1_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: gamma_set(3) = [0.01_dp, 0.05_dp, 0.1_dp]
real(kind=dp), parameter :: lambda_set(3) = [0.0_dp, 0.05_dp, 0.1_dp]
integer, parameter :: delta = 5

character(len=256) :: data_file
real(kind=dp), allocatable :: data(:, :), data_t(:, :)
real(kind=dp) :: test_error, train_error, best_test
integer, allocatable :: cpt_hat(:), best_cpt(:)
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

print *, "file               = ", trim(data_file)
print *, "p                  = ", size(data_t, 1)
print *, "n                  = ", size(data_t, 2)
print *, "gamma_set          = ", gamma_set
print *, "lambda_set         = ", lambda_set
print *, "delta              = ", delta

do j = 1, size(lambda_set)
    do i = 1, size(gamma_set)
        call solve_cv_dp_var1_1d(data_t, gamma_set(i), lambda_set(j), delta, cpt_hat, k_hat, test_error, train_error)
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

end program xsim_cv_dp_var1_file
