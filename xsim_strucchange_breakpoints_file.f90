program xsim_strucchange_breakpoints_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use strucchange_pkg_mod, only: solve_breakpoints_regression_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: h = 20
character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), rss_vec(:), bic_vec(:)
integer, allocatable :: cpt_hat(:)
integer(kind=long_int) :: t_start
integer :: best_m

call system_clock(t_start)
call get_data_file_arg("xstrucchange_breakpoints_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))

call solve_breakpoints_regression_1d(y, x, h, cpt_hat, rss_vec, bic_vec, best_m)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "p                  =", size(x, 2)
print *, "h                  =", h
print *, "max breaks         =", ubound(bic_vec, 1)
print *, "best m             =", best_m
if (size(cpt_hat) > 0) then
    write (*,'(A)', advance='no') " estimated          ="
    call print_int_list(cpt_hat)
else
    print *, "estimated          ="
end if
write (*,'(A)', advance='no') " RSS                ="
call print_real_list(rss_vec)
write (*,'(A)', advance='no') " BIC                ="
call print_real_list(bic_vec)

deallocate(dat, y, x, cpt_hat, rss_vec, bic_vec)
call print_wall_time(t_start)

contains

    subroutine print_int_list(v)
        integer, intent(in) :: v(:)
        integer :: i
        do i = 1, size(v)
            write (*,'(1X,I0)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_int_list

    subroutine print_real_list(v)
        real(kind=dp), intent(in) :: v(0:)
        integer :: i
        do i = lbound(v, 1), ubound(v, 1)
            write (*,'(1X,F0.12)', advance='no') v(i)
        end do
        write (*,*)
    end subroutine print_real_list

end program xsim_strucchange_breakpoints_file
