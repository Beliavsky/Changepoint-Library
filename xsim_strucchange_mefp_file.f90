program xsim_strucchange_mefp_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use strucchange_pkg_mod, only: solve_mefp_ols_cusum_1d
use util_mod, only: print_wall_time
implicit none

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), process(:), boundary(:)
real(kind=dp) :: alpha, critval
integer(kind=long_int) :: t_start
integer :: history_n, breakpoint, imax(1)

call system_clock(t_start)
call get_data_file_arg("xstrucchange_mefp_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))
history_n = 50
alpha = 0.05_dp

call solve_mefp_ols_cusum_1d(y(1:history_n), x(1:history_n, :), y, x, alpha, process, boundary, breakpoint, critval)
imax = maxloc(abs(process))

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "history n          =", history_n
print *, "p                  =", size(x, 2)
print *, "type               =", "mefp OLS-CUSUM"
print *, "alpha              =", alpha
print *, "critval            =", critval
print *, "breakpoint         =", breakpoint
print *, "argmax abs process =", imax(1) - 1
print *, "max abs process    =", maxval(abs(process))
print *, "process checksum   =", sum(process)
print *, "boundary checksum  =", sum(boundary)

deallocate(dat, y, x, process, boundary)
call print_wall_time(t_start)
end program xsim_strucchange_mefp_file
