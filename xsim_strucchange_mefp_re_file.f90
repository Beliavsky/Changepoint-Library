program xsim_strucchange_mefp_re_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use strucchange_pkg_mod, only: solve_mefp_re_1d
use util_mod, only: print_wall_time
implicit none

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), process(:, :), statistic(:), boundary(:)
real(kind=dp) :: alpha, critval
integer(kind=long_int) :: t_start
integer :: history_n, breakpoint

call system_clock(t_start)
call get_data_file_arg("xstrucchange_mefp_re_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))
history_n = 50
alpha = 0.05_dp

call solve_mefp_re_1d(y(1:history_n), x(1:history_n, :), y, x, alpha, process, statistic, boundary, breakpoint, critval)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "history n          =", history_n
print *, "p                  =", size(x, 2)
print *, "type               =", "mefp RE"
print *, "alpha              =", alpha
print *, "critval            =", critval
print *, "breakpoint         =", breakpoint
print *, "nrow process       =", size(process, 1)
print *, "ncol process       =", size(process, 2)
print *, "statistic checksum =", sum(statistic)
print *, "statistic max      =", maxval(statistic)
print *, "boundary checksum  =", sum(boundary)

deallocate(dat, y, x, process, statistic, boundary)
call print_wall_time(t_start)
end program xsim_strucchange_mefp_re_file
