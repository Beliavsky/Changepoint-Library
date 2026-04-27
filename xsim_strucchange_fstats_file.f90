program xsim_strucchange_fstats_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use strucchange_pkg_mod, only: solve_fstats_regression_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: from = 0.2_dp, to = 0.8_dp
character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), fstats(:)
integer, allocatable :: point_idx(:)
integer(kind=long_int) :: t_start
integer :: breakpoint
real(kind=dp) :: min_rss

call system_clock(t_start)
call get_data_file_arg("xstrucchange_fstats_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))

call solve_fstats_regression_1d(y, x, from, to, fstats, point_idx, breakpoint, min_rss)

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "p                  =", size(x, 2)
print *, "from               =", from
print *, "to                 =", to
print *, "breakpoint         =", breakpoint
print *, "max F              =", maxval(fstats)
print *, "min RSS            =", min_rss
print *, "Fstats checksum    =", sum(fstats)

deallocate(dat, y, x, fstats, point_idx)
call print_wall_time(t_start)
end program xsim_strucchange_fstats_file
