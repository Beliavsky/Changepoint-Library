program xsim_strucchange_me_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use strucchange_pkg_mod, only: solve_efp_me_1d
use util_mod, only: print_wall_time
implicit none

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), process(:, :)
real(kind=dp) :: h
integer(kind=long_int) :: t_start
integer :: imax(2)

call system_clock(t_start)
call get_data_file_arg("xstrucchange_me_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))
h = 0.25_dp

call solve_efp_me_1d(y, x, h, process)
imax = maxloc(abs(process))

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "p                  =", size(x, 2)
print *, "type               =", "ME"
print *, "h                  =", h
print *, "nrow process       =", size(process, 1)
print *, "ncol process       =", size(process, 2)
print *, "imax row           =", imax(1)
print *, "imax col           =", imax(2)
print *, "max abs process    =", maxval(abs(process))
print *, "process checksum   =", sum(process)

deallocate(dat, y, x, process)
call print_wall_time(t_start)
end program xsim_strucchange_me_file
