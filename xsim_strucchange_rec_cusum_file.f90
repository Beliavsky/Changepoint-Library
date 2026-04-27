program xsim_strucchange_rec_cusum_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use strucchange_pkg_mod, only: solve_efp_rec_cusum_1d
use util_mod, only: print_wall_time
implicit none

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), process(:)
integer(kind=long_int) :: t_start
integer :: imax(1)

call system_clock(t_start)
call get_data_file_arg("xstrucchange_rec_cusum_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))

call solve_efp_rec_cusum_1d(y, x, process)
imax = maxloc(abs(process))

print *, "file               =", trim(data_file)
print *, "n                  =", size(y)
print *, "p                  =", size(x, 2)
print *, "type               =", "Rec-CUSUM"
print *, "argmax abs process =", imax(1) - 1
print *, "max abs process    =", maxval(abs(process))
print *, "process checksum   =", sum(process)

deallocate(dat, y, x, process)
call print_wall_time(t_start)
end program xsim_strucchange_rec_cusum_file
