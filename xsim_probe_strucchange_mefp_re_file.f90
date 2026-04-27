program xsim_probe_strucchange_mefp_re_file
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use strucchange_pkg_mod, only: solve_mefp_re_1d
implicit none

character(len=256) :: data_file
real(kind=dp), allocatable :: dat(:, :), y(:), x(:, :), process(:, :), statistic(:), boundary(:)
real(kind=dp) :: alpha, critval
integer :: history_n, breakpoint, idx(6), j

call get_data_file_arg("xstrucchange_mefp_re_data.txt", data_file)
call read_matrix_file(data_file, dat)
allocate(y(size(dat, 1)), x(size(dat, 1), size(dat, 2) - 1))
y = dat(:, 1)
x = dat(:, 2:size(dat, 2))
history_n = 50
alpha = 0.05_dp

call solve_mefp_re_1d(y(1:history_n), x(1:history_n, :), y, x, alpha, process, statistic, boundary, breakpoint, critval)

idx = [1, 2, 10, 20, 40, 70]
print *, "critval =", critval
print *, "breakpoint =", breakpoint
do j = 1, size(idx)
    print *, "row =", idx(j), "proc1 =", process(idx(j), 1), "proc2 =", process(idx(j), 2), "stat =", statistic(idx(j))
end do

deallocate(dat, y, x, process, statistic, boundary)
end program xsim_probe_strucchange_mefp_re_file
