program xsim_ecp_eagglo_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use ecp_pkg_mod, only: solve_e_agglo_default_2d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: alpha = 1.0_dp
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:, :), fit(:)
integer, allocatable :: estimates(:), cluster(:), progression(:, :), merged(:, :)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xecp_eagglo_data.txt", data_file)
call read_matrix_file(data_file, x)

call solve_e_agglo_default_2d(x, alpha, estimates, cluster, fit, progression, merged)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(x, 1)
print *, "p                  = ", size(x, 2)
print *, "alpha              = ", alpha
print *, "estimates          = ", join_ints(estimates)
print *, "cluster checksum   = ", real(sum(cluster), dp)
print *, "fit length         = ", size(fit)
print *, "fit checksum       = ", sum(fit)
print *, "progression dim    = ", size(progression, 1), size(progression, 2)
print *, "progression checksum = ", real(sum(pack(progression, progression > 0)), dp)
print *, "merged checksum    = ", real(sum(merged), dp)

call print_wall_time(t_start)

contains

function join_ints(vals) result(out)
integer, intent(in) :: vals(:)
character(len=256) :: out
integer :: i, pos
out = ""
pos = 1
do i = 1, size(vals)
    write(out(pos:), '(I0)') vals(i)
    pos = len_trim(out) + 1
    if (i < size(vals)) then
        out(pos:pos) = " "
        pos = pos + 1
    end if
end do
end function join_ints

end program xsim_ecp_eagglo_file
