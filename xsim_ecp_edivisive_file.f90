program xsim_ecp_edivisive_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use ecp_pkg_mod, only: solve_e_divisive_fixedk_2d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: k = 2, min_size = 10
real(kind=dp), parameter :: alpha = 1.0_dp
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:, :)
integer, allocatable :: estimates(:), cluster(:)
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xecp_edivisive_data.txt", data_file)
call read_matrix_file(data_file, x)

call solve_e_divisive_fixedk_2d(x, k, min_size, alpha, estimates, cluster)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(x, 1)
print *, "p                  = ", size(x, 2)
print *, "k                  = ", k
print *, "min.size           = ", min_size
print *, "alpha              = ", alpha
print *, "k.hat              = ", size(estimates) - 1
print *, "estimates          = ", join_ints(estimates)
print *, "cluster checksum   = ", real(sum(cluster), dp)

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

end program xsim_ecp_edivisive_file
