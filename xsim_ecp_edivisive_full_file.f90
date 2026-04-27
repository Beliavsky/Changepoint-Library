program xsim_ecp_edivisive_full_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use ecp_pkg_mod, only: solve_e_divisive_permtest_2d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: sig_lvl = 0.1_dp, alpha = 1.0_dp
integer, parameter :: r = 19, min_size = 10
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:, :), perms_real(:, :)
integer, allocatable :: perms(:, :), order_found(:), estimates(:), permutations(:), cluster(:)
real(kind=dp), allocatable :: p_values(:)
integer(kind=long_int) :: t_start
integer :: considered_last

call system_clock(t_start)
call get_data_file_arg("xecp_edivisive_full_data.txt", data_file)
call read_matrix_file(data_file, x)
call read_matrix_file("xecp_edivisive_full_perms.txt", perms_real)
allocate(perms(size(perms_real, 1), size(perms_real, 2)))
perms = nint(perms_real)

call solve_e_divisive_permtest_2d(x, sig_lvl, r, min_size, alpha, perms, order_found, estimates, considered_last, p_values, permutations, cluster)

print *, "file               = ", trim(data_file)
print *, "perm_file          = ", "xecp_edivisive_full_perms.txt"
print *, "n                  = ", size(x, 1)
print *, "p                  = ", size(x, 2)
print *, "sig.lvl            = ", sig_lvl
print *, "R                  = ", r
print *, "min.size           = ", min_size
print *, "alpha              = ", alpha
print *, "k.hat              = ", size(estimates) - 1
print *, "estimates          = ", join_ints(estimates)
print *, "order.found        = ", join_ints(order_found)
print *, "considered.last    = ", considered_last
print *, "p.values           = ", join_reals(p_values)
print *, "permutations       = ", join_ints(permutations)
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

function join_reals(vals) result(out)
real(kind=dp), intent(in) :: vals(:)
character(len=256) :: out
integer :: i, pos
out = ""
pos = 1
do i = 1, size(vals)
    write(out(pos:), '(F0.12)') vals(i)
    pos = len_trim(out) + 1
    if (i < size(vals)) then
        out(pos:pos) = " "
        pos = pos + 1
    end if
end do
end function join_reals

end program xsim_ecp_edivisive_full_file
