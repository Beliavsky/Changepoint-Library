program xsim_ecp_cp3o_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use ecp_pkg_mod, only: solve_e_cp3o_2d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: kmax = 3, minsize = 10
real(kind=dp), parameter :: alpha = 1.0_dp
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:, :), gofm(:)
integer, allocatable :: estimates(:), cploc(:, :), cploc_lens(:)
integer :: number, i
integer(kind=long_int) :: t_start

call system_clock(t_start)
call get_data_file_arg("xecp_cp3o_data.txt", data_file)
call read_matrix_file(data_file, x)

call solve_e_cp3o_2d(x, kmax, minsize, alpha, number, estimates, gofm, cploc, cploc_lens)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(x, 1)
print *, "p                  = ", size(x, 2)
print *, "K                  = ", kmax
print *, "minsize            = ", minsize
print *, "alpha              = ", alpha
print *, "number             = ", number
print *, "estimates          = ", join_ints(estimates)
print *, "gof length         = ", size(gofm)
print *, "gof checksum       = ", sum(gofm)
print *, "cpLoc checksum     = ", real(sum(cploc), dp)
print *, "cpLoc lengths      = ", join_ints(cploc_lens)
do i = 1, size(cploc, 1)
    if (cploc_lens(i) > 0) print *, "cpLoc ", i, " = ", join_ints(cploc(i, 1:cploc_lens(i)))
end do

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

end program xsim_ecp_cp3o_file
