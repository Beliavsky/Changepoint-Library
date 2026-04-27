program xsim_changepointnp_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoint_np_pkg_mod, only: solve_cpt_np_pelt_1d, changepoint_np_penalty_value
use util_mod, only: print_wall_time
implicit none

character(len=*), parameter :: penalty = "MBIC"
integer, parameter :: minseglen = 2, nquantiles = 10
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:), lastchangelike(:)
integer, allocatable :: cpts(:), lastchangecpts(:), numchangecpts(:)
integer(kind=long_int) :: t_start
real(kind=dp) :: pen_value

call system_clock(t_start)
call get_data_file_arg("xchangepointnp_data.txt", data_file)
call read_series_file(data_file, x)

pen_value = changepoint_np_penalty_value(penalty, size(x), 1)
call solve_cpt_np_pelt_1d(x, penalty, minseglen, nquantiles, cpts, lastchangecpts, lastchangelike, numchangecpts)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(x)
print *, "penalty            = ", penalty
print *, "pen.value          = ", pen_value
print *, "minseglen          = ", minseglen
print *, "nquantiles         = ", nquantiles
print *, "cpts               = ", join_ints(cpts)
print *, "lastchangecpts checksum = ", real(sum(lastchangecpts), dp)
print *, "lastchangelike checksum = ", sum(lastchangelike)
print *, "numchangecpts checksum = ", real(sum(numchangecpts), dp)

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

end program xsim_changepointnp_file
