program xsim_capa_file
! read a univariate series from a text file and detect collective and point anomalies by CAPA
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use changepoint_mod, only: solve_capa_l2_1d
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: segment_penalty = 24.0_dp ! CAPA penalty for collective anomalies
real(kind=dp), parameter :: point_penalty = 14.0_dp ! CAPA penalty for point anomalies
integer, parameter :: min_segment_length = 10 ! minimum allowed collective anomaly length
integer, parameter :: max_segment_length = 500 ! maximum allowed collective anomaly length

character(len=256) :: data_file
real(kind=dp), allocatable :: x(:)
integer, allocatable :: segment_starts(:), segment_ends(:), point_locs(:)
integer(kind=long_int) :: t_start
integer :: i

call system_clock(t_start)
call get_data_file_arg("xcapa_data.txt", data_file)
call read_series_file(data_file, x)
call solve_capa_l2_1d(x, segment_penalty, point_penalty, segment_starts, segment_ends, point_locs, min_segment_length, max_segment_length)

print *, "file             =", trim(data_file)
print *, "n                =", size(x)
print *, "segment penalty  =", segment_penalty
print *, "point penalty    =", point_penalty
write(*, '(A)', advance='no') 'segment anomalies='
if (size(segment_starts) > 0) then
    do i = 1, size(segment_starts)
        write(*, '(" [", I0, ",", I0, ")")', advance='no') segment_starts(i), segment_ends(i)
    end do
else
    write(*, '(A)', advance='no') ' none'
end if
print *
write(*, '(A)', advance='no') 'point anomalies  ='
if (size(point_locs) > 0) then
    do i = 1, size(point_locs)
        write(*, '(" [", I0, ",", I0, ")")', advance='no') point_locs(i), point_locs(i) + 1
    end do
else
    write(*, '(A)', advance='no') ' none'
end if
print *

deallocate(x, segment_starts, segment_ends, point_locs)
call print_wall_time(t_start)

end program xsim_capa_file
