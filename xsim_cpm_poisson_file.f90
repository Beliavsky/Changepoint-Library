program xsim_cpm_poisson_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use cpm_pkg_mod, only: solve_cpm_detect_poisson_1d, solve_cpm_process_poisson_1d, solve_cpm_batch_poisson_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: startup = 20
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:), thresholds(:), batch_threshold_vec(:), ds_detect(:), ds_batch(:)
integer, allocatable :: cps(:), dts(:)
integer(kind=long_int) :: t_start
integer :: cp_detect, dt_detect, cp_batch

call system_clock(t_start)
call get_data_file_arg("xcpm_poisson_data.txt", data_file)
call read_series_file(data_file, x)
call read_series_file("xcpm_poisson_thresholds.txt", thresholds)
call read_series_file("xcpm_poisson_batch_threshold.txt", batch_threshold_vec)

call solve_cpm_detect_poisson_1d(x, thresholds, startup, ds_detect, cp_detect, dt_detect)
call solve_cpm_process_poisson_1d(x, thresholds, startup, cps, dts)
call solve_cpm_batch_poisson_1d(x, batch_threshold_vec(1), ds_batch, cp_batch)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(x)
print *, "type               = Poisson"
print *, "ARL0               = NA"
print *, "startup            = ", startup
print *, "alpha              = NA"
print *, "detect changePoint = ", cp_detect
print *, "detect detectionTime = ", dt_detect
print *, "detect Ds checksum = ", sum(ds_detect)
print *, "process changePoints = ", join_ints(cps)
print *, "process detectionTimes = ", join_ints(dts)
print *, "batch changePoint  = ", cp_batch
print *, "batch threshold    = ", batch_threshold_vec(1)
print *, "batch Ds checksum  = ", sum(ds_batch)

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

end program xsim_cpm_poisson_file
