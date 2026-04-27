program xsim_cpm_joint_hawkins_file
use kind_mod, only: dp, long_int
use compare_io_mod, only: read_series_file, get_data_file_arg
use cpm_pkg_mod, only: solve_cpm_detect_joint_hawkins_1d, solve_cpm_batch_joint_hawkins_1d
use util_mod, only: print_wall_time
implicit none

integer, parameter :: startup = 20
real(kind=dp), parameter :: arl0 = 500.0_dp
character(len=256) :: data_file
real(kind=dp), allocatable :: x(:), thresholds(:), batch_threshold_vec(:), ds_detect(:), ds_batch(:)
integer(kind=long_int) :: t_start
integer :: cp_detect, dt_detect, cp_batch

call system_clock(t_start)
call get_data_file_arg("xcpm_joint_hawkins_data.txt", data_file)
call read_series_file(data_file, x)
call read_series_file("xcpm_joint_hawkins_thresholds.txt", thresholds)
call read_series_file("xcpm_joint_hawkins_batch_threshold.txt", batch_threshold_vec)

call solve_cpm_detect_joint_hawkins_1d(x, thresholds, startup, ds_detect, cp_detect, dt_detect)
call solve_cpm_batch_joint_hawkins_1d(x, batch_threshold_vec(1), ds_batch, cp_batch)

print *, "file               = ", trim(data_file)
print *, "n                  = ", size(x)
print *, "type               = JointHawkins"
print *, "ARL0               = ", arl0
print *, "startup            = ", startup
print *, "alpha              = NA"
print *, "detect changePoint = ", cp_detect
print *, "detect detectionTime = ", dt_detect
print *, "detect Ds checksum = ", sum(ds_detect)
print *, "batch changePoint  = ", cp_batch
print *, "batch threshold    = ", batch_threshold_vec(1)
print *, "batch Ds checksum  = ", sum(ds_batch)

call print_wall_time(t_start)

end program xsim_cpm_joint_hawkins_file
