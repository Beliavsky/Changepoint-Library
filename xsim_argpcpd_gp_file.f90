program xsim_argpcpd_gp_file
! read a univariate series from a text file and estimate changepoints by a GP-based ArgpCpd approximation
use kind_mod, only: dp, long_int
use changepoint_mod, only: solve_argpcpd_gp_1d, solve_argpcpd_gp_rust_1d, reduce_changepoints_normal_1d
use compare_io_mod, only: read_series_file, get_data_file_arg
use util_mod, only: print_wall_time
implicit none

real(kind=dp), parameter :: scale = 3.0_dp ! scale of the constant kernel factor
real(kind=dp), parameter :: length_scale = 10.0_dp ! length scale of the RBF kernel
real(kind=dp), parameter :: noise_level = 0.01_dp ! white-noise kernel level
integer, parameter :: max_lag = 12 ! maximum autoregressive lag in the lag-feature representation
real(kind=dp), parameter :: alpha0 = 2.0_dp ! scale-gamma alpha parameter
real(kind=dp), parameter :: beta0 = 1.0_dp ! scale-gamma beta parameter
real(kind=dp), parameter :: logistic_hazard_h = -5.0_dp ! hazard scale in logit units
real(kind=dp), parameter :: logistic_hazard_a = 1.0_dp ! logistic hazard slope parameter
real(kind=dp), parameter :: logistic_hazard_b = 1.0_dp ! logistic hazard offset parameter
real(kind=dp), parameter :: epsilon = 1.0e-10_dp ! tail-probability threshold for truncating run-length states
integer, parameter :: max_train = 600 ! maximum GP training suffix length used in each predictive update

character(len=256) :: data_file
character(len=32) :: arg, mode
real(kind=dp), allocatable :: x(:)
integer, allocatable :: map_cps(:), reduced_cps(:), runlen_argmax(:)
integer(kind=long_int) :: t_start
integer :: ios, target_ncp, debug_start, debug_end, t

call system_clock(t_start)
call get_data_file_arg("xargpcpd_data.txt", data_file)
call read_series_file(data_file, x)

target_ncp = 0
mode = "dense"
debug_start = 0
debug_end = 0
call get_command_argument(2, arg, status=ios)
if (ios == 0 .and. len_trim(arg) > 0) then
    read(arg, *, iostat=ios) target_ncp
    if (ios /= 0) then
        mode = trim(adjustl(arg))
        call get_command_argument(3, arg, status=ios)
        if (ios == 0 .and. len_trim(arg) > 0) read(arg, *, iostat=ios) target_ncp
    end if
end if
call get_command_argument(4, arg, status=ios)
if (ios == 0 .and. len_trim(arg) > 0) read(arg, *, iostat=ios) debug_start
call get_command_argument(5, arg, status=ios)
if (ios == 0 .and. len_trim(arg) > 0) read(arg, *, iostat=ios) debug_end

select case (trim(mode))
case ("dense")
    call solve_argpcpd_gp_1d(x, map_cps, runlen_argmax, scale, length_scale, noise_level, max_lag, alpha0, beta0, &
                             logistic_hazard_h, logistic_hazard_a, logistic_hazard_b, epsilon, max_train)
case ("rust")
    call solve_argpcpd_gp_rust_1d(x, map_cps, runlen_argmax, scale, length_scale, noise_level, max_lag, alpha0, beta0, &
                                  logistic_hazard_h, logistic_hazard_a, logistic_hazard_b, epsilon)
case default
    error stop "mode must be 'dense' or 'rust'"
end select

if (target_ncp >= 1 .and. target_ncp < size(map_cps)) then
    call reduce_changepoints_normal_1d(x, map_cps, target_ncp, reduced_cps, min_seg_len=max_lag+1)
end if

print *, "file        =", trim(data_file)
print *, "n           =", size(x)
print *, "mode        =", trim(mode)
print *, "max_lag     =", max_lag
if (trim(mode) == "dense") print *, "max_train   =", max_train
print *, "hazard_h    =", logistic_hazard_h
if (size(map_cps) > 0) print *, "estimated   =", map_cps
if (allocated(reduced_cps)) print *, "reduced     =", reduced_cps
if (debug_start > 0) then
    if (debug_end <= 0) debug_end = min(size(runlen_argmax), debug_start + 20)
    debug_start = max(1, min(debug_start, size(runlen_argmax)))
    debug_end = max(debug_start, min(debug_end, size(runlen_argmax)))
    print *, "runlen_argmax window:"
    do t = debug_start, debug_end
        write(*, '(i6,1x,i6)') t, runlen_argmax(t)
    end do
end if

if (allocated(reduced_cps)) deallocate(reduced_cps)
deallocate(x, map_cps, runlen_argmax)
call print_wall_time(t_start)

end program xsim_argpcpd_gp_file
