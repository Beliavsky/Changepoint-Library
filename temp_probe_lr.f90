program probe_lr
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file
use changepoints_pkg_var1_mod, only: solve_local_refine_cv_var1_1d
implicit none
real(kind=dp), allocatable :: data(:,:), data_t(:,:)
integer, allocatable :: cpt_init(:), refined(:)
real(kind=dp) :: zeta_best
allocate(cpt_init(2))
cpt_init = [56, 80]
call read_matrix_file('xcv_dp_var1_data.txt', data)
data_t = transpose(data)
call solve_local_refine_cv_var1_1d(cpt_init, data_t, [0.01_dp,0.05_dp,0.1_dp,0.2_dp], 5, refined, zeta_best)
print *, 'zeta=', zeta_best
if (allocated(refined)) then
  write(*,'(A)',advance='no') 'refined='
  if (size(refined) > 0) then
    write(*,'(1X,I0,1X,I0)') refined(1), refined(2)
  else
    write(*,*)
  end if
end if
end program probe_lr
