program probe_obj
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file
use changepoints_pkg_var1_mod
implicit none
real(kind=dp), allocatable :: data(:,:), data_t(:,:), x_curr(:,:), x_futu(:,:)
real(kind=dp) :: v
integer :: eta
call read_matrix_file('xcv_dp_var1_data.txt', data)
data_t = transpose(data)
x_curr = data_t(:,1:size(data_t,2)-1)
x_futu = data_t(:,2:size(data_t,2))
do eta = 62, 110
  if (eta==70 .or. eta==71 .or. eta==72 .or. eta==106) then
    v = obj_lr_var1_1d(57, 120, eta, x_futu, x_curr, 0.05_dp)
    print *, 'eta=', eta, ' obj=', v
  end if
end do
end program probe_obj
