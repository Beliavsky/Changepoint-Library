program inspect_part
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file
use changepoints_pkg_var1_mod, only: solve_dp_var1_1d
implicit none
real(kind=dp), allocatable :: data(:,:), data_t(:,:), x_curr(:,:), x_futu(:,:)
integer, allocatable :: partition(:), cpt_hat(:)
integer :: i
call read_matrix_file('xdp_var1_data.txt', data)
data_t = transpose(data)
x_curr = data_t(:,1:size(data_t,2)-1)
x_futu = data_t(:,2:size(data_t,2))
call solve_dp_var1_1d(x_futu, x_curr, 0.2_dp, 0.0_dp, 5, partition, cpt_hat)
print *, 'k=', size(cpt_hat)
if (size(cpt_hat) > 0) then
  write(*,'(A)',advance='no') 'cpt='
  do i=1,size(cpt_hat)
    write(*,'(1X,I0)',advance='no') cpt_hat(i)
  end do
  write(*,*)
end if
print *, 'sum=', sum(partition)
open(unit=10,file='temp_var1_partition_f.txt',status='replace')
do i=1,size(partition)
  write(10,'(I0,1X)',advance='no') partition(i)
end do
write(10,*)
close(10)
end program inspect_part
