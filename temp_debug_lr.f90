program debug_lr
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file
use changepoints_pkg_var1_mod
implicit none
real(kind=dp), allocatable :: data(:,:), data_t(:,:), data_temp(:,:), x_curr(:,:), x_futu(:,:), x_curr_train(:,:), x_curr_test(:,:), x_futu_train(:,:), x_futu_test(:,:)
real(kind=dp) :: zb
integer :: est
call read_matrix_file('xcv_dp_var1_data.txt', data)
data_t = transpose(data)
data_temp = data_t
x_curr = data_temp(:,1:size(data_temp,2)-1)
x_futu = data_temp(:,2:size(data_temp,2))
x_curr_train = x_curr(:,1:size(x_curr,2):2)
x_curr_test = x_curr(:,2:size(x_curr,2):2)
x_futu_train = x_futu(:,1:size(x_futu,2):2)
x_futu_test = x_futu(:,2:size(x_futu,2):2)
call glasso_error_test_var1_1d(1, 40, x_futu_train, x_curr_train, x_futu_test, x_curr_test, 5, [0.01_dp,0.05_dp,0.1_dp,0.2_dp], zb)
print *, 'z1=', zb
call find_one_change_grouplasso_var1_1d(2, 50, x_futu, x_curr, 5, zb, est)
print *, 'est1full=', est
call glasso_error_test_var1_1d(20, 60, x_futu_train, x_curr_train, x_futu_test, x_curr_test, 5, [0.01_dp,0.05_dp,0.1_dp,0.2_dp], zb)
print *, 'z2=', zb
call find_one_change_grouplasso_var1_1d(40, 110, x_futu, x_curr, 5, zb, est)
print *, 'est2full=', est
end program debug_lr
