program xsim_probe_local_refine_cv_var1_file
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use changepoints_pkg_var1_mod, only: probe_local_refine_cv_var1_case_1d
implicit none

character(len=256) :: filename
real(kind=dp), allocatable :: data_in(:, :), data(:, :)

call get_data_file_arg('xcv_dp_var1_data.txt', filename)
call read_matrix_file(trim(filename), data_in)
allocate(data(size(data_in, 2), size(data_in, 1)))
data = transpose(data_in)

write(*, '(A,A)') 'file               = ', trim(filename)
write(*, '(A,I0)') 'p                  = ', size(data, 1)
write(*, '(A,I0)') 'n                  = ', size(data, 2)
call probe_local_refine_cv_var1_case_1d(data)

deallocate(data_in, data)
end program xsim_probe_local_refine_cv_var1_file
