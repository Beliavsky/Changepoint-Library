program xsim_segmented_multi_file
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use segmented_pkg_mod, only: eval_segmented_lm2_draws
implicit none

character(len=256) :: data_file
character(len=*), parameter :: draws_file = "xsegmented_multi_draws.txt"
real(kind=dp), allocatable :: dat(:, :), draws(:, :), fitted_mean(:), param_means(:)

call get_data_file_arg("xsegmented_multi_data.txt", data_file)
call read_matrix_file(data_file, dat)
call read_matrix_file(draws_file, draws)
call eval_segmented_lm2_draws(dat(:, 1), draws, fitted_mean, param_means)

print *, "file               =", trim(data_file)
print *, "draws_file         =", draws_file
print *, "n                  =", size(dat, 1)
print *, "psi checksum       =", sum(param_means(1:2))
print *, "coef checksum      =", sum(param_means(3:6))
print *, "fitted checksum    =", sum(fitted_mean)
print *, "fitted head checksum =", sum(fitted_mean(:min(10, size(fitted_mean))))
print *, "intercept          =", param_means(3)
print *, "x slope            =", param_means(4)
print *, "hinge1 slope       =", param_means(5)
print *, "hinge2 slope       =", param_means(6)
end program xsim_segmented_multi_file
