program xsim_mcp_ar_file
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use mcp_pkg_mod, only: eval_mcp_ar1_draws
implicit none

character(len=256) :: data_file
character(len=*), parameter :: draws_file = "xmcp_ar_draws.txt"
real(kind=dp), allocatable :: dat(:, :), draws(:, :), fitted_mean(:), param_means(:)

call get_data_file_arg("xmcp_ar_data.txt", data_file)
call read_matrix_file(data_file, dat)
call read_matrix_file(draws_file, draws)
call eval_mcp_ar1_draws(dat(:, 1), draws, fitted_mean, param_means)

print *, "file               =", trim(data_file)
print *, "draws_file         =", draws_file
print *, "n                  =", size(dat, 1)
print *, "ndraw              =", size(draws, 1)
print *, "param mean checksum =", sum(param_means)
print *, "fitted checksum    =", sum(fitted_mean)
print *, "fitted head checksum =", sum(fitted_mean(:min(10, size(fitted_mean))))
print *, "ar mean            =", param_means(1)
print *, "int mean           =", param_means(2)
print *, "sigma mean         =", param_means(3)
end program xsim_mcp_ar_file
