program xsim_mcp_arsigma_file
use kind_mod, only: dp
use compare_io_mod, only: read_matrix_file, get_data_file_arg
use mcp_pkg_mod, only: eval_mcp_arsigma_draws
implicit none

character(len=256) :: data_file
character(len=*), parameter :: draws_file = "xmcp_arsigma_draws.txt"
real(kind=dp), allocatable :: dat(:, :), draws(:, :), fitted_mean(:), sigma_mean(:), param_means(:)
real(kind=dp) :: cp_mean

call get_data_file_arg("xmcp_arsigma_data.txt", data_file)
call read_matrix_file(data_file, dat)
call read_matrix_file(draws_file, draws)
call eval_mcp_arsigma_draws(dat(:, 1), draws, fitted_mean, sigma_mean, cp_mean, param_means)

print *, "file               =", trim(data_file)
print *, "draws_file         =", draws_file
print *, "n                  =", size(dat, 1)
print *, "ndraw              =", size(draws, 1)
print *, "cp mean            =", cp_mean
print *, "param mean checksum =", sum(param_means)
print *, "fitted checksum    =", sum(fitted_mean)
print *, "fitted head checksum =", sum(fitted_mean(:min(10, size(fitted_mean))))
print *, "sigma fitted checksum =", sum(sigma_mean)
print *, "sigma head checksum =", sum(sigma_mean(:min(10, size(sigma_mean))))
print *, "ar checksum        =", sum(param_means(1:5))
print *, "sigma param checksum =", param_means(8) + param_means(9) + param_means(10)
end program xsim_mcp_arsigma_file
