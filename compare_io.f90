module compare_io_mod
use kind_mod, only: dp
implicit none
private
public :: read_series_file, read_time_series_file, read_matrix_file, read_series_labels_file, read_true_bkps_file, get_data_file_arg

contains

subroutine get_data_file_arg(default_name, filename)
character(len=*), intent(in) :: default_name
character(len=*), intent(out) :: filename
integer :: status

call get_command_argument(1, filename, status=status)
if (status /= 0 .or. len_trim(filename) == 0) filename = default_name
end subroutine get_data_file_arg

subroutine read_series_file(filename, x)
character(len=*), intent(in) :: filename
real(kind=dp), allocatable, intent(out) :: x(:)

integer :: unit_in, ios, nrow, i
character(len=1024) :: line

open(newunit=unit_in, file=filename, status="old", action="read", iostat=ios)
if (ios /= 0) error stop "could not open input file"

nrow = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (is_data_line(line)) nrow = nrow + 1
end do

if (nrow <= 0) error stop "input file contains no numeric data"

allocate(x(nrow))
rewind(unit_in)

i = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (.not. is_data_line(line)) cycle
    i = i + 1
    read(line, *, iostat=ios) x(i)
    if (ios /= 0) error stop "could not parse numeric data line"
end do

close(unit_in)
end subroutine read_series_file

subroutine read_time_series_file(filename, t, x)
character(len=*), intent(in) :: filename
real(kind=dp), allocatable, intent(out) :: t(:), x(:)

integer :: unit_in, ios, nrow, i
character(len=1024) :: line

open(newunit=unit_in, file=filename, status="old", action="read", iostat=ios)
if (ios /= 0) error stop "could not open input file"

nrow = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (is_data_line(line)) nrow = nrow + 1
end do

if (nrow <= 0) error stop "input file contains no numeric data"

allocate(t(nrow), x(nrow))
rewind(unit_in)

i = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (.not. is_data_line(line)) cycle
    i = i + 1
    read(line, *, iostat=ios) t(i), x(i)
    if (ios /= 0) error stop "could not parse time-series data line"
end do

close(unit_in)
end subroutine read_time_series_file

subroutine read_matrix_file(filename, x)
character(len=*), intent(in) :: filename
real(kind=dp), allocatable, intent(out) :: x(:, :)

integer :: unit_in, ios, nrow, ncol, i
character(len=1024) :: line

open(newunit=unit_in, file=filename, status="old", action="read", iostat=ios)
if (ios /= 0) error stop "could not open input file"

nrow = 0
ncol = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (.not. is_data_line(line)) cycle
    nrow = nrow + 1
    if (ncol == 0) ncol = count_tokens(line)
end do

if (nrow <= 0 .or. ncol <= 0) error stop "input file contains no numeric matrix data"

allocate(x(nrow, ncol))
rewind(unit_in)

i = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (.not. is_data_line(line)) cycle
    i = i + 1
    read(line, *, iostat=ios) x(i, :)
    if (ios /= 0) error stop "could not parse matrix data line"
end do

close(unit_in)
end subroutine read_matrix_file

subroutine read_series_labels_file(filename, x, labels)
character(len=*), intent(in) :: filename
real(kind=dp), allocatable, intent(out) :: x(:)
integer, allocatable, intent(out) :: labels(:)

integer :: unit_in, ios, nrow, i
character(len=1024) :: line

open(newunit=unit_in, file=filename, status="old", action="read", iostat=ios)
if (ios /= 0) error stop "could not open input file"

nrow = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (is_data_line(line)) nrow = nrow + 1
end do

if (nrow <= 0) error stop "input file contains no numeric data"

allocate(x(nrow), labels(nrow))
rewind(unit_in)

i = 0
do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (.not. is_data_line(line)) cycle
    i = i + 1
    read(line, *, iostat=ios) x(i), labels(i)
    if (ios /= 0) error stop "could not parse labeled numeric data line"
end do

close(unit_in)
end subroutine read_series_labels_file

subroutine read_true_bkps_file(filename, truth)
character(len=*), intent(in) :: filename
integer, allocatable, intent(out) :: truth(:)

integer :: unit_in, ios
character(len=1024) :: line

open(newunit=unit_in, file=filename, status="old", action="read", iostat=ios)
if (ios /= 0) error stop "could not open input file"

do
    read(unit_in, "(A)", iostat=ios) line
    if (ios /= 0) exit
    if (index(adjustl(line), "# true_bkps =") == 1) then
        call parse_true_bkps(line, truth)
        exit
    end if
end do

close(unit_in)
if (.not. allocated(truth)) error stop "input file is missing '# true_bkps =' metadata"
end subroutine read_true_bkps_file

subroutine parse_true_bkps(line, truth)
character(len=*), intent(in) :: line
integer, allocatable, intent(out) :: truth(:)
character(len=1024) :: payload
integer :: i, count_vals

payload = adjustl(line(index(line, "=") + 1:))
count_vals = 0
do i = 1, len_trim(payload)
    if (payload(i:i) /= " " .and. (i == 1 .or. payload(i-1:i-1) == " ")) count_vals = count_vals + 1
end do
allocate(truth(count_vals))
read(payload, *) truth
end subroutine parse_true_bkps

logical function is_data_line(line)
character(len=*), intent(in) :: line
character(len=:), allocatable :: trimmed

trimmed = adjustl(trim(line))
is_data_line = len(trimmed) > 0
if (.not. is_data_line) return
if (trimmed(1:1) == "#") is_data_line = .false.
end function is_data_line

integer function count_tokens(line)
character(len=*), intent(in) :: line
integer :: i
character(len=:), allocatable :: trimmed

trimmed = adjustl(trim(line))
count_tokens = 0
do i = 1, len_trim(trimmed)
    if (trimmed(i:i) /= " " .and. (i == 1 .or. trimmed(i-1:i-1) == " ")) count_tokens = count_tokens + 1
end do
end function count_tokens

end module compare_io_mod
