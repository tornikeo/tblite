
module tblite_debug
  use mctc_env, only : wp
  implicit none
  public :: print_cpp_array_3d, print_cpp_array_2d, print_cpp_array_1d

contains
subroutine print_cpp_array_3d(name, A)
  implicit none
  character(len=*) :: name
  integer :: nx, ny, nz
  integer :: i, j, k
  real(wp), intent(in) :: A(:,:,:)
  nx = size(A,1)
  ny = size(A,2)
  nz = size(A,3)
  
  write(*, '(A)', advance="no") name
  write(*, '(A)', advance="no") "["
  write(*, '(I0)', advance="no") nx
  write(*, '(A)', advance="no") "]["
  write(*, '(I0)', advance="no") ny
  write(*, '(A)', advance="no") "]["
  write(*, '(I0)', advance="no") nz
  write(*, '(A)', advance="no") "]"
  write(*, '(A)', advance="no") "   {"  ! open the outermost brace
  do i = 1, nx
     write(*, '(A)', advance="no") "  {"
     do j = 1, ny
        write(*, '(A)', advance="no") "    {"
        do k = 1, nz
           write(*, '(F13.8)', advance="no") A(i, j, k)
           if (k < nz) then
              write(*, '(A)', advance="no") ", "
           endif
        end do
        write(*, '(A)', advance="no") "}"
        if (j < ny) then
           write(*, '(A)', advance="no") ","
          !  write(*, '(A)') ""  ! newline
        endif
     end do
     write(*, '(A)', advance="no") "  }"
     if (i < nx) then
        write(*, '(A)', advance="no") ","
     else
        ! write(*, '(A)') ""
     endif
  end do
  write(*, '(A)') "};"  ! close the outermost brace and add semicolon
  print*, ""
end subroutine print_cpp_array_3d

subroutine print_cpp_array_2d(name, A)
  implicit none
  character(len=*), intent(in) :: name
  integer :: nx, ny, i, j
  real(wp), intent(in) :: A(:, :)
  nx = size(A,1)
  ny = size(A,2)
  
  write(*, '(A)', advance="no") name
  write(*, '(A)', advance="no") "["
  write(*, '(I0)', advance="no") nx
  write(*, '(A)', advance="no") "]["
  write(*, '(I0)', advance="no") ny
  write(*, '(A)', advance="no") "] = {"
  do i = 1, nx
     write(*, '(A)', advance="no") " {"
     do j = 1, ny
        write(*, '(F13.8)', advance="no") A(i,j)
        if (j < ny) then
           write(*, '(A)', advance="no") ", "
        endif
     end do
     write(*, '(A)', advance="no") "}"
     if (i < nx) then
        write(*, '(A)', advance="no") ", "
     endif
  end do
  write(*, '(A)') "};"
end subroutine print_cpp_array_2d

subroutine print_cpp_array_1d(name, A)
  implicit none
  character(len=*), intent(in) :: name
  integer :: n, i
  real(wp), intent(in) :: A(:)
  n = size(A)
  
  write(*, '(A)', advance="no") name
  write(*, '(A)', advance="no") "["
  write(*, '(I0)', advance="no") n
  write(*, '(A)', advance="no") "] = {"
  do i = 1, n
     write(*, '(F13.8)', advance="no") A(i)
     if (i < n) then
        write(*, '(A)', advance="no") ", "
     endif
  end do
  write(*, '(A)') "};"
end subroutine print_cpp_array_1d


end module tblite_debug