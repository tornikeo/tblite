
module tblite_debug
  use mctc_env, only : wp
  implicit none
  public :: print_cpp_array_3d, print_cpp_array_2d, print_cpp_array_1d

contains

subroutine print_cpp_array_4d(name, A)
  implicit none
  character(len=*), intent(in) :: name
  integer :: nx, ny, nz, nw
  integer :: i, j, k, l
  real(wp), intent(in) :: A(:,:,:,:)
  
  nx = size(A,1)
  ny = size(A,2)
  nz = size(A,3)
  nw = size(A,4)
  
  write(*, '(A)', advance="no") name
  write(*, '(A)', advance="no") "["
  write(*, '(I0)', advance="no") nx
  write(*, '(A)', advance="no") "]["
  write(*, '(I0)', advance="no") ny
  write(*, '(A)', advance="no") "]["
  write(*, '(I0)', advance="no") nz
  write(*, '(A)', advance="no") "]["
  write(*, '(I0)', advance="no") nw
  write(*, '(A)', advance="no") "] = {"
  
  do i = 1, nx
     write(*, '(A)', advance="no") " {"
     do j = 1, ny
        write(*, '(A)', advance="no") "  {"
        do k = 1, nz
           write(*, '(A)', advance="no") "   {"
           do l = 1, nw
              write(*, '(F13.8)', advance="no") A(i,j,k,l)
              if (l < nw) then
                 write(*, '(A)', advance="no") ", "
              endif
           end do
           write(*, '(A)', advance="no") "   }"
           if (k < nz) then
              write(*, '(A)', advance="no") ", "
           endif
        end do
        write(*, '(A)', advance="no") "  }"
        if (j < ny) then
           write(*, '(A)', advance="no") ", "
        endif
     end do
     write(*, '(A)', advance="no") " }"
     if (i < nx) then
        write(*, '(A)', advance="no") ", "
     endif
  end do
  write(*, '(A)') "};"
  
end subroutine print_cpp_array_4d

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
        write(*, '(F20.15)', advance="no") A(i,j)
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