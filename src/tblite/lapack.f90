! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> Proxy module to reexport high-level linear algebra package wrappers
module tblite_lapack
   use tblite_lapack_getrf, only : getrf => wrap_getrf
   use tblite_lapack_getri, only : getri => wrap_getri
   use tblite_lapack_getrs, only : getrs => wrap_getrs
   implicit none
   public
end module tblite_lapack
