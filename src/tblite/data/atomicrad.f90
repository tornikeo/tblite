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

!> Atomic Radii of the Elements
!>
!> M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
!> in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!> edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!> corrected Nov. 17, 2010 for the 92nd edition.
module tblite_data_atomicrad
   use mctc_env, only : wp
   use mctc_io_convert, only : aatoau
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_atomic_rad, atomic_rad


   !> Get atomic radius for a species
   interface get_atomic_rad
      module procedure :: get_atomic_rad_symbol
      module procedure :: get_atomic_rad_number
   end interface get_atomic_rad


   integer, parameter :: max_elem = 118

   !> Atomic radii
   real(wp), parameter :: atomic_rad(max_elem) = aatoau * [ &
      & 0.32_wp, 0.37_wp, 1.30_wp, 0.99_wp, 0.84_wp, 0.75_wp, 0.71_wp, 0.64_wp, &
      & 0.60_wp, 0.62_wp, 1.60_wp, 1.40_wp, 1.24_wp, 1.14_wp, 1.09_wp, 1.04_wp, &
      & 1.00_wp, 1.01_wp, 2.00_wp, 1.74_wp, 1.59_wp, 1.48_wp, 1.44_wp, 1.30_wp, &
      & 1.29_wp, 1.24_wp, 1.18_wp, 1.17_wp, 1.22_wp, 1.20_wp, 1.23_wp, 1.20_wp, &
      & 1.20_wp, 1.18_wp, 1.17_wp, 1.16_wp, 2.15_wp, 1.90_wp, 1.76_wp, 1.64_wp, &
      & 1.56_wp, 1.46_wp, 1.38_wp, 1.36_wp, 1.34_wp, 1.30_wp, 1.36_wp, 1.40_wp, &
      & 1.42_wp, 1.40_wp, 1.40_wp, 1.37_wp, 1.36_wp, 1.36_wp, 2.38_wp, 2.06_wp, &
      & 1.94_wp, 1.84_wp, 1.90_wp, 1.88_wp, 1.86_wp, 1.85_wp, 1.83_wp, 1.82_wp, &
      & 1.81_wp, 1.80_wp, 1.79_wp, 1.77_wp, 1.77_wp, 1.78_wp, 1.74_wp, 1.64_wp, &
      & 1.58_wp, 1.50_wp, 1.41_wp, 1.36_wp, 1.32_wp, 1.30_wp, 1.30_wp, 1.32_wp, &
      & 1.44_wp, 1.45_wp, 1.50_wp, 1.42_wp, 1.48_wp, 1.46_wp, 2.42_wp, 2.11_wp, &
      & 2.01_wp, 1.90_wp, 1.84_wp, 1.83_wp, 1.80_wp, 1.80_wp, 1.73_wp, 1.68_wp, &
      & 1.68_wp, 1.68_wp, 1.65_wp, 1.67_wp, 1.73_wp, 1.76_wp, 1.61_wp, 1.57_wp, &
      & 1.49_wp, 1.43_wp, 1.41_wp, 1.34_wp, 1.29_wp, 1.28_wp, 1.21_wp, 1.22_wp, &
      & 1.36_wp, 1.43_wp, 1.62_wp, 1.75_wp, 1.65_wp, 1.57_wp]


contains


!> Get atomic radius for species with a given symbol
elemental function get_atomic_rad_symbol(symbol) result(radius)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> atomic radius
   real(wp) :: radius

   radius = get_atomic_rad(to_number(symbol))

end function get_atomic_rad_symbol


!> Get atomic radius for species with a given atomic number
elemental function get_atomic_rad_number(number) result(radius)

   !> Atomic number
   integer, intent(in) :: number

   !> atomic radius
   real(wp) :: radius

   if (number > 0 .and. number <= size(atomic_rad, dim=1)) then
      radius = atomic_rad(number)
   else
      radius = -1.0_wp
   end if

end function get_atomic_rad_number


end module tblite_data_atomicrad
