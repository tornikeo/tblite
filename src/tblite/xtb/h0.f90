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

!> Implementation of the effective core Hamiltonian used in the extended tight binding.
module tblite_xtb_h0
   use, intrinsic :: iso_fortran_env
   use iso_c_binding
   use iso_c_binding, only : c_int, c_double
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list
   use tblite_basis_type, only : basis_type, cgto_type
   use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto, maxl, msao
   use tblite_scf_potential, only : potential_type
   use tblite_xtb_spec, only : tb_h0spec
   implicit none

   interface
      subroutine get_vec_(xyz_iat, xyz_jat, trans, vec) bind(C, name="get_vec_")
         use mctc_env, only : wp
         implicit none
         real(wp), intent(in) :: xyz_iat(3)  ! (3, n_atoms)
         real(wp), intent(in) :: xyz_jat(3)  ! (3, n_atoms)
         real(wp), intent(in) :: trans(3) ! (3, n_trans)
         real(wp), intent(out) :: vec(3)
      end subroutine get_vec_
   end interface

   interface 
      !> Notice: Dimensions are passed in C order. 
      !> On Fortran side, dimensions are passed in reverse (arr, size(arr, 2), size(arr, 1)).
      !> On C side, dimension are correct, without change. No transpose is done.
      subroutine cuda_get_hamiltonian_kernel( nao, nelem, &
        mol_nat, mol_nid, mol_nbd, & !> structure_type
        mol_id, mol_id_dim1, &
        mol_num, mol_num_dim1, &
        mol_xyz, mol_xyz_dim1, mol_xyz_dim2, &
        mol_uhf, &
        mol_charge, &
        mol_lattice, mol_lattice_dim1, mol_lattice_dim2, &
        mol_periodic, mol_periodic_dim1, &
        mol_bond, mol_bond_dim1, mol_bond_dim2, &
        trans, trans_dim1, trans_dim2, & !> trans
        alist_inl, alist_inl_dim1, & !> adjacency_list
        alist_nnl, alist_nnl_dim1, &
        alist_nlat, alist_nlat_dim1, &
        alist_nltr, alist_nltr_dim1, &
        bas_maxl, bas_nsh, bas_nao, bas_intcut, bas_min_alpha, & !> basis_type
        bas_nsh_id, bas_nsh_id_dim1, &
        bas_nsh_at, bas_nsh_at_dim1, &
        bas_nao_sh, bas_nao_sh_dim1, &
        bas_iao_sh, bas_iao_sh_dim1, &
        bas_ish_at, bas_ish_at_dim1, &
        bas_ao2at, bas_ao2at_dim1, &
        bas_ao2sh, bas_ao2sh_dim1, &
        bas_sh2at, bas_sh2at_dim1, &
        cgto, cgto_dim1, cgto_dim2, & 
        h0_selfenergy, h0_selfenergy_dim1, h0_selfenergy_dim2, & !> tb_hamiltonian
        h0_kcn, h0_kcn_dim1, h0_kcn_dim2, &
        h0_kq1, h0_kq1_dim1, h0_kq1_dim2, &
        h0_kq2, h0_kq2_dim1, h0_kq2_dim2, &
        h0_hscale, h0_hscale_dim1, h0_hscale_dim2, h0_hscale_dim3, h0_hscale_dim4, &
        h0_shpoly, h0_shpoly_dim1, h0_shpoly_dim2, &
        h0_rad, h0_rad_dim1, &
        h0_refocc, h0_refocc_dim1, h0_refocc_dim2, & !> other vars
        selfenergy, overlap, dpint, qpint, hamiltonian &
      ) bind(C, name="cuda_get_hamiltonian_kernel_")
        use iso_c_binding
        use tblite_basis_type, only : basis_type, cgto_type
        use tblite_adjlist, only : adjacency_list

        implicit none
        integer(c_int), value :: nao, nelem

        !> structure_type
        integer(c_int), value :: mol_nat
        integer(c_int), value :: mol_nid
        integer(c_int), value :: mol_nbd
        integer(c_int), intent(in) :: mol_id(*)
        integer(c_int), value :: mol_id_dim1
        integer(c_int), intent(in) :: mol_num(*)
        integer(c_int), value :: mol_num_dim1
        real(c_double), intent(in) :: mol_xyz(*)
        integer(c_int), value :: mol_xyz_dim1, mol_xyz_dim2
        integer(c_int), value :: mol_uhf
        real(c_double), value :: mol_charge
        real(c_double), intent(in) :: mol_lattice(*)
        integer(c_int), value :: mol_lattice_dim1, mol_lattice_dim2
        integer(c_int), intent(in) :: mol_periodic(*)
        integer(c_int), value :: mol_periodic_dim1
        integer(c_int), intent(in) :: mol_bond(*)
        integer(c_int), value :: mol_bond_dim1, mol_bond_dim2

        !> trans
        real(c_double), intent(in) :: trans(*)
        integer(c_int), value :: trans_dim1
        integer(c_int), value :: trans_dim2

        !> adjacency_list
        integer(c_int), intent(in) :: alist_inl(*)
        integer(c_int), value :: alist_inl_dim1
        integer(c_int), intent(in) :: alist_nnl(*)
        integer(c_int), value :: alist_nnl_dim1
        integer(c_int), intent(in) :: alist_nlat(*)
        integer(c_int), value :: alist_nlat_dim1
        integer(c_int), intent(in) :: alist_nltr(*)
        integer(c_int), value :: alist_nltr_dim1

        !> basis_type
        integer(c_int), value :: bas_maxl
        integer(c_int), value :: bas_nsh
        integer(c_int), value :: bas_nao
        real(c_double), value :: bas_intcut
        real(c_double), value :: bas_min_alpha
        integer(c_int), intent(in) :: bas_nsh_id(*)
        integer(c_int), value :: bas_nsh_id_dim1
        integer(c_int), intent(in) :: bas_nsh_at(*)
        integer(c_int), value :: bas_nsh_at_dim1
        integer(c_int), intent(in) :: bas_nao_sh(*)
        integer(c_int), value :: bas_nao_sh_dim1
        integer(c_int), intent(in) :: bas_iao_sh(*)
        integer(c_int), value :: bas_iao_sh_dim1
        integer(c_int), intent(in) :: bas_ish_at(*)
        integer(c_int), value :: bas_ish_at_dim1
        integer(c_int), intent(in) :: bas_ao2at(*)
        integer(c_int), value :: bas_ao2at_dim1
        integer(c_int), intent(in) :: bas_ao2sh(*)
        integer(c_int), value :: bas_ao2sh_dim1
        integer(c_int), intent(in) :: bas_sh2at(*)
        integer(c_int), value :: bas_sh2at_dim1
        type(cgto_type), intent(in) :: cgto(*)
        integer(c_int), value :: cgto_dim1, cgto_dim2

        !> tb_hamiltonian
        real(c_double), intent(in) :: h0_selfenergy(*)
        integer(c_int), value :: h0_selfenergy_dim1, h0_selfenergy_dim2
        real(c_double), intent(in) :: h0_kcn(*)
        integer(c_int), value :: h0_kcn_dim1, h0_kcn_dim2
        real(c_double), intent(in) :: h0_kq1(*)
        integer(c_int), value :: h0_kq1_dim1, h0_kq1_dim2
        real(c_double), intent(in) :: h0_kq2(*)
        integer(c_int), value :: h0_kq2_dim1, h0_kq2_dim2
        real(c_double), intent(in) :: h0_hscale(*)
        integer(c_int), value :: h0_hscale_dim1, h0_hscale_dim2, h0_hscale_dim3, h0_hscale_dim4
        real(c_double), intent(in) :: h0_shpoly(*)
        integer(c_int), value :: h0_shpoly_dim1, h0_shpoly_dim2
        real(c_double), intent(in) :: h0_rad(*)
        integer(c_int), value :: h0_rad_dim1
        real(c_double), intent(in) :: h0_refocc(*)
        integer(c_int), value :: h0_refocc_dim1, h0_refocc_dim2

        real(c_double), intent(in) :: selfenergy(*) !(:)
        real(c_double), intent(out) :: overlap(*) !(:, :)
        real(c_double), intent(out) :: dpint(*) !(:, :, :)
        real(c_double), intent(out) :: qpint(*) !(:, :, :)
        real(c_double), intent(out) :: hamiltonian(*) !(:,:)
      end subroutine
   end interface

   private
   public :: tb_hamiltonian, new_hamiltonian
   public :: get_selfenergy, get_hamiltonian, cuda_get_hamiltonian, get_occupation, get_hamiltonian_gradient

   type :: tb_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: selfenergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq1(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq2(:, :)
      !> Enhancement factor to scale the Hamiltonian elements
      real(wp), allocatable :: hscale(:, :, :, :)
      !> Polynomial coefficients for distance dependent enhancement factor
      real(wp), allocatable :: shpoly(:, :)
      !> Atomic radius for polynomial enhancement
      real(wp), allocatable :: rad(:)
      !> Reference occupation numbers
      real(wp), allocatable :: refocc(:, :)
   end type tb_hamiltonian
      
contains


!> Constructor for a new Hamiltonian object, consumes a Hamiltonian specification
   subroutine new_hamiltonian(self, mol, bas, spec)
      type(tb_hamiltonian), intent(out) :: self
      type(structure_type), intent(in) :: mol
      type(basis_type), intent(in) :: bas
      class(tb_h0spec), intent(in) :: spec

      integer :: mshell

      mshell = maxval(bas%nsh_id)
      allocate(self%selfenergy(mshell, mol%nid), self%kcn(mshell, mol%nid), &
      & self%kq1(mshell, mol%nid), self%kq2(mshell, mol%nid))
      call spec%get_selfenergy(mol, bas, self%selfenergy)
      call spec%get_cnshift(mol, bas, self%kcn)
      call spec%get_q1shift(mol, bas, self%kq1)
      call spec%get_q2shift(mol, bas, self%kq2)

      allocate(self%hscale(mshell, mshell, mol%nid, mol%nid))
      call spec%get_hscale(mol, bas, self%hscale)

      allocate(self%shpoly(mshell, mol%nid), self%rad(mol%nid))
      call spec%get_rad(mol, bas, self%rad)
      call spec%get_shpoly(mol, bas, self%shpoly)

      allocate(self%refocc(mshell, mol%nid))
      call spec%get_reference_occ(mol, bas, self%refocc)
   end subroutine new_hamiltonian


   subroutine get_selfenergy(h0, id, ish_at, nshell, cn, qat, selfenergy, dsedcn, dsedq)
      type(tb_hamiltonian), intent(in) :: h0
      integer, intent(in) :: id(:)
      integer, intent(in) :: ish_at(:)
      integer, intent(in) :: nshell(:)
      real(wp), intent(in), optional :: cn(:)
      real(wp), intent(in), optional :: qat(:)
      real(wp), intent(out) :: selfenergy(:)
      real(wp), intent(out), optional :: dsedcn(:)
      real(wp), intent(out), optional :: dsedq(:)

      integer :: iat, izp, ish, ii

      selfenergy(:) = 0.0_wp
      if (present(dsedcn)) dsedcn(:) = 0.0_wp
      if (present(dsedq)) dsedq(:) = 0.0_wp
      do iat = 1, size(id)
         izp = id(iat)
         ii = ish_at(iat)
         do ish = 1, nshell(izp)
            selfenergy(ii+ish) = h0%selfenergy(ish, izp)
         end do
      end do
      if (present(cn)) then
         if (present(dsedcn)) then
            do iat = 1, size(id)
               izp = id(iat)
               ii = ish_at(iat)
               do ish = 1, nshell(izp)
                  selfenergy(ii+ish) = selfenergy(ii+ish) - h0%kcn(ish, izp) * cn(iat)
                  dsedcn(ii+ish) = -h0%kcn(ish, izp)
               end do
            end do
         else
            do iat = 1, size(id)
               izp = id(iat)
               ii = ish_at(iat)
               do ish = 1, nshell(izp)
                  selfenergy(ii+ish) = selfenergy(ii+ish) - h0%kcn(ish, izp) * cn(iat)
               end do
            end do
         end if
      end if
      if (present(qat)) then
         if (present(dsedq)) then
            do iat = 1, size(id)
               izp = id(iat)
               ii = ish_at(iat)
               do ish = 1, nshell(izp)
                  selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
                  dsedq(ii+ish) = -h0%kq1(ish, izp) - h0%kq2(ish, izp)*2*qat(iat)
               end do
            end do
         else
            do iat = 1, size(id)
               izp = id(iat)
               ii = ish_at(iat)
               do ish = 1, nshell(izp)
                  selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
               end do
            end do
         end if
      end if

   end subroutine get_selfenergy

   subroutine print_cgtos(bas)
    type(basis_type), intent(in) :: bas
    integer :: i, j, k
    do i = 1, size(bas%cgto, 1)
      do j = 1, size(bas%cgto, 2)
        ! write in C style, with %f format (limit number of digits )
        print*, " "
        write(*, "(A, I1, A, I1, A)") "cgto(", i, ",", j, ")"
        write(*, "(A, I3)") "  ang: ", bas%cgto(i, j)%ang
        write(*, "(A, I3)") "  nprim: ", bas%cgto(i, j)%nprim
        write(*, "(A, I3)") "  alpha: "
        do k = 1, bas%cgto(i, j)%nprim
          write(*, "(A, F12.8)", advance="no") " ", bas%cgto(i, j)%alpha(k)
        end do
        print*, " "
        write(*, "(A, I3)") "  coeff: "
        do k = 1, bas%cgto(i, j)%nprim
          write(*, "(A, F12.8)", advance="no") " ", bas%cgto(i, j)%coeff(k)
        end do
        print*, " "
      end do
    end do
   end subroutine print_cgtos

   subroutine print_adjlist(alist)
      type(adjacency_list), intent(in) :: alist
      integer :: i, j
      write(*, "(A)", advance="no") "inl: "
      do i = 1, size(alist%inl, 1)
        write(*, "(I3)", advance="no") alist%inl(i)
      end do
      print*,""  ! Add a newline after printing the inl values
      write(*, "(A)", advance="no") "nnl: "
      do i = 1, size(alist%nnl, 1)
        write(*, "(I3)", advance="no") alist%nnl(i)
      end do
      print*,""  ! Add a newline after printing the nnl values
      write(*, "(A)", advance="no") "nlat: "
      do i = 1, size(alist%nlat, 1)
        write(*, "(I3)", advance="no") alist%nlat(i)
      end do
      print*,""  ! Add a newline after printing the nlat values
      write(*, "(A)", advance="no") "nltr: "
      do i = 1, size(alist%nltr, 1)
        write(*, "(I3)", advance="no") alist%nltr(i)
      end do
      print*,""  ! Add a newline after printing the nltr values
   end subroutine print_adjlist

  subroutine print_tb_hamiltonian(h0)
    type(tb_hamiltonian), intent(in) :: h0
    integer :: i, j, k, l, m
    write(*, "(A)", advance="no") "selfenergy: "
    do i = 1, size(h0%selfenergy, 1)
      do j = 1, size(h0%selfenergy, 2)
        write(*, "(F12.8)", advance="no") h0%selfenergy(i, j)
      end do
    end do
    print*,""  ! Add a newline after printing the selfenergy values
    write(*, "(A)", advance="no") "kcn: "
    do i = 1, size(h0%kcn, 1)
      do j = 1, size(h0%kcn, 2)
        write(*, "(F12.8)", advance="no") h0%kcn(i, j)
      end do
    end do
    print*,""  ! Add a newline after printing the kcn values
    write(*, "(A)", advance="no") "kq1: "
    do i = 1, size(h0%kq1, 1)
      do j = 1, size(h0%kq1, 2)
        write(*, "(F12.8)", advance="no") h0%kq1(i, j)
      end do
    end do
    print*,""  ! Add a newline after printing the kq1 values
    write(*, "(A)", advance="no") "kq2: "
    do i = 1, size(h0%kq2, 1)
      do j = 1, size(h0%kq2, 2)
        write(*, "(F12.8)", advance="no") h0%kq2(i, j)
      end do
    end do
    print*,""  ! Add a newline after printing the kq2 values
    write(*, "(A)", advance="no") "hscale: "
    do i = 1, size(h0%hscale, 1)
      do j = 1, size(h0%hscale, 2)
        do k = 1, size(h0%hscale, 3)
          do l = 1, size(h0%hscale, 4)
            write(*, "(F12.8)", advance="no") h0%hscale(i, j, k, l)
          end do
        end do
      end do
    end do
    print*,""  ! Add a newline after printing the hscale values
    write(*, "(A)", advance="no") "shpoly: "
    do i = 1, size(h0%shpoly, 1)
      do j = 1, size(h0%shpoly, 2)
        write(*, "(F12.8)", advance="no") h0%shpoly(i, j)
      end do
    end do
    print*,""  ! Add a newline after printing the shpoly values
    write(*, "(A)", advance="no") "rad: "
    do i = 1, size(h0%rad, 1)
      write(*, "(F12.8)", advance="no") h0%rad(i)
    end do
    print*,""  ! Add a newline after printing the rad values
    write(*, "(A)", advance="no") "refocc: "
    do i = 1, size(h0%refocc, 1)
      do j = 1, size(h0%refocc, 2)
        write(*, "(F12.8)", advance="no") h0%refocc(i, j)
      end do
    end do
    print*,""  ! Add a newline after printing the refocc values
  end subroutine print_tb_hamiltonian

  subroutine print_basis_type(bas)
    type(basis_type), intent(in) :: bas
    integer :: i, j, k
    write(*, "(A, I3, A)") "maxl: ", bas%maxl, " "
    write(*, "(A, I3, A)") "nsh: ", bas%nsh, " "
    write(*, "(A, I3, A)") "nao: ", bas%nao, " "
    write(*, "(A, F12.8, A)") "intcut: ", bas%intcut, " "
    write(*, "(A, F12.8, A)") "min_alpha: ", bas%min_alpha, " "

    write(*, "(A)", advance="no") "nsh_id: "
    do i = 1, size(bas%nsh_id, 1)
      write(*, "(I3)", advance="no") bas%nsh_id(i)
    end do
    print*,""  ! Add a newline after printing the nsh_id values
    write(*, "(A)", advance="no") "nsh_at: "
    do i = 1, size(bas%nsh_at, 1)
      write(*, "(I3)", advance="no") bas%nsh_at(i)
    end do
    print*,""  ! Add a newline after printing the nsh_at values
    write(*, "(A)", advance="no") "nao_sh: "
    do i = 1, size(bas%nao_sh, 1)
      write(*, "(I3)", advance="no") bas%nao_sh(i)
    end do
    print*,""  ! Add a newline after printing the nao_sh values
    write(*, "(A)", advance="no") "iao_sh: "
    do i = 1, size(bas%iao_sh, 1)
      write(*, "(I3)", advance="no") bas%iao_sh(i)
    end do
    print*,""  ! Add a newline after printing the iao_sh values
    write(*, "(A)", advance="no") "ish_at: "
    do i = 1, size(bas%ish_at, 1)
      write(*, "(I3)", advance="no") bas%ish_at(i)
    end do
    print*,""  ! Add a newline after printing the ish_at values
    write(*, "(A)", advance="no") "ao2at: "
    do i = 1, size(bas%ao2at, 1)
      write(*, "(I3)", advance="no") bas%ao2at(i)
    end do
    print*,""  ! Add a newline after printing the ao2at values
    write(*, "(A)", advance="no") "ao2sh: "
    do i = 1, size(bas%ao2sh, 1)
      write(*, "(I3)", advance="no") bas%ao2sh(i)
    end do
    print*,""  ! Add a newline after printing the ao2sh values
    write(*, "(A)", advance="no") "sh2at: "
    do i = 1, size(bas%sh2at, 1)
      write(*, "(I3)", advance="no") bas%sh2at(i)
    end do
    print*,""  ! Add a newline after printing the sh2at values
  end subroutine print_basis_type

  subroutine cuda_get_hamiltonian(mol, trans, alist, bas, h0, selfenergy, overlap, dpint, qpint, &
  & hamiltonian)
    use iso_c_binding
    implicit none
    !> Molecular structure data
    type(structure_type), intent(in) :: mol
    !> Lattice points within a given realspace cutoff
    real(wp), intent(in) :: trans(:, :)
    !> Neighbour list
    type(adjacency_list), intent(in) :: alist
    !> Basis set information
    type(basis_type), intent(in) :: bas
    !> Hamiltonian interaction data
    type(tb_hamiltonian), intent(in) :: h0
    !> Diagonal elememts of the Hamiltonian
    real(wp), intent(in) :: selfenergy(:)
    !> Overlap integral matrix
    real(wp), intent(out) :: overlap(:, :)
    !> Dipole moment integral matrix
    real(wp), intent(out) :: dpint(:, :, :)
    !> Quadrupole moment integral matrix
    real(wp), intent(out) :: qpint(:, :, :)
    !> Effective Hamiltonian
    real(wp), intent(out) :: hamiltonian(:, :)

    integer(kind=c_int) :: nao 
    integer(kind=c_int) :: nelem
    integer(c_int), allocatable :: periodic_as_integers(:)
    integer :: i, j, k
    nao = size(hamiltonian, 1)
    nelem = size(selfenergy, 1)
    !> transform logical periodic to integer

    !> 0 = false, 1 = true
    allocate(periodic_as_integers(size(mol%periodic, 1)))
    periodic_as_integers = 0
    do i = 1, size(mol%periodic, 1)
      if (mol%periodic(i)) periodic_as_integers(i) = 1
    end do

    !  overlap = 1
    !  dpint = 2
    !  qpint = 3
    !  hamiltonian = 4
    ! since we don't support lattices yet, if trans is not zero, error
    ! if (any(trans /= 0.0_wp)) then
    !    print*, "Error: Non-zero translation vector provided."
    !    stop
    ! end if
    ! if (any(mol%periodic)) then
    !    print*, "Error: Periodic boundary conditions not supported in CUDA yet."
    !    stop
    ! end if

    !> Print all cgtos in order
    ! call print_cgtos(bas)
    print*, "================= FORTRAN ================="
    ! call print_adjlist(alist)
    ! call print_tb_hamiltonian(h0)
    call print_basis_type(bas)
    
    call cuda_get_hamiltonian_kernel( nao, nelem, &
      !> structure_type
      mol%nat,&
      mol%nid, &
      mol%nbd, &
      mol%id, size(mol%id, 1), &
      mol%num, size(mol%num, 1), &
      mol%xyz, size(mol%xyz, 2), size(mol%xyz, 1), &
      mol%uhf, &
      mol%charge, &
      mol%lattice, size(mol%lattice, 2), size(mol%lattice, 1), &
      periodic_as_integers, size(periodic_as_integers, 1), &
      mol%bond, size(mol%bond, 2), size(mol%bond, 1), &
      !> trans
      trans, size(trans, 2), size(trans, 1), &
      !> adjacency_list
      alist%inl, size(alist%inl, 1), &
      alist%nnl, size(alist%nnl, 1), &
      alist%nlat, size(alist%nlat, 1), &
      alist%nltr, size(alist%nltr, 1), &
        !> basis_type
      bas%maxl, bas%nsh, bas%nao, bas%intcut, bas%min_alpha, &
      bas%nsh_id, size(bas%nsh_id, 1), &
      bas%nsh_at, size(bas%nsh_at, 1), &
      bas%nao_sh, size(bas%nao_sh, 1), &
      bas%iao_sh, size(bas%iao_sh, 1), &
      bas%ish_at, size(bas%ish_at, 1), &
      bas%ao2at, size(bas%ao2at, 1), &
      bas%ao2sh, size(bas%ao2sh, 1), &
      bas%sh2at, size(bas%sh2at, 1), &
      bas%cgto, size(bas%cgto, 2), size(bas%cgto, 1), &
      !> tb_hamiltonian
      h0%selfenergy, size(h0%selfenergy,2), size(h0%selfenergy,1), &
      h0%kcn, size(h0%kcn,2), size(h0%kcn,1), &
      h0%kq1, size(h0%kq1,2), size(h0%kq1,1), &
      h0%kq2, size(h0%kq2,2), size(h0%kq2,1), &
      h0%hscale, size(h0%hscale,4), size(h0%hscale,3), size(h0%hscale,2), size(h0%hscale,1), &
      h0%shpoly, size(h0%shpoly,2), size(h0%shpoly,1), &
      h0%rad, size(h0%rad,1), &
      h0%refocc, size(h0%refocc,2), size(h0%refocc,1), &
      !> selfenergy, overlap, dpint, qpint, hamiltonian
      selfenergy, overlap, dpint, qpint, hamiltonian)

    print*,"";
  end subroutine cuda_get_hamiltonian

   subroutine get_hamiltonian(mol, trans, alist, bas, h0, selfenergy, overlap, dpint, qpint, &
   & hamiltonian)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: trans(:, :)
      !> Neighbour list
      type(adjacency_list), intent(in) :: alist
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: h0
      !> Diagonal elememts of the Hamiltonian
      real(wp), intent(in) :: selfenergy(:)
      !> Overlap integral matrix
      real(wp), intent(out) :: overlap(:, :)
      !> Dipole moment integral matrix
      real(wp), intent(out) :: dpint(:, :, :)
      !> Quadrupole moment integral matrix
      real(wp), intent(out) :: qpint(:, :, :)
      !> Effective Hamiltonian
      real(wp), intent(out) :: hamiltonian(:, :)

      integer :: i,j,l;
      integer :: iat, jat, izp, jzp, itr, k, img, inl
      integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
      real(wp) :: rr, r2, vec(3), cutoff2, hij, shpoly, dtmpj(3), qtmpj(6)
      real(wp), allocatable :: stmp(:), dtmpi(:, :), qtmpi(:, :)

      overlap(:, :) = 0.0_wp
      dpint(:, :, :) = 0.0_wp
      qpint(:, :, :) = 0.0_wp
      hamiltonian(:, :) = 0.0_wp

      allocate(stmp(msao(bas%maxl)**2), dtmpi(3, msao(bas%maxl)**2), qtmpi(6, msao(bas%maxl)**2))
 
      ! $omp parallel do schedule(runtime) default(none) &
      ! $omp shared(mol, bas, trans, alist, overlap, dpint, qpint, hamiltonian, h0, selfenergy) &
      ! $omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, k) &
      ! $omp private(r2, vec, stmp, dtmpi, qtmpi, dtmpj, qtmpj, hij, shpoly, rr, inl, img)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         inl = alist%inl(iat)
         do img = 1, alist%nnl(iat)
            jat = alist%nlat(img+inl)
            itr = alist%nltr(img+inl)
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp, dtmpi, qtmpi)

                  shpoly = (1.0_wp + h0%shpoly(ish, izp)*rr) &
                     * (1.0_wp + h0%shpoly(jsh, jzp)*rr)

                  hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                     * h0%hscale(jsh, ish, jzp, izp) * shpoly

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        call shift_operator(vec, stmp(ij), dtmpi(:, ij), qtmpi(:, ij), &
                        & dtmpj, qtmpj)
                        ! $omp atomic
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           + stmp(ij)

                        do k = 1, 3
                           ! $omp atomic
                           dpint(k, jj+jao, ii+iao) = dpint(k, jj+jao, ii+iao) &
                              + dtmpi(k, ij)
                        end do

                        do k = 1, 6
                           ! $omp atomic
                           qpint(k, jj+jao, ii+iao) = qpint(k, jj+jao, ii+iao) &
                              + qtmpi(k, ij)
                        end do

                        ! $omp atomic
                        hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                           + stmp(ij) * hij

                        if (iat /= jat) then
                           ! $omp atomic
                           overlap(ii+iao, jj+jao) = overlap(ii+iao, jj+jao) &
                              + stmp(ij)
                           do k = 1, 3
                              ! $omp atomic
                              dpint(k, ii+iao, jj+jao) = dpint(k, ii+iao, jj+jao) &
                                 + dtmpj(k)
                           end do

                           do k = 1, 6
                              ! $omp atomic
                              qpint(k, ii+iao, jj+jao) = qpint(k, ii+iao, jj+jao) &
                                 + qtmpj(k)
                           end do
                           ! $omp atomic
                           hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                              + stmp(ij) * hij
                        end if
                     end do
                  end do

               end do
            end do

         end do
      end do

      ! $omp parallel do schedule(runtime) default(none) &
      ! $omp shared(mol, bas, trans, cutoff2, overlap, dpint, qpint, hamiltonian, h0, selfenergy) &
      ! $omp private(iat, izp, itr, is, ish, jsh, ii, jj, iao, jao, nao, ij) &
      ! $omp private(r2, vec, stmp, dtmpi, qtmpi, hij, shpoly, rr)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         vec(:) = 0.0_wp
         r2 = 0.0_wp
         rr = sqrt(sqrt(r2) / (h0%rad(izp) + h0%rad(izp)))
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is+ish)
            do jsh = 1, bas%nsh_id(izp)
               jj = bas%iao_sh(is+jsh)
               call multipole_cgto(bas%cgto(jsh, izp), bas%cgto(ish, izp), &
               & r2, vec, bas%intcut, stmp, dtmpi, qtmpi)

               shpoly = (1.0_wp + h0%shpoly(ish, izp)*rr) &
                  * (1.0_wp + h0%shpoly(jsh, izp)*rr)

               hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(is+jsh)) &
                  * shpoly

               nao = msao(bas%cgto(jsh, izp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)
                     overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                        + stmp(ij)

                     dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                        + dtmpi(:, ij)

                     qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                        + qtmpi(:, ij)

                     hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                        + stmp(ij) * hij
                  end do
               end do

            end do
         end do

      end do

   end subroutine get_hamiltonian

  subroutine print_4d_array_literal(arr)
    implicit none
    integer, parameter :: wd = kind(1.0d0)
    real(wd), intent(in) :: arr(:,:,:,:)
    integer :: i, j, k, l
    integer :: d1, d2, d3, d4
    character(len=1000) :: line
    character(len=32) :: num_str
  
    d1 = size(arr, 1)
    d2 = size(arr, 2)
    d3 = size(arr, 3)
    d4 = size(arr, 4)
  
    print *, "["
    do i = 1, d1
      print *, "  ["
      do j = 1, d2
        print *, "    ["
        do k = 1, d3
          line = "      ["
          do l = 1, d4
            write(num_str, '(F0.6)') arr(i,j,k,l)
            line = trim(line) // trim(adjustl(num_str))
            if (l < d4) line = trim(line) // ", "
          end do
          line = trim(line) // "]"
          if (k < d3) then
            print *, trim(line) // ","
          else
            print *, trim(line)
          end if
        end do
        if (j < d2) then
          print *, "    ],"
        else
          print *, "    ]"
        end if
      end do
      if (i < d1) then
        print *, "  ],"
      else
        print *, "  ]"
      end if
    end do
    print *, "]"
  end subroutine print_4d_array_literal

  subroutine print_int_matrix(mat)
    implicit none
    integer, intent(in) :: mat(:,:)
    integer :: i, j
    character(len=100) :: line
    character(len=16) :: int_str
  
    do i = 1, size(mat, 1)
      line = ""
      do j = 1, size(mat, 2)
        write(int_str, '(I0)') mat(i, j)
        line = trim(line) // " " // trim(adjustl(int_str))
      end do
      print *, trim(line)
    end do
  end subroutine print_int_matrix

  subroutine print_3d_array(arr)
    implicit none
    integer, parameter :: wd = kind(1.0d0)
    real(wd), intent(in) :: arr(:,:,:)
    integer :: i, j, k
    integer :: d1, d2, d3
    character(len=100) :: line
    character(len=32) :: num_str

    d1 = size(arr, 1)
    d2 = size(arr, 2)
    d3 = size(arr, 3)

    print *, "["
    do i = 1, d1
      print *, "  ["
      do j = 1, d2
        line = "    ["
        do k = 1, d3
          write(num_str, '(F0.6)') arr(i,j,k)
          line = trim(line) // trim(adjustl(num_str))
          if (k < d3) line = trim(line) // ", "
        end do
        line = trim(line) // "]"
        if (j < d2) then
          print *, trim(line) // ","
        else
          print *, trim(line)
        end if
      end do
      if (i < d1) then
        print *, "  ],"
      else
        print *, "  ]"
      end if
    end do
    print *, "]"
  end subroutine print_3d_array

  subroutine print_cgto_alpha(bas)
    integer :: i,j,k
    type(basis_type), intent(in) :: bas

    write(*,'(A)',advance='no') '['
    do i = 1,3
      write(*,'(A)',advance='no') '['
      do j = 1,4
        write(*,'(A)',advance='no') '['
        do k = 1,12
          write(*,'(F10.6,",")',advance='no') bas%cgto(i,j)%alpha(k)
        end do
        write(*,'(A)', advance='yes') '],'
      end do
      write(*,'(A)', advance='yes') '],'
    end do
    write(*,'(A)', advance='yes') ']'
  end subroutine print_cgto_alpha

  subroutine print_cgto_coeff(bas)
    integer :: i,j,k
    type(basis_type), intent(in) :: bas

    write(*,'(A)',advance='no') '['
    do i = 1,3
      write(*,'(A)',advance='no') '['
      do j = 1,4
        write(*,'(A)',advance='no') '['
        do k = 1,12
          write(*,'(F10.6,",")',advance='no') bas%cgto(i,j)%coeff(k)
        end do
        write(*,'(A)', advance='yes') '],'
      end do
      write(*,'(A)', advance='yes') '],'
    end do
    write(*,'(A)', advance='yes') ']'
  end subroutine print_cgto_coeff


  subroutine print_cgto_nprim(bas)
    integer :: i,j
    type(basis_type), intent(in) :: bas

    write(*,'(A)',advance='no') '['
    do i = 1,3
      write(*,'(A)',advance='no') '['
      do j = 1,4
        write(*,'(I4,",")',advance='no') bas%cgto(i,j)%nprim
      end do
      write(*,'(A)', advance='yes') '],'
    end do
    write(*,'(A)', advance='yes') ']'
  end subroutine print_cgto_nprim


  subroutine get_hamiltonian_gradient(mol, trans, list_, bas, h0, selfenergy, dsedcn, &
   & pot, pmat, xmat, dEdcn, gradient, sigma)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: trans(:, :)
      !> Neighbour list
      type(adjacency_list), intent(in) :: list_
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: h0
      !> Diagonal elememts of the Hamiltonian
      real(wp), intent(in) :: selfenergy(:)
      !> Derivative of the diagonal elements of the Hamiltonian w.r.t. the coordination number
      real(wp), intent(in) :: dsedcn(:)
      !> Density dependent potential shifts on the Hamiltonian
      type(potential_type), intent(in) :: pot
      !> Density matrix
      real(wp), intent(in) :: pmat(:, :, :)
      !> Energy weighted density matrix
      real(wp), intent(in) :: xmat(:, :, :)

      !> Derivative of the electronic energy w.r.t. the coordination number
      real(wp), intent(inout) :: dEdcn(:)
      !> Derivative of the electronic energy w.r.t. coordinate displacements
      real(wp), intent(inout) :: gradient(:, :)
      !> Derivative of the electronic energy w.r.t. strain deformations
      real(wp), intent(inout) :: sigma(:, :)

      integer :: iat, jat, izp, jzp, itr, img, inl, spin, nspin
      integer :: ish, jsh, is_, js, ii, jj, iao, jao, nao, ij
      real(wp) :: rr, r2, vec(3), cutoff2, hij, shpoly, dshpoly, dG(3), hscale
      real(wp) :: sval, dcni, dcnj, dhdcni, dhdcnj, hpij, pij
      real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :), tmp(:)
      real(wp), allocatable :: dstmp(:, :), ddtmpi(:, :, :), dqtmpi(:, :, :)
      real(wp), allocatable :: ddtmpj(:, :, :), dqtmpj(:, :, :)
      integer :: i,j,k
      real(wp) :: row(12)
      
      nspin = size(pmat, 3)

      allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & dtmp(3, msao(bas%maxl)**2), ddtmpi(3, 3, msao(bas%maxl)**2), &
      & qtmp(6, msao(bas%maxl)**2), dqtmpi(3, 6, msao(bas%maxl)**2), &
      & ddtmpj(3, 3, msao(bas%maxl)**2), dqtmpj(3, 6, msao(bas%maxl)**2))
      allocate(tmp(3))

      ! call print_cgto_alpha(bas)
      ! write(*,*) '===='
      ! call print_cgto_coeff(bas)
      ! write(*,*) '===='
      ! call print_cgto_nprim(bas)
      ! write(*,*) '===='
      
      ! tmp(:) = 0.0_wp;
      ! $omp parallel do schedule(runtime) default(none) reduction(+:dEdcn, gradient, sigma) &
      ! $omp shared(mol, bas, trans, h0, selfenergy, dsedcn, pot, pmat, xmat, list, nspin) &
      ! $omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, spin, &
      ! $omp& r2, vec, stmp, dtmp, qtmp, dstmp, ddtmpi, dqtmpi, ddtmpj, dqtmpj, hij, shpoly, &
      ! $omp& dshpoly, dG, dcni, dcnj, dhdcni, dhdcnj, hpij, rr, sval, hscale, pij, inl, img)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is_ = bas%ish_at(iat)
         inl = list_%inl(iat)
         do img = 1, list_%nnl(iat)
            jat = list_%nlat(img+inl)
            itr = list_%nltr(img+inl)
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            if (iat == jat) cycle
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is_+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call multipole_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                  & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddtmpj, dqtmpj, &
                  & ddtmpi, dqtmpi)

                  shpoly = (1.0_wp + h0%shpoly(ish, izp)*rr) &
                  & * (1.0_wp + h0%shpoly(jsh, jzp)*rr)
                  dshpoly = ((1.0_wp + h0%shpoly(ish, izp)*rr)*h0%shpoly(jsh, jzp)*rr &
                  & + (1.0_wp + h0%shpoly(jsh, jzp)*rr)*h0%shpoly(ish, izp)*rr) &
                  & * 0.5_wp / r2

                  hscale = h0%hscale(jsh, ish, jzp, izp)
                  hij = 0.5_wp * (selfenergy(is_+ish) + selfenergy(js+jsh)) * hscale
                  dhdcni = dsedcn(is_+ish) * shpoly * hscale
                  dhdcnj = dsedcn(js+jsh) * shpoly * hscale

                  dG(:) = 0.0_wp
                  dcni = 0.0_wp
                  dcnj = 0.0_wp
                  ! print*, jsh, jzp, ish, izp
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        do spin = 1, nspin
                           pij = pmat(jj+jao, ii+iao, spin)
                           hpij = pij * hij * shpoly
                           sval = 2*hpij - 2*xmat(jj+jao, ii+iao, spin) &
                              - pij * (pot%vao(jj+jao, spin) + pot%vao(ii+iao, spin))

                           tmp = matmul(ddtmpi(:, :, ij), pot%vdp(:, iat, spin))
                           dG(:) = dG + sval * dstmp(:, ij) &
                              + 2*hpij*stmp(ij) * dshpoly / shpoly * vec &
                              - pij * matmul(ddtmpi(:, :, ij), pot%vdp(:, iat, spin)) &
                              - pij * matmul(ddtmpj(:, :, ij), pot%vdp(:, jat, spin)) &
                              - pij * matmul(dqtmpi(:, :, ij), pot%vqp(:, iat, spin)) &
                              - pij * matmul(dqtmpj(:, :, ij), pot%vqp(:, jat, spin))

                           dcni = dcni + dhdcni * pmat(jj+jao, ii+iao, spin) * stmp(ij)
                           dcnj = dcnj + dhdcnj * pmat(jj+jao, ii+iao, spin) * stmp(ij)
                        end do
                     end do
                  end do
                  dEdcn(iat) = dEdcn(iat) + dcni
                  dEdcn(jat) = dEdcn(jat) + dcnj
                  gradient(:, iat) = gradient(:, iat) + dG
                  gradient(:, jat) = gradient(:, jat) - dG
                  sigma(:, :) = sigma + 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
                     + spread(dG, 1, 3) * spread(vec, 2, 3))

               end do
            end do

         end do
      end do

      ! $omp parallel do schedule(runtime) default(none) reduction(+:dEdcn) &
      ! $omp shared(mol, bas, dsedcn, pmat, nspin) &
      ! $omp private(iat, izp, jzp, is, ish, ii, iao, dcni, dhdcni, spin)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is_ = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is_+ish)
            dhdcni = dsedcn(is_+ish)
            dcni = 0.0_wp
            do iao = 1, msao(bas%cgto(ish, izp)%ang)
               do spin = 1, nspin
                  dcni = dcni + dhdcni * pmat(ii+iao, ii+iao, spin)
               end do
            end do
            dEdcn(iat) = dEdcn(iat) + dcni
         end do
      end do

   end subroutine get_hamiltonian_gradient


   subroutine get_occupation(mol, bas, h0, nocc, n0at, n0sh)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Hamiltonian interaction data
      type(tb_hamiltonian), intent(in) :: h0
      !> Occupation number
      real(wp), intent(out) :: nocc
      !> Reference occupation for each atom
      real(wp), intent(out) :: n0at(:)
      !> Reference occupation for each shell
      real(wp), intent(out) :: n0sh(:)

      integer :: iat, ish, izp, ii

      nocc = -mol%charge
      n0at(:) = 0.0_wp
      n0sh(:) = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            nocc = nocc + h0%refocc(ish, izp)
            n0at(iat) = n0at(iat) + h0%refocc(ish, izp)
            n0sh(ii+ish) = n0sh(ii+ish) + h0%refocc(ish, izp)
         end do
      end do

   end subroutine get_occupation


!> Shift multipole operator from Ket function (center i) to Bra function (center j),
!> the multipole operator on the Bra function can be assembled from the lower moments
!> on the Ket function and the displacement vector using horizontal shift rules.
   pure subroutine shift_operator(vec, s, di, qi, dj, qj)
      !> Displacement vector of center i and j
      real(wp),intent(in) :: vec(:)
      !> Overlap integral between basis functions
      real(wp),intent(in) :: s
      !> Dipole integral with operator on Ket function (center i)
      real(wp),intent(in) :: di(:)
      !> Quadrupole integral with operator on Ket function (center i)
      real(wp),intent(in) :: qi(:)
      !> Dipole integral with operator on Bra function (center j)
      real(wp),intent(out) :: dj(:)
      !> Quadrupole integral with operator on Bra function (center j)
      real(wp),intent(out) :: qj(:)

      real(wp) :: tr

      ! Create dipole operator on Bra function from Ket function and shift contribution
      ! due to monopol displacement
      dj(1) = di(1) + vec(1)*s
      dj(2) = di(2) + vec(2)*s
      dj(3) = di(3) + vec(3)*s

      ! For the quadrupole operator on the Bra function we first construct the shift
      ! contribution from the dipole and monopol displacement, since we have to remove
      ! the trace contribution from the shift and the moment integral on the Ket function
      ! is already traceless
      qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
      qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
      qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
      qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
      qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
      qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
      ! Now collect the trace of the shift contribution
      tr = 0.5_wp * (qj(1) + qj(3) + qj(6))

      ! Finally, assemble the quadrupole operator on the Bra function from the operator
      ! on the Ket function and the traceless shift contribution
      qj(1) = qi(1) + 1.5_wp * qj(1) - tr
      qj(2) = qi(2) + 1.5_wp * qj(2)
      qj(3) = qi(3) + 1.5_wp * qj(3) - tr
      qj(4) = qi(4) + 1.5_wp * qj(4)
      qj(5) = qi(5) + 1.5_wp * qj(5)
      qj(6) = qi(6) + 1.5_wp * qj(6) - tr
   end subroutine shift_operator


end module tblite_xtb_h0
