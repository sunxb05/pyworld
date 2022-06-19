#!/usr/bin/env python
#
#  gen1int: compute one-electron integrals using rotational London atomic-orbitals
#  Copyright 2009-2011 Bin Gao, and Andreas Thorvaldsen
#
#  gen1int is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  gen1int is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with gen1int. If not, see <http://www.gnu.org/licenses/>.

__author__ = "Bin Gao, and Andreas J. Thorvaldsen"
__copyright__ = "Copyright 2009-2011"
__credits__ = ["Radovan Bast", "Kenneth Ruud"]
__license__ = "LGPLv3"
__version__ = "0.2.1"
__maintainer__ = "Bin Gao"
__email__ = "bin.gao@uit.no"
__status__ = "Development"

#C  2009-11-21, Bin Gao:
#C  * modified from Gen1Int Version 0.1.0
#
#C> \brief module of property
#C> \details contains the definition of property and related subroutines
#C> \author Bin Gao
#C> \date 2009-11-21
#module g1int_property
#  ! precision
#  use xprecision
#  ! tracking errors
#  use xtrack
#  ! character processing
#  use xcharacter
#  ! molecule
#  use xmolecule
#  ! basis set
#  use xbasis
#  ! generalized one-electron integral
#  use g1int_integral
#
#  implicit none
#  !> length of the name of the property
#  integer, private, parameter :: len_name_prop = 8
#  !> definition of property
#  type, public :: g1prop_t
#    private
#    !> name of the property
#    character*(len_name_prop) :: prop_name = ''
#    !> range of generalized one-electron integrals to compute,
#    !> it would be useful when you just want to calculate, such as,
#    !> the dipole length integrals along X and Y directions,
#    !> you could then simply set range_int1 = (/1,2/)
#    integer, allocatable :: range_int1(:)
#    !> generalized one-electron integrals belong to this property
#    type(g1int1_t), allocatable :: int1(:)
#  end type g1prop_t
#
#  public :: g1int_property_set
#  public :: g1int_property_dump
#  public :: g1int_property_size
#  public :: g1int_property_calc
#  public :: g1int_property_clean
#
#  private :: g1int_property_alloc_int1
#
#  contains
#
#  ! (1) check the validity of the given name of property, or
#  ! (2) generate the labels of the integrals belong to this property, or
#  ! (3) create the integrals (and their operators) of this property for calculations
#  subroutine g1int_property_set( prop_name, info_deriv,
#                                 this_prop, prop_int_param, set_opt1, this_mol, this_bas, 
#                                 prop_rfun, prop_aux_order, )
#    implicit none
#  !> \param prop_int_param contains some additional information of the property,
#  !>        such as the order of Cartesian and spherical multipole moment integrals
#    integer, optional, intent(in) :: prop_int_param(:)
#
#
#    ! given name to check if it is a defined property
#    character*(*), optional, intent(in) :: prop_name
#    ! type of r_C function, could be '1', '2', 'D' or 'X' mean 1/r_C, 1/r_C^2,
#    ! Delta(r_C) or no r_C function
#    character, optional, intent(out) :: prop_rfun
#    ! order of auxiliary functions from the property
#    ! should be the sum of
#    ! (1) powers of multipole moments,
#    ! (2) orders of derivatives of 1/r_C function,
#    ! (3) orders of derivatives with respect to electrons,
#    ! of the operators of property
#    ! -1 means this property does not need auxiliary function
#    !
#    ! NB: since the geometric derivatives is defined outside the operators,
#    !     its order will be taken care by other subroutines, such as those
#    !     in gen1int_input.F90
#    integer, optional, intent(out) :: prop_aux_order
#    ! number of centers belong to this property (+2 by consider the centers of basis sets)
#    integer, optional, intent(out) :: prop_num_cent
#    ! property
#    type(g1prop), optional, intent(inout) :: this_prop
#    ! if create the operators for the integrals of the property
#    logical, optional, intent(in) :: set_opt1
#    !////////////////////////////////////////////////////////////////
#    !                       local variables
#    !////////////////////////////////////////////////////////////////
#    ! name to check
#    character*(x_max_lenc) chk_name
#    ! if we are only asked to check the name of property
#    logical chk_prop
#    ! if create the operators for the integrals of the property
#    logical new_op1
#    ! prefactor
#    real(xp) prefact1
#    ! incremental recorder for integrals
#    integer int1
#    ! incremental recorder for atoms
#    integer iatom
#    ! incremental recorders for XYZ directions
#    integer ixyz, jxyz
#    ! stuff for generating labels of operators
#    integer iord1
#    integer nx, ny, nz
#    integer mx, my, mz
#    integer lx, ly, lz
#    integer kx, ky, kz
#    character*1, parameter :: xpower1(3)  = (/'X','Y','Z'/)
#    character*2, parameter :: xpower2(6)  = (/'XX','XY','YY','XZ','YZ','ZZ'/)
#    character*3, parameter :: xpower3(10) = (/'XXX','XXY','XYY','YYY','XXZ', &
#                                              'XYZ','YYZ','XZZ','YZZ','ZZZ'/)
#    character*7, parameter :: xderiv2(6)  = (/'dd/dxdx','dd/dxdy','dd/dydy', &
#                                              'dd/dxdz','dd/dydz','dd/dzdz'/)
#    ! for derivatives of London atomic orbitals
#    type(g1lao), pointer :: currt_bderv_lao
#    ! push current subroutine into the stack
#    call xsub_enter('g1int_property_set')
#    ! it is not just check the name
#    if ( present(this_prop) ) then
#      ! check if GEN1INT initialized
#      if (xlib_uninit) call xstop(trim(xlib_name)//' is not initialized!')
#      chk_prop = .false.
#      ! if create the operators for the integrals of the property
#      if ( present(set_opt1) ) then
#        new_op1 = set_opt1
#      else
#        new_op1 = .true.
#      end if
#      ! name to check
#      chk_name = this_prop%prop_name
#    ! if we are only asked to check the name
#    else
#      chk_prop = present(prop_rfun) .and. present(prop_aux_order)
#      if ( chk_prop ) then
#        ! name to check
#        chk_name = prop_name
#        ! do not create the operators
#        new_op1 = .false.
#      ! if there is not enough arguments
#      else
#        if ( xlib_lprint.ge.5 ) &
#          write(x_lupri,1000)'Warning! Noting will be done in g1int_property_set!'
#        ! pop the stack
#        call xsub_exit
#        return
#      end if
#    end if
#    ! check the given name of a property
#    select case( trim(chk_name) )
#    !========================================
#    ! One-electron potential energy integral
#    !========================================
#!fixme: maybe this is not efficient for calculations ...
#    case('*1ELPOT')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 0
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,1/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! label of the computed integrals
#        this_prop%g1prop_int1(1)%label = 'POTENERG_S'
#        ! create operators for the computed integrals
#        if ( new_op1 ) then
#          ! the number of operators
#          this_prop%g1prop_int1(1)%num_g1opt = num_atoms
#          ! allocate operators
#          allocate( this_prop%g1prop_int1(1)%g1int_opt(num_atoms), stat=x_ierr )
#          if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#          ! loop over atoms
#          do iatom = 1, num_atoms
#            ! create operators
#            call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(1)% &
#                                                       g1int_opt(iatom),        &
#!fixme: we should take minus sigh, right?
#                                          prefactor  = -charge_atoms(iatom),    &
#                                          r_func     = '1',                     &
#                                          atom_rcent = iatom )
#          end do ! loop over atoms
#        end if
#      end if
#    !===============================================
#    ! Electronic angular momentum around the nuclei
#    !===============================================
#    case('*ANGLON')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3/) )
#        ! loop over angular momentum
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#          this_prop%g1prop_int1(int1)%label = xpower1(int1)//'ANGLON__N'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = 2
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            select case(int1)
#            ! operators -(y_N*d/dz - z_N*d/dy)
#            case(1)
#              nx=0; ny=1; nz=0
#              mx=0; my=0; mz=1
#            ! operators -(z_N*d/dx - x_N*d/dz)
#            case(2)
#              nx=0; ny=0; nz=1
#              mx=1; my=0; mz=0
#            ! operators -(x_N*d/dy - y_N*d/dx)
#            case(3)
#              nx=1; ny=0; nz=0
#              mx=0; my=1; mz=0
#            end select
#            ! create operators
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(1),               &
#                                          prefactor   = -one,                       &
#                                          mult_power  = (/nx,ny,nz/),               &
#                                          ! the center of multipole moments is always the center of ket
#                                          atom_mcent  = -1,                         &
#                                          ederv_order = (/mx,my,mz/) )
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(2),               &
#                                          mult_power  = (/mx,my,mz/),               &
#                                          ! the center of multipole moments is always the center of ket
#                                          atom_mcent  = -1,                         &
#                                          ederv_order = (/nx,ny,nz/) )
#          end if
#        end do ! loop over angular momentum
#      end if
#    !===============================================
#    ! Electronic angular momentum around the origin
#    !===============================================
#    case('*ANGMOM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3/) )
#        ! loop over angular momentum
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#          this_prop%g1prop_int1(int1)%label = xpower1(int1)//'ANGMOM__A'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = 2
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            select case(int1)
#            ! operators -(y_O*d/dz - z_O*d/dy)
#            case(1)
#              nx=0; ny=1; nz=0
#              mx=0; my=0; mz=1
#            ! operators -(z_O*d/dx - x_O*d/dz)
#            case(2)
#              nx=0; ny=0; nz=1
#              mx=1; my=0; mz=0
#            ! operators -(x_O*d/dy - y_O*d/dx)
#            case(3)
#              nx=1; ny=0; nz=0
#              mx=0; my=1; mz=0
#            end select
#            ! create operators
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(1),               &
#                                          prefactor   = -one,                       &
#                                          mult_power  = (/nx,ny,nz/),               &
#                                          coord_mcent = gaug_orig,                  &
#                                          ederv_order = (/mx,my,mz/) )
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(2),               &
#                                          mult_power  = (/mx,my,mz/),               &
#                                          coord_mcent = gaug_orig,                  &
#                                          ederv_order = (/nx,ny,nz/) )
#          end if
#        end do ! loop over angular momentum
#      end if
#    !===========================================================
#    ! Cartesian multipole moment integrals (of order g1_iorcar)
#    !===========================================================
#    case('*CARMOM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(g1_iorcar+1)*(g1_iorcar+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian powers
#        do nz = 0, g1_iorcar
#          do ny = 0, g1_iorcar-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian power of X direction
#              nx = g1_iorcar - (ny+nz)
#              ! label of the computed integrals
#!fixme: if we consider more than 100th orders? do we?
#              this_prop%g1prop_int1(int1)%label = 'CM'//xfigure(nx/10)//xfigure(mod(nx,10))  &
#                                                     //xfigure(ny/10)//xfigure(mod(ny,10))  &
#                                                     //xfigure(nz/10)//xfigure(mod(nz,10))//'_S'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                ! create operators
#                call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                           g1int_opt(1),               &
#                                              mult_power = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian powers
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian powers
#      end if
#    !=======================
#    ! Darwin type integrals
#    !=======================
#    case('*DARWIN')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'D'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,1/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! label of the computed integrals
#        this_prop%g1prop_int1(1)%label = 'DARWIN___S'
#        ! create operators for the computed integrals
#        if ( new_op1 ) then
#          ! the number of operators
#          this_prop%g1prop_int1(1)%num_g1opt = num_atoms
#          ! allocate operators
#          allocate( this_prop%g1prop_int1(1)%g1int_opt(num_atoms), stat=x_ierr )
#          if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#          ! prefactor
#          prefact1 = pi*alpha2_const*half
#          ! loop over atoms
#          do iatom = 1, num_atoms
#            ! create operators
#            call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(1)%      &
#                                                       g1int_opt(iatom),             &
#                                          prefactor  = prefact1*charge_atoms(iatom), &
#                                          r_func     = 'D',                          &
#                                          atom_rcent = iatom )
#          end do ! loop over atoms
#        end if
#      end if
#    !==================================
#    ! Dipole moment (length) integrals
#    !==================================
#    case('*DIPLEN')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! order of multipole moments
#        iord1 = 1
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(iord1+1)*(iord1+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian powers
#        do nz = 0, iord1
#          do ny = 0, iord1-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian power of X direction
#              nx = iord1 - (ny+nz)
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = xpower1(int1)//'DIPLEN__S'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                ! create operators
#                call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                           g1int_opt(1),               &
#                                              mult_power = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian powers
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian powers
#      end if
#    !===========================
#    ! Dipole velocity integrals
#    !===========================
#    case('*DIPVEL')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! order of total derivatives
#        iord1 = 1
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(iord1+1)*(iord1+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian derivatives
#        do nz = 0, iord1
#          do ny = 0, iord1-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian derivative of X direction
#              nx = iord1 - (ny+nz)
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = xpower1(int1)//'DIPVEL__A'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                ! create operators
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(1),               &
#                                              ederv_order = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian derivatives
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian derivatives
#      end if
#    !================================================================
#    ! Overlap integrals for the small component in DPT ( <opG|opG> )
#    !================================================================
#    case('*DPTOVL')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! order of total derivatives
#        iord1 = 2
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(iord1+1)*(iord1+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian derivatives
#        do nz = 0, iord1
#          do ny = 0, iord1-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian derivative of X direction
#              nx = iord1 - (ny+nz)
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = xderiv2(int1)//'__S'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                ! create operators
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(1),               &
#                                              prefactor   = -fourth,                    &
#                                              ederv_order = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian derivatives
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian derivatives
#      end if
#    !===========================================================
#    ! London orbital contribution to diamagnetic susceptibility
#    !===========================================================
#    case('*DSUSLH')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 2
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,6/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! loop over integrals
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#!fixme: if there are more than 100 atoms
#          this_prop%g1prop_int1(int1)%label = xpower2(int1)//'DSUSLH_N'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            select case(int1)
#            ! XX
#            case(1)
#              kx=2; ky=0; kz=0
#            ! XY
#            case(2)
#              kx=1; ky=1; kz=0
#            ! YY
#            case(3)
#              kx=0; ky=2; kz=0
#            ! XZ
#            case(4)
#              kx=1; ky=0; kz=1
#            ! YZ
#            case(5)
#              kx=0; ky=1; kz=1
#            ! ZZ
#            case(6)
#              kx=0; ky=0; kz=2
#            end select
#            ! magnetic derivatives from London atomic orbitals
#            call gen1int_bderv_london( (/kx,ky,kz/) )
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = ( num_atoms + 3 )*num_bderv_lao
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt &
#                      (this_prop%g1prop_int1(int1)%num_g1opt), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            !-------------------------------
#            ! one-electron potential energy
#            !-------------------------------
#            ! initialize the index of operator
#            ixyz = 0
#            ! loop over atoms
#            do iatom = 1, num_atoms
#              ! create operators
#              currt_bderv_lao => head_bderv_lao
#              ! loop over generated expansion terms
#              do jxyz = 1, num_bderv_lao
#                ! index of operator
#                ixyz = ixyz + 1
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)%  &
#                                                            g1int_opt(ixyz),             &
#                                              prefactor   = -charge_atoms(iatom)         &
#                                                          * currt_bderv_lao%prefactor,   &
#                                              ndisp_power = currt_bderv_lao%ndisp_power, &
#                                              mult_power  = currt_bderv_lao%mult_power,  &
#                                              r_func      = '1',                         &
#                                              atom_rcent  = iatom )
#                ! next value
#                currt_bderv_lao => currt_bderv_lao%next_lao
#              end do ! loop over generated expansion terms
#            end do ! loop over atoms
#            !---------------------------
#            ! electronic kinetic energy
#            !---------------------------
#            ! loop over Cartesian derivatives
#            do nz = 0, 1
#              do ny = 0, 1-nz
#                ! generate Cartesian derivatives of X direction
#                nx = 1 - (ny+nz)
#                ! create operators
#                currt_bderv_lao => head_bderv_lao
#                ! loop over generated expansion terms
#                do jxyz = 1, num_bderv_lao
#                  ! index of operator
#                  ixyz = ixyz + 1
#                  ! since we create first order derivatives, we need to multiply by 2
#                  call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)%      &
#                                                              g1int_opt(ixyz),                 &
#                                                prefactor   = -half*currt_bderv_lao%prefactor, &
#                                                ndisp_power = currt_bderv_lao%ndisp_power,     &
#                                                mult_power  = currt_bderv_lao%mult_power,      &
#                                                ederv_order = (/2*nx,2*ny,2*nz/) )
#                  ! next value
#                  currt_bderv_lao => currt_bderv_lao%next_lao
#                end do ! loop over generated expansion terms
#              end do ! loop over Cartesian derivatives
#            end do ! loop over Cartesian derivatives
#            ! clean linked list
#            call gen1int_clean_laolist
#          end if
#        end do ! loop over integrals
#      end if
#    !===============================================================
#    ! Diamagnetic susceptiblity without London orbital contribution
#    !===============================================================
#    case('*DSUSNL')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals, the integrals are,
#      ! XX: y^2 + z^2, YY: z^2 + x^2, ZZ: x^2 + y^2, 
#      ! XY: -x*y, YZ: -y*z, and XZ: -x*z
#      else
#        ! order of multipole moments
#        iord1 = 2
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(iord1+1)*(iord1+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian powers
#        do nz = 0, iord1
#          do ny = 0, iord1-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian power of X direction
#              nx = iord1 - (ny+nz)
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = xpower2(int1)//'DSUSNL_N'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! XX, YY, and ZZ components
#                if ( nx.eq.2 .or. ny.eq.2 .or. nz.eq.2 ) then
#                  ! the number of operators
#                  this_prop%g1prop_int1(int1)%num_g1opt = 2
#                  ! allocate operators
#                  allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#                  if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                  ! create operators
#                  call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                             g1int_opt(1),               &
#                                                prefactor  = fourth,                     &
#                                                mult_power = (/nz,nx,ny/),               &
#                                                ! the center of multipole moments is always the center of ket
#                                                atom_mcent = -1 )
#                  call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                             g1int_opt(2),               &
#                                                prefactor  = fourth,                     &
#                                                mult_power = (/ny,nz,nx/),               &
#                                                ! the center of multipole moments is always the center of ket
#                                                atom_mcent = -1 )
#                ! XY, XZ, and YZ components
#                else
#                  ! the number of operators
#                  this_prop%g1prop_int1(int1)%num_g1opt = 1
#                  ! allocate operators
#                  allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                  if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                  ! create operators, XY: -x*y, YZ: -y*z, and XZ: -x*z
#                  call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                             g1int_opt(1),               &
#                                                prefactor  = -fourth,                    &
#                                                mult_power = (/nx,ny,nz/),               &
#                                                ! the center of multipole moments is always the center of ket
#                                                atom_mcent = -1 )
#                end if
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian powers
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian powers
#      end if
#    !==============================================================
#    ! Electric field gradient at the individual nuclei (Cartesian)
#    !==============================================================
#    case('*EFGCAR')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 2
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals, the integrals are,
#      ! XX: ( 3x^2 - r^2 )/r^5,
#      ! YY: ( 3y^2 - r^2 )/r^5,
#      ! ZZ: ( 3z^2 - r^2 )/r^5,
#      ! XY: 3xy/r^5, YZ: 3yz/r^5, and XZ: 3xz/r^5
#      else
#        ! order of multipole moments
#        iord1 = 2
#        ! allocate computed integrals
#!fixme: considering molecular symmetry
#        call g1int_property_alloc_int1( this_prop, (/1,num_atoms*(iord1+1)*(iord1+2)/2/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! index for generating label
#          ixyz = 0
#          ! loop over Cartesian powers
#          do nz = 0, iord1
#            do ny = 0, iord1-nz
#              ! index of integral
#              int1 = int1 + 1
#              ! index for generating label
#              ixyz = ixyz + 1
#              ! this is not the integral we want
#              if ( int1.lt.this_prop%range_g1int(1) ) then
#                cycle
#              ! we get the required integral
#              else if ( int1.le.this_prop%range_g1int(2) ) then
#                ! generate Cartesian power of X direction
#                nx = iord1 - (ny+nz)
#                ! label of the computed integrals
#                this_prop%g1prop_int1(int1)%label = xpower2(ixyz)//'EFG'// &
#                                                   xfigure(iatom/10)//    &
#                                                   xfigure(mod(iatom,10))//'1_S'
#                ! create operators for the computed integrals
#                if ( new_op1 ) then
#                  ! XX, YY, and ZZ components
#                  if ( nx.eq.2 .or. ny.eq.2 .or. nz.eq.2 ) then
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 3
#                    ! prefactor
#                    prefact1 = third2
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(3), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                  ! XY, XZ, and YZ components
#                  else
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 1
#                    ! prefactor
#                    prefact1 = one
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                  end if
#                  ! create operator 2x^2/r^5, 2z^2/r^5, 2z^2/r^5, 3xy/r^5, 3yz/r^5, or 3xz/r^5
#                  call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                             g1int_opt(1),               &
#                                                prefactor  = prefact1,                   &
#                                                rderv_order = (/nx,ny,nz/),              &
#                                                r_func      = '1',                       &
#                                                atom_rcent  = iatom )
#                  ! extra operators for XX, YY, and ZZ components
#                  if ( nx.eq.2 .or. ny.eq.2 .or. nz.eq.2 ) then
#                    ! create operator -x^2/r^5, -y^2/r^5, or -z^2/r^5
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(2),               &
#                                                  prefactor   = -third,                     &
#                                                  rderv_order = (/ny,nz,nx/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                    ! create operator -x^2/r^5, -y^2/r^5, or -z^2/r^5
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(3),               &
#                                                  prefactor   = -third,                     &
#                                                  rderv_order = (/nz,nx,ny/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                  end if
#                end if
#              ! we have got all the required integrals, exit
#              else
#                exit
#              end if
#            end do ! loop over Cartesian powers
#            ! we have got all the required integrals, exit
#            if ( int1.ge.this_prop%range_g1int(2) ) exit
#          end do ! loop over Cartesian powers
#        end do ! loop over atoms
#      end if
#    !======================================
#    ! One-electron Fermi contact integrals
#    !======================================
#    case('*FC')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'D'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! prefactor
#        prefact1 = four*pi*ge_factor/three
#        ! loop over atoms
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#!fixme: if there are more than 100 atoms
#          this_prop%g1prop_int1(int1)%label = 'FC_'//name_atoms(int1)(1:3)  &
#                                                  //xfigure(int1/10)       &
#                                                  //xfigure(mod(int1,10))//'_S'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = 1
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                       g1int_opt(1),               &
#                                          prefactor  = prefact1,                   &
#                                          r_func     = 'D',                        &
#                                          atom_rcent = int1 )
#          end if
#        end do ! loop over atoms
#      end if
#    !=====================================
#    ! Electronic kinetic energy integrals
#    !=====================================
#    case('*KINENE')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,1/) )
#        ! label of the computed integrals
#        this_prop%g1prop_int1(1)%label = 'KINENERG_S'
#        ! create operators for the computed integrals
#        if ( new_op1 ) then
#          ! the number of operators
#          this_prop%g1prop_int1(1)%num_g1opt = 3
#          ! allocate operators
#          allocate( this_prop%g1prop_int1(1)%g1int_opt(3), stat=x_ierr )
#          if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#          ! initialize the index of operator
#          ixyz = 0
#          ! loop over Cartesian derivatives
#          do nz = 0, 1
#            do ny = 0, 1-nz
#              ! index of operator
#              ixyz = ixyz + 1
#              ! generate Cartesian derivatives of X direction
#              nx = 1 - (ny+nz)
#              ! create operators
#              ! since we create first order derivatives, we need to multiply by 2
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(1)% &
#                                                          g1int_opt(ixyz),         &
#                                            prefactor   = -half,                   &
#                                            ederv_order = (/2*nx,2*ny,2*nz/) )
#            end do ! loop over Cartesian derivatives
#          end do ! loop over Cartesian derivatives
#        end if
#      end if
#    !================================================
#    ! London orbital contribution to magnetic moment
#    !================================================
#    case('*LONMOM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! loop over integrals
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#!fixme: if there are more than 100 atoms
#          this_prop%g1prop_int1(int1)%label = xpower1(int1)//'LONMOM__N'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = num_atoms + num_atoms + 6
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt &
#                      (this_prop%g1prop_int1(int1)%num_g1opt), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            select case(int1)
#            ! Y_{MN}*z-Z_{MN}*y
#            case(1)
#              kx=0; ky=1; kz=0
#              lx=0; ly=0; lz=1
#            ! Z_{MN}*x-X_{MN}*z
#            case(2)
#              kx=0; ky=0; kz=1
#              lx=1; ly=0; lz=0
#            ! X_{MN}*y-Y_{MN}*x
#            case(3)
#              kx=1; ky=0; kz=0
#              lx=0; ly=1; lz=0
#            end select
#            !-------------------------------
#            ! one-electron potential energy
#            !-------------------------------
#            ! initialize the index of operator
#            ixyz = 0
#            ! loop over atoms
#            do iatom = 1, num_atoms
#              ! index of operator
#              ixyz = ixyz + 1
#              ! create operators
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)%  &
#                                                          g1int_opt(ixyz),             &
#                                            prefactor   = -half*charge_atoms(iatom),   &
#                                            ndisp_power = (/kx,ky,kz/),                &
#                                            mult_power  = (/lx,ly,lz/),                &
#                                            r_func      = '1',                         &
#                                            atom_rcent  = iatom )
#              ! index of operator
#              ixyz = ixyz + 1
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(ixyz),            &
#                                            prefactor   = half*charge_atoms(iatom),   &
#                                            ndisp_power = (/lx,ly,lz/),               &
#                                            mult_power  = (/kx,ky,kz/),               &
#                                            r_func      = '1',                        &
#                                            atom_rcent  = iatom )
#            end do ! loop over atoms
#            !---------------------------
#            ! electronic kinetic energy
#            !---------------------------
#            ! loop over Cartesian derivatives
#            do nz = 0, 1
#              do ny = 0, 1-nz
#                ! generate Cartesian derivatives of X direction
#                nx = 1 - (ny+nz)
#                ! index of operator
#                ixyz = ixyz + 1
#                ! create operators
#                ! since we create first order derivatives, we need to multiply by 2
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(ixyz),            &
#                                              prefactor   = -fourth,                    &
#                                              ndisp_power = (/kx,ky,kz/),               &
#                                              mult_power  = (/lx,ly,lz/),               &
#                                              ederv_order = (/2*nx,2*ny,2*nz/) )
#                ! index of operator
#                ixyz = ixyz + 1
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(ixyz),            &
#                                              prefactor   = fourth,                     &
#                                              ndisp_power = (/lx,ly,lz/),               &
#                                              mult_power  = (/kx,ky,kz/),               &
#                                              ederv_order = (/2*nx,2*ny,2*nz/) )
#              end do ! loop over Cartesian derivatives
#            end do ! loop over Cartesian derivatives
#          end if
#        end do ! loop over integrals
#      end if
#    !==============================================
#    ! One-electron contribution to magnetic moment
#    !==============================================
#    case('*MAGMOM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! loop over dipole moments
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#!fixme: if there are more than 100 atoms
#          this_prop%g1prop_int1(int1)%label = 'dh/dB'//xpower1(int1)//'___A'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = num_atoms + num_atoms + 8
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt &
#                      (this_prop%g1prop_int1(int1)%num_g1opt), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            select case(int1)
#            ! Y_{MN}*z-Z_{MN}*y
#            case(1)
#              kx=0; ky=1; kz=0
#              lx=0; ly=0; lz=1
#            ! Z_{MN}*x-X_{MN}*z
#            case(2)
#              kx=0; ky=0; kz=1
#              lx=1; ly=0; lz=0
#            ! X_{MN}*y-Y_{MN}*x
#            case(3)
#              kx=1; ky=0; kz=0
#              lx=0; ly=1; lz=0
#            end select
#            !-------------------------------
#            ! one-electron potential energy
#            !-------------------------------
#            ! initialize the index of operator
#            ixyz = 0
#            ! loop over atoms
#            do iatom = 1, num_atoms
#              ! index of operator
#              ixyz = ixyz + 1
#              ! create operators
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)%  &
#                                                          g1int_opt(ixyz),             &
#                                            prefactor   = -half*charge_atoms(iatom),   &
#                                            ndisp_power = (/kx,ky,kz/),                &
#                                            mult_power  = (/lx,ly,lz/),                &
#                                            r_func      = '1',                         &
#                                            atom_rcent  = iatom )
#              ! index of operator
#              ixyz = ixyz + 1
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(ixyz),            &
#                                            prefactor   = half*charge_atoms(iatom),   &
#                                            ndisp_power = (/lx,ly,lz/),               &
#                                            mult_power  = (/kx,ky,kz/),               &
#                                            r_func      = '1',                        &
#                                            atom_rcent  = iatom )
#            end do ! loop over atoms
#            !---------------------------
#            ! electronic kinetic energy
#            !---------------------------
#            ! loop over Cartesian derivatives
#            do nz = 0, 1
#              do ny = 0, 1-nz
#                ! generate Cartesian derivatives of X direction
#                nx = 1 - (ny+nz)
#                ! index of operator
#                ixyz = ixyz + 1
#                ! create operators
#                ! since we create first order derivatives, we need to multiply by 2
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(ixyz),            &
#                                              prefactor   = -fourth,                    &
#                                              ndisp_power = (/kx,ky,kz/),               &
#                                              mult_power  = (/lx,ly,lz/),               &
#                                              ederv_order = (/2*nx,2*ny,2*nz/) )
#                ! index of operator
#                ixyz = ixyz + 1
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(ixyz),            &
#                                              prefactor   = fourth,                     &
#                                              ndisp_power = (/lx,ly,lz/),               &
#                                              mult_power  = (/kx,ky,kz/),               &
#                                              ederv_order = (/2*nx,2*ny,2*nz/) )
#              end do ! loop over Cartesian derivatives
#            end do ! loop over Cartesian derivatives
#            !-----------------------------------------------
#            ! electronic angular momentum around the nuclei
#            !-----------------------------------------------
#            ! index of operator
#            ixyz = ixyz + 1
#            ! create operators
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(ixyz),            &
#                                          prefactor   = -half,                      &
#                                          mult_power  = (/kx,ky,kz/),               &
#                                          ! the center of multipole moments is always the center of ket
#                                          atom_mcent  = -1,                         &
#                                          ederv_order = (/lx,ly,lz/) )
#            ! index of operator
#            ixyz = ixyz + 1
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(ixyz),            &
#                                          prefactor   = half,                       &
#                                          mult_power  = (/lx,ly,lz/),               &
#                                          ! the center of multipole moments is always the center of ket
#                                          atom_mcent  = -1,                         &
#                                          ederv_order = (/kx,ky,kz/) )
#          end if
#        end do ! loop over integrals
#      end if
#    !=========================
#    ! Mass velocity integrals
#    !=========================
#    case('*MASSVE')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,1/) )
#        ! label of the computed integrals
#        this_prop%g1prop_int1(1)%label = 'MASSVELO_S'
#        ! create operators for the computed integrals
#        if ( new_op1 ) then
#          ! the number of operators
#          this_prop%g1prop_int1(1)%num_g1opt = 6
#          ! allocate operators
#          allocate( this_prop%g1prop_int1(1)%g1int_opt(6), stat=x_ierr )
#          if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#          ! initialize the index of operator
#          ixyz = 0
#          ! loop over Cartesian derivatives
#          do nz = 0, 2
#            do ny = 0, 2-nz
#              ! index of operator
#              ixyz = ixyz + 1
#              ! generate Cartesian derivatives of X direction
#              nx = 2 - (ny+nz)
#              ! for d^4/dx^4, d^4/dy^4, and d^4/dz^4, we use alpha^2/8 (alpha, fine structure constant)
#              if ( nx.eq.2 .or. ny.eq.2 .or. nz.eq.2 ) then
#                prefact1 = -alpha2_const*eighth
#              ! for d^4/dx^2dy^2, d^4/dx^2dz^2, and d^4/dy^2dz^2, we use alpha^2/4
#              else
#                prefact1 = -alpha2_const*fourth
#              end if
#              ! create operators
#              ! since we create first order derivatives, we need to multiply by 2
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(1)% &
#                                                          g1int_opt(ixyz),         &
#                                            prefactor   = prefact1,                &
#                                            ederv_order = (/2*nx,2*ny,2*nz/) )
#            end do ! loop over Cartesian derivatives
#          end do ! loop over Cartesian derivatives
#        end if
#      end if
#    !=========================================
#    ! Electric field at the individual nuclei
#    !=========================================
#    case('*NELFLD')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3*num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over derivative directions: X, Y, Z
#          do ixyz = 1, 3
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = 'NEF_'//xfigure(int1/100)         &
#                                                       //xfigure(mod(int1,100)/10) &
#                                                       //xfigure(mod(mod(int1,100),10))//'__S'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                select case(ixyz)
#                case(1)
#                  nx=1; ny=0; nz=0
#                case(2)
#                  nx=0; ny=1; nz=0
#                case(3)
#                  nx=0; ny=0; nz=1
#                end select
#                ! create operators
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(1),               &
#                                              rderv_order = (/nx,ny,nz/),               &
#                                              r_func      = '1',                        &
#                                              atom_rcent  = iatom )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over derivative directions
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !====================================
#    ! Nuclear shielding tensor integrals
#    !====================================
#    case('*NST')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 3
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals, the integrals are,
#      !
#      ! y_N*y_K/r_K^3+z_N*z_K/r_K^3 - (Y_{MN}*z-Z_{MN}*y)*(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy), (001 x)
#      ! -x_N*y_K/r_K^3              - (Y_{MN}*z-Z_{MN}*y)*(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz), (001 y)
#      ! -x_N*z_K/r_K^3              - (Y_{MN}*y-Z_{MN}*x)*(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx), (001 z)
#      ! -y_N*x_K/r_K^3              - (Z_{MN}*x-X_{MN}*z)*(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy), (002 x)
#      ! x_N*x_K/r_K^3+z_N*z_K/r_K^3 - (Z_{MN}*x-X_{MN}*z)*(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz), (002 y)
#      ! -y_N*z_K/r_K^3              - (Z_{MN}*x-X_{MN}*z)*(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx), (002 z)
#      ! -z_N*x_K/r_K^3              - (X_{MN}*y-Y_{MN}*x)*(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy), (003 x)
#      ! -z_N*y_K/r_K^3              - (X_{MN}*y-Y_{MN}*x)*(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz), (003 y)
#      ! x_N*x_K/r_K^3+y_N*y_K/r_K^3 - (X_{MN}*y-Y_{MN}*x)*(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx), (003 z)
#      !
#      ! where K runs for all atoms
#      else   
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,9*num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! initialize the index of integral
#        int1 = 0
#        ! initialize the index of the nuclear magnetic moment coordinate
#        iord1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over the nuclear magnetic moment coordinate
#          do jxyz = 1, 3
#            ! index of the nuclear magnetic moment coordinate
#            iord1 = iord1 + 1
#            ! loop over the x, y, z components of the magnetic field
#            do ixyz = 1, 3
#              ! index of integral
#              int1 = int1 + 1
#              ! this is not the integral we want
#              if ( int1.lt.this_prop%range_g1int(1) ) then
#                cycle
#              ! we get the required integral
#              else if ( int1.le.this_prop%range_g1int(2) ) then
#                ! label of the computed integrals
#                this_prop%g1prop_int1(int1)%label = xfigure(iord1/100)//              &
#                                                   xfigure(mod(iord1,100)/10)//      &
#                                                   xfigure(mod(mod(iord1,100),10))// &
#                                                   '_NST'//xpower1(ixyz)//'_S'
#                ! create operators for the computed integrals
#                if ( new_op1 ) then
#                  ! XX, YY, ZZ components
#                  if ( ixyz.eq.jxyz ) then
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 6
#                    select case(ixyz)
#                    ! (y_N*y_K + z_N*z_K)/r_K^3 - (Y_{MN}*z-Z_{MN}*y)*(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy)
#                    case(1)
#                      nx=0; ny=1; nz=0
#                    ! (z_N*z_K + x_N*x_K)/r_K^3 - (Z_{MN}*x-X_{MN}*z)*(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz)
#                    case(2)
#                      nx=0; ny=0; nz=1
#                    ! (x_N*x_K + y_N*y_K)/r_K^3 - (X_{MN}*y-Y_{MN}*x)*(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx)
#                    case(3)
#                      nx=1; ny=0; nz=0
#                    end select
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(6), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(1),               &
#                                                  prefactor   = half,                       &
#                                                  mult_power  = (/nx,ny,nz/),               &
#                                                  atom_mcent  = -1,                         &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(2),               &
#                                                  prefactor   = half,                       &
#                                                  mult_power  = (/nz,nx,ny/),               &
#                                                  atom_mcent  = -1,                         &
#                                                  rderv_order = (/nz,nx,ny/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(3),               &
#                                                  prefactor   = half,                       &
#                                                  ndisp_power = (/nx,ny,nz/),               &
#                                                  mult_power  = (/nz,nx,ny/),               &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/nz,nx,ny/) )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(4),               &
#                                                  prefactor   = -half,                      &
#                                                  ndisp_power = (/nx,ny,nz/),               &
#                                                  mult_power  = (/nz,nx,ny/),               &
#                                                  rderv_order = (/nz,nx,ny/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/nx,ny,nz/) )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(5),               &
#                                                  prefactor   = -half,                      &
#                                                  ndisp_power = (/nz,nx,ny/),               &
#                                                  mult_power  = (/nx,ny,nz/),               &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/nz,nx,ny/) )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(6),               &
#                                                  prefactor   = half,                       &
#                                                  ndisp_power = (/nz,nx,ny/),               &
#                                                  mult_power  = (/nx,ny,nz/),               &
#                                                  rderv_order = (/nz,nx,ny/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/nx,ny,nz/) )
#                  else
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 5
#                    ! nuclear magnetic moment coordinate
#                    select case(jxyz)
#                    ! operators -x_N*y_K/r_K^3, and -x_N*z_K/r_K^3
#                    ! operators -(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy)
#                    case(1)
#                      mx=1; my=0; mz=0
#                    ! operators -y_N*x_K/r_K^3, and -y_N*z_K/r_K^3
#                    ! operators -(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz)
#                    case(2)
#                      mx=0; my=1; mz=0
#                    ! operators -z_N*x_K/r_K^3, and -z_N*y_K/r_K^3
#                    ! operators -(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx)
#                    case(3)
#                      mx=0; my=0; mz=1
#                    end select
#                    ! x, y, z components of the magnetic field
#                    select case(ixyz)
#                    ! operators -y_N*x_K/r_K^3, and -z_N*x_K/r_K^3
#                    ! operators Y_{MN}*z-Z_{MN}*y
#                    case(1)
#                      nx=1; ny=0; nz=0
#                    ! operators -x_N*y_K/r_K^3, and -z_N*y_K/r_K^3
#                    ! operators Z_{MN}*x-X_{MN}*z
#                    case(2)
#                      nx=0; ny=1; nz=0
#                    ! operators -x_N*z_K/r_K^3, and -y_N*z_K/r_K^3
#                    ! operators X_{MN}*y-Y_{MN}*x
#                    case(3)
#                      nx=0; ny=0; nz=1
#                    end select
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(5), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(1),               &
#                                                  prefactor   = -half,                      &
#                                                  mult_power  = (/mx,my,mz/),               &
#                                                  atom_mcent  = -1,                         &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(2),               &
#                                                  prefactor   = half,                       &
#                                                  ndisp_power = (/nz,nx,ny/),               &
#                                                  mult_power  = (/ny,nz,nx/),               &
#                                                  rderv_order = (/mz,mx,my/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/my,mz,mx/) )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(3),               &
#                                                  prefactor   = -half,                      &
#                                                  ndisp_power = (/nz,nx,ny/),               &
#                                                  mult_power  = (/ny,nz,nx/),               &
#                                                  rderv_order = (/my,mz,mx/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/mz,mx,my/) )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(4),               &
#                                                  prefactor   = -half,                      &
#                                                  ndisp_power = (/ny,nz,nx/),               &
#                                                  mult_power  = (/nz,nx,ny/),               &
#                                                  rderv_order = (/mz,mx,my/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/my,mz,mx/) )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(5),               &
#                                                  prefactor   = half,                       &
#                                                  ndisp_power = (/ny,nz,nx/),               &
#                                                  mult_power  = (/nz,nx,ny/),               &
#                                                  rderv_order = (/my,mz,mx/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom,                      &
#                                                  ederv_order = (/mz,mx,my/) )
#
#                  end if ! if XX, YY, ZZ components
#                end if ! if create operators
#              ! we have got all the required integrals, exit
#              else
#                exit
#              end if
#            end do ! loop over the x, y, z components of the magnetic field
#            ! we have got all the required integrals, exit
#            if ( int1.ge.this_prop%range_g1int(2) ) exit
#          end do ! loop over the nuclear magnetic moment coordinate
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !======================================================
#    ! Nuclear shielding integrals with common gauge origin
#    !======================================================
#    case('*NSTCGO')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 2
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals, the integrals are,
#      !
#      ! y_O*y_K/r_K^3+z_O*z_K/r_K^3, (001 x)
#      ! -x_O*y_K/r_K^3,              (001 y)
#      ! -x_O*z_K/r_K^3,              (001 z)
#      ! -y_O*x_K/r_K^3,              (002 x)
#      ! x_O*x_K/r_K^3+z_O*z_K/r_K^3, (002 y)
#      ! -y_O*z_K/r_K^3,              (002 z)
#      ! -z_O*x_K/r_K^3,              (003 x)
#      ! -z_O*y_K/r_K^3,              (003 y)
#      ! x_O*x_K/r_K^3+y_O*y_K/r_K^3, (003 z)
#      !
#      ! where K runs for all atoms
#      !
#      ! NOTE: DALTON uses X,Y,Z for K, 001, 002, ..., for O
#      else   
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,9*num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! initialize the index of integral
#        int1 = 0
#        ! initialize the index of the nuclear magnetic moment coordinate
#        iord1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over the nuclear magnetic moment coordinate
#          do jxyz = 1, 3
#            ! index of the nuclear magnetic moment coordinate
#            iord1 = iord1 + 1
#            ! loop over the x, y, z components of the magnetic field
#            do ixyz = 1, 3
#              ! index of integral
#              int1 = int1 + 1
#              ! this is not the integral we want
#              if ( int1.lt.this_prop%range_g1int(1) ) then
#                cycle
#              ! we get the required integral
#              else if ( int1.le.this_prop%range_g1int(2) ) then
#                ! label of the computed integrals
#                this_prop%g1prop_int1(int1)%label = xfigure(iord1/100)//              &
#                                                   xfigure(mod(iord1,100)/10)//      &
#                                                   xfigure(mod(mod(iord1,100),10))// &
#                                                   'NSCO'//xpower1(ixyz)//'_A'
#                ! create operators for the computed integrals
#                if ( new_op1 ) then
#                  ! XX, YY, ZZ components
#                  if ( ixyz.eq.jxyz ) then
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 2
#                    select case(ixyz)
#                    ! operators (y_O*y_K + z_O*z_K)/r_K^3
#                    case(1)
#                      nx=0; ny=1; nz=0
#                    ! operators (z_O*z_K + x_O*x_K)/r_K^3
#                    case(2)
#                      nx=0; ny=0; nz=1
#                    ! operators (x_O*x_K + y_O*y_K)/r_K^3
#                    case(3)
#                      nx=1; ny=0; nz=0
#                    end select
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(1),               &
#                                                  prefactor   = half,                       &
#                                                  mult_power  = (/nx,ny,nz/),               &
#                                                  coord_mcent = gaug_orig,                  &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(2),               &
#                                                  prefactor   = half,                       &
#                                                  mult_power  = (/nz,nx,ny/),               &
#                                                  coord_mcent = gaug_orig,                  &
#                                                  rderv_order = (/nz,nx,ny/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                  else
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 1
#                    ! nuclear magnetic moment coordinate
#                    select case(jxyz)
#                    ! operators -x_O*y_K/r_K^3, and -x_O*z_K/r_K^3
#                    case(1)
#                      mx=1; my=0; mz=0
#                    ! operators -y_O*x_K/r_K^3, and -y_O*z_K/r_K^3
#                    case(2)
#                      mx=0; my=1; mz=0
#                    ! operators -z_O*x_K/r_K^3, and -z_O*y_K/r_K^3
#                    case(3)
#                      mx=0; my=0; mz=1
#                    end select
#                    ! x, y, z components of the magnetic field
#                    select case(ixyz)
#                    ! operators -y_O*x_K/r_K^3, and -z_O*x_K/r_K^3
#                    case(1)
#                      nx=1; ny=0; nz=0
#                    ! operators -x_O*y_K/r_K^3, and -z_O*y_K/r_K^3
#                    case(2)
#                      nx=0; ny=1; nz=0
#                    ! operators -x_O*z_K/r_K^3, and -y_O*z_K/r_K^3
#                    case(3)
#                      nx=0; ny=0; nz=1
#                    end select
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(1),               &
#                                                  prefactor   = -half,                      &
#                                                  mult_power  = (/mx,my,mz/),               &
#                                                  coord_mcent = gaug_orig,                  &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                  end if ! if XX, YY, ZZ components
#                end if ! if create operators
#              ! we have got all the required integrals, exit
#              else
#                exit
#              end if
#            end do ! loop over the x, y, z components of the magnetic field
#            ! we have got all the required integrals, exit
#            if ( int1.ge.this_prop%range_g1int(2) ) exit
#          end do ! loop over the nuclear magnetic moment coordinate
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !===================================================================
#    ! London orbital contribution to nuclear shielding tensor integrals
#    !===================================================================
#    case('*NSTLON')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 3
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals, the integrals are,
#      !
#      ! -(Y_{MN}*z-Z_{MN}*y)*(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy), (001 x)
#      ! -(Z_{MN}*x-X_{MN}*z)*(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy), (001 y)
#      ! -(X_{MN}*y-Y_{MN}*x)*(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy), (001 z)
#      ! -(Y_{MN}*z-Z_{MN}*y)*(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz), (002 x)
#      ! -(Z_{MN}*x-X_{MN}*z)*(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz), (002 y)
#      ! -(X_{MN}*y-Y_{MN}*x)*(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz), (002 z)
#      ! -(Y_{MN}*z-Z_{MN}*y)*(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx), (003 x)
#      ! -(Z_{MN}*x-X_{MN}*z)*(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx), (003 y)
#      ! -(X_{MN}*y-Y_{MN}*x)*(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx), (003 z)
#      !
#      ! where K runs for all atoms
#      !
#      ! NOTE: DALTON uses X,Y,Z for N, 001, 002, ..., for K
#      else   
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,9*num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! initialize the index of integral
#        int1 = 0
#        ! initialize the index of the nuclear magnetic moment coordinate
#        iord1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over the nuclear magnetic moment coordinate
#          do jxyz = 1, 3
#            select case(jxyz)
#            ! operators -(y_K/r_K^3*d/dz - z_K/r_K^3*d/dy)
#            case(1)
#              mx=0; my=1; mz=0
#              nx=0; ny=0; nz=1
#            ! operators -(z_K/r_K^3*d/dx - x_K/r_K^3*d/dz)
#            case(2)
#              mx=0; my=0; mz=1
#              nx=1; ny=0; nz=0
#            ! operators -(x_K/r_K^3*d/dy - y_K/r_K^3*d/dx)
#            case(3)
#              mx=1; my=0; mz=0
#              nx=0; ny=1; nz=0
#            end select
#            ! index of the nuclear magnetic moment coordinate
#            iord1 = iord1 + 1
#            ! loop over the x, y, z components of the magnetic field
#            do ixyz = 1, 3
#              ! index of integral
#              int1 = int1 + 1
#              ! this is not the integral we want
#              if ( int1.lt.this_prop%range_g1int(1) ) then
#                cycle
#              ! we get the required integral
#              else if ( int1.le.this_prop%range_g1int(2) ) then
#                ! label of the computed integrals
#                this_prop%g1prop_int1(int1)%label = xfigure(iord1/100)//              &
#                                                   xfigure(mod(iord1,100)/10)//      &
#                                                   xfigure(mod(mod(iord1,100),10))// &
#                                                   'NSLO'//xpower1(ixyz)//'_N'
#                ! create operators for the computed integrals
#                if ( new_op1 ) then
#                  ! the number of operators
#                  this_prop%g1prop_int1(int1)%num_g1opt = 4
#                  select case(ixyz)
#                  ! operators Y_{MN}*z-Z_{MN}*y
#                  case(1)
#                    kx=0; ky=1; kz=0
#                    lx=0; ly=0; lz=1
#                  ! operators Z_{MN}*x-X_{MN}*z
#                  case(2)
#                    kx=0; ky=0; kz=1
#                    lx=1; ly=0; lz=0
#                  ! operators X_{MN}*y-Y_{MN}*x
#                  case(3)
#                    kx=1; ky=0; kz=0
#                    lx=0; ly=1; lz=0
#                  end select
#                  ! allocate operators
#                  allocate( this_prop%g1prop_int1(int1)%g1int_opt(4), stat=x_ierr )
#                  if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                  call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                              g1int_opt(1),               &
#                                                prefactor   = half,                       &
#                                                ndisp_power = (/kx,ky,kz/),               &
#                                                mult_power  = (/lx,ly,lz/),               &
#                                                rderv_order = (/mx,my,mz/),               &
#                                                r_func      = '1',                        &
#                                                atom_rcent  = iatom,                      &
#                                                ederv_order = (/nx,ny,nz/) )
#                  call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                              g1int_opt(2),               &
#                                                prefactor   = -half,                      &
#                                                ndisp_power = (/kx,ky,kz/),               &
#                                                mult_power  = (/lx,ly,lz/),               &
#                                                rderv_order = (/nx,ny,nz/),               &
#                                                r_func      = '1',                        &
#                                                atom_rcent  = iatom,                      &
#                                                ederv_order = (/mx,my,mz/) )
#                  call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                              g1int_opt(3),               &
#                                                prefactor   = -half,                      &
#                                                ndisp_power = (/lx,ly,lz/),               &
#                                                mult_power  = (/kx,ky,kz/),               &
#                                                rderv_order = (/mx,my,mz/),               &
#                                                r_func      = '1',                        &
#                                                atom_rcent  = iatom,                      &
#                                                ederv_order = (/nx,ny,nz/) )
#                  call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                              g1int_opt(4),               &
#                                                prefactor   = half,                       &
#                                                ndisp_power = (/lx,ly,lz/),               &
#                                                mult_power  = (/kx,ky,kz/),               &
#                                                rderv_order = (/nx,ny,nz/),               &
#                                                r_func      = '1',                        &
#                                                atom_rcent  = iatom,                      &
#                                                ederv_order = (/mx,my,mz/) )
#                end if ! if create operators
#              ! we have got all the required integrals, exit
#              else
#                exit
#              end if
#            end do ! loop over the x, y, z components of the magnetic field
#            ! we have got all the required integrals, exit
#            if ( int1.ge.this_prop%range_g1int(2) ) exit
#          end do ! loop over the nuclear magnetic moment coordinate
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !=================================================================
#    ! Nuclear shielding integrals without London orbital contribution
#    !=================================================================
#    case('*NSTNOL')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 2
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals, the integrals are,
#      !
#      ! y_N*y_K/r_K^3+z_N*z_K/r_K^3, (001 x)
#      ! -x_N*y_K/r_K^3,              (001 y)
#      ! -x_N*z_K/r_K^3,              (001 z)
#      ! -y_N*x_K/r_K^3,              (002 x)
#      ! x_N*x_K/r_K^3+z_N*z_K/r_K^3, (002 y)
#      ! -y_N*z_K/r_K^3,              (002 z)
#      ! -z_N*x_K/r_K^3,              (003 x)
#      ! -z_N*y_K/r_K^3,              (003 y)
#      ! x_N*x_K/r_K^3+y_N*y_K/r_K^3, (003 z)
#      !
#      ! where K runs for all atoms
#      !
#      ! NOTE: DALTON uses X,Y,Z for K, 001, 002, ..., for N
#      else   
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,9*num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! initialize the index of integral
#        int1 = 0
#        ! initialize the index of the nuclear magnetic moment coordinate
#        iord1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over the nuclear magnetic moment coordinate
#          do jxyz = 1, 3
#            ! index of the nuclear magnetic moment coordinate
#            iord1 = iord1 + 1
#            ! loop over the x, y, z components of the magnetic field
#            do ixyz = 1, 3
#              ! index of integral
#              int1 = int1 + 1
#              ! this is not the integral we want
#              if ( int1.lt.this_prop%range_g1int(1) ) then
#                cycle
#              ! we get the required integral
#              else if ( int1.le.this_prop%range_g1int(2) ) then
#                ! label of the computed integrals
#                this_prop%g1prop_int1(int1)%label = xfigure(iord1/100)//              &
#                                                   xfigure(mod(iord1,100)/10)//      &
#                                                   xfigure(mod(mod(iord1,100),10))// &
#                                                   'NSNL'//xpower1(ixyz)//'_N'
#                ! create operators for the computed integrals
#                if ( new_op1 ) then
#                  ! XX, YY, ZZ components
#                  if ( ixyz.eq.jxyz ) then
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 2
#                    select case(ixyz)
#                    ! operators (y_N*y_K + z_N*z_K)/r_K^3
#                    case(1)
#                      nx=0; ny=1; nz=0
#                    ! operators (z_N*z_K + x_N*x_K)/r_K^3
#                    case(2)
#                      nx=0; ny=0; nz=1
#                    ! operators (x_N*x_K + y_N*y_K)/r_K^3
#                    case(3)
#                      nx=1; ny=0; nz=0
#                    end select
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(1),               &
#                                                  prefactor   = half,                       &
#                                                  mult_power  = (/nx,ny,nz/),               &
#                                                  atom_mcent  = -1,                         &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(2),               &
#                                                  prefactor   = half,                       &
#                                                  mult_power  = (/nz,nx,ny/),               &
#                                                  atom_mcent  = -1,                         &
#                                                  rderv_order = (/nz,nx,ny/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                  else
#                    ! the number of operators
#                    this_prop%g1prop_int1(int1)%num_g1opt = 1
#                    ! nuclear magnetic moment coordinate
#                    select case(jxyz)
#                    ! operators -x_N*y_K/r_K^3, and -x_N*z_K/r_K^3
#                    case(1)
#                      mx=1; my=0; mz=0
#                    ! operators -y_N*x_K/r_K^3, and -y_N*z_K/r_K^3
#                    case(2)
#                      mx=0; my=1; mz=0
#                    ! operators -z_N*x_K/r_K^3, and -z_N*y_K/r_K^3
#                    case(3)
#                      mx=0; my=0; mz=1
#                    end select
#                    ! x, y, z components of the magnetic field
#                    select case(ixyz)
#                    ! operators -y_N*x_K/r_K^3, and -z_N*x_K/r_K^3
#                    case(1)
#                      nx=1; ny=0; nz=0
#                    ! operators -x_N*y_K/r_K^3, and -z_N*y_K/r_K^3
#                    case(2)
#                      nx=0; ny=1; nz=0
#                    ! operators -x_N*z_K/r_K^3, and -y_N*z_K/r_K^3
#                    case(3)
#                      nx=0; ny=0; nz=1
#                    end select
#                    ! allocate operators
#                    allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                    if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                    call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                                g1int_opt(1),               &
#                                                  prefactor   = -half,                      &
#                                                  mult_power  = (/mx,my,mz/),               &
#                                                  atom_mcent  = -1,                         &
#                                                  rderv_order = (/nx,ny,nz/),               &
#                                                  r_func      = '1',                        &
#                                                  atom_rcent  = iatom )
#                  end if ! if XX, YY, ZZ components
#                end if ! if create operators
#              ! we have got all the required integrals, exit
#              else
#                exit
#              end if
#            end do ! loop over the x, y, z components of the magnetic field
#            ! we have got all the required integrals, exit
#            if ( int1.ge.this_prop%range_g1int(2) ) exit
#          end do ! loop over the nuclear magnetic moment coordinate
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !================================
#    ! Potential energy at the nuclei
#    !================================
#    case('*NUCPOT')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 0
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! loop over atoms
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#!fixme: if there are more than 100 atoms
#          this_prop%g1prop_int1(int1)%label = 'POT.E_'//name_atoms(int1)(1:2)//'_S'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = 1
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                       g1int_opt(1),               &
#!fixme: without charge??
#                                          prefactor  = one,                        &
#                                          r_func     = '1',                        &
#                                          atom_rcent = int1 )
#          end if
#        end do ! loop over atoms
#      end if
#    !==========================
#    ! One-electron Hamiltonian
#    !==========================
#!fixme: maybe this operator is not efficient for calculations ...
#    case('*ONEHAMIL')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 0
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,1/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! label of the computed integrals
#        this_prop%g1prop_int1(1)%label = 'ONEHAMIL_S'
#        ! create operators for the computed integrals
#        if ( new_op1 ) then
#          ! the number of operators
#          this_prop%g1prop_int1(1)%num_g1opt = num_atoms + 3
#          ! allocate operators
#          allocate( this_prop%g1prop_int1(1)%g1int_opt &
#                    (this_prop%g1prop_int1(1)%num_g1opt), stat=x_ierr )
#          if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#          !-------------------------------
#          ! one-electron potential energy
#          !-------------------------------
#          ! loop over atoms
#          do iatom = 1, num_atoms
#            ! create operators
#            call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(1)% &
#                                                       g1int_opt(iatom),        &
#                                          prefactor  = -charge_atoms(iatom),    &
#                                          r_func     = '1',                     &
#                                          atom_rcent = iatom )
#          end do ! loop over atoms
#          !---------------------------
#          ! electronic kinetic energy
#          !---------------------------
#          ! initialize the index of operator
#          ixyz = num_atoms
#          ! loop over Cartesian derivatives
#          do nz = 0, 1
#            do ny = 0, 1-nz
#              ! index of operator
#              ixyz = ixyz + 1
#              ! generate Cartesian derivatives of X direction
#              nx = 1 - (ny+nz)
#              ! create operators
#              ! since we create first order derivatives, we need to multiply by 2
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(1)% &
#                                                          g1int_opt(ixyz),         &
#                                            prefactor   = -half,                   &
#                                            ederv_order = (/2*nx,2*ny,2*nz/) )
#            end do ! loop over Cartesian derivatives
#          end do ! loop over Cartesian derivatives
#        end if
#      end if
#    !===================
#    ! Overlap integrals
#    !===================
#    case('*OVERLAP')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,1/) )
#        ! label of the computed integrals
#        this_prop%g1prop_int1(1)%label = 'OVERLAP__S'
#        ! create operators for the computed integrals
#        if ( new_op1 ) then
#          ! the number of operators
#          this_prop%g1prop_int1(1)%num_g1opt = 1
#          ! allocate operators
#          allocate( this_prop%g1prop_int1(1)%g1int_opt(1), stat=x_ierr )
#          if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#          ! create operators
#          call gen1int_operator_create( new_op1  = this_prop%g1prop_int1(1)%g1int_opt(1) )
#        end if
#      end if
#    !===================================
#    ! Paramagnetic spin-orbit integrals
#    !===================================
#    case('*PSO')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '1'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = 2
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 3
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3*num_atoms/) )
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        this_prop%num_int1_cent = 3
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over angular momentum
#          do ixyz = 1, 3
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = 'PSO_'//xfigure(int1/100)         &
#                                                       //xfigure(mod(int1,100)/10) &
#                                                       //xfigure(mod(mod(int1,100),10))//'__A'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 2
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#!fixme: not -(y_k*d/dz - z_k*d/dy)/r_k^3??
#                select case(ixyz)
#                ! operators (y_k*d/dz - z_k*d/dy)/r_k^3
#                case(1)
#                  nx=0; ny=1; nz=0
#                  mx=0; my=0; mz=1
#                ! operators (z_k*d/dx - x_k*d/dz)/r_k^3
#                case(2)
#                  nx=0; ny=0; nz=1
#                  mx=1; my=0; mz=0
#                ! operators (x_k*d/dy - y_k*d/dx)/r_k^3
#                case(3)
#                  nx=1; ny=0; nz=0
#                  mx=0; my=1; mz=0
#                end select
#                ! create operators
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(1),               &
#                                              rderv_order = (/nx,ny,nz/),               &
#                                              r_func      = '1',                        &
#                                              atom_rcent  = iatom,                      &
#                                              ederv_order = (/mx,my,mz/) )
#                call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                            g1int_opt(2),               &
#                                              prefactor   = -one,                       &
#                                              rderv_order = (/mx,my,mz/),               &
#                                              r_func      = '1',                        & 
#                                              atom_rcent  = iatom,                      &
#                                              ederv_order = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over angular momentum
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !=============================================
#    ! First magnetic derivative of overlap matrix
#    !=============================================
#    case('*S1MAG')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3/) )
#        ! loop over integrals
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#          this_prop%g1prop_int1(int1)%label = 'dS/dB'//xpower1(int1)//'___A'
#          ! create operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = 2
#            ! allocate operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            select case(int1)
#            ! Y_{MN}*z-Z_{MN}*y
#            case(1)
#              kx=0; ky=1; kz=0
#              lx=0; ly=0; lz=1
#            ! Z_{MN}*x-X_{MN}*z
#            case(2)
#              kx=0; ky=0; kz=1
#              lx=1; ly=0; lz=0
#            ! X_{MN}*y-Y_{MN}*x
#            case(3)
#              kx=1; ky=0; kz=0
#              lx=0; ly=1; lz=0
#            end select
#            ! create operators
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(1),               &
#                                          prefactor   = half,                       &
#                                          ndisp_power = (/kx,ky,kz/),               &
#                                          mult_power  = (/lx,ly,lz/) )
#            call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                        g1int_opt(2),               &
#                                          prefactor   = -half,                      &
#                                          ndisp_power = (/lx,ly,lz/),               &
#                                          mult_power  = (/kx,ky,kz/) )
#          end if
#        end do ! loop over integrals
#      end if
#    !=========================
#    ! Second moment integrals
#    !=========================
#    case('*SECMOM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! order of multipole moments
#        iord1 = 2
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(iord1+1)*(iord1+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian powers
#        do nz = 0, iord1
#          do ny = 0, iord1-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian power of X direction
#              nx = iord1 - (ny+nz)
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = xpower2(int1)//'SECMOM_S'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                ! create operators
#                call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                           g1int_opt(1),               &
#                                              mult_power = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian powers
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian powers
#      end if
#    !=================================================================
#    ! Bra-half-derivative overlap integrals not to be antisymmetrized
#    !=================================================================
#!fixme: consider molecular symmetry
#    case('*SQHDOL')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3*num_atoms/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over Cartesian derivatives
#          do ixyz = 1, 3
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = 'SQHDL'//xfigure(int1/100)               &
#                                                        //xfigure(mod(int1,100)/10)       &
#                                                        //xfigure(mod(mod(int1,100),10))//'_N'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                select case (ixyz)
#                case (1)
#                  nx=1; ny=0; nz=0
#                case (2)
#                  nx=0; ny=1; nz=0
#                case (3)
#                  nx=0; ny=0; nz=1
#                end select
#                ! create operators
#                call gen1int_operator_create( new_op1   = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(1),               &
#                                              prefactor = -one,                       &
#                                              atom_bra  = iatom,                      &
#                                              gderv_bra = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian derivatives
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !=================================================================
#    ! Ket-half-derivative overlap integrals not to be antisymmetrized
#    !=================================================================
#!fixme: consider molecular symmetry
#    case('*SQHDOR')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,3*num_atoms/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over atoms
#        do iatom = 1, num_atoms
#          ! loop over Cartesian derivatives
#          do ixyz = 1, 3
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = 'SQHDR'//xfigure(int1/100)               &
#                                                        //xfigure(mod(int1,100)/10)       &
#                                                        //xfigure(mod(mod(int1,100),10))//'_N'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                select case (ixyz)
#                case (1)
#                  nx=1; ny=0; nz=0
#                case (2)
#                  nx=0; ny=1; nz=0
#                case (3)
#                  nx=0; ny=0; nz=1
#                end select
#                ! create operators
#                call gen1int_operator_create( new_op1   = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(1),               &
#                                              prefactor = -one,                       &
#                                              atom_ket  = iatom,                      &
#                                              gderv_ket = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian derivatives
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over atoms
#      end if
#    !======================================
#    ! Traceless theta quadrupole integrals
#    !======================================
#    case('*THETA')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals, the integrals are,
#      ! XX: ( 2x^2 - y^2 - z^2 )/2,
#      ! YY: ( 2y^2 - x^2 - z^2 )/2,
#      ! ZZ: ( 2z^2 - x^2 - y^2 )/2,
#      ! XY: 3xy/2, YZ: 3yz/2, and XZ: 3xz/2
#      else
#        ! order of multipole moments
#        iord1 = 2
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(iord1+1)*(iord1+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian powers
#        do nz = 0, iord1
#          do ny = 0, iord1-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian power of X direction
#              nx = iord1 - (ny+nz)
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = xpower2(int1)//'THETA__S'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! XX, YY, and ZZ components
#                if ( nx.eq.2 .or. ny.eq.2 .or. nz.eq.2 ) then
#                  ! the number of operators
#                  this_prop%g1prop_int1(int1)%num_g1opt = 3
#                  ! prefactor
#                  prefact1 = one
#                  ! allocate operators
#                  allocate( this_prop%g1prop_int1(int1)%g1int_opt(3), stat=x_ierr )
#                  if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                ! XY, XZ, and YZ components
#                else
#                  ! the number of operators
#                  this_prop%g1prop_int1(int1)%num_g1opt = 1
#                  ! prefactor
#                  prefact1 = second3
#                  ! allocate operators
#                  allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                  if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                end if
#                ! create operator x^2, y^2, z^2, 3xy/2, 3yz/2, or 3xz/2
#                call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                           g1int_opt(1),               &
#                                              prefactor  = prefact1,                   &
#                                              mult_power = (/nx,ny,nz/) )
#                ! extra operators for XX, YY, and ZZ components
#                if ( nx.eq.2 .or. ny.eq.2 .or. nz.eq.2 ) then
#                  ! create operator -0.5x^2, -0.5y^2, or -0.5z^2
#                  call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                             g1int_opt(2),               &
#                                                prefactor  = -half,                      &
#                                                mult_power = (/ny,nz,nx/) )
#                  ! create operator -0.5x^2, -0.5y^2, or -0.5z^2
#                  call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                             g1int_opt(3),               &
#                                                prefactor  = -half,                      &
#                                                mult_power = (/nz,nx,ny/) )
#                end if
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian powers
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian powers
#      end if
#    !=========================
#    ! Third moments integrals
#    !=========================
#    case('*THIRDM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      ! we are going to create the computed integrals
#      else
#        ! order of multipole moments
#        iord1 = 3
#        ! allocate computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,(iord1+1)*(iord1+2)/2/) )
#        ! initialize the index of integral
#        int1 = 0
#        ! loop over Cartesian powers
#        do nz = 0, iord1
#          do ny = 0, iord1-nz
#            ! index of integral
#            int1 = int1 + 1
#            ! this is not the integral we want
#            if ( int1.lt.this_prop%range_g1int(1) ) then
#              cycle
#            ! we get the required integral
#            else if ( int1.le.this_prop%range_g1int(2) ) then
#              ! generate Cartesian power of X direction
#              nx = iord1 - (ny+nz)
#              ! label of the computed integrals
#              this_prop%g1prop_int1(int1)%label = xpower3(int1)//'_3MOM_S'
#              ! create operators for the computed integrals
#              if ( new_op1 ) then
#                ! the number of operators
#                this_prop%g1prop_int1(int1)%num_g1opt = 1
#                ! allocate operators
#                allocate( this_prop%g1prop_int1(int1)%g1int_opt(1), stat=x_ierr )
#                if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#                ! create operators
#                call gen1int_operator_create( new_op1    = this_prop%g1prop_int1(int1)% &
#                                                           g1int_opt(1),               &
#                                              mult_power = (/nx,ny,nz/) )
#              end if
#            ! we have got all the required integrals, exit
#            else
#              exit
#            end if
#          end do ! loop over Cartesian powers
#          ! we have got all the required integrals, exit
#          if ( int1.ge.this_prop%range_g1int(2) ) exit
#        end do ! loop over Cartesian powers
#      end if
#    !===================================================
#    !> 2NDMM: Second moment integrals in momentum space
#    !===================================================
#    case('*2NDMM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      !> There are 6 integrals in total,
#      !> \f{eqnarray}{
#      !>   x\mathrm{d}x + \mathrm{d}x x & = & 2x\mathrm{d}x + 1,
#      !>   \\x\mathrm{d}y + \mathrm{d}y x,
#      !>   \\x\mathrm{d}z + \mathrm{d}x z,
#      !>   \\y\mathrm{d}y + \mathrm{d}y y & = & 2y\mathrm{d}y + 1,
#      !>   \\y\mathrm{d}z + \mathrm{d}y z,
#      !>   \\z\mathrm{d}z + \mathrm{d}z z & = & 2z\mathrm{d}z + 1.
#      !> \f}
#      else
#        ! allocates computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,6/) )
#        ! loops over integrals
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#          this_prop%g1prop_int1(int1)%label = xpower2(int1)//'2NDMM__A'
#          ! creates operators for the computed integrals
#          if ( new_op1 ) then
#            ! the number of operators
#            this_prop%g1prop_int1(int1)%num_g1opt = 2
#            ! allocates operators
#            allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#            if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#            if ( int1 == 1 .or. int1 == 3 .or. int1 == 6 ) then
#              select case (int1)
#              !> operators \f$2x\mathrm{d}x + 1\f$
#              case (1)
#                nx=1; ny=0; nz=0
#              !> operators \f$2y\mathrm{d}y + 1\f$
#              case (3)
#                nx=0; ny=1; nz=0
#              !> operators \f$2z\mathrm{d}z + 1\f$
#              case (6)
#                nx=0; ny=0; nz=1
#              end select
#              ! creates operators
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(1),               &
#                                            prefactor   = two,                        &
#                                            mult_power  = (/nx,ny,nz/),               &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/nx,ny,nz/) )
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(2) )
#            else
#              select case (int1)
#              !> operators \f$x\mathrm{d}y + \mathrm{d}y x\f$
#              case (2)
#                nx=1; ny=0; nz=0
#                mx=0; my=1; mz=0
#              !> operators \f$x\mathrm{d}z + \mathrm{d}x z\f$
#              case (4)
#                nx=1; ny=0; nz=0
#                mx=0; my=0; mz=1
#              !> operators \f$y\mathrm{d}z + \mathrm{d}y z\f$
#              case (5)
#                nx=0; ny=1; nz=0
#                mx=0; my=0; mz=1
#              end select
#              ! creates operators
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(1),               &
#                                            mult_power  = (/nx,ny,nz/),               &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/mx,my,mz/) )
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(2),               &
#                                            mult_power  = (/mx,my,mz/),               &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/nx,ny,nz/) )
#            end if
#          end if
#        end do ! loops over integrals
#      end if
#    !===========================================
#    ! Third moments integrals in momentum space
#    !===========================================
#    case('*3RDMM')
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = 'X'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#        ! number of centers belong to this property (+2 by consider the centers of basis sets)
#        if ( present(prop_num_cent) ) prop_num_cent = 2
#      !> there are 10 integrals in total,
#      !> \f{eqnarray}{
#      !>   x^2\mathrm{d}x + x\mathrm{d}x x + \mathrm{d}x x^2
#      !>   & = & 3x^2\mathrm{d}x + 3x,
#      !>   \\x^2\mathrm{d}y + x\mathrm{d}x y + \mathrm{d}x xy
#      !>   & = & x^2\mathrm{d}y + 2xy\mathrm{d}x + y,
#      !>   \\x^2\mathrm{d}z + x\mathrm{d}x z + \mathrm{d}x xz
#      !>   & = & x^2\mathrm{d}z + 2xz\mathrm{d}x + z,
#      !>   \\xy\mathrm{d}y + x\mathrm{d}y y + \mathrm{d}x y^2
#      !>   & = & y^2\mathrm{d}x + 2xy\mathrm{d}y + x,
#      !>   \\xy\mathrm{d}z + x\mathrm{d}y z + \mathrm{d}x yz
#      !>   & = & xy\mathrm{d}z + xz\mathrm{d}y + yz\mathrm{d}x,
#      !>   \\xz\mathrm{d}z + x\mathrm{d}z z + \mathrm{d}x z^2
#      !>   & = & z^2\mathrm{d}x + 2xz\mathrm{d}z + x,
#      !>   \\y^2\mathrm{d}y + y\mathrm{d}y y + \mathrm{d}y y^2
#      !>   & = & 3y^2\mathrm{d}y + 3y,
#      !>   \\y^2\mathrm{d}z + y\mathrm{d}y z + \mathrm{d}y yz
#      !>   & = & y^2\mathrm{d}z + 2yz\mathrm{d}y + z,
#      !>   \\yz\mathrm{d}z + y\mathrm{d}z z + \mathrm{d}y z^2
#      !>   & = & z^2\mathrm{d}y + 2yz\mathrm{d}z + y,
#      !>   \\z^2\mathrm{d}z + z\mathrm{d}z z + \mathrm{d}z z^2
#      !>   & = & 3z^2\mathrm{d}z + 3z.
#      !> \f}
#      else
#        ! allocates computed integrals
#        call g1int_property_alloc_int1( this_prop, (/1,10/) )
#        ! loops over integrals
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          ! label of the computed integrals
#          this_prop%g1prop_int1(int1)%label = xpower3(int1)//'3RDMM_A'
#          ! creates operators for the computed integrals
#          if ( new_op1 ) then
#            !> operators \f$xy\mathrm{d}z + xz\mathrm{d}y + yz\mathrm{d}x\f$
#            if ( int1 == 6 ) then
#              ! the number of operators
#              this_prop%g1prop_int1(int1)%num_g1opt = 3
#              ! allocates operators
#              allocate( this_prop%g1prop_int1(int1)%g1int_opt(3), stat=x_ierr )
#              if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#              ! creates operators
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(1),               &
#                                            mult_power  = (/1,1,0/),                  &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/0,0,1/) )
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(2),               &
#                                            mult_power  = (/1,0,1/),                  &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/0,1,0/) )
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(3),               &
#                                            mult_power  = (/0,1,1/),                  &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/1,0,0/) )
#            else if ( int1 == 1 .or. int1 == 4 .or. int1 == 10 ) then
#              ! the number of operators
#              this_prop%g1prop_int1(int1)%num_g1opt = 2
#              ! allocates operators
#              allocate( this_prop%g1prop_int1(int1)%g1int_opt(2), stat=x_ierr )
#              if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#              select case (int1)
#              !> operators \f$3x^2\mathrm{d}x + 3x\f$
#              case (1)
#                nx=2; ny=0; nz=0
#                mx=1; my=0; mz=0
#              !> operators \f$3y^2\mathrm{d}y + 3y\f$
#              case (4)
#                nx=0; ny=2; nz=0
#                mx=0; my=1; mz=0
#              !> operators \f$3z^2\mathrm{d}z + 3z\f$
#              case (10)
#                nx=0; ny=0; nz=2
#                mx=0; my=0; mz=1
#              end select
#              ! creates operators
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(1),               &
#                                            prefactor   = three,                      &
#                                            mult_power  = (/nx,ny,nz/),               &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/mx,my,mz/) )
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(2),               &
#                                            prefactor   = three,                      &
#                                            mult_power  = (/mx,my,mz/),               &
#                                            coord_mcent = elec_orig )
#            else
#              ! the number of operators
#              this_prop%g1prop_int1(int1)%num_g1opt = 3
#              ! allocates operators
#              allocate( this_prop%g1prop_int1(int1)%g1int_opt(3), stat=x_ierr )
#              if ( x_ierr.ne.0 ) call xstop('allocate g1int_opt!')
#              select case (int1) 
#              !> operators \f$x^2\mathrm{d}y + 2xy\mathrm{d}x + y\f$
#              case (2)
#                nx=1; ny=0; nz=0
#                mx=0; my=1; mz=0
#              !> operators \f$y^2\mathrm{d}x + 2xy\mathrm{d}y + x\f$
#              case (3)
#                nx=0; ny=1; nz=0
#                mx=1; my=0; mz=0
#              !> operators \f$x^2\mathrm{d}z + 2xz\mathrm{d}x + z\f$
#              case (5)
#                nx=1; ny=0; nz=0
#                mx=0; my=0; mz=1
#              !> operators \f$y^2\mathrm{d}z + 2yz\mathrm{d}y + z\f$
#              case (7)
#                nx=0; ny=1; nz=0
#                mx=0; my=0; mz=1
#              !> operators \f$z^2\mathrm{d}x + 2xz\mathrm{d}z + x\f$
#              case (8)
#                nx=0; ny=0; nz=1
#                mx=1; my=0; mz=0
#              !> operators \f$z^2\mathrm{d}y + 2yz\mathrm{d}z + y\f$
#              case (9)
#                nx=0; ny=0; nz=1
#                mx=0; my=1; mz=0
#              end select
#              ! creates operators
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(1),               &
#                                            mult_power  = (/2*nx,2*ny,2*nz/),         &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/mx,my,mz/) )
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(2),               &
#                                            prefactor   = two,                        &
#                                            mult_power  = (/nx+mx,ny+my,nz+mz/),      &
#                                            coord_mcent = elec_orig,                  &
#                                            ederv_order = (/nx,ny,nz/)  )
#              call gen1int_operator_create( new_op1     = this_prop%g1prop_int1(int1)% &
#                                                          g1int_opt(3),               &
#                                            mult_power  = (/mx,my,mz/),               &
#                                            coord_mcent = elec_orig )
#            end if
#          end if
#        end do ! loops over integrals
#      end if
#    !================================
#    ! error, property is not defined
#    !================================
#    case default
#      if ( chk_prop ) then
#        ! the type of r_C function
#        prop_rfun = '?'
#        ! order of auxilary function needed by this property
#        ! -1 means this property does not need auxiliary function
#        prop_aux_order = -1
#      else
#        write(x_lupri,1000)'Oops! I have no idea about the property: '//trim(chk_name)
#        call xstop(trim(chk_name)//' is not recognized!')
#      end if
#    end select
#    ! if we are asked to create the operators for the integrals of the property
#    if ( new_op1 ) then
#      ! operators have been created
#      this_prop%set_opt1 = .false.
#      ! dump the operators in readable format
#      if ( xlib_lprint.ge.5 ) then
#        write(x_lupri,'()')
#        write(x_lupri,1000)'Property: '//trim(chk_name)
#        write(x_lupri,1000)'Computed integrals: ',this_prop%range_g1int
#        ! loop over the computed integrals
#        do int1 = this_prop%range_g1int(1), this_prop%range_g1int(2)
#          write(x_lupri,1000)'Label of integral: '//trim(this_prop%g1prop_int1(int1)%label)
#          call gen1int_operator_dump( this_prop%g1prop_int1(int1)%g1int_opt )
#        end do
#      end if
#    end if
#    ! pop the stack
#    call xsub_exit
#    return
#1000  format('PRCR ',A,2I6)
#  end subroutine g1int_property_set
#
#  ! allocate one-electron integrals for a property
#  subroutine g1int_property_alloc_int1( this_prop, prop_range_int1 )
#    implicit none
#    ! property
#    type(g1prop), intent(inout) :: this_prop
#    ! the range of integrals for this property
#    integer, intent(in) :: prop_range_int1(2)
#    !////////////////////////////////////////////////////////////////
#    !                       local variables
#    !////////////////////////////////////////////////////////////////
#    ! incremental recorder of integrals
#    integer int1
#    ! push current subroutine into the stack
#    call xsub_enter('g1int_property_alloc_int1')
#    ! the range of computed integrals
#    if ( allocated(this_prop%range_g1int) ) then
#      this_prop%range_g1int(1) = max( this_prop%range_g1int(1), prop_range_int1(1) )
#      ! if this_prop%range_g1int(2) <= 0, we compute until the last integral
#      if ( this_prop%range_g1int(2).le.0 ) then
#        this_prop%range_g1int(2) = prop_range_int1(2)
#      else
#        this_prop%range_g1int(2) = min( this_prop%range_g1int(2), prop_range_int1(2) )
#      end if
#    else
#      allocate( this_prop%range_g1int(2), stat=x_ierr )
#      if ( x_ierr.ne.0 ) call xstop('allocate range_g1int!')
#      this_prop%range_g1int = prop_range_int1
#    end if
#    ! check the validity of the range of computed integrals
#    if ( this_prop%range_g1int(1).gt.this_prop%range_g1int(2) ) then
#      write(x_lupri,1000)'Range of computed integrals: ',this_prop%range_g1int
#      call xstop('Incorrect range of computed integrals!')
#    end if
#    ! allocate computed integrals now ...
#    if ( allocated(this_prop%g1prop_int1) ) deallocate(this_prop%g1prop_int1)
#    allocate( this_prop%g1prop_int1(this_prop%range_g1int(1):this_prop%range_g1int(2)), stat=x_ierr )
#    if ( x_ierr.ne.0 ) call xstop('Faild to allocate computed integrals!')
#    ! pop the stack
#    call xsub_exit
#    return
#1000  format('PRI1 ',A,2I8)
#  end subroutine g1int_property_alloc_int1
#
#end module g1int_property
#
