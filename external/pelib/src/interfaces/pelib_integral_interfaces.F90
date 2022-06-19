!
!   Polarizable Embedding (PE) library
!   Copyright (C) 2013-2018, 2020 The PE library developers. See the CONTRIBUTORS file
!                                 in the top-level directory of this distribution.
!
!   This file is part of the PE library.
!
!   The PE library is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License as
!   published by the Free Software Foundation, either version 3 of the
!   License, or (at your option) any later version.
!
!   The PE library is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
!
!   Contact information:
!
!   Jogvan Magnus Haugaard Olsen
!   E-mail: foeroyingur@gmail.com
!

module pelib_integral_interfaces

    use pelib_precision

    implicit none

    private

    public :: init_pelib_integral_interfaces, finalize_pelib_integral_interfaces
    public :: Tk_integrals, multipole_integrals

    integer(ip), save :: nnucs
    real(rp), dimension(:), allocatable, save :: core_com
    real(rp), dimension(3), parameter :: zero_origin = [0.0, 0.0, 0.0]
    real(rp), dimension(:,:), allocatable, save :: core_coords

contains

#if defined(DUMMY_INTEGRALS)

subroutine init_pelib_integral_interfaces(com, coords)

    real(rp), dimension(:), intent(in) :: com
    real(rp), dimension(:,:), intent(in) :: coords

    print *, com
    print *, coords
    error stop

end subroutine init_pelib_integral_interfaces

subroutine finalize_pelib_integral_interfaces()

    error stop

end subroutine finalize_pelib_integral_interfaces

subroutine Tk_integrals(inttype, coord, Tk_ints)

    character(*), intent(in) :: inttype
    real(rp), dimension(:), intent(in) :: coord
    real(rp), dimension(:,:,:), intent(out) :: Tk_ints

    print *, inttype
    print *, coord
    print *, Tk_ints
    error stop

end subroutine Tk_integrals

subroutine multipole_integrals(multipole_origin, ints)

    real(rp), dimension(:), intent(in) :: multipole_origin
    real(rp), dimension(:,:,:), intent(out) :: ints

    print *, multipole_origin
    print *, ints
    error stop

end subroutine multipole_integrals

#elif defined(GEN1INT_INTEGRALS)

subroutine init_pelib_integral_interfaces(com, coords)

    real(rp), dimension(:), intent(in) :: com
    real(rp), dimension(:,:), intent(in) :: coords

    if ((size(com) /= 3) .or. (size(coords, 1) /= 3)) then
        error stop 'wrong input size in init_pelib_integral_interfaces'
    end if

    nnucs = size(coords, 2)

    if (.not. allocated(core_com)) allocate(core_com(3))
    if (.not. allocated(core_coords)) allocate(core_coords(3,nnucs))

    core_com = com
    core_coords = coords

end subroutine init_pelib_integral_interfaces

!-----------------------------------------------------------------------------

subroutine finalize_pelib_integral_interfaces()

    deallocate(core_com, core_coords)

end subroutine finalize_pelib_integral_interfaces

!-----------------------------------------------------------------------------

subroutine Tk_integrals(inttype, coord, Tk_ints)

    use gen1int
    use gen1int_api
    use gen1int_matrix

    ! in-/output variables
    character(*), intent(in) :: inttype
    real(rp), dimension(:), intent(in) :: coord
    real(rp), dimension(:,:,:), intent(out) :: Tk_ints

    ! variables related to Gen1Int
    integer(ip) :: io = 6
    integer(ip) :: printlvl = 0
    integer(ip) :: num_ao, num_prop, prop_sym, ierr
    integer(ip) :: num_geo_total, num_geo_bra, num_geo_ket
    ! non-zero components for the operator, the first dimension is for bra and
    ! ket sub-shells, the last is the number of non-zero components, which should
    ! be 1 for non-relativistic calculations
#if defined (PRG_DIRAC)
    integer(ip) nnz_comp(2,2)
#else
    integer(ip) nnz_comp(2,1)
#endif
    type(one_prop_t) :: prop_operator
    type(nary_tree_t) :: nary_tree_total, nary_tree_bra, nary_tree_ket
    type(matrix), dimension(:), allocatable :: intmats
    logical(ip) :: triangular, symmetric

    ! variables related to in-/output from PElib
    integer(ip) :: nbas, nnbas, ncomps, nprops
    integer(ip) :: i, j, k, l, x, y, z
    integer(ip), dimension(1) :: idx
    integer(ip), dimension(:,:), allocatable :: row2col
    real(rp), dimension(1) :: sign_factor
    real(rp), dimension(3,1) :: site_coord

    Tk_ints = 0.0

    site_coord(:,1) = coord

    nnbas = size(Tk_ints, 1)
    nbas = int(0.5 * (sqrt(1.0 + 8.0 * nnbas) - 1.0))
    ncomps = size(Tk_ints, 2)
    k = int(0.5 * (sqrt(1.0 + 8.0 * ncomps) - 1.0)) - 1
    nprops = size(Tk_ints, 3)

    if (inttype == 'potential_derivative' .or. inttype == 'gaussian_potential') then
        if (nprops /= 1) error stop 'array size inconsistency 1'
    else if (inttype == 'molecular_gradient') then
        if (nprops /= 3 * nnucs) error stop 'array size inconsistency 2'
    else if (inttype == 'magnetic_gradient') then
        if (nprops /= 3) error stop 'array size inconsistency 3'
    end if

    if (mod(k, 2) == 0) then
        sign_factor = - 1.0
    else if (mod(k, 2) /= 0) then
        sign_factor = 1.0
    end if

    if (inttype == 'magnetic_gradient') then
        sign_factor = - sign_factor
    end if

    idx = -1
    if (inttype == 'potential_derivative' .or.&
       &inttype == 'molecular_gradient' .or.&
       &inttype == 'magnetic_gradient') then
        call OnePropCreate(prop_name=INT_POT_ENERGY, &
                           one_prop=prop_operator, &
                           info_prop=ierr, &
                           idx_nuclei=idx, &
                           coord_nuclei=site_coord, &
                           charge_nuclei=sign_factor, &
                           order_geo_pot=k)
!    else if (inttype == 'gaussian') then
!        call OnePropCreate(prop_name=INT_GAUSSIAN_POT, &
!                           one_prop=prop_operator, &
!                           info_prop=ierr, &
!                           idx_gauorg=idx, &
!                           gaupot_origin=site_coord, &
!                           gaupot_charge=sign_factor, &
!                           gaupot_expt=gauexp, &
!                           order_geo_pot=k)
    end if
    if (ierr /= 0) then
        error stop 'failed to create property operator in Tk_integrals'
    end if

    if (inttype == 'magnetic_gradient') then
        call OnePropSetMag(one_prop=prop_operator, &
                           order_mag=1, &
                           order_mag_bra=0, &
                           order_mag_ket=0)
        call OnePropSetLAO(one_prop=prop_operator, &
                           gauge_origin=core_com, &
                           origin_London_PF=zero_origin)
    end if

    ! creates N-ary tree for geometric derivatives
    idx = 0
    call Gen1IntAPINaryTreeCreate(max_num_cent=0, &
                                  order_geo=0, &
                                  num_geo_atoms=0, &
                                  idx_geo_atoms=idx, &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0, &
                                  order_geo=0, &
                                  num_geo_atoms=0, &
                                  idx_geo_atoms=idx, &
                                  nary_tree=nary_tree_ket)
    if (inttype == 'molecular_gradient') then
        call Gen1IntAPINaryTreeCreate(max_num_cent=1, &
                                      order_geo=1, &
                                      num_geo_atoms=0, &
                                      idx_geo_atoms=idx, &
                                      nary_tree=nary_tree_total)
    else
        call Gen1IntAPINaryTreeCreate(max_num_cent=0, &
                                      order_geo=0, &
                                      num_geo_atoms=0,&
                                      idx_geo_atoms=idx, &
                                      nary_tree=nary_tree_total)
    end if

    ! gets the number of property integrals and their symmetry
    call OnePropGetNumProp(one_prop=prop_operator, num_prop=num_prop)
    call OnePropGetSymmetry(one_prop=prop_operator, prop_sym=prop_sym)

    ! gets the number of geometric derivatives
    call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)

    ! updates the number and symmetry of property integrals
    num_prop = num_prop * num_geo_bra * num_geo_ket * num_geo_total

    if (inttype == 'potential_derivative' .or. inttype == 'gaussian_potential') then
        nprops = 1
    else if (inttype == 'molecular_gradient') then
        nprops = 3 * nnucs
    else if (inttype == 'magnetic_gradient') then
        nprops = 3
    end if
    if (num_prop /= nprops * ncomps) error stop 'wrong number of components'


    call Gen1IntAPIGetNumAO(num_ao=num_ao)
    if (num_ao /= nbas) error stop 'array size inconsistency 4'

    select case(prop_sym)
        case(SYMM_INT_MAT)
            triangular = .true.
            symmetric = (prop_sym == SYMM_INT_MAT)
        case(ANTI_INT_MAT)
            triangular = .true.
            symmetric = (prop_sym == ANTI_INT_MAT)
        case(SQUARE_INT_MAT)
            triangular = .false.
            symmetric = .false.
        case default
            triangular = .false.
            symmetric = .false.
    end select
    if ((.not. symmetric) .or. (.not. triangular)) then
        error stop 'array property inconsistency'
    end if

    allocate(intmats(num_prop), stat=ierr)
    if (ierr /= 0) error stop 'failed to allocate matrices'

    allocate(row2col(3,ncomps), stat=ierr)
    if (ierr /= 0) error stop 'failed to allocate matrices'

    ! help array for conversion between Gen1Int (xx xy yy ...) and PElib (xx xy xz ..)
    i = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z > k .or. x+y+z < k) cycle
                row2col(1,i) = x
                row2col(2,i) = y
                row2col(3,i) = z
                i = i + 1
            end do
        end do
    end do

    if (inttype == 'magnetic_gradient') then
        i = 1
        do z = 0, k
            do y = 0, k
                do x = 0, k
                    if (x+y+z /= k) cycle
                    do j = 1, ncomps
                        if (x == row2col(1,j) .and.&
                            y == row2col(2,j) .and.&
                            z == row2col(3,j)) then
                            do l = 1, nprops
                                call MatAssociate(work_alpha=Tk_ints(:,j,l),&
                                                  num_row=nbas, &
                                                  A=intmats(i), &
                                                  info_mat=ierr, &
                                                  triangular=triangular, &
                                                  symmetric=symmetric)
                                if (ierr /= 0) error stop 'failed to associate matrices'
                                i = i + 1
                            end do
                        end if
                    end do
                end do
            end do
        end do
    else
        i = 1
        do l = 1, nprops
            do z = 0, k
                do y = 0, k
                    do x = 0, k
                        if (x+y+z /= k) cycle
                        do j = 1, ncomps
                            if (x == row2col(1,j) .and.&
                                y == row2col(2,j) .and.&
                                z == row2col(3,j)) then
                                call MatAssociate(work_alpha=Tk_ints(:,j,l), &
                                                  num_row=nbas, &
                                                  A=intmats(i), &
                                                  info_mat=ierr, &
                                                  triangular=triangular, &
                                                  symmetric=symmetric)
                                if (ierr /= 0) error stop 'failed to associate matrices'
                                i = i + 1
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end if

    ! sets the non-zero components for the one-electron operator
    ! (large,large)
    nnz_comp(1,1) = 1
    nnz_comp(2,1) = 1
#if defined (PRG_DIRAC)
    ! (small,small)
    nnz_comp(1,2) = 2
    nnz_comp(2,2) = 2
#endif

    ! calculates the integrals, please refer to the comments in subroutine
    ! \fn(Gen1IntOnePropGetIntExpt) in gen1int_api.F90
    call Gen1IntOnePropGetIntExpt(nnz_comp=nnz_comp, &
                                  one_prop=prop_operator, &
                                  nary_tree_bra=nary_tree_bra, &
                                  nary_tree_ket=nary_tree_ket, &
                                  nary_tree_total=nary_tree_total, &
                                  num_ints=num_prop, &
                                  val_ints=intmats, &
                                  num_dens=1, &
                                  io_viewer=io, &
                                  level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)

    do i = 1, num_prop
        call MatNullify(A=intmats(i))
    end do

    deallocate(intmats)

end subroutine Tk_integrals

subroutine multipole_integrals(multipole_origin, ints)

    use gen1int
    use gen1int_api
    use gen1int_matrix

    real(rp), dimension(3), intent(in) :: multipole_origin
    real(rp), dimension(:,:,:), intent(out) :: ints

    ! variables related to Gen1Int
    integer(ip) :: io = 6
    integer(ip) :: printlvl = 0
    integer(ip) :: num_ao, num_prop, prop_sym, ierr
    integer(ip) :: num_geo_total, num_geo_bra, num_geo_ket
    ! non-zero components for the operator, the first dimension is for bra and
    ! ket sub-shells, the last is the number of non-zero components, which should
    ! be 1 for non-relativistic calculations
#if defined (PRG_DIRAC)
    integer(ip) nnz_comp(2,2)
#else
    integer(ip) nnz_comp(2,1)
#endif
    type(one_prop_t) :: prop_operator
    type(nary_tree_t) :: nary_tree_total, nary_tree_bra, nary_tree_ket
    type(matrix), dimension(:), allocatable :: intmats
    logical(ip) :: triangular, symmetric

    ! variables related to in-/output from PElib
    integer(ip) :: nbas, nnbas, ncomps, nprops
    integer(ip) :: i, j, k, l, x, y, z
    integer(ip), dimension(1) :: idx
    integer(ip), dimension(:,:), allocatable :: row2col

    ints = 0.0

    nnbas = size(ints, 1)
    nbas = int(0.5 * (sqrt(1.0 + 8.0 * nnbas) - 1.0))
    ncomps = size(ints, 2)
    k = int(0.5 * (sqrt(1.0 + 8.0 * ncomps) - 1.0)) - 1
    nprops = size(ints, 3)

    if (nprops /= 1) then
        error stop "Invalid dimension in multipole integral"
    end if

    call OnePropCreate(prop_name=INT_CART_MULTIPOLE, &
                       one_prop=prop_operator, &
                       info_prop=ierr, &
                       dipole_origin=multipole_origin, &
                       order_mom=k)

    ! creates N-ary tree for geometric derivatives
    idx = 0
    call Gen1IntAPINaryTreeCreate(max_num_cent=0, &
                                  order_geo=0, &
                                  num_geo_atoms=0, &
                                  idx_geo_atoms=idx, &
                                  nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0, &
                                  order_geo=0, &
                                  num_geo_atoms=0, &
                                  idx_geo_atoms=idx, &
                                  nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeCreate(max_num_cent=0, &
                                  order_geo=0, &
                                  num_geo_atoms=0,&
                                  idx_geo_atoms=idx, &
                                  nary_tree=nary_tree_total)

    ! gets the number of property integrals and their symmetry
    call OnePropGetNumProp(one_prop=prop_operator, num_prop=num_prop)
    call OnePropGetSymmetry(one_prop=prop_operator, prop_sym=prop_sym)

    ! gets the number of geometric derivatives
    call NaryTreeGetNumGeo(nary_tree=nary_tree_bra, num_unique_geo=num_geo_bra)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_ket, num_unique_geo=num_geo_ket)
    call NaryTreeGetNumGeo(nary_tree=nary_tree_total, num_unique_geo=num_geo_total)

    ! updates the number and symmetry of property integrals
    num_prop = num_prop * num_geo_bra * num_geo_ket * num_geo_total

    if (num_prop /= nprops * ncomps) error stop 'wrong number of components'


    call Gen1IntAPIGetNumAO(num_ao=num_ao)
    if (num_ao /= nbas) error stop 'array size inconsistency 4'

    select case(prop_sym)
        case(SYMM_INT_MAT)
            triangular = .true.
            symmetric = (prop_sym == SYMM_INT_MAT)
        case(ANTI_INT_MAT)
            triangular = .true.
            symmetric = (prop_sym == ANTI_INT_MAT)
        case(SQUARE_INT_MAT)
            triangular = .false.
            symmetric = .false.
        case default
            triangular = .false.
            symmetric = .false.
    end select

    if ((.not. symmetric) .or. (.not. triangular)) then
        error stop 'array property inconsistency'
    end if

    allocate(intmats(num_prop), stat=ierr)
    if (ierr /= 0) error stop 'failed to allocate matrices'

    allocate(row2col(3,ncomps), stat=ierr)
    if (ierr /= 0) error stop 'failed to allocate matrices'

    ! help array for conversion between Gen1Int (xx xy yy ...) and PElib (xx xy xz ..)
    i = 1
    do x = k, 0, -1
        do y = k, 0, -1
            do z = k, 0, -1
                if (x+y+z > k .or. x+y+z < k) cycle
                row2col(1,i) = x
                row2col(2,i) = y
                row2col(3,i) = z
                i = i + 1
            end do
        end do
    end do

    i = 1
    do l = 1, nprops
        do z = 0, k
            do y = 0, k
                do x = 0, k
                    if (x+y+z /= k) cycle
                    do j = 1, ncomps
                        if (x == row2col(1,j) .and.&
                            y == row2col(2,j) .and.&
                            z == row2col(3,j)) then
                            call MatAssociate(work_alpha=ints(:,j,l), &
                                              num_row=nbas, &
                                              A=intmats(i), &
                                              info_mat=ierr, &
                                              triangular=triangular, &
                                              symmetric=symmetric)
                            if (ierr /= 0) error stop 'failed to associate matrices'
                            i = i + 1
                        end if
                    end do
                end do
            end do
        end do
    end do

    ! sets the non-zero components for the one-electron operator
    ! (large,large)
    nnz_comp(1,1) = 1
    nnz_comp(2,1) = 1
#if defined (PRG_DIRAC)
    ! (small,small)
    nnz_comp(1,2) = 2
    nnz_comp(2,2) = 2
#endif

    ! calculates the integrals, please refer to the comments in subroutine
    ! \fn(Gen1IntOnePropGetIntExpt) in gen1int_api.F90
    call Gen1IntOnePropGetIntExpt(nnz_comp=nnz_comp, &
                                  one_prop=prop_operator, &
                                  nary_tree_bra=nary_tree_bra, &
                                  nary_tree_ket=nary_tree_ket, &
                                  nary_tree_total=nary_tree_total, &
                                  num_ints=num_prop, &
                                  val_ints=intmats, &
                                  num_dens=1, &
                                  io_viewer=io, &
                                  level_print=printlvl)

    call OnePropDestroy(one_prop=prop_operator)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_bra)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_ket)
    call Gen1IntAPINaryTreeDestroy(nary_tree=nary_tree_total)

    do i = 1, num_prop
        call MatNullify(A=intmats(i))
    end do

    deallocate(intmats)

end subroutine multipole_integrals

#endif

end module pelib_integral_interfaces
