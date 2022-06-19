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

module pde_utils

#ifdef VAR_PDE
    use hdf5
#endif

    use pelib_precision

    implicit none

    private

    ! public subroutines/functions
    public :: pde_save_density, pde_twoints
    public :: pde_get_num_core_nuclei, pde_get_fragment_density

contains

#ifdef VAR_PDE

subroutine pde_save_density(ao_denmat, ew_denmat, num_frag_bas)

    use pelib_mpi
    use pelib_options
    use pelib_constants
    use pelib_potential_derivatives
    use pelib_integral_interfaces
    use pelib_blas_interfaces

    real(rp), dimension(:), intent(in) :: ao_denmat
    real(rp), dimension(:,:), intent(in) :: ew_denmat
    integer(ip), intent(in) :: num_frag_bas

    integer(ip) :: i
    integer(ip) :: num_core_nuclei
    real(rp) :: nucel_energy
    real(rp), dimension(:), allocatable :: Zcore_ints
    real(rp), dimension(:,:), allocatable :: Rcore, Zcore
    real(rp), dimension(:,:,:), allocatable :: temp_ints
    real(rp), dimension(:), allocatable :: fd_static_field
    real(rp), dimension(:,:), allocatable :: Ftmp

    !hdf5 IO variables
    integer(4) :: error
    integer(HSIZE_T), dimension(1) :: dim_1d
    integer(HSIZE_T), dimension(2) :: dim_2d
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: group_id
    integer(HID_T) :: space_id

    ! get root id and number of MPI processes
#if defined(VAR_MPI)
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)
    master = 0
#else
    myid = 0
    nprocs = 1
    master = 0
#endif

    ! initialize variables
    ndens = 1
    nbas = int(num_frag_bas, ip)
    nnbas = nbas * (nbas + 1) / 2
    n2bas = nbas * nbas

    if (nprocs == 1) then
        site_start = 1
        site_finish = nsites
#if defined(VAR_MPI)
    else
        call mpi_bcast(ndens, 1, impi, master, comm, ierr)
        call mpi_bcast(nbas, 1, impi, master, comm, ierr)
        call mpi_bcast(nnbas, 1, impi, master, comm, ierr)
        call mpi_bcast(n2bas, 1, impi, master, comm, ierr)
        call mpi_bcast(ao_denmat(1), nnbas, rmpi, master, comm, ierr)
        call mpi_sync()
#endif
    end if

    ! get electric field from fragment density at polarizable sites
    ! TODO: better solution for neglecting polarization
    ! TODO: potential from fragment density at surface points
    if (pelib_polar .and. npols > 0) then
        allocate(fd_static_field(3*npols))
        allocate(Ftmp(3*npols,1))
        Ftmp = 0.0
        call electron_fields(Ftmp, ao_denmat)
        fd_static_field = Ftmp(:,1)
        Ftmp = 0.0
        call nuclear_fields(Ftmp(:,1))
        fd_static_field = fd_static_field + Ftmp(:,1)
        deallocate(Ftmp)
    end if

    if (myid == master) then
        call h5open_f(error)
        call h5fopen_f(trim(h5pdefile), H5F_ACC_RDWR_F, file_id, error)
        call h5gopen_f(file_id, 'core_fragment', group_id, error)
        call h5dopen_f(group_id, 'num_nuclei', dset_id, error)
        dim_1d = 1
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_core_nuclei, dim_1d, error)
        call h5dclose_f(dset_id, error)
        allocate(Zcore(1,num_core_nuclei), Rcore(3,num_core_nuclei))
        call h5dopen_f(group_id, 'coordinates', dset_id, error)
        dim_2d = [3, num_core_nuclei]
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Rcore, dim_2d, error)
        call h5dclose_f(dset_id, error)
        call h5dopen_f(group_id, 'charges', dset_id, error)
        dim_2d = [1, num_core_nuclei]
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Zcore, dim_2d, error)
        call h5dclose_f(dset_id, error)
        ! remember that we assume everything is in bohr!

        ! calculate nuclear - electron energy contribution
        allocate(Zcore_ints(nbas*(nbas+1)/2))
        allocate(temp_ints(nbas*(nbas+1)/2,1,1))
        Zcore_ints = 0.0
        do i = 1, num_core_nuclei
            call Tk_integrals('potential_derivative', Rcore(:,i), temp_ints)
            Zcore_ints = Zcore_ints + Zcore(1,i) * temp_ints(:,1,1)
        end do
        nucel_energy = dot(ao_denmat, Zcore_ints)
        deallocate(Rcore, Zcore, Zcore_ints, temp_ints)

        dim_1d = 1
        call h5screate_simple_f(1_4, dim_1d, space_id, error)
        call h5dcreate_f(group_id, 'nuclear-electron energy', H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, nucel_energy, dim_1d, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)
        call h5gclose_f(group_id, error)

        call h5gopen_f(file_id, 'fragment', group_id, error)

        dim_1d = 1
        call h5screate_simple_f(1_4, dim_1d, space_id, error)
        call h5dcreate_f(group_id, 'num_bas', H5T_NATIVE_INTEGER, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nbas, dim_1d, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

        dim_1d = nbas*(nbas+1)/2
        call h5screate_simple_f(1_4, dim_1d, space_id, error)
        call h5dcreate_f(group_id, 'density matrix', H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ao_denmat, dim_1d, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

        dim_2d = [nbas, nbas]
        call h5screate_simple_f(2_4, dim_2d, space_id, error)
        call h5dcreate_f(group_id, 'energy-weighted density matrix', H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ew_denmat, dim_2d, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

        dim_1d = 1
        call h5screate_simple_f(1_4, dim_1d, space_id, error)
        call h5dcreate_f(group_id, 'num_pols', H5T_NATIVE_INTEGER, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, npols, dim_1d, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

        dim_1d = 3*npols
        call h5screate_simple_f(1_4, dim_1d, space_id, error)
        call h5dcreate_f(group_id, 'electric fields', H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fd_static_field, dim_1d, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

        call h5gclose_f(group_id, error)
        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end if

    if (allocated(fd_static_field)) deallocate(fd_static_field)

end subroutine pde_save_density

subroutine pde_twoints(core_fckmat, intmol_overlap, num_dimer_bas)

    use pelib_options

    real(rp), intent(in), dimension(:) :: core_fckmat
    real(rp), intent(in), dimension(:,:) :: intmol_overlap
    integer(ip), intent(in) :: num_dimer_bas

    integer(ip) :: i, j, k
    integer(ip) :: num_frag_bas
    integer(ip) :: num_core_bas
    real(rp), dimension(:), allocatable :: packed_repmat
    real(rp), dimension(:,:), allocatable :: full_repmat
    real(rp), dimension(:,:), allocatable :: ew_denmat

    !hdf5 IO variables
    integer(4) :: error
    integer(HSIZE_T), dimension(1) :: dim_1d
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: group_id
    integer(HID_T) :: space_id

    call h5open_f(error)
    call h5fopen_f(trim(h5pdefile), H5F_ACC_RDWR_F, file_id, error)
    call h5gopen_f(file_id, 'fragment', group_id, error)
    call h5dopen_f(group_id, 'num_bas', dset_id, error)
    dim_1d = 1
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_frag_bas, dim_1d, error)
    call h5dclose_f(dset_id, error)
    allocate(ew_denmat(num_frag_bas,num_frag_bas))
    call h5dopen_f(group_id, 'energy-weighted density matrix', dset_id, error)
    dim_1d = num_frag_bas*(num_frag_bas+1)/2
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, ew_denmat, dim_1d, error)
    call h5dclose_f(dset_id, error)
    call h5gclose_f(group_id, error)

    num_core_bas = int(num_dimer_bas, ip) - num_frag_bas

    !open interface, file and get group id
    call h5gopen_f(file_id, 'core_fragment', group_id, error)

    ! store number of basis functions
    dim_1d = 1
    call h5screate_simple_f(1_4, dim_1d, space_id, error)
    call h5dcreate_f(group_id, 'num_bas', H5T_NATIVE_INTEGER, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, num_core_bas, dim_1d, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    ! store fckmat
    dim_1d = num_core_bas*(num_core_bas+1)/2
    call h5screate_simple_f(1_4, dim_1d, space_id, error)
    call h5dcreate_f(group_id, 'electrostatic matrix', H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, core_fckmat, dim_1d, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    allocate(full_repmat(num_core_bas,num_core_bas))
    full_repmat = - matmul(matmul(intmol_overlap, ew_denmat), transpose(intmol_overlap))
    deallocate(ew_denmat)
    allocate(packed_repmat(num_core_bas*(num_core_bas+1)/2))
    k = 1
    do i = 1, num_core_bas
        do j = 1, i
            packed_repmat(k) = full_repmat(j,i)
            k = k + 1
        end do
    end do
    deallocate(full_repmat)

    ! store exchange-repulsion matrix
    dim_1d = num_core_bas*(num_core_bas+1)/2
    call h5screate_simple_f(1_4, dim_1d, space_id, error)
    call h5dcreate_f(group_id, 'exchange-repulsion matrix', H5T_NATIVE_DOUBLE, space_id, dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, packed_repmat, dim_1d, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(space_id, error)

    deallocate(packed_repmat)

    ! close group, file, and interface
    call h5gclose_f(group_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)

end subroutine pde_twoints

subroutine pde_get_fragment_density(denmat, num_bas)

    use pelib_options

    real(rp), dimension(:), allocatable, intent(out) :: denmat
    integer(ip), intent(out) :: num_bas

    !hdf5 IO variables
    integer(4) :: error
    integer(HSIZE_T), dimension(1) :: dim_1d
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: group_id

    call h5open_f(error)
    call h5fopen_f(trim(h5pdefile), H5F_ACC_RDWR_F, file_id, error)
    call h5gopen_f(file_id, 'fragment', group_id, error)
    call h5dopen_f(group_id, 'num_bas', dset_id, error)
    dim_1d = 1
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_bas, dim_1d, error)
    call h5dclose_f(dset_id, error)
    allocate(denmat(num_bas*(num_bas+1)/2))
    call h5dopen_f(group_id, 'density matrix', dset_id, error)
    dim_1d = num_bas*(num_bas+1)/2
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, denmat, dim_1d, error)
    call h5dclose_f(dset_id, error)
    call h5gclose_f(group_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)

end subroutine pde_get_fragment_density

function pde_get_num_core_nuclei() result(num_nuclei)

    use pelib_options

    integer(ip) :: num_nuclei

    !hdf5 IO variables
    integer(4) :: error
    integer(HSIZE_T), dimension(1) :: dim_1d
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: group_id

    call h5open_f(error)
    call h5fopen_f(trim(h5pdefile), H5F_ACC_RDWR_F, file_id, error)
    call h5gopen_f(file_id, 'core_fragment', group_id, error)
    call h5dopen_f(group_id, 'num_nuclei', dset_id, error)
    dim_1d = 1
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_nuclei, dim_1d, error)
    call h5dclose_f(dset_id, error)
    call h5gclose_f(group_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)

end function pde_get_num_core_nuclei

#else

subroutine pde_save_density(ao_denmat, ew_denmat, num_frag_bas)

    use pelib_mpi
    use pelib_options
    use pelib_constants
    use pelib_potential_derivatives
    use pelib_integral_interfaces
    use pelib_blas_interfaces

    real(rp), dimension(:), intent(in) :: ao_denmat
    real(rp), dimension(:,:), intent(in) :: ew_denmat
    integer(ip), intent(in) :: num_frag_bas

    error stop 'Not compiled with PDE support (recompile with -DENABLE_PDE=ON)'

end subroutine pde_save_density

subroutine pde_twoints(core_fckmat, intmol_overlap, num_dimer_bas)

    use pelib_options

    real(rp), intent(in), dimension(:) :: core_fckmat
    real(rp), intent(in), dimension(:,:) :: intmol_overlap
    integer(ip), intent(in) :: num_dimer_bas

    error stop 'Not compiled with PDE support (recompile with -DENABLE_PDE=ON)'

end subroutine pde_twoints

subroutine pde_get_fragment_density(denmat, num_bas)

    use pelib_options

    real(rp), dimension(:), allocatable, intent(out) :: denmat
    integer(ip), intent(out) :: num_bas

    error stop 'Not compiled with PDE support (recompile with -DENABLE_PDE=ON)'

end subroutine pde_get_fragment_density

function pde_get_num_core_nuclei() result(num_nuclei)

    use pelib_options

    integer(ip) :: num_nuclei

    error stop 'Not compiled with PDE support (recompile with -DENABLE_PDE=ON)'

end function pde_get_num_core_nuclei
#endif

end module pde_utils
