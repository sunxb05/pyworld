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
module pelib_potential_derivatives

    use pelib_precision

    implicit none

    private

    public :: electron_potentials, electron_fields, electron_fields_fmm, electron_fields_multipole
    public :: nuclear_potentials, nuclear_fields
    public :: multipole_potentials, multipole_fields
    public :: fragment_density_fields

contains

!------------------------------------------------------------------------------

subroutine electron_potentials(Vels, denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces, only: Tk_integrals
    use pelib_blas_interfaces, only: dot

    real(rp), dimension(:,:), intent(out) :: Vels
    real(rp), dimension(:), intent(in) :: denmats

    integer(ip) :: site
    integer(ip) :: i, j, l, m
    real(rp), dimension(:,:,:), allocatable :: Vel_ints

    allocate(Vel_ints(nnbas,1,1))

    i = 1
    do site = surp_start, surp_finish
        call Tk_integrals('potential_derivative', Rsp(:,site), Vel_ints)
        do j = 1, ndens
            l = (j - 1) * nnbas + 1
            m = j * nnbas
            Vels(i,j) = dot(denmats(l:m), Vel_ints(:,1,1))
        end do
        i = i + 1
    end do

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + surpdists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vels(1,i), surpdists, displs, rmpi, master, comm, ierr)
        end do
    else if (myid /= master) then
        do i = 1, ndens
            call mpi_gatherv(Vels(1,i), surpdists(myid), rmpi, dummy_real, surpdists, displs, rmpi, master, comm, ierr)
        end do
    end if
#endif

    if (pelib_debug) then
         do i = 1, nsurp
             write (luout,*) 'i, Vels(i)' , i, Vels(i,:)
         end do
    end if

    deallocate(Vel_ints)

end subroutine electron_potentials

!------------------------------------------------------------------------------

subroutine electron_fields(Fels, denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces
    use pelib_blas_interfaces, only: nrm2, dot

    real(rp), dimension(:,:), intent(out) :: Fels
    real(rp), dimension(:), intent(in) :: denmats

    logical(lp) :: skip
    integer(ip) :: site
    integer(ip) :: i, j, k, l, m
    real(rp) :: factor
    real(rp), dimension(3) :: Rms
    real(rp), dimension(:,:,:), allocatable :: Fel_ints

    allocate(Fel_ints(nnbas,3,1))

    i = 0
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        if (pelib_savden) then
            skip = .false.
            do j = 1, nnucs
                Rms = Rs(:,site) - Rm(:,j)
                if (nrm2(Rms) <= 1.0) skip = .true.
            end do
            if (skip) then
                i = i + 3
                cycle
            end if
        end if
        call Tk_integrals('potential_derivative', Rs(:,site), Fel_ints)
        do j = 1, 3
            do k = 1, ndens
                l = (k - 1) * nnbas + 1
                m = k * nnbas
                Fels(i+j,k) = dot(denmats(l:m), Fel_ints(:,j,1))
            end do
        end do
        if (pelib_core_damp) then
            call get_damping_factor(Rs(:,site), P11s(:,site), factor)
            do k = 1, ndens
                Fels(i+1:i+3,k) = factor * Fels(i+1:i+3,k)
            end do
        endif
        i = i + 3
    end do

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + poldists(i-1)
        end do
        do i = 1, ndens
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fels(1,i), poldists, displs, rmpi, master, comm, ierr)
        end do
    else if (myid /= master) then
        do i = 1, ndens
            call mpi_gatherv(Fels(1,i), poldists(myid), rmpi, dummy_real, poldists, displs, rmpi, master, comm, ierr)
        end do
    end if
#endif

    deallocate(Fel_ints)
end subroutine electron_fields

subroutine electron_fields_fmm(Fels, denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces
    use pelib_blas_interfaces, only: nrm2, dot
    use fmm, only: get_boxes, interpolate_electron_fields

    real(rp), dimension(:,:), intent(out) :: Fels
    real(rp), dimension(:), intent(in) :: denmats

    logical(lp) :: skip
    integer(ip) :: site
    integer(ip) :: i, j, k, l, m, idx, field_order, multipole_size, ncomp, component, p
    real(rp) :: factor
    real(rp), dimension(3) :: Rms
    real(rp), dimension(:,:,:), allocatable :: Fel_ints
    real(rp), dimension(:,:), allocatable :: fmm_coordinates, box_coordinates
    real(rp), dimension(:,:,:), allocatable :: box_local_expansion, fmm_field
    integer(ip) :: iwork, work_start, work_stop, work_size
    integer(ip), dimension(:), allocatable :: box_indices

    allocate(fmm_coordinates(npols, 3), source=0.0_rp)
    allocate(fmm_field(npols, 4, ndens))
    l = 1
    do i = 1, nsites
        if (zeroalphas(i)) cycle
        fmm_coordinates(l, :) = Rs(:, i)
        l = l + 1
    end do

    call get_boxes(fmm_coordinates, transpose(Rm), qmmm_interaction_theta, fmm_ncrit, qmmm_interaction_order, idx_near_field, box_coordinates, box_indices)

    ! near field
    allocate(Fel_ints(nnbas,3,1))
    work_size = size(idx_near_field) / nprocs
    work_start = 1 + myid * work_size
    work_stop = (myid + 1) * work_size
    if (myid == nprocs - 1) work_stop = size(idx_near_field)

    if (incore_nearfield .and. .not. allocated(near_field_integrals)) then
        allocate(near_field_integrals(nnbas, 3, 1, work_start:work_stop))
        do iwork = work_start, work_stop
            idx = idx_near_field(iwork)
            call Tk_integrals('potential_derivative', fmm_coordinates(idx, :), near_field_integrals(:, 1:3, 1:1, iwork))
        end do
    end if

    do iwork = work_start, work_stop
        idx = idx_near_field(iwork)
        if (incore_nearfield) then
            Fel_ints(:,:,1) = near_field_integrals(:, :, 1, iwork)
        else
            call Tk_integrals('potential_derivative', fmm_coordinates(idx,:), Fel_ints)
        end if
        do j = 1, 3
            do k = 1, ndens
                l = (k - 1) * nnbas + 1
                m = k * nnbas
                Fels(3*(idx-1)+j,k) = dot(denmats(l:m), Fel_ints(:,j,1))
            end do
        end do
        if (pelib_core_damp) then
            error stop "Not implemented"
        endif
    end do
    deallocate(Fel_ints)
#if defined(VAR_MPI)
    if (nprocs > 1) then
        if (myid == master) then
            call mpi_reduce(MPI_IN_PLACE, Fels(1, 1), size(Fels), rmpi, MPI_SUM, 0, comm, ierr)
        else
            call mpi_reduce(Fels(1, 1), Fels(1, 1), size(Fels), rmpi, MPI_SUM, 0, comm, ierr)
        end if
    end if
#endif


    ! get local expansion from fields due to electrons
    multipole_size = (qmmm_interaction_order + 1) * (qmmm_interaction_order + 2) * (qmmm_interaction_order + 3) / 6
    work_size = size(box_coordinates, 1) / nprocs
    work_start = 1 + myid * work_size
    work_stop = (myid + 1) * work_size
    if (myid == nprocs - 1) work_stop = size(box_coordinates, 1)
    allocate(box_local_expansion(size(box_coordinates,1), multipole_size, ndens))
    if (.not. allocated(Fel_ints_fmm)) then 
        allocate(Fel_ints_fmm(nnbas,multipole_size,1, work_start:work_stop))
        do iwork = work_start, work_stop
            l = 1
            do p = 0, qmmm_interaction_order
                ncomp = (p+1)*(p+2)/2
                call Tk_integrals('potential_derivative', box_coordinates(iwork, :), Fel_ints_fmm(:, l:l+ncomp-1, 1:1, iwork))
                Fel_ints_fmm(:, l:l+ncomp-1, 1:1, iwork) = Fel_ints_fmm(:, l:l+ncomp-1, 1:1, iwork)
                l = l + ncomp
            end do
        end do
    end if
    box_local_expansion = 0.0d0
    do iwork = work_start, work_stop
        j = 1
        do p = 0, qmmm_interaction_order
            ncomp = (p+1)*(p+2)/2
            do component = 1, ncomp
                do k = 1, ndens
                    l = (k - 1) * nnbas + 1
                    m = k * nnbas
                    box_local_expansion(iwork,j,k) = dot(denmats(l:m), Fel_ints_fmm(:,j,1, iwork)) * (-1)**(p+1)
                end do
                j = j + 1
            end do
        end do
    end do
#ifdef VAR_MPI
    call mpi_allreduce(MPI_IN_PLACE, box_local_expansion(1, 1, 1), size(box_local_expansion), rmpi, MPI_SUM, comm, ierr)
#endif

    ! interpolate local expansion
    field_order = 1
    fmm_field = 0.0d0
    call interpolate_electron_fields(comm, fmm_coordinates, box_local_expansion, box_indices, qmmm_interaction_theta, fmm_ncrit, qmmm_interaction_order, 1, fmm_field)

    if (myid == master) then
        do k = 1, ndens
            do i = 0, npols - 1
                Fels(3*i+1, k) = Fels(3*i+1, k) + fmm_field(i+1, 2, k)
                Fels(3*i+2, k) = Fels(3*i+2, k) + fmm_field(i+1, 3, k)
                Fels(3*i+3, k) = Fels(3*i+3, k) + fmm_field(i+1, 4, k)
            end do
        end do
    end if
end subroutine electron_fields_fmm

subroutine electron_fields_multipole(Fels, denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces
    use pelib_blas_interfaces, only: nrm2, dot
    use pelib_interaction_tools

    real(rp), dimension(:,:), intent(out) :: Fels
    real(rp), dimension(:), intent(in) :: denmats

    logical(lp) :: skip
    integer(ip) :: site
    integer(ip) :: i, j, k, l, m, idx, field_order, multipole_size, ncomp, component, p
    real(rp) :: factor
    real(rp), dimension(:,:,:), allocatable :: Fel_ints
    real(rp), dimension(:,:), allocatable :: box_coordinates
    real(rp), dimension(:,:,:), allocatable :: multipoles, multipole_field ! [N x multipole_length x ndens]
    integer(ip) :: iwork, work_start, work_stop, work_size
    integer(ip) :: num_operators

    if (.not. allocated(multipole_interaction_operators)) then
        call build_multipole_interaction_operators
    end if

    ! near field
    allocate(Fel_ints(nnbas,3,1))
    work_size = size(idx_near_field) / nprocs
    work_start = 1 + myid * work_size
    work_stop = (myid + 1) * work_size
    if (myid == nprocs - 1) work_stop = size(idx_near_field)

    if (incore_nearfield .and. .not. allocated(near_field_integrals)) then
        allocate(near_field_integrals(nnbas, 3, 1, work_start:work_stop))
        do iwork = work_start, work_stop
            call Tk_integrals('potential_derivative', mm_coordinates_near(:, iwork), near_field_integrals(:, 1:3, 1:1, iwork))
        end do
    end if

    do iwork = work_start, work_stop
        if (incore_nearfield) then
            Fel_ints(:,:,1) = near_field_integrals(:, :, 1, iwork)
        else
            call Tk_integrals('potential_derivative', mm_coordinates_near(:, iwork), Fel_ints)
        end if
        do j = 1, 3
            do k = 1, ndens
                l = (k - 1) * nnbas + 1
                m = k * nnbas
                Fels(3*(idx_near_field(iwork)-1)+j,k) = dot(denmats(l:m), Fel_ints(:,j,1))
            end do
        end do
        if (pelib_core_damp) then
            error stop "Not implemented"
        endif
    end do
    deallocate(Fel_ints)
#if defined(VAR_MPI)
    if (nprocs > 1) then
        if (myid == master) then
            call mpi_reduce(MPI_IN_PLACE, Fels(1, 1), size(Fels), rmpi, MPI_SUM, 0, comm, ierr)
        else
            call mpi_reduce(Fels(1, 1), Fels(1, 1), size(Fels), rmpi, MPI_SUM, 0, comm, ierr)
        end if
    end if
#endif
    ! get multipoles from QM density
    if (size(idx_far_field) > 0) then 
        num_operators = size(multipole_interaction_operators, 3)
        multipole_size = size(multipole_interaction_operators, 2)
        allocate(multipoles(num_operators, multipole_size, ndens))
        do i = 1, num_operators
            do j = 1, multipole_size
                do k = 1, ndens
                    l = (k - 1) * nnbas + 1
                    m = k * nnbas
                    multipoles(i, j, k) = -dot(denmats(l:m), multipole_interaction_operators(:, j, i))
                end do
            end do
        end do
        l = 1
        ! evaluate field due to QM multipoles
        allocate(multipole_field(size(idx_far_field), 4, ndens))
        field_order = 1
        call get_multipole_field(multipole_field, multipole_operator_centers, mm_coordinates_far, multipoles, field_order)
        if (myid == master) then
            do k = 1, ndens
                do i = 1, size(idx_far_field)
                    Fels(3*(idx_far_field(i)-1)+1, k) = multipole_field(i, 2, k)
                    Fels(3*(idx_far_field(i)-1)+2, k) = multipole_field(i, 3, k)
                    Fels(3*(idx_far_field(i)-1)+3, k) = multipole_field(i, 4, k)
                end do
            end do
        end if
        deallocate(multipoles)
        deallocate(multipole_field)
    end if
end subroutine electron_fields_multipole

!------------------------------------------------------------------------------

subroutine nuclear_potentials(Vnucs)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: Tk_tensor

    real(rp), dimension(:), intent(out) :: Vnucs

    logical(lp) :: lexist
    integer(ip) :: lu, site
    integer(ip) :: i, j
    real(rp), dimension(3) :: Rmsp
    real(rp), dimension(1) :: Tmsp

    if (myid == master) then
        inquire(file='pelib_nuclear_potential.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist) then
        if (myid == master) then
            call openfile('pelib_nuclear_potential.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Vnucs
            close(lu)
        end if
    else
        Vnucs = 0.0
        i = 1
        do site = surp_start, surp_finish
            do j = 1, nnucs
                Rmsp = Rsp(:,site) - Rm(:,j)
                call Tk_tensor(Tmsp, Rmsp)
                Vnucs(i) = Vnucs(i) + Zm(1,j) * Tmsp(1)
            end do
            i = i + 1
        end do
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + surpdists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vnucs(1), surpdists, displs, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_gatherv(Vnucs(1), surpdists(myid), rmpi, dummy_real, surpdists, displs, rmpi, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            call openfile('pelib_nuclear_potential.bin', lu, 'new', 'unformatted')
            rewind(lu)
            write(lu) Vnucs
            close(lu)
        end if
    end if

    if (pelib_debug) then
        do i=1, nsurp
            write (luout,*) 'Vnucs(i)' ,i, Vnucs(i)
        end do
    end if

end subroutine nuclear_potentials

!------------------------------------------------------------------------------

subroutine nuclear_fields(Fnucs)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: Tk_tensor
    use pelib_blas_interfaces, only: nrm2

    real(rp), dimension(:), intent(out) :: Fnucs

    logical(lp) :: lexist, skip
    integer(ip) :: lu, site
    integer(ip) :: i, j, k
    real(rp) :: factor
    real(rp), dimension(3) :: Rms, Tms

    if (myid == master) then
        inquire(file='pelib_nuclear_field.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist .and. ((scf_cycle > 1) .or. pelib_restart)) then
        if (myid == master) then
            call openfile('pelib_nuclear_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fnucs
            close(lu)
        end if
    else
        Fnucs = 0.0
        i = 0
        do site = site_start, site_finish
            if (zeroalphas(site)) cycle
            if (pelib_savden) then
                skip = .false.
                do j = 1, nnucs
                    Rms = Rs(:,site) - Rm(:,j)
                    if (nrm2(Rms) <= 1.0) skip = .true.
                end do
                if (skip) then
                    i = i + 3
                    cycle
                end if
            end if
            do j = 1, nnucs
                Rms = Rs(:,site) - Rm(:,j)
                call Tk_tensor(Tms, Rms)
                do k = 1, 3
                    Fnucs(i+k) = Fnucs(i+k) - Zm(1,j) * Tms(k)
                end do
            end do
            if (pelib_core_damp) then
                call get_damping_factor(Rs(:,site), P11s(:,site), factor)
                Fnucs(i+1:i+3) = factor * Fnucs(i+1:i+3)
            end if
            i = i + 3
        end do
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + poldists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fnucs(1), poldists, displs, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_gatherv(Fnucs(1), poldists(myid), rmpi, dummy_real, poldists, displs, rmpi, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            call openfile('pelib_nuclear_field.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Fnucs
            close(lu)
        end if
    end if

end subroutine nuclear_fields

!------------------------------------------------------------------------------

#ifdef VAR_PDE

subroutine fragment_density_fields(Ffd)

    use hdf5

    use pelib_options

    real(rp), dimension(:), intent(out) :: Ffd

    integer(ip) :: num_fd_pols

    !hdf5 IO variables
    integer(4) :: error
    integer(HSIZE_T), dimension(1) :: dim_1d
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id

    Ffd = 0.0

    call h5open_f(error)
    call h5fopen_f(trim(h5pdefile), H5F_ACC_RDWR_F, file_id, error)

    call h5dopen_f(file_id, 'num_fields', dset_id, error)
    dim_1d = 1
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_fd_pols, dim_1d, error)
    call h5dclose_f(dset_id, error)

    if (num_fd_pols /= npols) then
        error stop 'array length inconsistency'
    end if

    call h5dopen_f(file_id, 'electric fields', dset_id, error)
    dim_1d = 3 * npols
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Ffd, dim_1d, error)
    call h5dclose_f(dset_id, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

end subroutine fragment_density_fields

#else

subroutine fragment_density_fields(Ffd)

    use pelib_options

    real(rp), dimension(:), intent(out) :: Ffd

    error stop 'Not compiled with PDE support (recompile with -DENABLE_PDE=ON)'

end subroutine fragment_density_fields

#endif

!------------------------------------------------------------------------------

subroutine multipole_potentials(Vmuls)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: multipole_derivative

    real(rp), dimension(:), intent(out) :: Vmuls

    logical(lp) :: lexist
    integer(ip) :: lu
    integer(ip) :: i, j, k
    real(rp), dimension(1) :: Fi
    real(rp), dimension(3) :: Rji

    if (myid == master) then
        inquire(file='pelib_multipole_potential.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist .and. ((scf_cycle > 1) .or. pelib_restart)) then
        if (myid == master) then
            call openfile('pelib_multipole_potential.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Vmuls
            close(lu)
        end if
    else
        Vmuls = 0.0
        k = 1
        do i = 1, nsurp
            do j = site_start, site_finish
                Fi = 0.0
                Rji = Rsp(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (maxval(abs(M0s(:,j))) >= zero) then
                        call multipole_derivative(Fi, Rji, M0s(:,j))
                    end if
                end if
                if (lmul(1)) then
                    if (maxval(abs(M1s(:,j))) >= zero) then
                        call multipole_derivative(Fi, Rji, M1s(:,j))
                    end if
                end if
                if (lmul(2)) then
                    if (maxval(abs(M2s(:,j))) >= zero) then
                        call multipole_derivative(Fi, Rji, M2s(:,j))
                    end if
                end if
                if (lmul(3)) then
                    if (maxval(abs(M3s(:,j))) >= zero) then
                        call multipole_derivative(Fi, Rji, M3s(:,j))
                    end if
                end if
                if (lmul(4)) then
                    if (maxval(abs(M4s(:,j))) >= zero) then
                        call multipole_derivative(Fi, Rji, M4s(:,j))
                    end if
                end if
                if (lmul(5)) then
                    if (maxval(abs(M5s(:,j))) >= zero) then
                        call multipole_derivative(Fi, Rji, M5s(:,j))
                    end if
                end if
                Vmuls(k) = Vmuls(k) + Fi(1)
            end do
            k = k + 1
        end do
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Vmuls(1), nsurp, rmpi, mpi_sum, master, comm, ierr)
        else if (myid /= master) then
            call mpi_reduce(Vmuls(1), dummy_real, nsurp, rmpi, mpi_sum, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            call openfile('pelib_multipole_potential.bin', lu, 'new', 'unformatted')
            rewind(lu)
            write(lu) Vmuls
            close(lu)
        end if
     end if

    if (pelib_debug) then
        do i=1, nsurp
            write (luout,*) 'Vmuls(i)' ,i, Vmuls(i)
        end do
    end if

end subroutine multipole_potentials

!------------------------------------------------------------------------------

subroutine multipole_fields(Fmuls)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: multipole_derivative, multipole_derivative_damped
    use fmm, only: field_direct, field_fmm

    real(rp), dimension(:), intent(out) :: Fmuls

    logical(lp) :: exclude, lexist
    integer(ip) :: lu
    integer(ip) :: i, j, k, l
    real(rp) :: alpha_i, alpha_j
    real(rp), dimension(3) :: Rji, Fi
    real(rp), parameter :: i3 = 1.0 / 3.0

    ! FMM related quantities
    real(rp), dimension(:,:), allocatable :: fmm_coordinates, fmm_multipoles, fmm_field
    real(rp), dimension(:), allocatable :: fmm_damping_factors
    integer(ip), dimension(:,:), allocatable :: fmm_exclusions
    character(len=6) :: fmm_damp_type = '      '

    if (myid == master) then
        inquire(file='pelib_multipole_field.bin', exist=lexist)
    end if

    if (.not. pelib_old_field) then
        allocate(fmm_coordinates(nsites, 3), source=0.0_rp)
        if (lmul(0)) k = 0
        if (lmul(1)) k = 1
        if (lmul(2)) k = 2
        if (lmul(3)) k = 3
        if (lmul(4)) k = 4
        if (lmul(5)) k = 5
        allocate(fmm_multipoles(nsites, (k+1)*(k+2)*(k+3)/6), source=0.0_rp)
        allocate(fmm_exclusions(nsites, mul_lexlst))
        allocate(fmm_field(nsites, 4))
        fmm_exclusions = transpose(mul_exclists)
        fmm_coordinates = transpose(Rs)
        if (lmul(0)) then
            fmm_multipoles(:, 1:1) = transpose(M0s)
        end if
        if (lmul(1)) then
            fmm_multipoles(:, 2:4) = transpose(M1s)
        end if
        if (lmul(2)) then
            fmm_multipoles(:, 5:10) = transpose(M2s)
        end if
        if (lmul(3)) then
            fmm_multipoles(:, 11:20) = transpose(M3s)
        end if
        if (lmul(4)) then
            fmm_multipoles(:, 21:35) = transpose(M4s)
        end if
        if (lmul(5)) then
            fmm_multipoles(:, 36:56) = transpose(M5s)
        end if
        if (pelib_mul_damp .or. pelib_amoeba_damp) then
            allocate(fmm_damping_factors(npols))
            if (pelib_mul_damp) then
                do i = 1, npols
                    alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
                    fmm_damping_factors(i) = sqrt(mul_damp) / alpha_i**(1.0 / 6.0)
                end do
            end if
            if (pelib_amoeba_damp) then
                do i = 1, npols
                    alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
                    fmm_damping_factors(i) = amoeba_damp**(1.0 / 6.0) / alpha_i**(1.0 / 6.0)
                end do
            end if
        end if
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif
    if (lexist .and. ((scf_cycle > 1) .or. pelib_restart)) then
        if (myid == master) then
            call openfile('pelib_multipole_field.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Fmuls
            close(lu)
        end if
    else
        Fmuls = 0.0
        if (.not. pelib_old_field) then
            fmm_field = 0.0
            if (pelib_mul_damp .or. pelib_amoeba_damp) then
                if (pelib_mul_damp) fmm_damp_type = 'THOLE '
                if (pelib_amoeba_damp) fmm_damp_type = 'AMOEBA'
                if (pelib_fmm) then
                    call field_fmm(comm, fmm_coordinates, fmm_multipoles, fmm_exclusions, fmm_theta, fmm_ncrit, fmm_expansion_order, 1, fmm_field, damp_type=fmm_damp_type, damping_factors=fmm_damping_factors)
                else
                    call field_direct(comm, fmm_coordinates, fmm_multipoles, fmm_exclusions, 1, fmm_field, damp_type=fmm_damp_type, damping_factors=fmm_damping_factors)
                end if
            else
                if (pelib_fmm) then
                    call field_fmm(comm, fmm_coordinates, fmm_multipoles, fmm_exclusions, fmm_theta, fmm_ncrit, fmm_expansion_order, 1, fmm_field)
                else
                    call field_direct(comm, fmm_coordinates, fmm_multipoles, fmm_exclusions, 1, fmm_field)
                end if
            end if
            l = 0
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                Fmuls(l+1) = fmm_field(i, 2)
                Fmuls(l+2) = fmm_field(i, 3)
                Fmuls(l+3) = fmm_field(i, 4)
                l = l + 3
            end do
        else
            l = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                if (pelib_mul_damp .or. pelib_amoeba_damp) then
                    alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
                end if
                do j = site_start, site_finish
                    if (i == j) then
                        cycle
                    end if
                    exclude = .false.
                    do k = 1, mul_lexlst
                        if (mul_exclists(k,i) == mul_exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (exclude) cycle
                    if (pelib_mul_damp .or. pelib_amoeba_damp) then
                        alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                    end if
                    Fi = 0.0
                    Rji = Rs(:,i) - Rs(:,j)
                    if (lmul(0)) then
                        if (maxval(abs(M0s(:,j))) >= zero) then
                            if (pelib_mul_damp) then
                                call multipole_derivative_damped(Fi, Rji, M0s(:,j), alpha_i, alpha_j, mul_damp)
                            else if (pelib_amoeba_damp) then
                                call multipole_derivative_damped(Fi, Rji, M0s(:,j), alpha_i, alpha_j, amoeba_damp)
                            else
                                call multipole_derivative(Fi, Rji, M0s(:,j))
                            end if
                        end if
                    end if
                    if (lmul(1)) then
                        if (maxval(abs(M1s(:,j))) >= zero) then
                            if (pelib_mul_damp) then
                                call multipole_derivative_damped(Fi, Rji, M1s(:,j), alpha_i, alpha_j, mul_damp)
                            else if (pelib_amoeba_damp) then
                                call multipole_derivative_damped(Fi, Rji, M1s(:,j), alpha_i, alpha_j, amoeba_damp)
                            else
                                call multipole_derivative(Fi, Rji, M1s(:,j))
                            end if
                        end if
                    end if
                    if (lmul(2)) then
                        if (maxval(abs(M2s(:,j))) >= zero) then
                            if (pelib_mul_damp) then
                                call multipole_derivative_damped(Fi, Rji, M2s(:,j), alpha_i, alpha_j, mul_damp)
                            else if (pelib_amoeba_damp) then
                                call multipole_derivative_damped(Fi, Rji, M2s(:,j), alpha_i, alpha_j, amoeba_damp)
                            else
                                call multipole_derivative(Fi, Rji, M2s(:,j))
                            end if
                        end if
                    end if
                    if (lmul(3)) then
                        if (maxval(abs(M3s(:,j))) >= zero) then
                            if (pelib_mul_damp) then
                                call multipole_derivative_damped(Fi, Rji, M3s(:,j), alpha_i, alpha_j, mul_damp)
                            else if (pelib_amoeba_damp) then
                                call multipole_derivative_damped(Fi, Rji, M3s(:,j), alpha_i, alpha_j, amoeba_damp)
                            else
                                call multipole_derivative(Fi, Rji, M3s(:,j))
                            end if
                        end if
                    end if
                    if (lmul(4)) then
                        if (maxval(abs(M4s(:,j))) >= zero) then
                            if (pelib_mul_damp .or. pelib_amoeba_damp) then
                                error stop 'damping not implemented for fourth-order multipoles'
                            else
                                call multipole_derivative(Fi, Rji, M4s(:,j))
                            end if
                        end if
                    end if
                    if (lmul(5)) then
                        if (maxval(abs(M5s(:,j))) >= zero) then
                            if (pelib_mul_damp .or. pelib_amoeba_damp) then
                                error stop 'damping not implemented for fifth-order multipoles'
                            else
                                call multipole_derivative(Fi, Rji, M5s(:,j))
                            end if
                        end if
                    end if
                    Fmuls(l:l+2) = Fmuls(l:l+2) + Fi
                end do
                l = l + 3
            end do
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                call mpi_reduce(mpi_in_place, Fmuls(1), 3 * npols, rmpi, mpi_sum, master, comm, ierr)
            else if (myid /= master) then
                call mpi_reduce(Fmuls(1), dummy_real, 3 * npols, rmpi, mpi_sum, master, comm, ierr)
            end if
#endif
        end if
        if (myid == master) then
            call openfile('pelib_multipole_field.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Fmuls
            close(lu)
        end if
     end if

end subroutine multipole_fields

!------------------------------------------------------------------------------

subroutine get_damping_factor(Ri, P11i, factor)

    use pelib_options
    use pelib_multipole_tools, only: damping_coefficient
    use pelib_blas_interfaces, only: nrm2

    real(rp), dimension(3), intent(in) :: Ri
    real(rp), dimension(6), intent(in) :: P11i
    real(rp), intent(out) :: factor

    integer(ip) :: j, jnuc
    real(rp) :: Rold, Rnew, alpha_i, alpha_m, v
    real(rp), dimension(3) :: Rmi
    real(rp), parameter :: i3 = 1.0 / 3.0

    ! Thole damping
    ! JPC A 102 (1998) 2399 & Mol. Sim. 32 (2006) 471
    ! v = a * u , where
    ! a = 2.1304 (default)
    ! u = R / (alpha_i * alpha_j)**(1/6)

    ! locate the atom closest to Ri
    jnuc = - 1
    Rold = 1.0e9
    do j = 1, nnucs
        Rmi = Ri - Rm(:,j)
        Rnew = nrm2(Rmi)
        if (Rnew <= Rold) then
            Rold = Rnew
            jnuc = j
        endif
    enddo

    if (jnuc == - 1) error stop 'Core damping failed because no atom was found'

    Rmi = Ri - Rm(:,jnuc)
    alpha_i = (P11i(1) + P11i(4) + P11i(6)) * i3
    alpha_m = (core_alphas(1,jnuc) + core_alphas(4,jnuc) + core_alphas(6,jnuc)) * i3

    call damping_coefficient(Rmi, alpha_i, alpha_m, core_damp, v)

    factor = 1.0 - (1.0 + v + 0.5 * v**2) * exp(- v)

end subroutine get_damping_factor

!------------------------------------------------------------------------------

end module pelib_potential_derivatives
