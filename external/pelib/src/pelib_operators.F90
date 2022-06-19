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
module pelib_operators

    use pelib_precision

    implicit none

    private

    public :: get_fock, get_magnetic_gradient, get_molecular_gradient
    public :: get_effdipole
    public :: es_fragment_densities

contains

!------------------------------------------------------------------------------

subroutine get_fock(denmats, fckmats, energies)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile

    real(rp), dimension(:), intent(in) :: denmats
    real(rp), dimension(:), intent(out), optional :: fckmats
    real(rp), dimension(:), intent(out), optional :: energies

    integer(ip) :: i
    logical(lp) :: es, pol

    es = .false.
    pol = .false.

    if ((any(lmul) .or. pelib_fd) .and. .not. response) es = .true.
    if (pelib_polar) pol = .true.

    if (fock .or. energy) then
        allocate(Epe(ndens))
        allocate(Ees(3,ndens))
        allocate(Efd(4,ndens))
        allocate(Epol(3,ndens))
        allocate(Esol(3,ndens))
        allocate(Elj(2))
        Epe = 0.0
        Ees = 0.0
        Efd = 0.0
        Epol = 0.0
        Esol = 0.0
        Elj = 0.0
    end if

    if (fock) then
        if (mixed) then
            if (es) call electrostatic_operator(denmats(1:nnbas), fckmats)
            if (pol) call polarization_operator(denmats(1:2*nnbas), fckmats)
        else
            if (es) call electrostatic_operator(denmats, fckmats)
            if (pol) call polarization_operator(denmats, fckmats)
        end if
    else if (energy) then
        if (es) call electrostatic_operator(denmats)
        if (pol) call polarization_operator(denmats)
    else if (response) then
        if (pol) call polarization_operator(denmats, fckmats)
    end if

    if (fock .or. energy) then
        if (lvdw) call lj_energy
        if (myid == master) then
            do i = 1, ndens
                Epe(i) = sum(Ees(1:2,i)) + sum(Efd(:,i)) + sum(Epol(:,i)) + sum(Esol(:,i)) + sum(Elj(:))
            end do
            energies = Epe
            if (energy) call print_info()
        end if
    end if

#if defined(VAR_MPI)
    if (fock .or. response) then
        if (myid == master .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, fckmats(1), ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
        else if (myid /= master) then
            call mpi_reduce(fckmats(1), dummy_real, ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
        end if
    end if
#endif

    if (fock .or. energy) then
        deallocate(Epe)
        deallocate(Ees)
        deallocate(Efd)
        deallocate(Epol)
        deallocate(Esol)
        deallocate(Elj)
    end if

end subroutine get_fock

!------------------------------------------------------------------------------

subroutine print_info()

    use pelib_options
    use pelib_utils, only: openfile
    use pelib_blas_interfaces, only: nrm2

    integer(ip) :: i, j
    integer(ip) :: lu, nwarn
    real(rp) :: indnrm
    real(rp), dimension(3) :: indtot
    real(rp), dimension(:,:), allocatable :: Mkinds
    logical(lp) :: lexist

    write(luout, '(/4x,a)') '.--------------------------------------------------.'
    write(luout, '(4x,a)')  '| Final results from polarizable embedding library |'
    write(luout, '(4x,a/)') '`--------------------------------------------------´'

    if (ndens > 1) then
        write(luout, '(a)') 'INFO: only printing information for first&
                            & input density matrix'
    end if

    if (lpol(1,1)) then
        allocate(Mkinds(3*npols,1))
        inquire(file='pelib_induced_moments.bin', exist=lexist)

        if (lexist) then
            call openfile('pelib_induced_moments.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if

        ! check induced dipoles
        nwarn = 0
        j = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            indnrm = nrm2(Mkinds(j:j+2,1))
            if (indnrm > 1.0 .and. (nwarn <= 10 .or. pelib_verbose)) then
                write(luout, '(a,i6,a,f8.4)') 'WARNING: induced dipole at site ', i, ' has a magnitude of: ', indnrm
                nwarn = nwarn + 1
            else if (nwarn > 10) then
                write(luout, '(a)') 'WARNING: there are more than ten induced dipoles with magnitude > 1.0 a.u.'
                write(luout, '(a)') '         (use .VERBOSE to print all of them)'
                exit
            end if
            j = j + 3
        end do
        if (nwarn > 0) then
            write(luout, '(a)') 'WARNING: there are induced dipoles with magnitude larger than 1.0 a.u.'
            write(luout, '(a)') '         This can be an indication that there are problems with the structure or'
            write(luout, '(a)') '         something else that can cause over-polarization. Some over-polarization'
            write(luout, '(a)') '         effects can be removed by damping the electric field on polarizable sites'
            write(luout, '(a/)') '         (please see the manual for damping options).'
        end if
    end if

 10 format(12x,a16,5x,f20.12)
 20 format(/7x,a)
    write(luout, '(8x,a)') 'Polarizable embedding energy contributions:'
    write(luout, '(7x,a)') '---------------------------------------------'
    if (any(lmul)) then
        write(luout, 20) 'Electrostatic contributions:'
        write(luout, 10) 'Electronic      ', Ees(1,1)
        write(luout, 10) 'Nuclear         ', Ees(2,1)
        write(luout, 10) 'Multipole       ', Ees(3,1)
        write(luout, 10) 'Total           ', sum(Ees(:,1))
    end if
    if (any(lpol)) then
        write(luout, 20) 'Polarization contributions:'
        write(luout, 10) 'Electronic      ', Epol(1,1)
        write(luout, 10) 'Nuclear         ', Epol(2,1)
        if (any(lmul)) then
            write(luout, 10) 'Multipole       ', Epol(3,1)
        end if
        write(luout, 10) 'Total           ', sum(Epol(:,1))
    end if
    if (lvdw) then
        write(luout, 20) 'LJ energy contribution:'
        write(luout, 10) 'Dispersion 6    ', Elj(1)
        write(luout, 10) 'Repulsion 12    ', Elj(2)
        write(luout, 10) 'Total           ', sum(Elj)
    end if
    if (pelib_sol) then
        write(luout, 20) 'Continuum solvation contributions:'
        write(luout, 10) 'Electronic      ', Esol(1,1)
        write(luout, 10) 'Nuclear         ', Esol(2,1)
        if (any(lmul)) then
            write(luout, 10) 'Multipole       ', Esol(3,1)
        end if
        write(luout, 10) 'Total           ', sum(Esol(:,1))
    end if
    if (pelib_fd) then
        write(luout, 20) 'Fragment density contributions:'
        write(luout, 10) 'Electronic      ', Efd(1,1)
        write(luout, 10) 'Nuclear         ', Efd(2,1)
        if (lpol(1,1)) then
            write(luout, 10) 'Polarization    ', Efd(3,1)
        end if
        if (pelib_repuls) then
            write(luout, 10) 'Repulsion       ', Efd(4,1)
        end if
        write(luout, 10) 'Total           ', sum(Efd(:,1))
    end if
    write(luout,'(/7x,a17,9x,f20.12)') 'Total PE energy: ', sum(Ees(1:2,1)) + &
                                                            sum(Efd(:,1)) + &
                                                            sum(Epol(:,1)) + &
                                                            sum(Esol(:,1)) + &
                                                            sum(Elj(:))
    write(luout,'(12x,a)') '(excl. multipole-multipole interaction energy)'
    if (lpol(1,1) .and. pelib_verbose) then
        write(luout, '(/15x,a)') 'Polarizable embedding information:'
        write(luout, '(14x,a/)') '------------------------------------'
        write(luout, '(23x,a)') 'Induced dipole moments'
        write(luout, '(7x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
        j = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            write(luout, '(5x,i6,3f15.8)') i, Mkinds(j:j+2,1)
            j = j + 3
        end do
        indtot = 0.0
        j = 1
        do i = 1, npols
            indtot(1) = indtot(1) + Mkinds(j,1)
            indtot(2) = indtot(2) + Mkinds(j+1,1)
            indtot(3) = indtot(3) + Mkinds(j+2,1)
            j = j + 3
        end do
        write(luout, '(/22x,a)') 'Total induced dipole moment'
        write(luout, '(21x,a,14x,a,14x,a)') 'X', 'Y', 'Z'
        write(luout, '(11x,3f15.8/)') indtot
    end if

    if (lpol(1,1)) deallocate(Mkinds)

end subroutine print_info

!-----------------------------------------------------------------------------

subroutine electrostatic_operator(denmats, fckmats)

    use pelib_mpi
    use pelib_options

    real(rp), dimension(:), intent(in) :: denmats
    real(rp), dimension(:), intent(inout), optional :: fckmats

    if (fock) then
        call es_multipoles(denmats, fckmats)
    else if (energy) then
        call es_multipoles(denmats)
    end if

    if (myid == master) then
        if (pelib_fd) then
            if (fock) then
                call es_fragment_densities(denmats, fckmats)
            else if (energy) then
                call es_fragment_densities(denmats)
            end if
        end if
    end if

end subroutine electrostatic_operator

!------------------------------------------------------------------------------

#ifdef VAR_PDE

subroutine es_fragment_densities(denmats, fckmats)

    use hdf5

    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: Tk_tensor
    use pelib_blas_interfaces, only: dot

    real(rp), dimension(:), intent(in) :: denmats
    real(rp), dimension(:), intent(inout), optional :: fckmats

    logical(lp) :: lexist
    integer(ip) :: lu
    integer(ip) :: i, j, k
    integer(ip) :: num_bas
    real(rp) :: Ene, Enn
    real(rp), dimension(ndens) :: Een, Eee
    real(rp), dimension(1) :: Tfm
    real(rp), dimension(3) :: Rfm
    real(rp), dimension(:), allocatable :: fd_esmat, fd_repmat
    real(rp), dimension(:,:), allocatable :: Zfd_ints
    real(rp), dimension(:), allocatable :: tmpfckmats

    !hdf5 IO variables
    integer(4) :: error
    integer(HSIZE_T), dimension(1) :: dim_1d
    integer(HSIZE_T), dimension(2) :: dim_2d
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id

    Efd = 0.0
    Eee = 0.0; Een = 0.0; Ene = 0.0; Enn = 0.0

    if (.not. pelib_repuls .and. .not. pelib_fdes) then
        return
    end if

    inquire(file='pde_electrostatics.bin', exist=lexist)
    if (lexist .and. fock .and. ((scf_cycle > 1) .or. pelib_restart)) then
        allocate(tmpfckmats(ndens * nnbas))
        call openfile('pde_electrostatics.bin', lu, 'old', 'unformatted')
        rewind(lu)
        read(lu) Efd(2,1), tmpfckmats
        close(lu)
        do i = 1, ndens
            j = (i - 1) * nnbas + 1
            k = i * nnbas
            Efd(1,i) = dot(denmats(j:k), tmpfckmats(j:k))
        end do
        fckmats = fckmats + tmpfckmats
        deallocate(tmpfckmats)
    else
        if (fock) then
            allocate(tmpfckmats(ndens * nnbas))
            tmpfckmats = 0.0
        end if
        call h5open_f(error)
        call h5fopen_f(trim(h5pdefile), H5F_ACC_RDONLY_F, file_id, error)

        call h5dopen_f(file_id, 'num_bas', dset_id, error)
        dim_1d = 1
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_bas, dim_1d, error)
        call h5dclose_f(dset_id, error)

        if (nnbas /= num_bas*(num_bas+1)/2) then
            error stop 'size inconsistency in es_fragment_densities'
        end if

        if (pelib_fdes) then
            call h5dopen_f(file_id, 'num_nuclei', dset_id, error)
            dim_1d = 1
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdnucs, dim_1d, error)
            call h5dclose_f(dset_id, error)

            allocate(Rfd(3,fdnucs), Zfd(1,fdnucs))

            call h5dopen_f(file_id, 'nuclear coordinates', dset_id, error)
            dim_2d = [3, fdnucs]
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Rfd, dim_2d, error)
            call h5dclose_f(dset_id, error)

            call h5dopen_f(file_id, 'nuclear charges', dset_id, error)
            dim_2d = [1, fdnucs]
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Zfd, dim_2d, error)
            call h5dclose_f(dset_id, error)

            call h5dopen_f(file_id, 'nuclear-electron energy', dset_id, error)
            dim_1d = 1
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Ene, dim_1d, error)
            call h5dclose_f(dset_id, error)

            allocate(fd_esmat(nnbas))

            call h5dopen_f(file_id, 'electrostatic matrix', dset_id, error)
            dim_1d = nnbas
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fd_esmat, dim_1d, error)
            call h5dclose_f(dset_id, error)

            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                if (fock) tmpfckmats(j:k) = tmpfckmats(j:k) + fd_esmat
                Eee(i) = dot(denmats(j:k), fd_esmat)
            end do
            deallocate(fd_esmat)

            allocate(Zfd_ints(nnbas,1))
            Zfd_ints = 0.0
            do i = 1, fdnucs
                do j = 1, nnucs
                    Rfm = Rm(:,j) - Rfd(:,i)
                    call Tk_tensor(Tfm, Rfm)
                    Enn = Enn + Zm(1,j) * Zfd(1,i) * Tfm(1)
                end do
                call Mk_integrals('potential_derivative', Zfd_ints, Rfd(:,i), Zfd(:,i))
            end do
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                if (fock) tmpfckmats(j:k) = tmpfckmats(j:k) + Zfd_ints(:,1)
                Een(i) = dot(denmats(j:k), Zfd_ints(:,1))
            end do

            do i = 1, ndens
                Efd(1,i) = Een(i) + Eee(i)
            end do
            Efd(2,1) = Ene + Enn

            deallocate(Zfd_ints)
            deallocate(Rfd, Zfd)
        end if

        if (pelib_repuls) then
            allocate(fd_repmat(nnbas))
            call h5dopen_f(file_id, 'exchange-repulsion matrix', dset_id, error)
            dim_1d = nnbas
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fd_repmat, dim_1d, error)
            call h5dclose_f(dset_id, error)
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                if (fock) tmpfckmats(j:k) = tmpfckmats(j:k) + fd_repmat
                Efd(4,i) = dot(denmats(j:k), fd_repmat)
            end do
            deallocate(fd_repmat)
        end if

        call h5fclose_f(file_id, error)
        call h5close_f(error)

        if (fock) then
            call openfile('pde_electrostatics.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Efd(2,1), tmpfckmats
            close(lu)
            fckmats = fckmats + tmpfckmats
            deallocate(tmpfckmats)
        end if
    end if

end subroutine es_fragment_densities

#else

subroutine es_fragment_densities(denmats, fckmats)

    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: Tk_tensor
    use pelib_blas_interfaces, only: dot

    real(rp), dimension(:), intent(in) :: denmats
    real(rp), dimension(:), intent(inout), optional :: fckmats

    error stop 'Not compiled with PDE support (recompile with -DENABLE_PDE=ON)'

end subroutine es_fragment_densities

#endif

!------------------------------------------------------------------------------

subroutine es_multipoles(denmats, fckmats)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_blas_interfaces, only: dot
    use pelib_multipole_tools, only: multipole_interaction_energy, &
                                  prefactors, Tk_tensor

    real(rp), dimension(:), intent(in) :: denmats
    real(rp), dimension(:), intent(inout), optional :: fckmats

    logical(lp) :: lexist, exclude
    integer(ip) :: lu
    integer(ip) :: site, ncomps
    integer(ip) :: i, j, k
    real(rp) :: tmp_energy
    real(rp), dimension(3) :: Rsm, Rij
    real(rp), dimension(:), allocatable :: Tsm, factors
    real(rp), dimension(:,:), allocatable :: Mk_ints
    real(rp), dimension(:), allocatable :: tmpfckmats

    Ees = 0.0

    if (myid == master) then
        inquire(file='pelib_electrostatics.bin', exist=lexist)
    end if
#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif
    if (lexist .and. fock .and. ((scf_cycle > 1) .or. pelib_restart)) then
        if (myid == master) then
            allocate(tmpfckmats(ndens * nnbas))
            call openfile('pelib_electrostatics.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Ees(2,1), Ees(3,1), tmpfckmats
            close(lu)
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(1,i) = dot(denmats(j:k), tmpfckmats(j:k))
            end do
            fckmats = fckmats + tmpfckmats
            deallocate(tmpfckmats)
        end if
    else
        allocate(Mk_ints(nnbas,1))
        do site = site_start, site_finish
            ! electron - multipole interaction energy
            Mk_ints = 0.0
            if (lmul(0)) then
                if (maxval(abs(M0s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M0s(:,site))
                end if
            end if
            if (lmul(1)) then
                if (maxval(abs(M1s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M1s(:,site))
                end if
            end if
            if (lmul(2)) then
                if (maxval(abs(M2s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M2s(:,site))
                end if
            end if
            if (lmul(3)) then
                if (maxval(abs(M3s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M3s(:,site))
                end if
            end if
            if (lmul(4)) then
                if (maxval(abs(M4s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M4s(:,site))
                end if
            end if
            if (lmul(5)) then
                if (maxval(abs(M5s(:,site))) >= zero) then
                    call Mk_integrals('potential_derivative', Mk_ints, Rs(:,site), M5s(:,site))
                end if
            end if
            do i = 1, ndens
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                Ees(1,i) = Ees(1,i) + dot(denmats(j:k), Mk_ints(:,1))
                if (fock) fckmats(j:k) = fckmats(j:k) + Mk_ints(:,1)
            end do

            ! nuclei - multipole interaction energy
            if (lmul(0)) then
                if (maxval(abs(M0s(:,site))) >= zero) then
                    ncomps = size(M0s(:,site), 1)
                    allocate(Tsm(ncomps))
                    allocate(factors(ncomps))
                    do i = 1, nnucs
                        Rsm = Rm(:,i) - Rs(:,site)
                        call Tk_tensor(Tsm, Rsm)
                        call prefactors(factors)
                        do j = 1, ncomps
                            Ees(2,1) = Ees(2,1) + factors(j) * M0s(j,site) * Tsm(j) * Zm(1,i)
                        end do
                    end do
                    deallocate(Tsm)
                    deallocate(factors)
                end if
            end if
            if (lmul(1)) then
                if (maxval(abs(M1s(:,site))) >= zero) then
                    ncomps = size(M1s(:,site), 1)
                    allocate(Tsm(ncomps))
                    allocate(factors(ncomps))
                    do i = 1, nnucs
                        Rsm = Rm(:,i) - Rs(:,site)
                        call Tk_tensor(Tsm, Rsm)
                        call prefactors(factors)
                        do j = 1, ncomps
                            Ees(2,1) = Ees(2,1) + factors(j) * M1s(j,site) * Tsm(j) * Zm(1,i)
                        end do
                    end do
                    deallocate(Tsm)
                    deallocate(factors)
                end if
            end if
            if (lmul(2)) then
                if (maxval(abs(M2s(:,site))) >= zero) then
                    ncomps = size(M2s(:,site), 1)
                    allocate(Tsm(ncomps))
                    allocate(factors(ncomps))
                    do i = 1, nnucs
                        Rsm = Rm(:,i) - Rs(:,site)
                        call Tk_tensor(Tsm, Rsm)
                        call prefactors(factors)
                        do j = 1, ncomps
                            Ees(2,1) = Ees(2,1) + factors(j) * M2s(j,site) * Tsm(j) * Zm(1,i)
                        end do
                    end do
                    deallocate(Tsm)
                    deallocate(factors)
                end if
            end if
            if (lmul(3)) then
                if (maxval(abs(M3s(:,site))) >= zero) then
                    ncomps = size(M3s(:,site), 1)
                    allocate(Tsm(ncomps))
                    allocate(factors(ncomps))
                    do i = 1, nnucs
                        Rsm = Rm(:,i) - Rs(:,site)
                        call Tk_tensor(Tsm, Rsm)
                        call prefactors(factors)
                        do j = 1, ncomps
                            Ees(2,1) = Ees(2,1) + factors(j) * M3s(j,site) * Tsm(j) * Zm(1,i)
                        end do
                    end do
                    deallocate(Tsm)
                    deallocate(factors)
                end if
            end if
            if (lmul(4)) then
                if (maxval(abs(M4s(:,site))) >= zero) then
                    ncomps = size(M4s(:,site), 1)
                    allocate(Tsm(ncomps))
                    allocate(factors(ncomps))
                    do i = 1, nnucs
                        Rsm = Rm(:,i) - Rs(:,site)
                        call Tk_tensor(Tsm, Rsm)
                        call prefactors(factors)
                        do j = 1, ncomps
                            Ees(2,1) = Ees(2,1) + factors(j) * M4s(j,site) * Tsm(j) * Zm(1,i)
                        end do
                    end do
                    deallocate(Tsm)
                    deallocate(factors)
                end if
            end if
            if (lmul(5)) then
                if (maxval(abs(M5s(:,site))) >= zero) then
                    ncomps = size(M5s(:,site), 1)
                    allocate(Tsm(ncomps))
                    allocate(factors(ncomps))
                    do i = 1, nnucs
                        Rsm = Rm(:,i) - Rs(:,site)
                        call Tk_tensor(Tsm, Rsm)
                        call prefactors(factors)
                        do j = 1, ncomps
                            Ees(2,1) = Ees(2,1) + factors(j) * M5s(j,site) * Tsm(j) * Zm(1,i)
                        end do
                    end do
                    deallocate(Tsm)
                    deallocate(factors)
                end if
            end if
        end do
        deallocate(Mk_ints)

        if ((myid == master) .and. (pelib_skipmul .eqv. .false.)) then
            ! multipole - multipole interaction energy
            do i = 1, nsites
                do j = i + 1, nsites
                    exclude = .false.
                    do k = 1, mul_lexlst
                        if (mul_exclists(k,i) == mul_exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (exclude) then
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    if (lmul(0)) then
                        if (maxval(abs(M0s(:,i))) >= zero .and. maxval(abs(M0s(:,j))) >= zero) then
                            call multipole_interaction_energy(tmp_energy, M0s(:,i), Rij, M0s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                        end if
                    end if
                    if (lmul(1)) then
                        if (maxval(abs(M1s(:,i))) >= zero .and. maxval(abs(M1s(:,j))) >= zero) then
                            call multipole_interaction_energy(tmp_energy, M0s(:,i), Rij, M1s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M1s(:,i), Rij, M0s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M1s(:,i), Rij, M1s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                        end if
                    end if
                    if (lmul(2)) then
                        if (maxval(abs(M2s(:,i))) >= zero .and. maxval(abs(M2s(:,j))) >= zero) then
                            call multipole_interaction_energy(tmp_energy, M0s(:,i), Rij, M2s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M2s(:,i), Rij, M0s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M1s(:,i), Rij, M2s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M2s(:,i), Rij, M1s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M2s(:,i), Rij, M2s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                        end if
                    end if
                    if (lmul(3)) then
                        if (maxval(abs(M3s(:,i))) >= zero .and. maxval(abs(M3s(:,j))) >= zero) then
                            call multipole_interaction_energy(tmp_energy, M0s(:,i), Rij, M3s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M3s(:,i), Rij, M0s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M1s(:,i), Rij, M3s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M3s(:,i), Rij, M1s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M2s(:,i), Rij, M3s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M3s(:,i), Rij, M2s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M3s(:,i), Rij, M3s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                        end if
                    end if
                    if (lmul(4)) then
                        if (maxval(abs(M4s(:,i))) >= zero .and. maxval(abs(M4s(:,j))) >= zero) then
                            call multipole_interaction_energy(tmp_energy, M0s(:,i), Rij, M4s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M4s(:,i), Rij, M0s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M1s(:,i), Rij, M4s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M4s(:,i), Rij, M1s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M2s(:,i), Rij, M4s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M4s(:,i), Rij, M2s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M3s(:,i), Rij, M4s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M4s(:,i), Rij, M3s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M4s(:,i), Rij, M4s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                        end if
                    end if
                    if (lmul(5)) then
                        if (maxval(abs(M5s(:,i))) >= zero .and. maxval(abs(M5s(:,j))) >= zero) then
                            call multipole_interaction_energy(tmp_energy, M0s(:,i), Rij, M5s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M5s(:,i), Rij, M0s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M1s(:,i), Rij, M5s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M5s(:,i), Rij, M1s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M2s(:,i), Rij, M5s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M5s(:,i), Rij, M2s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M3s(:,i), Rij, M5s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M5s(:,i), Rij, M3s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M4s(:,i), Rij, M5s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M5s(:,i), Rij, M4s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                            call multipole_interaction_energy(tmp_energy, M5s(:,i), Rij, M5s(:,j))
                            Ees(3,1) = Ees(3,1) + tmp_energy
                        end if
                    end if
                end do
            end do
        end if

#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, Ees(1,1), 3 * ndens, rmpi, mpi_sum, master, comm, ierr)
        else if (myid /= master) then
            call mpi_reduce(Ees(1,1), dummy_real, 3 * ndens, rmpi, mpi_sum, master, comm, ierr)
        end if
#endif

        if (fock) then
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                allocate(tmpfckmats(ndens * nnbas))
                tmpfckmats = fckmats
                call mpi_reduce(mpi_in_place, fckmats(1), ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
            else if (myid /= master) then
                call mpi_reduce(fckmats(1), dummy_real, ndens * nnbas, rmpi, mpi_sum, master, comm, ierr)
            end if
#endif
            if (myid == master) then
                call openfile('pelib_electrostatics.bin', lu, 'unknown', 'unformatted')
                rewind(lu)
                write(lu) Ees(2,1), Ees(3,1), fckmats
                close(lu)
            end if
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                fckmats = tmpfckmats
                deallocate(tmpfckmats)
            end if
#endif
        end if
    end if

end subroutine es_multipoles

!------------------------------------------------------------------------------

subroutine polarization_operator(denmats, fckmats)

    use pelib_mpi
    use pelib_options
    use pelib_potential_derivatives
    use pelib_integral_interfaces
    use pelib_utils, only: openfile
    use pelib_blas_interfaces, only: dot
    use pelib_induced_moments, only: induced_moments
    use fmm, only: get_multipoles
    use pelib_multipole_tools, only: prefactors
    use pelib_interaction_tools, only: build_multipole_interaction_operators, get_multipole_field

    real(rp), dimension(:), intent(in), optional :: denmats
    real(rp), dimension(:), intent(inout), optional :: fckmats

    integer(ip) :: site
    integer(ip) :: i, j, k, l, m, n, idx, ncomp
    integer(ip) :: multipole_size
    integer(ip) :: lu
    logical(lp) :: lexist
    real(rp), dimension(3) :: indtot
    real(rp), dimension(:,:), allocatable :: Vels, Vtots, Qinds
    real(rp), dimension(:), allocatable :: Vnucs, Vmuls, Vfds
    real(rp), dimension(:,:,:), allocatable :: Vel_ints
    real(rp), dimension(:,:), allocatable :: Fels, Ftots, Mkinds
    real(rp), dimension(:), allocatable :: Fnucs, Fmuls, Ffds, Fext
    real(rp), dimension(:,:,:), allocatable :: Fel_ints
    ! fmm stuff
    real(rp), dimension(:,:), allocatable :: fmm_coordinates, box_coordinates
    real(rp), dimension(:,:,:), allocatable :: fmm_multipoles, box_multipoles, multipoles, mm_field_derivatives
    real(rp), dimension(:), allocatable :: factors
    integer(ip) :: iwork, work_start, work_stop, work_size
    integer(ip) :: start_idx, stop_idx
    integer(ip) :: num_operators

    if (lpol(1,1)) then
        allocate(Mkinds(3*npols,ndens), Fels(3*npols,ndens), Ftots(3*npols,ndens))
        Fels = 0.0; Ftots = 0.0
        allocate(Fnucs(3*npols), Fmuls(3*npols), Ffds(3*npols), Fext(3*npols))
        Fnucs = 0.0; Fmuls = 0.0; Ffds = 0.0; Fext= 0.0
    end if
    if (pelib_sol) then
        allocate(Qinds(nsurp,ndens), Vels(nsurp,ndens), Vtots(nsurp,ndens))
        Vels = 0.0; Vtots = 0.0
        allocate(Vnucs(nsurp), Vmuls(nsurp), Vfds(nsurp))
        Vnucs = 0.0; Vmuls = 0.0; Vfds = 0.0
    end if

    if (response) then
        if (lpol(1,1)) then
            if (trim(qmmm_interaction_scheme) == 'EXACT') then
                call electron_fields(Fels, denmats)
            else if (trim(qmmm_interaction_scheme) == 'FMM_FOCK') then
                call electron_fields_fmm(Fels, denmats)
            else
                ! ESPF, CESPF, SINGLE_CENTER
                call electron_fields_multipole(Fels, denmats)
            end if
            if (myid == master) then
                do i = 1, ndens
                    Ftots(:,i) = Fels(:,i)
                end do
            end if
        end if
        if (pelib_sol) then
            call electron_potentials(Vels, denmats)
            if (myid == master) then
                do i = 1, ndens
                    if (pelib_diis) then
                        Vtots(:,i) = Vels(:,i)
                    else
                        Vtots(:,i) = - Vels(:,i)
                    end if
                end do
            end if
        end if
    else
        if (lpol(1,1)) then
            if (pelib_field) then
                do site = site_start, site_finish
                    if (zeroalphas(site)) cycle
                    do i = 1, 3
                        Fext(i+3*(site-1)) = efields(i)
                    end do
                end do
            end if
            if (mixed) then
                start_idx = 1+nnbas
                stop_idx  = 2*nnbas
            else
                start_idx = 1
                stop_idx = nnbas
            end if
            if (trim(qmmm_interaction_scheme) == 'EXACT') then
                call electron_fields(Fels, denmats(start_idx:stop_idx))
            else if (trim(qmmm_interaction_scheme) == 'FMM_FOCK') then
                call electron_fields_fmm(Fels, denmats(start_idx:stop_idx))
            else
                ! ESPF, CESPF, SINGLE_CENTER
                call electron_fields_multipole(Fels, denmats(start_idx:stop_idx))
            end if
            call nuclear_fields(Fnucs)
            call multipole_fields(Fmuls)
            if (myid == master) then
                if (pelib_fd .and. pelib_fdes) call fragment_density_fields(Ffds)
                do i = 1, ndens
                    Ftots(:,i) = Fels(:,i) + Fnucs + Fmuls + Ffds + Fext
                end do
            end if
        end if
        if (pelib_sol) then
            if (mixed) then
                call electron_potentials(Vels, denmats(1+nnbas:2*nnbas))
            else
                call electron_potentials(Vels, denmats(1:nnbas))
            end if
            call nuclear_potentials(Vnucs)
            if (.not. pelib_novmu) then
                call multipole_potentials(Vmuls)
            end if
            if (myid == master) then
                do i = 1, ndens
                    if (pelib_sol) then
                        if (pelib_diis) then
                            Vtots(:,i) = Vels(:,i) + Vnucs + Vmuls
                        else
                            Vtots(:,i) = - (Vels(:,i) + Vnucs + Vmuls)
                        end if
                    end if
                end do
            end if
        end if
    end if
    if (energy .and. (any(lpol) .or. pelib_sol)) then
        if (myid == master) then
            if (any(lpol)) then
                inquire(file='pelib_induced_moments.bin', exist=lexist)
                if (lexist) then
                    call openfile('pelib_induced_moments.bin', lu, 'old', 'unformatted')
                    rewind(lu)
                    read(lu) Mkinds
                    close(lu)
                end if
            end if
            if (pelib_sol) then
                inquire(file='pelib_induced_charges.bin', exist=lexist)
                if (lexist) then
                    call openfile('pelib_induced_charges.bin', lu, 'old', 'unformatted')
                    rewind(lu)
                    read(lu) Qinds
                    close(lu)
                end if
            end if
        end if
    else if (any(lpol) .and. pelib_sol) then
        call induced_moments(Mkinds=Mkinds, Fks=Ftots, Qinds=Qinds, Vs=Vtots)
    else if (any(lpol) .and. .not. pelib_sol) then
        call induced_moments(Mkinds=Mkinds, Fks=Ftots)
    else if (.not. any(lpol) .and. pelib_sol) then
        call induced_moments(Qinds=Qinds, Vs=Vtots)
    end if
    if (.not. (response) .and. myid == master) then
        do i = 1, ndens
            if (any(lpol)) then
                if (lpol(1,1)) then
                    if (mixed) then
                        Fels = 0.0
                        if (trim(qmmm_interaction_scheme) == 'EXACT') then
                            call electron_fields(Fels, denmats(1:nnbas))
                        else if (trim(qmmm_interaction_scheme) == 'FMM_FOCK') then
                            call electron_fields_fmm(Fels, denmats(1:nnbas))
                        else
                            ! ESPF, CESPF, SINGLE_CENTER
                            call electron_fields_multipole(Fels, denmats(1:nnbas))
                        end if
                    end if
                    Epol(1,i) = - 0.5 * dot(Mkinds(:,i), Fels(:,i))
                    Epol(2,i) = - 0.5 * dot(Mkinds(:,i), Fnucs)
                    if (any(lmul)) Epol(3,i) = - 0.5 * dot(Mkinds(:,i), Fmuls)
                end if
                if (pelib_fd .and. pelib_fdes) then
                    Efd(3,i) = - 0.5 * dot(Mkinds(:,i), Ffds)
                end if
            end if
            if (pelib_sol) then
                if (mixed) then
                    Vels = 0.0
                    call electron_potentials(Vels, denmats(1:nnbas))
                end if
                Esol(1,i) = 0.5 * dot(Qinds(:,i), Vels(:,i))
                Esol(2,i) = 0.5 * dot(Qinds(:,i), Vnucs)
                Esol(3,i) = 0.5 * dot(Qinds(:,i), Vmuls)
            end if
        end do
    end if
    if (myid == master .and. (pelib_debug .or. pelib_verbose) .and. .not. energy) then
        write(luout, '(/6x,a)') 'Info from polarizable embedding library:'
        write(luout, '(5x,a/)') '------------------------------------------'
        do i = 1, ndens
            write(luout, '(7x,a,i3)') 'Input density no.: ', i
            if (lpol(1,1) .and. .not. pelib_debug) then
                indtot = 0.0
                k = 1
                do j = 1, npols
                    indtot(1) = indtot(1) + Mkinds(k,i)
                    indtot(2) = indtot(2) + Mkinds(k+1,i)
                    indtot(3) = indtot(3) + Mkinds(k+2,i)
                    k = k + 3
                end do
                write(luout, '(/24x,a)') 'Total induced dipole moment'
                write(luout, '(23x,a,14x,a,14x,a)') 'X', 'Y', 'Z'
                write(luout, '(13x,3f15.8/)') indtot
            end if
            if (any(lpol) .and. pelib_debug) then
                write(luout, '(/15x,a)') 'Total electric field at polarizable sites:'
                write(luout, '(9x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
                k = 1
                do site = 1, nsites
                    if (zeroalphas(site)) cycle
                    write(luout, '(7x,i6,3f15.8)') site, Ftots(k:k+2,i)
                    k = k + 3
                end do
                write(luout, '(/25x,a)') 'Induced dipole moments'
                write(luout, '(9x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
                k = 1
                do site = 1, nsites
                    if (zeroalphas(site)) cycle
                    write(luout, '(7x,i6,3f15.8)') site, Mkinds(k:k+2,i)
                    k = k + 3
                end do
            end if
            if (pelib_sol .and. pelib_debug) then
                write(luout, '(28x,a)') 'Induced charges'
                do j = 1, nsurp
                    write(luout, '(25x,f15.8)') Qinds(j,i)
                end do
                write(luout, '(4x,a,f15.8/)') 'Sum of induced charges: ', sum(Qinds(1:nsurp,i))
            end if
            if (pelib_sol .and. .not. pelib_debug) then
                write(luout, '(4x,a,f15.8/)') 'Sum of induced charges: ', sum(Qinds(1:nsurp,i))
            end if
        end do
    end if
#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        do i = 1, ndens
            if (lpol(1,1)) then
                displs(0) = 0
                do j = 1, nprocs
                    displs(j) = displs(j-1) + poldists(j-1)
                end do
                call mpi_scatterv(Mkinds(1,i), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
            end if
            if (pelib_sol) then
                displs(0) = 0
                do j = 1, nprocs
                    displs(j) = displs(j-1) + surpdists(j-1)
                end do
                call mpi_scatterv(Qinds(1,i), surpdists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
            end if
        end do
    else if (myid /= master) then
        do i = 1, ndens
            if (lpol(1,1)) then
                call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,i), poldists(myid), rmpi, master, comm, ierr)
            end if
            if (pelib_sol) then
                call mpi_scatterv(dummy_real, surpdists, displs, rmpi, Qinds(1,i), surpdists(myid), rmpi, master, comm, ierr)
            end if
        end do
    end if
#endif
    if (fock .or. response) then
        if (any(lpol)) then
            ! approximate fock contributions
            if (trim(qmmm_interaction_scheme) == 'FMM_FOCK') then
#if defined(VAR_MPI)
                if (nprocs > 1) then
                    call mpi_bcast(Mkinds(1,1), 3*npols*ndens, rmpi, master, comm, ierr)
                end if
#endif
                allocate(fmm_coordinates(npols, 3))
                allocate(fmm_multipoles(npols, 4, ndens))
                l = 1
                do i = 1, nsites
                    if (zeroalphas(i)) cycle
                    fmm_coordinates(l, :) = Rs(:, i)
                    l = l + 1
                end do
                do i = 0, npols - 1
                    do n = 1, ndens
                        fmm_multipoles(i+1, 1, n) = 0.0
                        fmm_multipoles(i+1, 2, n) = Mkinds(3*i+1, n)
                        fmm_multipoles(i+1, 3, n) = Mkinds(3*i+2, n)
                        fmm_multipoles(i+1, 4, n) = Mkinds(3*i+3, n)
                    end do
                end do
                call get_multipoles(comm, fmm_coordinates, fmm_multipoles, transpose(Rm), qmmm_interaction_theta, fmm_ncrit, qmmm_interaction_order, idx_near_field, box_coordinates, box_multipoles)
                allocate(Fel_ints(nnbas,3,1))
                ! near field
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
                        call Tk_integrals('potential_derivative', fmm_coordinates(idx, :), Fel_ints)
                    end if
                    do j = 1, 3
                        do k = 1, ndens
                            l = (k - 1) * nnbas + 1
                            m = k * nnbas
                            if (mixed) then
                                fckmats(l:m) = fckmats(l:m) - 0.5 * fmm_multipoles(idx, 1+j, k) * Fel_ints(:,j,1)
                            else
                                fckmats(l:m) = fckmats(l:m) - fmm_multipoles(idx, 1+j, k) * Fel_ints(:,j,1)
                            end if
                        end do
                    end do
                end do
                ! far field
                deallocate(Fel_ints)
                multipole_size = (qmmm_interaction_order + 1) * (qmmm_interaction_order + 2) * (qmmm_interaction_order + 3) / 6
                work_size = size(box_coordinates, 1) / nprocs
                work_start = 1 + myid * work_size
                work_stop = (myid + 1) * work_size
                if (myid == nprocs - 1) work_stop = size(box_coordinates, 1)
                allocate(factors(multipole_size))
                l = 1
                do k = 0, qmmm_interaction_order
                    ncomp = (k+1)*(k+2)/2
                    call prefactors(factors(l:l+ncomp-1))
                    l = l + ncomp
                end do
                if (.not. allocated(Fel_ints_fmm)) then
                    allocate(Fel_ints_fmm(nnbas, multipole_size, 1, work_start:work_stop))
                    do iwork = work_start, work_stop
                        l = 1
                        do k = 0, qmmm_interaction_order
                            ncomp = (k+1)*(k+2)/2
                            call Tk_integrals('potential_derivative', box_coordinates(iwork, :), Fel_ints_fmm(:, l:l+ncomp-1, 1:1, iwork))
                            l = l + ncomp
                        end do
                    end do
                end if
                do iwork = work_start, work_stop
                    do j = 1, multipole_size
                        do k = 1, ndens
                            l = (k - 1) * nnbas + 1
                            m = k * nnbas
                            if (mixed) then
                                fckmats(l:m) = fckmats(l:m) + 0.5 * factors(j) * box_multipoles(iwork, j, k) * Fel_ints_fmm(:,j,1,iwork)
                            else
                                fckmats(l:m) = fckmats(l:m) + factors(j) * box_multipoles(iwork, j, k) * Fel_ints_fmm(:,j,1,iwork)
                            end if
                        end do
                    end do
                end do
                deallocate(factors)
            else if (trim(qmmm_interaction_scheme) == 'EXACT') then
                ! normal
                allocate(Fel_ints(nnbas,3,1))
                i = 0
                do site = site_start, site_finish
                    if (zeroalphas(site)) cycle
                    call Tk_integrals('potential_derivative', Rs(:,site), Fel_ints)
                    do j = 1, 3
                        do k = 1, ndens
                            l = (k - 1) * nnbas + 1
                            m = k * nnbas
                            if (mixed) then
                                fckmats(l:m) = fckmats(l:m) - 0.5 * Mkinds(i+j,k) * Fel_ints(:,j,1)
                            else
                                fckmats(l:m) = fckmats(l:m) - Mkinds(i+j,k) * Fel_ints(:,j,1)
                            end if
                        end do
                    end do
                    i = i + 3
                end do
                deallocate(Fel_ints)
            else
#if defined(VAR_MPI)
                if (nprocs > 1) then
                    call mpi_bcast(Mkinds(1,1), 3*npols*ndens, rmpi, master, comm, ierr)
                end if
#endif
                ! single-center, ESPF or CESPF
                ! ...or wrong input
                if (.not. allocated(multipole_interaction_operators)) then
                    call build_multipole_interaction_operators
                end if
                ! near-field
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
                allocate(Fel_ints(nnbas,3,1))
                do iwork = work_start, work_stop
                    idx = idx_near_field(iwork)
                    if (incore_nearfield) then
                        Fel_ints(:,:,1) = near_field_integrals(:, :, 1, iwork)
                    else
                        call Tk_integrals('potential_derivative', mm_coordinates_near(:, iwork), Fel_ints)
                    end if
                    do j = 1, 3
                        do k = 1, ndens
                            l = (k - 1) * nnbas + 1
                            m = k * nnbas
                            if (mixed) then
                                fckmats(l:m) = fckmats(l:m) - 0.5 * Mkinds(3*(idx-1)+j, k) * Fel_ints(:,j,1)
                            else
                                fckmats(l:m) = fckmats(l:m) - Mkinds(3*(idx-1)+j, k) * Fel_ints(:,j,1)
                            end if
                        end do
                    end do
                end do
                deallocate(Fel_ints)
                allocate(multipoles(size(idx_far_field), 4, ndens))
                multipoles(:, 1, :) = 0.0d0
                do i = 1, size(idx_far_field)
                    idx = idx_far_field(i)
                    do k = 1, ndens
                        multipoles(i, 2, k) = Mkinds(3*(idx-1)+1, k)
                        multipoles(i, 3, k) = Mkinds(3*(idx-1)+2, k)
                        multipoles(i, 4, k) = Mkinds(3*(idx-1)+3, k)
                    end do
                end do
                multipole_size = size(multipole_interaction_operators, 2) 
                num_operators = size(multipole_interaction_operators, 3)
                allocate(mm_field_derivatives(num_operators, multipole_size, ndens))
                call get_multipole_field(mm_field_derivatives, mm_coordinates_far, multipole_operator_centers, multipoles, qmmm_interaction_order)
                deallocate(multipoles)

#if defined(VAR_MPI)
                if (nprocs > 1) then
                    call mpi_bcast(mm_field_derivatives, size(mm_field_derivatives), rmpi, master, comm, ierr)
                end if
#endif

                work_size = multipole_size*num_operators / nprocs
                work_start = 1 + myid * work_size
                work_stop = (myid + 1) * work_size
                if (myid == nprocs - 1) work_stop = multipole_size*num_operators

                allocate(factors(multipole_size))
                l = 1
                do k = 0, qmmm_interaction_order
                    ncomp = (k+1)*(k+2)/2
                    call prefactors(factors(l:l+ncomp-1))
                    factors(l:l+ncomp-1) = factors(l:l+ncomp-1) * (-1)**k
                    l = l + ncomp
                end do

                do iwork = work_start, work_stop
                    i = 1 + ((iwork-1) / multipole_size)
                    j = 1 + mod(iwork - 1, multipole_size)
                    do k = 1, ndens
                        l = (k - 1) * nnbas + 1
                        m = k * nnbas
                        if (mixed) then
                            fckmats(l:m) = fckmats(l:m) + 0.5 * factors(j) * mm_field_derivatives(i, j, k) * multipole_interaction_operators(:, j, i)
                        else
                            fckmats(l:m) = fckmats(l:m) + factors(j) * mm_field_derivatives(i, j, k) * multipole_interaction_operators(:, j, i)
                        end if
                   end do
                end do

                deallocate(mm_field_derivatives)
                deallocate(factors)
            end if
        end if
        if (pelib_sol) then
            i = 1
            allocate(Vel_ints(nnbas,1, 1))
            do site = surp_start, surp_finish
                call Tk_integrals('potential_derivative', Rsp(:,site), Vel_ints)
                do k = 1, ndens
                    l = (k - 1) * nnbas + 1
                    m = k * nnbas
                    fckmats(l:m) = fckmats(l:m) + Qinds(i,k) * Vel_ints(:,1,1)
                end do
                i = i + 1
            end do
            deallocate(Vel_ints)
        end if
    end if

    if (any(lpol)) then
        deallocate(Mkinds, Fels, Ftots)
        deallocate(Fnucs, Fmuls, Ffds)
    end if
    if (pelib_sol) then
        deallocate(Qinds, Vels, Vtots)
        deallocate(Vnucs, Vmuls, Vfds)
    end if
end subroutine polarization_operator

!------------------------------------------------------------------------------

subroutine get_effdipole(fckmats)

    use pelib_mpi
    use pelib_options
    use pelib_constants
    use pelib_integral_interfaces
    use pelib_potential_derivatives
    use pelib_blas_interfaces, only: dot
    use pelib_utils, only: openfile
    use pelib_induced_moments, only: induced_moments, induced_moments_eef
    use pelib_multipole_tools, only: multipole_derivative

    real(rp), dimension(3) :: Rsc, Fp
    real(rp), dimension(3,3,nnucs+1+ncrds) :: lf_tensor
    real(rp), dimension(3,3) :: alpha_mol
    real(rp), dimension(:), intent(out) :: fckmats
    real(rp), dimension(:,:), allocatable :: Fext, Mkinds, Dqinds
    real(rp), dimension(:,:,:), allocatable :: Vel_ints, Fel_ints
    real(rp), dimension(:,:), allocatable :: nrmtmp
    real(rp) :: r_atom
    real(rp) :: rj_ts

    logical(lp), parameter :: local_debug = .false.
    integer(ip) :: i, l, m, j, k, site, atom
    real(rp) :: dqsum

    if (pelib_sol .and. .not. any(lpol)) then
        allocate(Dqinds(nsurp,3), Vel_ints(nnbas,1,1))
        Vel_ints = 0.0; Dqinds = 0.0
    else if (pelib_sol .and. any(lpol)) then
        error stop 'Coupling between fixsol and PE not implemented for EEF yet'
    else if (.not. pelib_sol .and. any(lpol)) then
        allocate(Fext(3*npols,3),Mkinds(3*npols,3),Fel_ints(nnbas,3,1))
        Fext = 0.0; Fel_ints = 0.0; Mkinds = 0.0
    else
        error stop 'You are doing something completely wrong. Back to the drawing board'
    end if

    if (pelib_sol) then
       if (myid == master) then
           allocate(nrmtmp(nsurp,3))
           nrmtmp = 0.0
           do i = 1, 3 ! Loop over x,y,z
               l = 1
               do j = 1, nsurp
                   if (Sa(i) <= zero) then
                        l = l + 1
                        cycle
                    end if
                   atom = idatm(l)
                   r_atom = all_coords(i,atom)
                   rj_ts   = Rsp(i,l)
                   nrmtmp(l,i) = (r_atom - rj_ts)/Radsp(atom)
                   l = l + 1
               end do
           end do
       end if
    end if

    if (any(lpol)) then
        l = 1
        do site = 1, nsites
            if (zeroalphas(site)) then
               cycle
            end if
            Fext(l,1) = 1.0
            Fext(l+1,2) = 1.0
            Fext(l+2,3) = 1.0
            l = l + 3
        end do
    end if

    if (pelib_sol .and. .not. any(lpol)) then
        if (myid == master) then
            call induced_moments_eef(Qinds=Dqinds, Vs=nrmtmp)
        end if
    else if (.not. pelib_sol .and. any(lpol)) then
        call induced_moments(Mkinds, Fext)
    else
        error stop 'You are doing something completely wrong. Back to the drawing board'
    end if

    if (pelib_debug .and. myid == master) then
       do j = 1, 3
           dqsum = 0.0
           do i = 1, nsurp
               dqsum = dqsum + Dqinds(i,j)
               write(luout,*) 'i,j, Dqinds(i,j), dqsum', i, j, Dqinds(i,j),dqsum
           end do
       end do
    end if

    if (local_debug .and. .not. pelib_sol) then
        write(luout, '(/15x,a)') 'Local field factor information:'
        write(luout, '(14x,a/)') '-------------------------------'
        do m = 1, 3
            j = 1
            write(luout, '(23x,a,i1)') 'Derivatives of induced dipole moments wrt. Fext_m: ', m
            write(luout, '(7x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                write(luout, '(5x,i6,3f15.8)') i, Mkinds(j:j+2,m)
                j = j + 3
            end do
        end do
    end if

    if (any(lpol)) then
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            do j = 1, 3
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + poldists(i-1)
                end do
                call mpi_scatterv(Mkinds(1,j), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
            end do
        else if (myid /= master) then
            do j = 1, 3
                call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,j), poldists(myid), rmpi, master, comm, ierr)
            end do
        end if
#endif

        l = 0
        do site = site_start, site_finish
            if (zeroalphas(site)) cycle
            call Tk_integrals('potential_derivative', Rs(:,site), Fel_ints)
            do m = 1, 3 !x, y, z dipole components
                do i = 1, 3 !x, y, z external field
                    j = (i - 1) * nnbas + 1
                    k = i * nnbas
                    fckmats(j:k) = fckmats(j:k) - Mkinds(l+m,i) * Fel_ints(:,m,1)
                end do
            end do
            l = l + 3
        end do

#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            call mpi_reduce(mpi_in_place, fckmats(1), 3 * nnbas , rmpi, mpi_sum, master, comm, ierr)
        else if (myid /= master) then
            call mpi_reduce(fckmats(1), dummy_real, 3 * nnbas, rmpi, mpi_sum, master, comm, ierr)
        end if
#endif
    end if

    if (pelib_sol) then
        j = 0
        if (myid == master) then
            do site = 1, nsurp
               call Tk_integrals('potential_derivative', Rsp(:,site), Vel_ints)
               do k = 1, 3
                  l = (k - 1) * nnbas + 1
                  m = k * nnbas
                  fckmats(l:m) = fckmats(l:m) + Dqinds(site,k) * Vel_ints(:,1,1)
               end do
               j = j + 1
            end do
            if (pelib_debug) then
                do i = 1, 3*nnbas
                    write(luout,*) 'eefmats from FixSol', fckmats(i),i
                end do
            end if
        end if
    end if

    if (.not. pelib_sol) then
        lf_tensor = 0.0
        alpha_mol = 0.0
        l = 1
        do site = site_start, site_finish
            if (zeroalphas(site)) cycle
!           lf_tensor at atomic sites
            do j = 1, nnucs
               Rsc = Rm(:,j) - Rs(:,site)
               do m = 1, 3
                  Fp = 0.0
                  call multipole_derivative(Fp, Rsc, Mkinds(l:l+2,m))
                  lf_tensor(:,m,j) = lf_tensor(:,m,j) + Fp
               end do
            end do
!           lf_tensor at COM
            Rsc = core_com(:) - Rs(:,site)
            do m = 1, 3
                Fp = 0.0
                call multipole_derivative(Fp, Rsc, Mkinds(l:l+2,m))
                lf_tensor(:,m,nnucs+1) = lf_tensor(:,m,nnucs+1) + Fp
            end do
!           lf_tensor at external sites given in input
            if (allocated(crds)) then
            write(luout,*) crds
                do j = 1, ncrds
                   Rsc = crds(:,j) - Rs(:,site)
                   do m = 1, 3
                      Fp = 0.0
                      call multipole_derivative(Fp, Rsc, Mkinds(l:l+2,m))
                      lf_tensor(:,m,nnucs+1+j) = lf_tensor(:,m,nnucs+1+j) + Fp
                   end do
                end do
            end if
            l = l + 3
        end do

#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            do m = 1, 3
                call mpi_reduce(mpi_in_place, alpha_mol(1,m), 3, rmpi, mpi_sum, master, comm, ierr)
                do j = 1, nnucs + 1 + ncrds
                    call mpi_reduce(mpi_in_place, lf_tensor(1,m,j), 3, rmpi, mpi_sum, master, comm, ierr)
                end do
            end do
        else if (myid /= master) then
            do m = 1, 3
                call mpi_reduce(alpha_mol(1,m), dummy_real, 3, rmpi, mpi_sum, master, comm, ierr)
                do j = 1, nnucs + 1 + ncrds
                    call mpi_reduce(lf_tensor(1,m,j), dummy_real, 3, rmpi, mpi_sum, master, comm, ierr)
                end do
            end do
        end if
#endif

        if (myid == master) then
            do j = 1, nnucs + 1 + ncrds
               do m = 1, 3
                   lf_tensor(m,m,j) = lf_tensor(m,m,j) + 1.0
               end do
            end do
        end if

        if (myid == master) then
            write(luout, '(/9x,a)') 'Effective external field (EEF) information:'
            write(luout, '(9x,a/)') '-------------------------------------------'

            do j = 1, nnucs + 1 + ncrds
               if (j > nnucs + 1) then
                   write(luout, '(20x,a,i5)') 'EEF tensor at site: ', j
               else if (j == nnucs + 1) then
                   write(luout, '(21x,a)') 'EEF tensor at COM'
               else
                   cycle
!                   write(luout, '(/20x,a,5i)') 'EEF tensor at nuclei: ', j
               end if
               do m = 1, 3
                   write(luout, '(5x,3f15.8)') (lf_tensor(i,m,j), i = 1,3)
               end do
            end do
            do j = 1, 3
                do m = 1, 3
                    alpha_mol(j,m) = dot(Fext(:,j),Mkinds(:,m))
                end do
            end do
            write(luout, '(/13x,a)') 'Environment polarizability (a.u.)'
            do m = 1, 3
                write(luout, '(5x,3f15.8)') (alpha_mol(m,i), i = 1,3)
            end do
        end if
    end if

    if (pelib_sol) then
        if (myid == master) then
           deallocate(Vel_ints)
           deallocate(Dqinds)
           deallocate(nrmtmp)
        end if
    end if
    if (any(lpol)) then
        deallocate(Fext)
        deallocate(Mkinds)
        deallocate(Fel_ints)
    end if

end subroutine get_effdipole

!------------------------------------------------------------------------------

subroutine get_molecular_gradient(denmats, fckmats, molgrads)

    use pelib_mpi
    use pelib_options

    real(rp), dimension(:), intent(in), optional :: denmats
    real(rp), dimension(:), intent(out), optional :: fckmats
    real(rp), dimension(:), intent(out), optional :: molgrads

    if (present(molgrads) .and. .not. present(denmats)) then
        error stop 'missing denmats in pelib_molecular_gradients'
    else if (present(denmats) .and. .not. present(molgrads)) then
        error stop 'missing molgrads in pelib_molecular_gradients'
    else if (.not. present(molgrads) .and. .not. present(fckmats)) then
        error stop 'nothing to do in pelib_molecular_gradients'
    end if

    if (any(lmul)) then
        if (present(molgrads) .and. present(fckmats)) then
            call es_gradients(denmats=denmats, fckmats=fckmats, molgrads=molgrads)
        else if (present(molgrads)) then
            call es_gradients(denmats=denmats, molgrads=molgrads)
        else if (present(fckmats)) then
            call es_gradients(fckmats=fckmats)
        end if
    end if
    ! true when having (iso)alphas
    if (lpol(1,1)) then
        if (present(molgrads) .and. present(fckmats)) then
            call pol_gradients(denmats=denmats, fckmats=fckmats, molgrads=molgrads)
        else if (present(molgrads)) then
            call pol_gradients(denmats=denmats, molgrads=molgrads)
        else if (present(fckmats)) then
            call pol_gradients(fckmats=fckmats)
        end if
    end if
    if (present(molgrads)) then
        if (lvdw) then
            call lj_gradients(molgrads)
        else
           write(luout, *) 'WARNING: No LJ contribution to geometrical derivatives is included'
        end if
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, molgrads(1), 3 * nnucs , rmpi, mpi_sum, master, comm, ierr)
    else if (myid /= master) then
        call mpi_reduce(molgrads(1), dummy_real, 3 * nnucs, rmpi, mpi_sum, master, comm, ierr)
    end if
#endif

end subroutine get_molecular_gradient

!------------------------------------------------------------------------------

subroutine pol_gradients(denmats, fckmats, molgrads)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_blas_interfaces, only: dot, spmv
    use pelib_multipole_tools, only: multipole_derivative
    use pelib_potential_derivatives
    use pelib_induced_moments, only: induced_moments

    real(rp), dimension(:), intent(in), optional :: denmats
    real(rp), dimension(:), intent(inout), optional :: molgrads
    real(rp), dimension(:), intent(inout), optional :: fckmats

    real(rp), dimension(3) :: Rms
    real(rp), dimension(6) :: grdFnucs
    real(rp), dimension(:), allocatable :: Fnucs, Fmuls !, Ffds
    real(rp), dimension(:,:), allocatable :: Mkinds, Fktots, Fels
    real(rp), dimension(:,:), allocatable :: Mk_ints

    integer(ip) :: i, j, k, l, site, lu
    logical(lp) :: lexist

    if (present(molgrads) .and. .not. present(denmats)) then
        error stop 'missing denmats in pol_gradients'
    else if (present(denmats) .and. .not. present(molgrads)) then
        error stop 'missing molgrads in pol_gradients'
    else if (.not. present(molgrads) .and. .not. present(fckmats)) then
        error stop 'nothing to do in pol_gradients'
    end if

    allocate(Mk_ints(nnbas,3*nnucs))
    allocate(Mkinds(3*npols,ndens))

    Mkinds = 0.0

!   Read in induced moments
    if (myid == master) then
        inquire(file='pelib_induced_moments.bin', exist=lexist)
        if (.not. lexist) then
            error stop 'pelib_induced_moments.bin does not exist in pol_gradients'
        end if
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist) then
        if (myid == master) then
            call openfile('pelib_induced_moments.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        end if
    else if (present(fckmats)) then
        error stop 'Updating induced dipoles according to a wrong density matrix'
    else
        allocate(Fels(3*npols,ndens), Fktots(3*npols,ndens))
        allocate(Fnucs(3*npols), Fmuls(3*npols))
        Fels = 0.0; Fktots = 0.0
        Fnucs = 0.0; Fmuls = 0.0
        if (trim(qmmm_interaction_scheme) == 'EXACT') then
            call electron_fields(Fels, denmats)
        else if (trim(qmmm_interaction_scheme) == 'FMM_FOCK') then
            call electron_fields_fmm(Fels, denmats)
        else
            ! ESPF, CESPF, SINGLE_CENTER
            call electron_fields_multipole(Fels, denmats)
        end if
        call nuclear_fields(Fnucs)
        call multipole_fields(Fmuls)
        do i = 1, ndens
            Fktots(1:3*npols,i) = Fels(:,i) + Fnucs + Fmuls !+ Ffds gradient of frozen densities not implemented yet.
        end do
        call induced_moments(Mkinds, Fktots)
        deallocate(Fels, Fktots, Fnucs, Fmuls)
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + poldists(i-1)
        end do
        call mpi_scatterv(Mkinds(1,1), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
    else if (myid /= master) then
        call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,1), poldists(myid), rmpi, master, comm, ierr)
    end if
#endif

    ! 1st derivatives of nuclear and electron fields multiplied with induced dipoles
    l = 1
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        if (present(molgrads)) then
            k = 1
            do i = 1, nnucs
                grdFnucs = 0.0
                Rms = Rs(:,site) - Rm(:,i)
                call multipole_derivative(grdFnucs, Rms, Zm(:,i))
                call spmv(grdFnucs, Mkinds(l:l+2,1), molgrads(k:k+2), 'L', -1.0_rp, 1.0_rp)
                k = k + 3
            end do
        end if

        Mk_ints = 0.0
        call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), Mkinds(l:l+2,1))
        l = l + 3

        if (present(molgrads)) then
            do i = 1, 3 * nnucs
                molgrads(i) = molgrads(i) + dot(denmats(1:nnbas), Mk_ints(:,i))
            end do
        end if
        if (present(fckmats)) then
            do i = 1, 3 * nnucs
                j = (i - 1) * nnbas + 1
                k = i * nnbas
                fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
            end do
        end if
    end do

    deallocate(Mkinds, Mk_ints)

end subroutine pol_gradients

!------------------------------------------------------------------------------

subroutine es_gradients(denmats, fckmats, molgrads)

!   Calculates gradient contribution from QM nuclear - multipole
!   and QM electron - multipole interaction energies

    use pelib_options
    use pelib_blas_interfaces, only: dot
    use pelib_multipole_tools, only: multipole_derivative

    real(rp), dimension(:), intent(in), optional :: denmats
    real(rp), dimension(:), intent(inout), optional :: fckmats
    real(rp), dimension(:), intent(inout), optional :: molgrads

    integer(ip) :: i, j, k, site
    real(rp), dimension(3) :: Rsm
    real(rp), dimension(3) :: mul_grad
    real(rp), dimension(:,:), allocatable :: Mk_ints

    if (present(molgrads) .and. .not. present(denmats)) then
        error stop 'missing denmats in es_gradients'
    else if (present(denmats) .and. .not. present(molgrads)) then
        error stop 'missing molgrads in es_gradients'
    else if (.not. present(molgrads) .and. .not. present(fckmats)) then
        error stop 'nothing to do in es_gradients'
    end if

    allocate(Mk_ints(nnbas,3*nnucs))

    Mk_ints = 0.0
    do site = site_start, site_finish
        ! 1st derivative of nuclei - multipole interaction energy
        if (present(molgrads)) then
            k = 1
            do i = 1, nnucs
                mul_grad = 0.0
                Rsm = Rm(:,i) - Rs(:,site)
                if (lmul(0)) then
                    if (maxval(abs(M0s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M0s(:,site))
                    end if
                end if
                if (lmul(1)) then
                    if (maxval(abs(M1s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M1s(:,site))
                    end if
                end if
                if (lmul(2)) then
                    if (maxval(abs(M2s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M2s(:,site))
                    end if
                end if
                if (lmul(3)) then
                    if (maxval(abs(M3s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M3s(:,site))
                    end if
                end if
                if (lmul(4)) then
                    if (maxval(abs(M4s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M4s(:,site))
                    end if
                end if
                if (lmul(5)) then
                    if (maxval(abs(M5s(:,site))) >= zero) then
                        call multipole_derivative(mul_grad, Rsm, M5s(:,site))
                    end if
                end if
                ! minus here because multipole_derivative gives opposite sign
                molgrads(k:k+2) = molgrads(k:k+2) - Zm(1,i) * mul_grad
                k = k + 3
            end do
        end if

        ! 1st derivative of electron - multipole interaction energy
        if (lmul(0)) then
            if (maxval(abs(M0s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M0s(:,site))
            end if
        end if
        if (lmul(1)) then
            if (maxval(abs(M1s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M1s(:,site))
            end if
        end if
        if (lmul(2)) then
            if (maxval(abs(M2s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M2s(:,site))
            end if
        end if
        if (lmul(3)) then
            if (maxval(abs(M3s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M3s(:,site))
            end if
        end if
        if (lmul(4)) then
            if (maxval(abs(M4s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M4s(:,site))
            end if
        end if
        if (lmul(5)) then
            if (maxval(abs(M5s(:,site))) >= zero) then
                call Mk_integrals('molecular_gradient', Mk_ints, Rs(:,site), M5s(:,site))
            end if
        end if
    end do
    if (present(molgrads)) then
        do i = 1, 3 * nnucs
            molgrads(i) = molgrads(i) + dot(denmats(1:nnbas), Mk_ints(:,i))
        end do
    end if
    if (present(fckmats)) then
        do i = 1, 3 * nnucs
            j = (i - 1) * nnbas + 1
            k = i * nnbas
            fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
        end do
    end if

    deallocate(Mk_ints)

end subroutine es_gradients

!------------------------------------------------------------------------------

subroutine lj_energy

    use pelib_mpi
    use pelib_options
    use pelib_constants
    use pelib_blas_interfaces, only: nrm2

    integer(ip) :: i, site
    real(rp) :: eps_sm, r_sm, term_6, norm
    real(rp), dimension(3) :: Rsm

! TODO: Also be able to handle bond midpoints
    Elj    = 0.0

    if (qmLJsites /= nnucs) error stop 'Number of LJ sites does not match the&
                                 & number of nuclei in the core region.'

    do site = site_start, site_finish
        if (Zs(1,site) <= zero) cycle
        do i = 1, nnucs
            r_sm = 0.5 * (LJs(1,site) + qmLJs(1,i)) * aa2bohr
            eps_sm = sqrt(LJs(2,site) * qmLJs(2,i)) * kcal2hartree
            Rsm = Rm(:,i) - Rs(:,site)
            norm = nrm2(Rsm)
            term_6 = (r_sm / norm)**(6.0)
            Elj(1) = Elj(1) - 2.0 * eps_sm * term_6
            Elj(2) = Elj(2) + eps_sm * term_6**(2.0)
        end do
    end do

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, Elj(1), 2, rmpi, mpi_sum, master, comm, ierr)
    else if (myid /= master) then
        call mpi_reduce(Elj(1), dummy_real, 2, rmpi, mpi_sum, master, comm, ierr)
    end if
#endif

end subroutine lj_energy

!------------------------------------------------------------------------------

subroutine lj_gradients(molgrads)

    use pelib_options
    use pelib_constants
    use pelib_blas_interfaces, only: nrm2

    real(rp), dimension(:), intent(inout) :: molgrads

    integer(ip) :: i, j, site
    real(rp) :: eps_sm, r_sm, term_6, norm
    real(rp), dimension(3) :: fac, Rsm

    do site = site_start, site_finish
        if (Zs(1,site) <= zero) cycle
        j = 1
        do i = 1, nnucs
            r_sm = 0.5 * (LJs(1,site) + qmLJs(1,i)) * aa2bohr
            eps_sm = (sqrt(LJs(2,site) * qmLJs(2,i))) * kcal2hartree
            Rsm = Rm(:,i) - Rs(:,site)
            norm = nrm2(Rsm)
            fac = - 12.0 * eps_sm * Rsm / (norm**(2.0))
            term_6 = (r_sm / norm)**(6.0)
            molgrads(j:j+2) = molgrads(j:j+2) + fac * (term_6**(2.0) - term_6)
            j = j + 3
        end do
    end do

end subroutine lj_gradients

!------------------------------------------------------------------------------

subroutine get_magnetic_gradient(fckmats)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile

    real(rp), dimension(:), intent(out) :: fckmats

    integer(ip) :: i, lu
    real(rp), dimension(:,:), allocatable :: Mkinds
    real(rp), dimension(:,:), allocatable :: Qinds
    logical(lp) :: lexist

    ! static multipole moment contribution
    if (lmul(0)) call lao_multipoles(M0s, fckmats)
    if (lmul(1)) call lao_multipoles(M1s, fckmats)
    if (lmul(2)) call lao_multipoles(M2s, fckmats)
    if (lmul(3)) call lao_multipoles(M3s, fckmats)
    if (lmul(4)) call lao_multipoles(M4s, fckmats)
    if (lmul(5)) call lao_multipoles(M5s, fckmats)

    ! cosmo contribution
    if(pelib_sol) then
        allocate(Qinds(1,nsurp))
        if (myid == master) then
            inquire(file='pelib_induced_charges.bin', exist=lexist)
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
        end if
#endif
        if (lexist) then
            if (myid == master) then
                call openfile('pelib_induced_charges.bin', lu, 'old', 'unformatted')
                rewind(lu)
                read(lu) Qinds
                close(lu)
            end if
        else
            error stop 'pelib_induced_charges.bin does not exist'
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(Qinds(1,1), nsurp, rmpi, master, comm, ierr)
        end if
#endif
        call lao_induced_charges(Qinds, fckmats)
        deallocate(Qinds, stat=ierr)
    endif

    ! polarization contribution
    if (lpol(1,1)) then
        allocate(Mkinds(3,npols))
        if (myid == master) then
            inquire(file='pelib_induced_moments.bin', exist=lexist)
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
        end if
#endif
        if (lexist) then
            if (myid == master) then
                call openfile('pelib_induced_moments.bin', lu, 'old', 'unformatted')
                rewind(lu)
                read(lu) Mkinds
                close(lu)
            end if
        else
            error stop 'pelib_induced_moments.bin does not exist'
        end if
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + poldists(i-1)
            end do
            call mpi_scatterv(Mkinds(1,1), poldists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_scatterv(dummy_real, poldists, displs, rmpi, Mkinds(1,1), poldists(myid), rmpi, master, comm, ierr)
        end if
#endif
        call lao_induced_multipoles(Mkinds, fckmats)
        deallocate(Mkinds)
    endif

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        call mpi_reduce(mpi_in_place, fckmats(1), 3 * nnbas, rmpi, mpi_sum, master, comm, ierr)
    else if (myid /= master) then
        call mpi_reduce(fckmats(1), dummy_real, 3 * nnbas, rmpi, mpi_sum, master, comm, ierr)
    end if
#endif

end subroutine get_magnetic_gradient

!------------------------------------------------------------------------------

subroutine lao_induced_charges(Qs, fckmats)

    use pelib_options

    real(rp), dimension(:,:), intent(in) :: Qs
    real(rp), dimension(:), intent(inout) :: fckmats

    integer(ip) :: site
    integer(ip) :: i, j, k
    real(rp), dimension(:,:), allocatable :: Mk_ints

    allocate(Mk_ints(nnbas,3))

    do site = surp_start, surp_finish
        Mk_ints = 0.0
        call Mk_integrals('magnetic_gradient', Mk_ints, Rsp(:,site), Qs(:,site))
        ! update for all directions of the magnetic field
        do i = 1, 3
            j = (i - 1) * nnbas + 1
            k = i * nnbas
            fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
        end do
    end do

    deallocate(Mk_ints)

end subroutine lao_induced_charges

!------------------------------------------------------------------------------

subroutine lao_multipoles(Mks, fckmats)

    use pelib_options

    real(rp), dimension(:,:), intent(in) :: Mks
    real(rp), dimension(:), intent(inout) :: fckmats

    integer(ip) :: i, j, k
    integer(ip) :: site
    real(rp), dimension(:,:), allocatable :: Mk_ints

    allocate(Mk_ints(nnbas,3))

    Mk_ints = 0.0
    do site = site_start, site_finish
        if (maxval(abs(Mks(:,site))) < zero) cycle
        call Mk_integrals('magnetic_gradient', Mk_ints, Rs(:,site), Mks(:,site))
    end do

    ! update for all directions of the magnetic field
    do i = 1, 3
        j = (i - 1) * nnbas + 1
        k = i * nnbas
        fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
    end do

    deallocate(Mk_ints)

end subroutine lao_multipoles

!------------------------------------------------------------------------------

subroutine lao_induced_multipoles(Mks, fckmats)

    use pelib_options

    real(rp), dimension(:,:), intent(in) :: Mks
    real(rp), dimension(:), intent(inout) :: fckmats

    integer(ip) :: i, j, k
    integer(ip) :: site
    real(rp), dimension(:,:), allocatable :: Mk_ints

    allocate(Mk_ints(nnbas,3))

    Mk_ints = 0.0
    i = 1
    do site = site_start, site_finish
        if (zeroalphas(site)) cycle
        call Mk_integrals('magnetic_gradient', Mk_ints, Rs(:,site), Mks(:,i))
        i = i + 1
    end do

    ! update for all directions of the magnetic field
    do i = 1, 3
        j = (i - 1) * nnbas + 1
        k = i * nnbas
        fckmats(j:k) = fckmats(j:k) + Mk_ints(:,i)
    end do

    deallocate(Mk_ints)

end subroutine lao_induced_multipoles

!------------------------------------------------------------------------------

subroutine Mk_integrals(inttype, Mk_ints, Rij, Mk)

    use pelib_options
    use pelib_integral_interfaces
    use pelib_multipole_tools, only: prefactors

    character(*), intent(in) :: inttype
    real(rp), dimension(:,:), intent(inout) :: Mk_ints
    real(rp), dimension(3), intent(in) :: Rij
    real(rp), dimension(:), intent(in) :: Mk

    integer(ip) :: i, j
    integer(ip) :: ncomps, nprops
    real(rp), dimension(:), allocatable :: factors
    real(rp), dimension(:,:,:), allocatable :: Tk_ints

    ncomps = size(Mk)

    if (inttype == 'potential_derivative') then
        nprops = 1
    else if (inttype == 'magnetic_gradient') then
        nprops = 3
    else if (inttype == 'molecular_gradient') then
        nprops = 3 * nnucs
    end if

    allocate(Tk_ints(nnbas,ncomps,nprops))

    call Tk_integrals(inttype, Rij, Tk_ints)

    allocate(factors(ncomps))
    call prefactors(factors)

    ! multiply T^(k) integrals with multipole to get M^(k) integrals
    do j = 1, nprops
        do i = 1, ncomps
            Mk_ints(:,j) = Mk_ints(:,j) + factors(i) * Mk(i) * Tk_ints(:,i,j)
        end do
    end do

    deallocate(factors)
    deallocate(Tk_ints)

end subroutine Mk_integrals

!------------------------------------------------------------------------------

end module pelib_operators
