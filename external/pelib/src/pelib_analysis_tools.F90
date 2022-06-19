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
module pelib_analysis_tools

#ifdef VAR_PDE
    use hdf5
#endif

    use pelib_precision

    implicit none

    private

    public :: compute_potential_cube, compute_field_cube
    public :: compute_mep_cube, compute_mef_cube, compute_mep_input

contains

!------------------------------------------------------------------------------

subroutine compute_potential_cube(denmats)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: multipole_derivative
    use pelib_induced_moments, only: induced_moments
    use pelib_potential_derivatives

    real(rp), dimension(:), intent(in), optional :: denmats
    integer(ip) :: i, j, k, l
    integer(ip) :: lu
    logical(lp) :: lexist
    character(len=3) :: idx_str
    real(rp), dimension(1) :: Vp
    real(rp), dimension(3) :: Rsg
    real(rp), dimension(:), allocatable :: Vpe
    real(rp), dimension(:,:), allocatable :: M1inds
    real(rp), dimension(:), allocatable :: Fnucs, Fmuls
    real(rp), dimension(:,:), allocatable :: Fels

    if (myid == master) then
        allocate(Vpe(npoints))
    else if (myid /= master) then
        allocate(Vpe(cubedists(myid)))
    end if

    Vpe = 0.0

    if (present(denmats)) then
        if (any(lpol)) then
            allocate(M1inds(3*npols,1))
            allocate(Fels(3*npols,1))
            Fels = 0.0
            call electron_fields(Fels, denmats)
            if (state_idx==0) then
                allocate(Fnucs(3*npols),Fmuls(3*npols))
                Fnucs = 0.0; Fmuls = 0.0
                call nuclear_fields(Fnucs)
                call multipole_fields(Fmuls)
                Fels(:,1) = Fels(:,1) + Fnucs + Fmuls
                deallocate(Fnucs,Fmuls)
            end if
            call induced_moments(Mkinds=M1inds, Fks=Fels)
            deallocate(Fels)
            write(luout, '(/15x,a)') 'Polarizable embedding information:'
            write(luout, '(14x,a/)') '------------------------------------'
            write(luout, '(23x,a)') 'Induced dipole moments'
            write(luout, '(7x,a,10x,a,14x,a,14x,a)') 'site', 'X', 'Y', 'Z'
            j = 1
            do i = 1, nsites
                if (zeroalphas(i)) cycle
                write(luout, '(5x,i6,3f15.8)') i, M1inds(j:j+2,1)
                j = j + 3
            end do
        end if
    end if

    if (any(lmul) .and. .not. (state_idx>=1)) then
        k = 1
        do i = cube_start, cube_finish
            Vp = 0.0
            do j = 1, nsites
                Rsg = Rg(:,i) - Rs(:,j)
                if (lmul(0)) then
                    call multipole_derivative(Vp, Rsg, M0s(:,j))
                end if
                if (lmul(1)) then
                    call multipole_derivative(Vp, Rsg, M1s(:,j))
                end if
                if (lmul(2)) then
                    call multipole_derivative(Vp, Rsg, M2s(:,j))
                end if
                if (lmul(3)) then
                    call multipole_derivative(Vp, Rsg, M3s(:,j))
                end if
                if (lmul(4)) then
                    call multipole_derivative(Vp, Rsg, M4s(:,j))
                end if
                if (lmul(5)) then
                    call multipole_derivative(Vp, Rsg, M5s(:,j))
                end if
            end do
            Vpe(k) = Vpe(k) + Vp(1)
            k = k + 1
        end do
    end if

    if (lpol(1,1)) then
        if (state_idx==0) then
            allocate(M1inds(3*npols,1))
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
                    read(lu) M1inds
                    close(lu)
                end if
            else
                write(luout, *) 'WARNING: cannot create cube file because pelib_induced_moments.bin not found'
                return
            end if
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds(1,1), 3 * npols, rmpi, master, comm, ierr)
        end if
#endif
        k = 1
        do i = cube_start, cube_finish
            l = 1
            Vp = 0.0
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsg = Rg(:,i) - Rs(:,j)
                call multipole_derivative(Vp, Rsg, M1inds(l:l+2,1))
                l = l + 3
            end do
            Vpe(k) = Vpe(k) + Vp(1)
            k = k + 1
        end do
        deallocate(M1inds)
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + cubedists(i-1)
        end do
        call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vpe, cubedists, displs, rmpi, master, comm, ierr)
    else if (myid /= master) then
        call mpi_gatherv(Vpe(1), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
    end if
#endif
    if (myid == master) then
        write(idx_str, '(I0.3)') state_idx
        call openfile('embedding_potential_'//idx_str//'.cube', lu, 'unknown', 'formatted')
        write(lu, '(a)') 'Embedding potential'
        write(lu, '(a)') 'Generated by the PE library'
        write(lu, '(i5,3f12.6)') nnucs, origin
        write(lu, '(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
        write(lu, '(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
        write(lu, '(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
        do j = 1, nnucs
            write(lu, '(i5,4f12.6)') int(Zm(1,j)), Zm(1,j), Rm(:,j)
        end do
        do i = 1, xsteps * ysteps
            j = (i - 1) * zsteps + 1
            k = j - 1 + zsteps
            write(lu, '(6e13.5)') Vpe(j:k)
        end do
        close(lu)
    end if

    deallocate(Vpe)

end subroutine compute_potential_cube

!------------------------------------------------------------------------------

subroutine compute_field_cube()

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_multipole_tools, only: multipole_derivative

    character(len=1) :: tcl
    character(len=99) :: cl
    integer(ip) :: i, j, k, l
    integer(ip) :: lu
    logical(lp) :: lexist
    real(rp), dimension(3) :: Rsg, Fp
    real(rp), dimension(:,:), allocatable :: Fpe, M1inds


    if (myid == master) then
        allocate(Fpe(npoints,3))
    else if (myid /= master) then
        allocate(Fpe(cubedists(myid),3))
    end if

    Fpe = 0.0

    if (any(lmul)) then
        k = 1
        do i = cube_start, cube_finish
            do j = 1, nsites
                Rsg = Rg(:,i) - Rs(:,j)
                if (lmul(0)) then
                    Fp = 0.0
                    call multipole_derivative(Fp, Rsg, M0s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(1)) then
                    Fp = 0.0
                    call multipole_derivative(Fp, Rsg, M1s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(2)) then
                    Fp = 0.0
                    call multipole_derivative(Fp, Rsg, M2s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(3)) then
                    Fp = 0.0
                    call multipole_derivative(Fp, Rsg, M3s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(4)) then
                    Fp = 0.0
                    call multipole_derivative(Fp, Rsg, M4s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
                if (lmul(5)) then
                    Fp = 0.0
                    call multipole_derivative(Fp, Rsg, M5s(:,j))
                    Fpe(k,:) = Fpe(k,:) + Fp
                end if
            end do
            k = k + 1
        end do
    end if

    if (lpol(1,1)) then
        allocate(M1inds(3*npols,1))
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
                read(lu) M1inds
                close(lu)
            end if
        else
            write(luout, *) 'WARNING: cannot create cube file because pelib_induced_moments.bin not found'
            return
        end if
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds(1,1), 3 * npols, rmpi, master, comm, ierr)
        end if
#endif
        k = 1
        do i = cube_start, cube_finish
            l = 1
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsg = Rg(:,i) - Rs(:,j)
                Fp = 0.0
                call multipole_derivative(Fp, Rsg, M1inds(l:l+2,1))
                Fpe(k,:) = Fpe(k,:) + Fp
                l = l + 3
            end do
            k = k + 1
        end do
        deallocate(M1inds)
    end if

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + cubedists(i-1)
        end do
        do i = 1, 3
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fpe(1,i), cubedists, displs, rmpi, master, comm, ierr)
        end do
    else if (myid /= master) then
        do i = 1, 3
            call mpi_gatherv(Fpe(1,i), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
        end do
    end if
#endif
    if (myid == master) then
        do l = 1, 3
            write(cl, *) l
            tcl = trim(adjustl(cl))
            call openfile('embedding_field_'//tcl//'.cube', lu, 'unknown', 'formatted')
            write(lu, '(a)') 'Embedding potential electric field component '//tcl
            write(lu, '(a)') 'Generated by the PE library'
            write(lu, '(i5,3f12.6)') nnucs, origin
            write(lu, '(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
            write(lu, '(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
            write(lu, '(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
            do j = 1, nnucs
                write(lu, '(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                write(lu, '(6e13.5)') Fpe(j:k,l)
            end do
            close(lu)
        end do
    end if

    deallocate(Fpe)

end subroutine compute_field_cube

!------------------------------------------------------------------------------

subroutine compute_mep_cube(denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces
    use pelib_blas_interfaces
    use pelib_utils, only: openfile
    use pelib_induced_moments, only: induced_moments
    use pelib_multipole_tools, only: multipole_derivative, multipole_derivative_damped

    real(rp), dimension(:), intent(in), optional :: denmats

    integer(ip) :: mulorder
    integer(ip) :: i, j, k, l
    integer(ip) :: lu, lum0, lum1, lum2, lum3, lum4, lum5
    logical(lp) :: exclude
    real(rp) :: alpha_i, alpha_j
    real(rp), parameter :: i3 = 1.0 / 3.0
    real(rp), dimension(1) :: Vs
    real(rp), dimension(3) :: Rji, Rsg, Rmp
    real(rp), dimension(:,:), allocatable :: Vqm, Vpe, Vind, Fmuls, M1inds
    real(rp), dimension(:,:,:), allocatable :: Tk_ints

    if (mep_qmcube) then
        if (.not. present(denmats)) then
            error stop 'Core-region MEP requested but no density matrix'
        end if
        if (myid == master) then
            allocate(Vqm(npoints,1))
        else if (myid /= master) then
            allocate(Vqm(cubedists(myid),1))
        end if
        Vqm = 0.0
        allocate(Tk_ints(nnbas,1,1))
        k = 1
        do i = cube_start, cube_finish
            call Tk_integrals('potential_derivative', Rg(:,i), Tk_ints)
            Vqm(k,1) = dot(denmats, Tk_ints(:,1,1))
            Vs = 0.0
            do j = 1, nnucs
                Rmp = Rg(:,i) - Rm(:,j)
                call multipole_derivative(Vs, Rmp, Zm(:,j))
            end do
            Vqm(k,1) = Vqm(k,1) + Vs(1)
            k = k + 1
        end do
        deallocate(Tk_ints)
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + cubedists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vqm(1,1), cubedists, displs, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_gatherv(Vqm(1,1), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            call openfile('qm_mep.cube', lu, 'new', 'formatted')
            write(lu,'(a)') 'Quantum molecular electrostatic potential'
            write(lu,'(a)') 'Generated by the PE library'
            write(lu,'(i5,3f12.6)') nnucs, origin
            write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
            write(lu,'(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
            write(lu,'(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
            do j = 1, nnucs
                write(lu,'(i5,4f12.6)') int(Zm(1,j)), Zm(1,j), Rm(:,j)
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                write(lu,'(6e13.5)') Vqm(j:k,1)
            end do
            close(lu)
        end if
        deallocate(Vqm)
    end if

    if (any(lmul) .and. mep_mulcube) then
        if (lmul(0)) mulorder = 0
        if (lmul(1)) mulorder = 1
        if (lmul(2)) mulorder = 2
        if (lmul(3)) mulorder = 3
        if (lmul(4)) mulorder = 4
        if (lmul(5)) mulorder = 5
        if (myid == master) then
            allocate(Vpe(npoints,0:mulorder))
        else if (myid /= master) then
            allocate(Vpe(cubedists(myid),0:mulorder))
        end if
        Vpe = 0.0
        k = 1
        do i = cube_start, cube_finish
            do j = 1, nsites
                Rsg = Rg(:,i) - Rs(:,j)
                if (lmul(0)) then
                    Vs = 0.0
                    call multipole_derivative(Vs, Rsg, M0s(:,j))
                    Vpe(k,0) = Vpe(k,0) + Vs(1)
                end if
                if (lmul(1)) then
                    Vs = 0.0
                    call multipole_derivative(Vs, Rsg, M1s(:,j))
                    Vpe(k,1) = Vpe(k,1) + Vs(1)
                end if
                if (lmul(2)) then
                    Vs = 0.0
                    call multipole_derivative(Vs, Rsg, M2s(:,j))
                    Vpe(k,2) = Vpe(k,2) + Vs(1)
                end if
                if (lmul(3)) then
                    Vs = 0.0
                    call multipole_derivative(Vs, Rsg, M3s(:,j))
                    Vpe(k,3) = Vpe(k,3) + Vs(1)
                end if
                if (lmul(4)) then
                    Vs = 0.0
                    call multipole_derivative(Vs, Rsg, M4s(:,j))
                    Vpe(k,4) = Vpe(k,4) + Vs(1)
                end if
                if (lmul(5)) then
                    Vs = 0.0
                    call multipole_derivative(Vs, Rsg, M5s(:,j))
                    Vpe(k,5) = Vpe(k,5) + Vs(1)
                end if
            end do
            k = k + 1
        end do
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + cubedists(i-1)
            end do
            do i = 0, mulorder
                call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vpe(1,i), cubedists, displs, rmpi, master, comm, ierr)
            end do
        else if (myid /= master) then
            do i = 0, mulorder
                call mpi_gatherv(Vpe(1,i), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
            end do
        end if
#endif
        if (myid == master) then
            if (mulorder >= 0) then
                call openfile('m0_mep.cube', lum0, 'new', 'formatted')
            end if
            if (mulorder >= 1) then
                call openfile('m1_mep.cube', lum1, 'new', 'formatted')
            end if
            if (mulorder >= 2) then
                call openfile('m2_mep.cube', lum2, 'new', 'formatted')
            end if
            if (mulorder >= 3) then
                call openfile('m3_mep.cube', lum3, 'new', 'formatted')
            end if
            if (mulorder >= 4) then
                call openfile('m4_mep.cube', lum4, 'new', 'formatted')
            end if
            if (mulorder >= 5) then
                call openfile('m5_mep.cube', lum5, 'new', 'formatted')
            end if
            do i = 1, mulorder + 1
                if (i == 1) then
                    lu = lum0
                else if (i == 2) then
                    lu = lum1
                else if (i == 3) then
                    lu = lum2
                else if (i == 4) then
                    lu = lum3
                else if (i == 5) then
                    lu = lum4
                else if (i == 6) then
                    lu = lum5
                end if
                write(lu,'(a)') 'Classical molecular electrostatic potential'
                write(lu,'(a)') 'Generated by the PE library'
                write(lu,'(i5,3f12.6)') nnucs, origin
                write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
                write(lu,'(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
                write(lu,'(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
                do j = 1, nnucs
                    write(lu,'(i5,4f12.6)') int(Zm(1,j)), Zm(1,j), Rm(:,j)
                end do
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                if (mulorder >= 0) then
                    write(lum0,'(6e13.5)') Vpe(j:k,0)
                end if
                if (mulorder >= 1) then
                    Vpe(j:k,0) = Vpe(j:k,0) + Vpe(j:k,1)
                    write(lum1,'(6e13.5)') Vpe(j:k,0)
                end if
                if (mulorder >= 2) then
                    Vpe(j:k,0) = Vpe(j:k,0) + Vpe(j:k,2)
                    write(lum2,'(6e13.5)') Vpe(j:k,0)
                end if
                if (mulorder >= 3) then
                    Vpe(j:k,0) = Vpe(j:k,0) + Vpe(j:k,3)
                    write(lum3,'(6e13.5)') Vpe(j:k,0)
                end if
                if (mulorder >= 4) then
                    Vpe(j:k,0) = Vpe(j:k,0) + Vpe(j:k,4)
                    write(lum4,'(6e13.5)') Vpe(j:k,0)
                end if
                if (mulorder >= 5) then
                    Vpe(j:k,0) = Vpe(j:k,0) + Vpe(j:k,5)
                    write(lum5,'(6e13.5)') Vpe(j:k,0)
                end if
            end do
            if (mulorder >= 0) then
                close(lum0)
            end if
            if (mulorder >= 1) then
                close(lum1)
            end if
            if (mulorder >= 2) then
                close(lum2)
            end if
            if (mulorder >= 3) then
                close(lum3)
            end if
            if (mulorder >= 4) then
                close(lum4)
            end if
            if (mulorder >= 5) then
                close(lum5)
            end if
        endif
        deallocate(Vpe)
    end if

    if (lpol(1,1) .and. mep_mulcube) then
        allocate(Fmuls(3*npols,1))
        Fmuls = 0.0
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            if (pelib_mul_damp .or. pelib_amoeba_damp) then
                alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
            end if
            do j = 1, nsites
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
                    if (alpha_j <= zero) then
                        error stop 'damping requires polarizabilities on all sites with multipoles'
                    end if
                end if
                Rji = Rs(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (maxval(abs(M0s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M0s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M0s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M0s(:,j))
                        end if
                    end if
                end if
                if (lmul(1)) then
                    if (maxval(abs(M1s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M1s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M1s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M1s(:,j))
                        end if
                    end if
                end if
                if (lmul(2)) then
                    if (maxval(abs(M2s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M2s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M2s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M2s(:,j))
                        end if
                    end if
                end if
                if (lmul(3)) then
                    if (maxval(abs(M3s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M3s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M3s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M3s(:,j))
                        end if
                    end if
                end if
                if (lmul(4)) then
                    if (maxval(abs(M4s(:,j))) >= zero) then
                        if (pelib_mul_damp .or. pelib_amoeba_damp) then
                            error stop 'damping not implemented for fourth order multipoles'
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M4s(:,j))
                        end if
                    end if
                end if
                if (lmul(5)) then
                    if (maxval(abs(M5s(:,j))) >= zero) then
                        if (pelib_mul_damp .or. pelib_amoeba_damp) then
                            error stop 'damping not implemented for fifth order multipoles'
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M5s(:,j))
                        end if
                    end if
                end if
            end do
            l = l + 3
        end do
        allocate(M1inds(3*npols,1))
        if (mep_extfld) then
            j = 1
            do i = 1, npols
                Fmuls(j:j+2,1) = Fmuls(j:j+2,1) + extfld
                j = j + 3
            end do
        end if
        call induced_moments(M1inds, Fmuls)
        deallocate(Fmuls)
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds(1,1), 3 * npols, rmpi, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            allocate(Vind(npoints,1))
        else if (myid /= master) then
            allocate(Vind(cubedists(myid),1))
        end if
        Vind = 0.0
        k = 1
        do i = cube_start, cube_finish
            l = 1
            Vs = 0.0
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsg = Rg(:,i) - Rs(:,j)
                call multipole_derivative(Vs, Rsg, M1inds(l:l+2,1))
                l = l + 3
            end do
            Vind(k,1) = Vind(k,1) + Vs(1)
            k = k + 1
        end do
        deallocate(M1inds)
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + cubedists(i-1)
            end do
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vind(1,1), cubedists, displs, rmpi, master, comm, ierr)
        else if (myid /= master) then
            call mpi_gatherv(Vind(1,1), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            call openfile('ind_mep.cube', lu, 'new', 'formatted')
            write(lu,'(a)') 'Classical induced molecular electrostatic potential'
            write(lu,'(a)') 'Generated by the PE library'
            write(lu,'(i5,3f12.6)') nnucs, origin
            write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
            write(lu,'(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
            write(lu,'(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
            do j = 1, nnucs
                write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
            end do
            do i = 1, xsteps * ysteps
                j = (i - 1) * zsteps + 1
                k = j - 1 + zsteps
                write(lu,'(6e13.5)') Vind(j:k,1)
            end do
            close(lu)
        end if
        deallocate(Vind)
    end if

end subroutine compute_mep_cube

!------------------------------------------------------------------------------

subroutine compute_mef_cube(denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces
    use pelib_blas_interfaces
    use pelib_utils, only: openfile
    use pelib_potential_derivatives, only: multipole_fields
    use pelib_induced_moments, only: induced_moments
    use pelib_multipole_tools, only: multipole_derivative, multipole_derivative_damped

    real(rp), dimension(:), intent(in), optional :: denmats

    character(len=1) :: tcl
    character(len=99) :: cl
    integer(ip) :: mulorder
    integer(ip) :: i, j, k, l
    integer(ip) :: lu, lum0, lum1, lum2, lum3, lum4, lum5
    logical(lp) :: exclude
    real(rp) :: alpha_i, alpha_j
    real(rp), parameter :: i3 = 1.0 / 3.0
    real(rp), dimension(3) :: Rji, Rsg, Rmg, Fs
    real(rp), dimension(:,:), allocatable :: Fqm, Find, Fmuls, M1inds
    real(rp), dimension(:,:,:), allocatable :: Fpe, Tk_ints

    if (mep_qmcube) then
        if (myid == master) then
            allocate(Fqm(npoints,3))
        else if (myid /= master) then
            allocate(Fqm(cubedists(myid),3))
        end if
        Fqm = 0.0
        allocate(Tk_ints(nnbas,3,1))
        k = 1
        do i = cube_start, cube_finish
            call Tk_integrals('potential_derivative', Rg(:,i), Tk_ints)
            do j = 1, 3
                Fqm(k,j) = dot(denmats, Tk_ints(:,j,1))
            end do
            do j = 1, nnucs
                Rmg = Rg(:,i) - Rm(:,j)
                Fs = 0.0
                call multipole_derivative(Fs, Rmg, Zm(:,j))
                Fqm(k,:) = Fqm(k,:) + Fs
            end do
            k = k + 1
        end do
        deallocate(Tk_ints)
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + cubedists(i-1)
            end do
            do i = 1, 3
            call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fqm(1,i), cubedists, displs, rmpi, master, comm, ierr)
            end do
        else if (myid /= master) then
            do i = 1, 3
                call mpi_gatherv(Fqm(1,i), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
            end do
        end if
#endif
        if (myid == master) then
            do l = 1, 3
                write(cl,*) l
                tcl = trim(adjustl(cl))
                call openfile('qm_field_'//tcl//'.cube', lu, 'new', 'formatted')
                write(lu,'(a)') 'Quantum electric field component '//tcl
                write(lu,'(a)') 'Generated by the PE library'
                write(lu,'(i5,3f12.6)') nnucs, origin
                write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
                write(lu,'(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
                write(lu,'(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
                do j = 1, nnucs
                    write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
                end do
                do i = 1, xsteps * ysteps
                    j = (i - 1) * zsteps + 1
                    k = j - 1 + zsteps
                    write(lu,'(6e13.5)') Fqm(j:k,l)
                end do
                close(lu)
            end do
        end if
        deallocate(Fqm)
    end if

    if (any(lmul) .and. mep_mulcube) then
        if (lmul(0)) mulorder = 0
        if (lmul(1)) mulorder = 1
        if (lmul(2)) mulorder = 2
        if (lmul(3)) mulorder = 3
        if (lmul(4)) mulorder = 4
        if (lmul(5)) mulorder = 5
        if (myid == master) then
            allocate(Fpe(npoints,3,0:mulorder))
        else if (myid /= master) then
            allocate(Fpe(cubedists(myid),3,0:mulorder))
        end if
        Fpe = 0.0
        k = 1
        do i = cube_start, cube_finish
            do j = 1, nsites
                Rsg = Rg(:,i) - Rs(:,j)
                if (lmul(0)) then
                    Fs = 0.0
                    call multipole_derivative(Fs, Rsg, M0s(:,j))
                    Fpe(k,:,0) = Fpe(k,:,0) + Fs
                end if
                if (lmul(1)) then
                    Fs = 0.0
                    call multipole_derivative(Fs, Rsg, M1s(:,j))
                    Fpe(k,:,1) = Fpe(k,:,1) + Fs
                end if
                if (lmul(2)) then
                    Fs = 0.0
                    call multipole_derivative(Fs, Rsg, M2s(:,j))
                    Fpe(k,:,2) = Fpe(k,:,2) + Fs
                end if
                if (lmul(3)) then
                    Fs = 0.0
                    call multipole_derivative(Fs, Rsg, M3s(:,j))
                    Fpe(k,:,3) = Fpe(k,:,3) + Fs
                end if
                if (lmul(4)) then
                    Fs = 0.0
                    call multipole_derivative(Fs, Rsg, M4s(:,j))
                    Fpe(k,:,4) = Fpe(k,:,4) + Fs
                end if
                if (lmul(5)) then
                    Fs = 0.0
                    call multipole_derivative(Fs, Rsg, M5s(:,j))
                    Fpe(k,:,5) = Fpe(k,:,5) + Fs
                end if
            end do
            k = k + 1
        end do
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + cubedists(i-1)
            end do
            do i = 1, 3
                do j = 0, mulorder
                    call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Fpe(1,i,j), cubedists, displs, rmpi, master, comm, ierr)
                end do
            end do
        else if (myid /= master) then
            do i = 1, 3
                do j = 0, mulorder
                    call mpi_gatherv(Fpe(1,i,j), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
                end do
            end do
        end if
#endif
        if (myid == master) then
            do l = 1, 3
                write(cl,*) l
                tcl = trim(adjustl(cl))
                if (mulorder >= 0) then
                    call openfile('m0_field_'//tcl//'.cube', lum0, 'new', 'formatted')
                end if
                if (mulorder >= 1) then
                    call openfile('m1_field_'//tcl//'.cube', lum1, 'new', 'formatted')
                end if
                if (mulorder >= 2) then
                    call openfile('m2_field_'//tcl//'.cube', lum2, 'new', 'formatted')
                end if
                if (mulorder >= 3) then
                    call openfile('m3_field_'//tcl//'.cube', lum3, 'new', 'formatted')
                end if
                if (mulorder >= 4) then
                    call openfile('m4_field_'//tcl//'.cube', lum4, 'new', 'formatted')
                end if
                if (mulorder >= 5) then
                    call openfile('m5_field_'//tcl//'.cube', lum5, 'new', 'formatted')
                end if
                do i = 1, mulorder + 1
                    if (i == 1) then
                        lu = lum0
                    else if (i == 2) then
                        lu = lum1
                    else if (i == 3) then
                        lu = lum2
                    else if (i == 4) then
                        lu = lum3
                    else if (i == 5) then
                        lu = lum4
                    else if (i == 6) then
                        lu = lum5
                    end if
                    write(lu,'(a)') 'Classical electric field component '//tcl
                    write(lu,'(a)') 'Generated by the PE library'
                    write(lu,'(i5,3f12.6)') nnucs, origin
                    write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
                    write(lu,'(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
                    write(lu,'(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
                    do j = 1, nnucs
                        write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
                    end do
                end do
                do i = 1, xsteps * ysteps
                    j = (i - 1) * zsteps + 1
                    k = j - 1 + zsteps
                    if (mulorder >= 0) then
                        write(lum0,'(6e13.5)') Fpe(j:k,l,0)
                    end if
                    if (mulorder >= 1) then
                        Fpe(j:k,l,0) = Fpe(j:k,l,0) + Fpe(j:k,l,1)
                        write(lum1,'(6e13.5)') Fpe(j:k,l,0)
                    end if
                    if (mulorder >= 2) then
                        Fpe(j:k,l,0) = Fpe(j:k,l,0) + Fpe(j:k,l,2)
                        write(lum2,'(6e13.5)') Fpe(j:k,l,0)
                    end if
                    if (mulorder >= 3) then
                        Fpe(j:k,l,0) = Fpe(j:k,l,0) + Fpe(j:k,l,3)
                        write(lum3,'(6e13.5)') Fpe(j:k,l,0)
                    end if
                    if (mulorder >= 4) then
                        Fpe(j:k,l,0) = Fpe(j:k,l,0) + Fpe(j:k,l,4)
                        write(lum4,'(6e13.5)') Fpe(j:k,l,0)
                    end if
                    if (mulorder >= 5) then
                        Fpe(j:k,l,0) = Fpe(j:k,l,0) + Fpe(j:k,l,5)
                        write(lum5,'(6e13.5)') Fpe(j:k,l,0)
                    end if
                end do
                if (mulorder >= 0) then
                    close(lum0)
                end if
                if (mulorder >= 1) then
                    close(lum1)
                end if
                if (mulorder >= 2) then
                    close(lum2)
                end if
                if (mulorder >= 3) then
                    close(lum3)
                end if
                if (mulorder >= 4) then
                    close(lum4)
                end if
                if (mulorder >= 5) then
                    close(lum5)
                end if
            end do
        end if
        deallocate(Fpe)
    end if

    if (lpol(1,1) .and. mep_mulcube) then
        allocate(Fmuls(3*npols,1))
        Fmuls = 0.0
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            if (pelib_mul_damp .or. pelib_amoeba_damp) then
                alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
            end if
            do j = 1, nsites
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
                    if (alpha_j <= zero) then
                        error stop 'damping requires polarizabilities on all sites with multipoles'
                    end if
                end if
                Rji = Rs(:,i) - Rs(:,j)
                if (lmul(0)) then
                    if (maxval(abs(M0s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M0s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M0s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M0s(:,j))
                        end if
                    end if
                end if
                if (lmul(1)) then
                    if (maxval(abs(M1s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M1s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M1s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M1s(:,j))
                        end if
                    end if
                end if
                if (lmul(2)) then
                    if (maxval(abs(M2s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M2s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M2s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M2s(:,j))
                        end if
                    end if
                end if
                if (lmul(3)) then
                    if (maxval(abs(M3s(:,j))) >= zero) then
                        if (pelib_mul_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M3s(:,j), alpha_i, alpha_j, mul_damp)
                        else if (pelib_amoeba_damp) then
                            call multipole_derivative_damped(Fmuls(l:l+2,1), Rji, M3s(:,j), alpha_i, alpha_j, amoeba_damp)
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M3s(:,j))
                        end if
                    end if
                end if
                if (lmul(4)) then
                    if (maxval(abs(M4s(:,j))) >= zero) then
                        if (pelib_mul_damp .or. pelib_amoeba_damp) then
                            error stop 'damping not implemented for fourth order multipoles'
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M4s(:,j))
                        end if
                    end if
                end if
                if (lmul(5)) then
                    if (maxval(abs(M5s(:,j))) >= zero) then
                        if (pelib_mul_damp .or. pelib_amoeba_damp) then
                            error stop 'damping not implemented for fifth order multipoles'
                        else
                            call multipole_derivative(Fmuls(l:l+2,1), Rji, M5s(:,j))
                        end if
                    end if
                end if
            end do
            l = l + 3
        end do
        allocate(M1inds(3*npols,1))
        if (mep_extfld) then
            j = 1
            do i = 1, npols
                Fmuls(j:j+2,1) = Fmuls(j:j+2,1) + extfld
                j = j + 3
            end do
        end if
        if (mep_lf) then
            Fmuls = 0.0
            do i = 1, npols
                Fmuls(lf_component+3*(i-1),1) = 1.0
            end do
        end if
        call induced_moments(M1inds, Fmuls)
        deallocate(Fmuls)
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(M1inds(1,1), 3 * npols, rmpi, master, comm, ierr)
        end if
#endif
        if (myid == master) then
            allocate(Find(npoints,3))
        else if (myid /= master) then
            allocate(Find(cubedists(myid),3))
        end if
        Find = 0.0
        k = 1
        do i = cube_start, cube_finish
            l = 1
            do j = 1, nsites
                if (zeroalphas(j)) cycle
                Rsg = Rg(:,i) - Rs(:,j)
                Fs = 0.0
                call multipole_derivative(Fs, Rsg, M1inds(l:l+2,1))
                Find(k,:) = Find(k,:) + Fs
                l = l + 3
            end do
            k = k + 1
        end do
        deallocate(M1inds)
#if defined(VAR_MPI)
        if (myid == master .and. nprocs > 1) then
            displs(0) = 0
            do i = 1, nprocs
                displs(i) = displs(i-1) + cubedists(i-1)
            end do
            do i = 1, 3
                call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Find(1,i), cubedists, displs, rmpi, master, comm, ierr)
            end do
        else if (myid /= master) then
            do i = 1, 3
                call mpi_gatherv(Find(1,i), cubedists(myid), rmpi, dummy_real, cubedists, displs, rmpi, master, comm, ierr)
            end do
        end if
#endif
        if (myid == master) then
            do l = 1, 3
                write(cl,*) l
                tcl = trim(adjustl(cl))
                call openfile('ind_field_'//tcl//'.cube', lu, 'new', 'formatted')
                write(lu,'(a)') 'Classical induced electric field component '//tcl
                write(lu,'(a)') 'Generated by the PE library'
                write(lu,'(i5,3f12.6)') nnucs, origin
                write(lu,'(i5,3f12.6)') xsteps, step(1), 0.0, 0.0
                write(lu,'(i5,3f12.6)') ysteps, 0.0, step(2), 0.0
                write(lu,'(i5,3f12.6)') zsteps, 0.0, 0.0, step(3)
                do j = 1, nnucs
                    write(lu,'(i5,4f12.6)') nint(Zm(1,j)), Zm(1,j), Rm(:,j)
                end do
                do i = 1, xsteps * ysteps
                    j = (i - 1) * zsteps + 1
                    k = j - 1 + zsteps
                    write(lu,'(6e13.5)') Find(j:k,l)
                end do
                close(lu)
            end do
        end if
        deallocate(Find)
    end if

end subroutine compute_mef_cube

!------------------------------------------------------------------------------

#ifdef VAR_PDE

subroutine compute_mep_input(denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces
    use pelib_blas_interfaces
    use pelib_multipole_tools, only: multipole_derivative

    real(rp), dimension(:), intent(in) :: denmats

    integer(ip) :: i, j, k
    integer(ip) :: num_coordinates
    integer(ip) :: quotient, remainder
    integer(ip) :: grid_start, grid_finish
    integer(ip), dimension(:), allocatable :: griddists, gridloops
    real(rp), dimension(1) :: Vtmp
    real(rp), dimension(3) :: Rmg
    real(rp), dimension(:), allocatable :: Vqm
    real(rp), dimension(:,:), allocatable :: Rgrid
    real(rp), dimension(:,:,:), allocatable :: Tk_ints

    !hdf5 IO variables
    integer(4) :: error
    integer(HSIZE_T), dimension(1) :: dim_1d
    integer(HSIZE_T), dimension(2) :: dim_2d
    integer(HID_T) :: file_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: space_id

#if defined(VAR_MPI)
    allocate(griddists(0:nprocs-1))
    allocate(gridloops(0:nprocs))
#endif

    if (myid == master) then
        ! read grid coordinates
        call h5open_f(error)
        call h5fopen_f(trim(h5gridfile), H5F_ACC_RDWR_F, file_id, error)
        call h5dopen_f(file_id, 'num_coordinates', dset_id, error)
        dim_1d = 1
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, num_coordinates, dim_1d, error)
        call h5dclose_f(dset_id, error)
        allocate(Rgrid(3,num_coordinates))
        call h5dopen_f(file_id, 'coordinates', dset_id, error)
        dim_2d = [3, num_coordinates]
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, Rgrid, dim_2d, error)
        call h5dclose_f(dset_id, error)

#if defined(VAR_MPI)
        quotient = num_coordinates / nprocs
        griddists = quotient
        if (nprocs * quotient < num_coordinates) then
            remainder = num_coordinates - nprocs * quotient
            do i = 1, remainder
                griddists(i-1) = griddists(i-1) + 1
            end do
        end if
        gridloops(0) = 0
        do i = 1, nprocs
            gridloops(i) = sum(griddists(0:i-1))
        end do
#endif
    end if

#if defined(VAR_MPI)
    call mpi_bcast(num_coordinates, 1, impi, master, comm, ierr)
    call mpi_bcast(griddists(0), nprocs, impi, master, comm, ierr)
    call mpi_bcast(gridloops(0), nprocs + 1, impi, master, comm, ierr)

    grid_start = gridloops(myid) + 1
    grid_finish = gridloops(myid+1)
#else
    grid_start = 1
    grid_finish = num_coordinates
#endif

#if defined(VAR_MPI)
    if (myid /= master) allocate(Rgrid(3,num_coordinates))
    call mpi_bcast(Rgrid(1,1), 3 * num_coordinates, rmpi, master, comm, ierr)
#endif

    ! compute potential on grid coordinates
    if (myid == master) then
        allocate(Vqm(num_coordinates))
    else if (myid /= master) then
        allocate(Vqm(griddists(myid)))
    end if
    Vqm = 0.0
    allocate(Tk_ints(nnbas,1,1))
    k = 1
    do i = grid_start, grid_finish
        call Tk_integrals('potential_derivative', Rgrid(:,i), Tk_ints)
        Vqm(k) = dot(denmats, Tk_ints(:,1,1))
        Vtmp = 0.0
        do j = 1, nnucs
            Rmg = Rgrid(:,i) - Rm(:,j)
            call multipole_derivative(Vtmp, Rmg, Zm(:,j))
        end do
        Vqm(k) = Vqm(k) + Vtmp(1)
        k = k + 1
    end do
    deallocate(Tk_ints)

#if defined(VAR_MPI)
    if (myid == master .and. nprocs > 1) then
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + griddists(i-1)
        end do
        call mpi_gatherv(mpi_in_place, dummy_int, rmpi, Vqm(1), griddists, displs, rmpi, master, comm, ierr)
    else if (myid /= master) then
        call mpi_gatherv(Vqm(1), griddists(myid), rmpi, dummy_real, griddists, displs, rmpi, master, comm, ierr)
    end if
#endif

    if (myid == master) then
        ! write potential to hdf5 file
        dim_1d = num_coordinates
        call h5screate_simple_f(1_4, dim_1d, space_id, error)
        call h5dcreate_f(file_id, 'potential', H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Vqm, dim_1d, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)
        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end if

    deallocate(Rgrid)
    deallocate(Vqm)

#if defined(VAR_MPI)
    deallocate(gridloops)
    deallocate(griddists)
#endif

end subroutine compute_mep_input

#else

subroutine compute_mep_input(denmats)

    use pelib_mpi
    use pelib_options
    use pelib_integral_interfaces
    use pelib_blas_interfaces
    use pelib_multipole_tools, only: multipole_derivative

    real(rp), dimension(:), intent(in) :: denmats

    error stop 'Not compiled with PDE support (recompile with -DENABLE_PDE=ON)'

end subroutine compute_mep_input

#endif

!------------------------------------------------------------------------------

end module pelib_analysis_tools
