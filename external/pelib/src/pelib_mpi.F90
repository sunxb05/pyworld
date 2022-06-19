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
module pelib_mpi

    use pelib_precision

#if defined(VAR_MPI)
#if defined(USE_MPI_MOD_F90)
    use mpi
#else
    include 'mpif.h'
#endif
    integer(ip), parameter :: comm = MPI_COMM_WORLD
    integer(ip), parameter :: impi = MPI_INTEGER
    integer(ip), parameter :: rmpi = MPI_REAL8
    integer(ip), parameter :: lmpi = MPI_LOGICAL
    integer(ip), parameter :: cmpi = MPI_CHARACTER
#else
    integer(ip), parameter :: comm = 666
#endif

    integer(ip), save :: master, myid, nprocs, ierr

contains

#if defined(VAR_MPI)
subroutine mpi_sync()

    use pelib_options
    use pelib_options
    use pelib_integral_interfaces

    implicit none

    integer(ip) :: i, j
    integer(ip) :: quotient, remainder

    if (.not. allocated(displs)) allocate(displs(0:nprocs))
    if (.not. allocated(siteloops)) allocate(siteloops(0:nprocs))
    if (.not. allocated(sitedists)) allocate(sitedists(0:nprocs-1))
    if (myid == master) then
        quotient = nsites / nprocs
        sitedists = quotient
        if (nprocs * quotient < nsites) then
            remainder = nsites - nprocs * quotient
            do i = 1, remainder
                sitedists(i-1) = sitedists(i-1) + 1
            end do
        end if
        siteloops(0) = 0
        do i = 1, nprocs
            siteloops(i) = sum(sitedists(0:i-1))
        end do
        displs(0) = 0
        do i = 1, nprocs
            displs(i) = displs(i-1) + 3 * sitedists(i-1)
        end do
    end if

    call mpi_bcast(nsites, 1, impi, master, comm, ierr)
    call mpi_bcast(siteloops(0), nprocs + 1, impi, master, comm, ierr)
    call mpi_bcast(sitedists(0), nprocs, impi, master, comm, ierr)

    site_start = siteloops(myid) + 1
    site_finish = siteloops(myid+1)

    call mpi_bcast(nnucs, 1, impi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(Zm)) allocate(Zm(1,nnucs))
    call mpi_bcast(Zm(1,1), nnucs, rmpi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(Rm)) allocate(Rm(3,nnucs))
    call mpi_bcast(Rm(1,1), 3 * nnucs, rmpi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(core_com)) allocate(core_com(3))
    call mpi_bcast(core_com(1), 3, rmpi, master, comm, ierr)
    if (myid /= master) call init_pelib_integral_interfaces(core_com, Rm)

    if (myid /= master .and. .not. allocated(Rs)) allocate(Rs(3,nsites))
    call mpi_bcast(Rs(1,1), 3 * nsites, rmpi, master, comm, ierr)

    if (myid /= master .and. .not. allocated(Zs)) allocate(Zs(1,nsites))
    call mpi_bcast(Zs(1,1), nsites, rmpi, master, comm, ierr)

    call mpi_bcast(pelib_polar, 1, lmpi, master, comm, ierr)
    call mpi_bcast(lpol(0,0), 9, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_sol, 1, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_diis, 1, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_nomuq, 1, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_novmu, 1, lmpi, master, comm, ierr)

    if (pelib_polar) then
        if (any(lpol)) then
            if (.not. allocated(poldists)) allocate(poldists(0:nprocs-1))
            if (myid == master) then
                poldists = 0
                do i = 1, nprocs
                    do j = siteloops(i-1) + 1, siteloops(i)
                        if (zeroalphas(j)) then
                            continue
                        else
                            poldists(i-1) = poldists(i-1) + 1
                        end if
                    end do
                end do
            end if
            poldists = 3 * poldists
            call mpi_bcast(poldists(0), nprocs, impi, master, comm, ierr)
            call mpi_bcast(npols, 1, impi, master, comm, ierr)
            call mpi_bcast(mul_lexlst, 1, impi, master, comm, ierr)
            if (myid /= master .and. .not. allocated(mul_exclists)) then
                allocate(mul_exclists(mul_lexlst,nsites))
            end if
            call mpi_bcast(mul_exclists(1,1), mul_lexlst * nsites, impi, master, comm, ierr)
            call mpi_bcast(pol_lexlst, 1, impi, master, comm, ierr)
            if (myid /= master .and. .not. allocated(pol_exclists)) then
                allocate(pol_exclists(pol_lexlst,nsites))
            end if
            call mpi_bcast(pol_exclists(1,1), pol_lexlst * nsites, impi, master, comm, ierr)
            if (myid /= master .and. .not. allocated(zeroalphas)) then
                allocate(zeroalphas(nsites))
            end if
            call mpi_bcast(zeroalphas(1), nsites, lmpi, master, comm, ierr)
        end if
        if (pelib_sol) then
            allocate(surploops(0:nprocs))
            allocate(surpdists(0:nprocs-1))
            if (myid == master) then
                quotient = nsurp / nprocs
                surpdists = quotient
                if (nprocs * quotient < nsurp) then
                    remainder = nsurp - nprocs * quotient
                    do i = 1, remainder
                        surpdists(i-1) = surpdists(i-1) + 1
                    end do
                end if
                surploops(0) = 0
                do i = 1, nprocs
                    surploops(i) = sum(surpdists(0:i-1))
                end do
            end if
            call mpi_bcast(nsurp, 1, impi, master, comm, ierr)
            call mpi_bcast(surpdists(0), nprocs, impi, master, comm, ierr)
            call mpi_bcast(surploops(0), nprocs + 1, impi, master, comm, ierr)
            surp_start = surploops(myid) + 1
            surp_finish = surploops(myid+1)
            if (myid /= master) allocate(Rsp(3,nsurp))
            call mpi_bcast(Rsp(1,1), 3 * nsurp, rmpi, master, comm, ierr)
            if (myid /= master) allocate(Sa(nsurp))
            call mpi_bcast(Sa(1), nsurp, rmpi, master, comm, ierr)
            call mpi_bcast(eps, 1, rmpi, master, comm, ierr)
            call mpi_bcast(epsinf, 1, rmpi, master, comm, ierr)
        end if
        call mpi_bcast(pelib_iter, 1, lmpi, master, comm, ierr)
        if (pelib_iter .and. any(lpol)) then
            if (myid /= master .and. .not. allocated(P11s)) allocate(P11s(6,nsites))
            call mpi_bcast(P11s(1,1), 6 * nsites, rmpi, master, comm, ierr)
        end if
        call mpi_bcast(pelib_nomb, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pelib_ind_damp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pelib_mul_damp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pelib_core_damp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pelib_amoeba_damp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(ind_damp, 1, rmpi, master, comm, ierr)
        call mpi_bcast(mul_damp, 1, rmpi, master, comm, ierr)
        call mpi_bcast(core_damp, 1, rmpi, master, comm, ierr)
        call mpi_bcast(amoeba_damp, 1, rmpi, master, comm, ierr)
        if (pelib_core_damp) then
            if (myid /= master .and. .not. allocated(core_alphas)) allocate(core_alphas(6,nnucs))
            call mpi_bcast(core_alphas, 6 * nnucs, rmpi, master, comm, ierr)
        endif
    end if

    call mpi_bcast(pelib_savden, 1, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_twoints, 1, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_fd, 1, lmpi, master, comm, ierr)
    if (pelib_fd) then
        call mpi_bcast(pelib_fdes, 1, lmpi, master, comm, ierr)
        call mpi_bcast(pelib_repuls, 1, lmpi, master, comm, ierr)
    end if

    call mpi_bcast(pelib_restart, 1, lmpi, master, comm, ierr)

    call mpi_bcast(pelib_lf, 1, lmpi, master, comm, ierr)
    if (pelib_lf) then
        call mpi_bcast(ncrds, 1, impi, master, comm, ierr)
        if (myid /= master .and. (ncrds > 0) .and. .not. allocated(crds)) allocate(crds(3,ncrds))
        call mpi_bcast(crds(1,1), 3 * ncrds, rmpi, master, comm, ierr)
    end if

    call mpi_bcast(lmul(0), 6, lmpi, master, comm, ierr)

    if (lmul(0)) then
        if (myid /= master .and. .not. allocated(M0s)) allocate(M0s(1,nsites))
        call mpi_bcast(M0s(1,1), nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(1)) then
        if (myid /= master .and. .not. allocated(M1s)) allocate(M1s(3,nsites))
        call mpi_bcast(M1s(1,1), 3 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(2)) then
        if (myid /= master .and. .not. allocated(M2s)) allocate(M2s(6,nsites))
        call mpi_bcast(M2s(1,1), 6 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(3)) then
        if (myid /= master .and. .not. allocated(M3s)) allocate(M3s(10,nsites))
        call mpi_bcast(M3s(1,1), 10 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(4)) then
        if (myid /= master .and. .not. allocated(M4s)) allocate(M4s(15,nsites))
        call mpi_bcast(M4s(1,1), 15 * nsites, rmpi, master, comm, ierr)
    end if
    if (lmul(5)) then
        if (myid /= master .and. .not. allocated(M5s)) allocate(M5s(21,nsites))
        call mpi_bcast(M5s(1,1), 21 * nsites, rmpi, master, comm, ierr)
    end if

    call mpi_bcast(lvdw, 1, lmpi, master, comm, ierr)
    if (lvdw) then
        if (myid /= master .and. .not. allocated(LJs)) allocate(LJs(2,nsites))
        call mpi_bcast(LJs(1,1), 2 * nsites, rmpi, master, comm, ierr)
        call mpi_bcast(qmLJsites, 1, impi, master, comm, ierr)
        if (myid /= master .and. .not. allocated(qmLJs)) allocate(qmLJs(2,qmLJsites))
        call mpi_bcast(qmLJs(1,1), 2 * qmLJsites, rmpi, master, comm, ierr)
    end if

    call mpi_bcast(pelib_mep, 1, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_cube, 1, lmpi, master, comm, ierr)

    if (pelib_mep .or. pelib_cube) then
        call mpi_bcast(cube_field, 1, lmpi, master, comm, ierr)
        call mpi_bcast(mep_field, 1, lmpi, master, comm, ierr)
        call mpi_bcast(mep_extfld, 1, lmpi, master, comm, ierr)
        call mpi_bcast(extfld(1), 3, rmpi, master, comm, ierr)
        call mpi_bcast(mep_qmcube, 1, lmpi, master, comm, ierr)
        call mpi_bcast(mep_mulcube, 1, lmpi, master, comm, ierr)
        call mpi_bcast(mep_input, 1, lmpi, master, comm, ierr)
        if (.not. allocated(cubeloops)) allocate(cubeloops(0:nprocs))
        if (.not. allocated(cubedists)) allocate(cubedists(0:nprocs-1))
        if (myid == master) then
            quotient = npoints / nprocs
            cubedists = quotient
            if (nprocs * quotient < npoints) then
                remainder = npoints - nprocs * quotient
                do i = 1, remainder
                    cubedists(i-1) = cubedists(i-1) + 1
                end do
            end if
            cubeloops(0) = 0
            do i = 1, nprocs
                cubeloops(i) = sum(cubedists(0:i-1))
            end do
        end if
        call mpi_bcast(npoints, 1, impi, master, comm, ierr)
        call mpi_bcast(cubedists(0), nprocs, impi, master, comm, ierr)
        call mpi_bcast(cubeloops(0), nprocs + 1, impi, master, comm, ierr)
        cube_start = cubeloops(myid) + 1
        cube_finish = cubeloops(myid+1)
        if (myid /= master .and. .not. allocated(Rg)) allocate(Rg(3,npoints))
        call mpi_bcast(Rg(1,1), 3 * npoints, rmpi, master, comm, ierr)
    end if

    call mpi_bcast(pelib_old_field, 1, lmpi, master, comm, ierr)
    call mpi_bcast(pelib_fmm, 1, lmpi, master, comm, ierr)
    call mpi_bcast(ESPF_bond_midpoints, 1, lmpi, master, comm, ierr)
    call mpi_bcast(incore_nearfield, 1, lmpi, master, comm, ierr)
    call mpi_bcast(fmm_expansion_order, 1, impi, master, comm, ierr)
    call mpi_bcast(fmm_ncrit, 1, impi, master, comm, ierr)
    call mpi_bcast(fmm_theta, 1, rmpi, master, comm, ierr)
    call mpi_bcast(qmmm_interaction_scheme, len(qmmm_interaction_scheme), cmpi, master, comm, ierr)
    call mpi_bcast(ESPF_solver, len(ESPF_solver), cmpi, master, comm, ierr)
    call mpi_bcast(ESPF_grid_type, len(ESPF_grid_type), cmpi, master, comm, ierr)
    call mpi_bcast(qmmm_interaction_rcut, 1, rmpi, master, comm, ierr)
    call mpi_bcast(qmmm_interaction_rfit, 1, rmpi, master, comm, ierr)
    call mpi_bcast(qmmm_interaction_theta, 1, rmpi, master, comm, ierr)
    call mpi_bcast(svdtol, 1, rmpi, master, comm, ierr)
    call mpi_bcast(qmmm_interaction_order, 1, impi, master, comm, ierr)

    synced = .true.

end subroutine mpi_sync
#endif

end module pelib_mpi
