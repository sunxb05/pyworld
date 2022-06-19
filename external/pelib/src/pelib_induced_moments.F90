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
module pelib_induced_moments

    use pelib_precision

    implicit none

    private

    public :: induced_moments, response_matrix, d_local, induced_moments_eef

contains

!------------------------------------------------------------------------------

subroutine induced_moments_eef(Mkinds, Fks, Qinds, Vs)

    use pelib_mpi
    use pelib_options

    real(rp), dimension(:,:), intent(out), optional :: Mkinds
    real(rp), dimension(:,:), intent(in), optional :: Fks
    real(rp), dimension(:,:), intent(out), optional :: Qinds
    real(rp), dimension(:,:), intent(in), optional :: Vs

    if (pelib_sol .and. .not. any(lpol)) then
        if (.not. present(Qinds) .and. .not. present(Vs)) then
            error stop 'missing in-/output pelib_options'
        end if
    else if (pelib_sol .and. any(lpol)) then
        if (.not. present(Qinds) .and. .not. present(Vs) .and.&
           & .not. present(Mkinds) .and. .not. present(Fks)) then
            error stop 'missing in-/output pelib_options'
        end if
    else if (.not. pelib_sol .and. any(lpol)) then
        if (.not. present(Mkinds) .and. .not. present(Fks)) then
            error stop 'missing in-/output pelib_options'
        end if
    else
        error stop 'pelib_polar should be .false. here'
    end if

    if (myid == master) then
        if (pelib_sol .and. .not. any(lpol)) then
            call direct_solver_eef(Qinds=Qinds, Vs=Vs)
        else if (pelib_sol .and. any(lpol)) then
            call direct_solver_eef(Mkinds=Mkinds, Fks=Fks, Qinds=Qinds, Vs=Vs)
        else if (.not. pelib_sol .and. any(lpol)) then
            call direct_solver_eef(Mkinds=Mkinds, Fks=Fks)
        end if
    end if

end subroutine induced_moments_eef
!------------------------------------------------------------------------------

subroutine induced_moments(Mkinds, Fks, Qinds, Vs)

    use pelib_mpi
    use pelib_options

    real(rp), dimension(:,:), intent(out), optional :: Mkinds
    real(rp), dimension(:,:), intent(in), optional :: Fks
    real(rp), dimension(:,:), intent(out), optional :: Qinds
    real(rp), dimension(:,:), intent(in), optional :: Vs

    if (pelib_sol .and. .not. any(lpol)) then
        if (.not. present(Qinds) .and. .not. present(Vs)) then
            error stop 'missing in-/output pelib_options'
        end if
    else if (pelib_sol .and. any(lpol)) then
        if (.not. present(Qinds) .and. .not. present(Vs) .and.&
           & .not. present(Mkinds) .and. .not. present(Fks)) then
            error stop 'missing in-/output pelib_options'
        end if
    else if (.not. pelib_sol .and. any(lpol)) then
        if (.not. present(Mkinds) .and. .not. present(Fks)) then
            error stop 'missing in-/output pelib_options'
        end if
    else
        error stop 'pelib_polar should be .false. here'
    end if

    if (pelib_iter) then
        if (pelib_diis) then
            if (pelib_sol .and. .not. any(lpol)) then
                call diis_solver(Qinds=Qinds, Vs=Vs)
            else if (pelib_sol .and. any(lpol)) then
                call diis_solver(Qinds=Qinds, Vs=Vs, Mkinds=Mkinds, Fks=Fks)
            else if (.not. pelib_sol .and. any(lpol)) then
                error stop 'DIIS solver can only be used with FIXSOL'
            end if
        else
            if (pelib_sol) then
                error stop 'standard solver is not compatible with fixsol'
            else if (any(lpol)) then
                call iterative_solver(Mkinds, Fks)
            end if
        end if
    else
        if (myid == master) then
            if (pelib_sol .and. .not. any(lpol)) then
                call direct_solver(Qinds=Qinds, Vs=Vs)
            else if (pelib_sol .and. any(lpol)) then
                call direct_solver(Mkinds=Mkinds, Fks=Fks, Qinds=Qinds, Vs=Vs)
            else if (.not. pelib_sol .and. any(lpol)) then
                call direct_solver(Mkinds=Mkinds, Fks=Fks)
            end if
        end if
    end if

end subroutine induced_moments

!------------------------------------------------------------------------------

subroutine iterative_solver(Mkinds, Fks)

    use pelib_mpi
    use pelib_options
    use pelib_utils, only: openfile
    use pelib_blas_interfaces, only: nrm2, spmv, dot
    use pelib_multipole_tools, only: Tk_tensor, Tk_damped_tensor
    use fmm, only: field_direct, field_fmm

    real(rp), dimension(:,:), intent(out) :: Mkinds
    real(rp), dimension(:,:), intent(in) :: Fks

    integer(ip) :: lu, iter
    integer(ip) :: i, j, k, l, m, n
    integer(ip), dimension(:), allocatable :: displacements
    logical(lp) :: exclude, lexist, restart, converged
    real(rp) :: alpha_i, alpha_j, norm, redthr
    real(rp), parameter :: i3 = 1.0 / 3.0
    real(rp), dimension(:), allocatable :: T2, Rij, Ftmp, M1tmp

    ! DIIS related quantities
    integer(ip) :: DIIS_size, DIIS_info
    integer(ip), parameter :: DIIS_max_size = 25
    integer(ip), dimension(1) :: DIIS_iter
    real(rp), dimension(:,:), allocatable :: DIIS_previous_induced_moments, DIIS_residuals
    real(rp), dimension(:,:), allocatable :: DIIS_Bmat
    real(rp), dimension(:), allocatable :: DIIS_weights
    integer(ip), dimension(:), allocatable  :: DIIS_ipiv

    ! FMM related quantities
    real(rp), dimension(:,:), allocatable :: fmm_coordinates, fmm_multipoles, fmm_field
    real(rp), dimension(:), allocatable :: fmm_damping_factors
    integer(ip), dimension(:,:), allocatable :: fmm_exclusions
    integer(ip), dimension(:), allocatable :: fmm_renumber
    character(len=6) :: fmm_damp_type = '      '

    if (myid == master) then
        if (pelib_redthr .and. fock .and. .not. pelib_restart) then
            redthr = 10.0**(- nint(log10(thriter)) - 2 * scf_cycle - redlvl)
            if (redthr * thriter > thriter) then
                write(luout,'(4x,a,es8.1)') 'INFO: using reduced threshold to&
                                            & determine induced dipole moments: ', redthr * thriter
            else
                redthr = 1.0
            end if
        else
            redthr = 1.0
        end if
    end if

    if (myid == master) then
        inquire(file='pelib_induced_moments.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist .and. (fock .or. energy) .and. (scf_cycle > 1 .or. pelib_restart) .and. .not. pelib_nomb) then
        restart = .true.
    else
        restart = .false.
    end if

    if (restart) then
        if (myid == master) then
            call openfile('pelib_induced_moments.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) Mkinds
            close(lu)
        else
            Mkinds = 0.0
        end if
    else
        Mkinds = 0.0
    end if

    allocate(displacements(0:nprocs))
    if (nprocs > 1) then
        displacements(0) = 0
#if defined(VAR_MPI)
        do i = 1, nprocs
            displacements(i) = displacements(i-1) + poldists(i-1)
        end do
#endif
    else
        displacements = 0
    end if

    allocate(T2(6), Rij(3), Ftmp(3*npols), M1tmp(3*npols))

    ! Allocate memory for the DIIS-related quantities
    allocate(DIIS_previous_induced_moments(3*npols,DIIS_max_size))
    allocate(DIIS_residuals(3*npols,DIIS_max_size))
    allocate(DIIS_Bmat(DIIS_max_size,DIIS_max_size))
    allocate(DIIS_weights(DIIS_max_size))
    allocate(DIIS_ipiv(DIIS_max_size))
    
    if (.not. pelib_old_field) then
        allocate(fmm_coordinates(npols, 3), source=0.0_rp)
        allocate(fmm_multipoles(npols, 4), source=0.0_rp)
        allocate(fmm_exclusions(npols, pol_lexlst), source=0_ip)
        allocate(fmm_field(npols, 4))
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            fmm_coordinates(l, :) = Rs(:, i)
            l = l + 1
        end do
        ! fmm_renumber the exclusion list from (1..nsites) to (1..npols)
        allocate(fmm_renumber(nsites), source=0_ip)
        l = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            fmm_renumber(i) = l
            l = l + 1
        end do
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            do j = 1, pol_lexlst
                if (pol_exclists(j, i) == 0) cycle
                fmm_exclusions(fmm_renumber(i), j) = fmm_renumber(pol_exclists(j, i))
            end do
        end do

        l = 1
        if (pelib_ind_damp .or. pelib_amoeba_damp) then
            allocate(fmm_damping_factors(npols), source=0.0_rp)
            if (pelib_ind_damp) then
                do i = 1, nsites
                    if (zeroalphas(i)) cycle
                    alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
                    fmm_damping_factors(l) = sqrt(ind_damp) / alpha_i**(1.0 / 6.0)
                    l = l + 1
                end do
            else if (pelib_amoeba_damp) then
                do i = 1, nsites
                    if (zeroalphas(i)) cycle
                    alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
                    fmm_damping_factors(l) = amoeba_damp**(1.0 / 6.0)  / alpha_i**(1.0 / 6.0) 
                    l = l + 1
                end do
            end if
        end if
    end if

    do n = 1, size(Mkinds, 2)
#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(Fks(1,n), 3*npols, rmpi, master, comm, ierr)
        end if
#endif
        if (.not. restart) then
            l = 1 + displacements(myid)
            do i = site_start, site_finish
                if (zeroalphas(i)) cycle
                call spmv(P11s(:,i), Fks(l:l+2,n), Mkinds(l:l+2,n), 'L')
                l = l + 3
            end do
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                call mpi_reduce(mpi_in_place, Mkinds(1,n), 3*npols, rmpi, mpi_sum, master, comm, ierr)
            else if (myid /= master) then
                call mpi_reduce(Mkinds(1,n), dummy_real, 3*npols, rmpi, mpi_sum, master, comm, ierr)
            end if
#endif
        end if

#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(Mkinds(1,n), 3*npols, rmpi, master, comm, ierr)
        end if
#endif

        if (pelib_nomb) cycle

        DIIS_previous_induced_moments = 0.0
        DIIS_residuals = 0.0
        DIIS_iter = 1

        iter = 1
        do
            ! Construct DIIS guess for induced dipoles
            if (iter > 2) then
                DIIS_size = min(iter, DIIS_max_size)
                DIIS_Bmat(1:DIIS_size,1:DIIS_size) = 0.0
                DIIS_Bmat(1:DIIS_size,DIIS_size) = -1.0
                DIIS_Bmat(DIIS_size,1:DIIS_size) = -1.0
                DIIS_Bmat(DIIS_size,DIIS_size) = 0.0
                do i = 1, DIIS_size-1
                    do j = 1, DIIS_size-1
                        DIIS_Bmat(j,i) = dot(DIIS_residuals(:,i), DIIS_residuals(:,j))
                    end do
                end do
                ! Allocate and fill right hand side
                DIIS_weights(1:DIIS_size-1) = 0.0
                DIIS_weights(DIIS_size) = -1.0
                ! Solve B x mu = rhs
                call dgesv(DIIS_size, 1, DIIS_Bmat(1:DIIS_size,1:DIIS_size), DIIS_size, &
                           DIIS_ipiv(1:DIIS_size), DIIS_weights(1:DIIS_size), DIIS_size, DIIS_info)
                if (DIIS_info /= 0) then
                    ! something went wrong so we continue without extrapolation
                    write(luout,*) 'INFO: something went wrong when solving DIIS equations for induced dipole solver'
                    write(luout,*) 'INFO: continuing with regular Jacobi iteration'
                else
                    ! Create DIIS guess for induced dipoles
                    Mkinds(:,n) = 0.0
                    do i = 1, DIIS_size-1
                        Mkinds(:,n) = Mkinds(:,n) + DIIS_weights(i) * DIIS_previous_induced_moments(:,i)
                    end do
                end if
            end if

            ! Compute field from induced dipoles on all sites
            Ftmp = 0.0
            if (.not. pelib_old_field) then
                fmm_field = 0.0
                do i = 0, npols - 1
                    fmm_multipoles(i+1, 1) = 0.0
                    fmm_multipoles(i+1, 2) = Mkinds(3*i+1, n)
                    fmm_multipoles(i+1, 3) = Mkinds(3*i+2, n)
                    fmm_multipoles(i+1, 4) = Mkinds(3*i+3, n)
                end do

                if (pelib_ind_damp .or. pelib_amoeba_damp) then
                    if (pelib_ind_damp) fmm_damp_type = 'THOLE '
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
                do i = 0, npols - 1
                    Ftmp(3*i+1) = fmm_field(i+1, 2)
                    Ftmp(3*i+2) = fmm_field(i+1, 3)
                    Ftmp(3*i+3) = fmm_field(i+1, 4)
                end do
#ifdef VAR_MPI
                call mpi_bcast(Ftmp(1), 3 * npols, rmpi, master, comm, ierr)
#endif
            else
                l = 1 + displacements(myid)
                do i = site_start, site_finish
                    if (zeroalphas(i)) cycle
                    if (pelib_ind_damp .or. pelib_amoeba_damp) then
                        alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
                    end if
                    m = 1
                    do j = 1, nsites
                        if (zeroalphas(j)) cycle
                        exclude = .false.
                        do k = 1, pol_lexlst
                            if (pol_exclists(k,i) == pol_exclists(1,j)) then
                                exclude = .true.
                                exit
                            end if
                        end do
                        if (i == j .or. exclude) then
                            m = m + 3
                            cycle
                        end if
                        Rij = Rs(:,j) - Rs(:,i)
                        if (pelib_amoeba_damp) then
                            alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                            call Tk_damped_tensor(T2, Rij, alpha_i, alpha_j, amoeba_damp)
                        else if (pelib_ind_damp) then
                            alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                            call Tk_damped_tensor(T2, Rij, alpha_i, alpha_j, ind_damp)
                        else
                            call Tk_tensor(T2, Rij)
                        end if
                        call spmv(T2, Mkinds(m:m+2,n), Ftmp(l:l+2), 'L', 1.0_rp, 1.0_rp)
                        m = m + 3
                    end do
                    l = l + 3
                end do
#if defined(VAR_MPI)
                if (nprocs > 1) then
                    call mpi_allreduce(mpi_in_place, Ftmp(1), 3*npols, rmpi, mpi_sum, comm, ierr)
                end if
#endif
            end if
            Ftmp = Ftmp + Fks(:,n)

            ! compute induced dipoles
            M1tmp = 0.0
            l = 1 + displacements(myid)
            do i = site_start, site_finish
                if (zeroalphas(i)) cycle
                call spmv(P11s(:,i), Ftmp(l:l+2), M1tmp(l:l+2), 'L')
                l = l + 3
            end do
#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_allreduce(mpi_in_place, M1tmp(1), 3*npols, rmpi, mpi_sum, comm, ierr)
            end if
#endif

            if (iter > DIIS_max_size) then
                DIIS_iter = minloc(abs(DIIS_weights(1:DIIS_size-1)))
            else
                DIIS_iter = iter
            end if
            DIIS_previous_induced_moments(:,DIIS_iter(1)) = M1tmp
            DIIS_residuals(:,DIIS_iter(1)) = M1tmp - Mkinds(:,n)

            norm = nrm2(M1tmp - Mkinds(:,n))
            Mkinds(:,n) = M1tmp

            if (myid == master) then
                if (norm < redthr * thriter) then
                    if (pelib_verbose .and. .not. energy) then
                        write(luout, '(6x,a,i2,a)') 'Induced dipole moments converged in ', iter, ' iterations.'
                    end if
                    converged = .true.
                else if (iter > 50) then
                    write(luout, *) 'could not converge induced dipole moments'
                    write(luout, *) 'Current norm: ', norm
                    write(luout, *) 'Target: ', redthr * thriter
                    error stop 'could not converge induced dipole moments'
                else
                    converged = .false.
                end if
            end if

#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_bcast(converged, 1, lmpi, master, comm, ierr)
            end if
#endif
            if (converged) then
                exit
            else
                iter = iter + 1
            end if
        end do
    end do

    if (fock) then
        if (myid == master) then
            call openfile('pelib_induced_moments.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Mkinds
            close(lu)
        end if
    end if

    deallocate(T2, Rij, Ftmp, M1tmp)
    deallocate(DIIS_previous_induced_moments)
    deallocate(DIIS_residuals)
    deallocate(DIIS_Bmat)
    deallocate(DIIS_weights)
    deallocate(DIIS_ipiv)

end subroutine iterative_solver

!------------------------------------------------------------------------------

subroutine direct_solver_eef(Mkinds, Fks, Qinds, Vs)

    use pelib_options
    use pelib_constants
    use pelib_blas_interfaces, only: spmv
    use pelib_multipole_tools, only: Tk_tensor


    real(rp), dimension(:,:), intent(out), optional :: Mkinds
    real(rp), dimension(:,:), intent(in), optional :: Fks
    real(rp), dimension(:,:), intent(out), optional :: Qinds
    real(rp), dimension(:,:), intent(in), optional :: Vs

    real(rp), dimension(:,:), allocatable :: D_full, Dinv, D
    real(rp), dimension(:,:), allocatable :: Minds_int, Ftot
    real(rp), dimension(:), allocatable :: work, B, Btmp, Rij, T2, nrmtmp, Ftmp
    real(rp), parameter :: fourpi = 4.0 * pi

    integer(ip), dimension(:,:), allocatable :: ipiv
    integer(ip) :: i, l, m, j, k, info, atom, nelems, id, ioff
    real(rp) :: eps_fac, dqsum

    if (pelib_noneq) then
        eps_fac = (epsinf - 1.0) / epsinf
    else
        eps_fac = (eps - 1.0) / eps
    end if

    if (present(Mkinds) .and. present(Qinds)) then
        nelems = 3 * npols + nsurp
    else if (present(Mkinds) .and. .not. present(Qinds)) then
        nelems = 3 * npols
    else if (.not. present(Mkinds) .and. present(Qinds)) then
        nelems = nsurp
    end if

    allocate(D_full(nelems,nelems))
    D_full = 0.0

    if (any(lpol)) then
        allocate(B(nelems*(nelems+1)/2))
        allocate(Btmp(nelems*(nelems+1)/2))

        call response_matrix(B)
        Btmp = 0.0
        l = 0
        j = 1
        do i = 1, nelems
            ioff = i + l
            l = ioff
            do m = i, nelems
                Btmp(ioff) = B(j)
                ioff = ioff + m
                j = j + 1
            end do
        end do

        call dsptsi(nelems,Btmp,D_full)
        deallocate(B)
        deallocate(Btmp)
        if (pelib_debug) then
            write(luout,*) 'D_full 1'
            call output(D_full,1,nelems,1,nelems,nelems,nelems,1,luout)
            id = 1
            do i = 1, nelems
                do j = 1, nelems
                    write(luout,*) 'D_full element id', D_full(i,j),id
                    id = id + 1
                end do
            end do
        end if


        allocate(D(nsurp,nsurp))
        call d_local(D)

        D_full(3*npols+1:,3*npols+1:) = D(:,:)
        deallocate(D)
        if (pelib_debug) then
            write(luout,*) 'D_full 2'
            call output(D_full,1,nelems,1,nelems,nelems,nelems,1,luout)
        end if
!   Replace T(1)*mu with T2*mu
!        D_full(3*npols+1:,1:3*npols) = 0.0
        if (pelib_debug) then
            write(luout,*) 'D_full 3'
            call output(D_full,1,nelems,1,nelems,nelems,nelems,1,luout)
        end if
        allocate(D(nsurp,3*npols))
        allocate(T2(6), Rij(3),nrmtmp(3),Ftmp(3))
        do j = 1, nsurp
            do i = 1, nsites
                Rij = Rsp(:,j) - Rs(:,i)
                call Tk_tensor(T2, Rij)
                do k = 1, 3
                   atom = idatm(j)
                   nrmtmp(k) = (all_coords(k,atom) - Rsp(k,j))/Radsp(atom)
                end do
                call spmv(T2, nrmtmp, Ftmp, 'L')
                D_full(3*npols+j,i) =  - Ftmp(1)
                D_full(3*npols+j,i+1) = - Ftmp(2)
                D_full(3*npols+j,i+2) = - Ftmp(3)
            end do
        end do
!        D_full(3*npols+1:,1:3*npols) = - D(:,:)
        if (pelib_debug) then
            write(luout,*) 'D_full 4'
            call output(D_full,1,nelems,1,nelems,nelems,nelems,1,luout)
        end if
        deallocate(T2, Rij,nrmtmp,Ftmp,D)
    else
        call d_local(D_full)
    end if

    allocate(ipiv(1,nelems))
    allocate(Dinv(nelems,nelems))

    Dinv = D_full
    call dgetrf(nelems,nelems,Dinv,nelems,ipiv,info)
    if (info /= 0) then
       write(luout,*) 'Factorization failed of classical D matrix'
       ERROR STOP 'Factorization failed of classical D matrix'
    end if

    allocate(work(nelems))
    call dgetri(nelems,Dinv,nelems,ipiv,work,nelems,info)
    if (info /= 0) then
       write(luout,*) 'Inversion failed of classical D matrix'
       ERROR STOP 'Inversion failed of classical D matrix'
    end if

    allocate(Minds_int(nelems,3))
    allocate(Ftot(nelems,3))
    Ftot = 0.0
    Minds_int = 0.0

    if (pelib_sol .and. .not. any(lpol)) then
        Ftot(:,:) = - Vs(:,:)
    else if (pelib_sol .and. any(lpol)) then
        Ftot(1:3*npols,1:3) = Fks(1:3*npols,1:3)
        Ftot(3*npols+1:3*npols+nsurp,1:3) = - Vs(1:nsurp,1:3)
    else if (.not. pelib_sol .and. any(lpol)) then
        Ftot(:,:) = Fks(:,:)
    end if
    do i = 1, 3 ! Loop over x,y,z
        Minds_int(:,i) = matmul(Dinv,Ftot(:,i))
        if (any(lpol)) then
            Mkinds(1:3*npols,i) = Minds_int(1:3*npols,i)
            do l = 1, nsurp
                Qinds(l,i) = Minds_int(3*npols+l,i)*Sa(l)*eps_fac
            end do
        else
            do l = 1, nsurp
                Qinds(l,i) = Minds_int(l,i)*Sa(l)*eps_fac
            end do
        end if
    end do

    if (pelib_debug) then
       do j = 1, 3
           dqsum = 0.0
           do i = 1, nelems
               dqsum = dqsum + Qinds(i,j)
               write(luout,*) 'dqsum', i, j, Qinds(i,j),dqsum
           end do
       end do
    end if

end subroutine direct_solver_eef
!------------------------------------------------------------------------------
subroutine direct_solver(Mkinds, Fks, Qinds, Vs)

    use pelib_options
    use pelib_utils, only: openfile
    use pelib_lapack_interfaces, only: lansp, pptrf, sptrf, pptrs, sptrs, ppcon, spcon

    real(rp), dimension(:,:), intent(out), optional :: Mkinds
    real(rp), dimension(:,:), intent(in), optional :: Fks
    real(rp), dimension(:,:), intent(out), optional :: Qinds
    real(rp), dimension(:,:), intent(in), optional :: Vs

    logical(lp) :: lexist
    integer(ip) :: nelems, lu, info
    integer(ip), dimension(:), allocatable :: ipiv
    real(rp) :: anorm, rcond
    real(rp), dimension(:), allocatable :: B
    real(rp), dimension(:,:), allocatable :: F
    real(rp) :: eps_fac

    if (present(Mkinds) .and. present(Qinds)) then
        nelems = 3 * npols + nsurp
    else if (present(Mkinds) .and. .not. present(Qinds)) then
        nelems = 3 * npols
    else if (.not. present(Mkinds) .and. present(Qinds)) then
        nelems = nsurp
    end if

    allocate(B(nelems*(nelems+1)/2))
    allocate(F(nelems,ndens))

    inquire(file='pelib_response_matrix.bin', exist=lexist)
    if (lexist .and. ((scf_cycle > 1) .or. pelib_restart)) then
        call openfile('pelib_response_matrix.bin', lu, 'old', 'unformatted')
        rewind(lu)
        if (chol) then
            read(lu) B
        else
            allocate(ipiv(nelems))
            read(lu) B, ipiv
        end if
        close(lu)
    else
        call response_matrix(B)
        if (pelib_debug) then
            anorm = lansp('1', B, 'L')
            write(luout, '(/4x,a,f15.8)') '1-norm of response matrix B: ', anorm
        end if

        if (chol) then
            call pptrf(B, 'L', info)
            if (info /= 0) then
                write(luout,*) 'INFO: Cholesky factorization of classical&
                               & response matrix failed.'
                write(luout,*) 'INFO: Cholesky factorization can be disabled with&
                               & NOCHOL under .DIRECT.'
                error stop 'cannot create classical response matrix.'
            end if
        else
            allocate(ipiv(nelems))
            call sptrf(B, 'L', ipiv, info)
            if (info /= 0) then
                error stop 'cannot create classical response matrix.'
            end if
        end if

        if (pelib_debug) then
            if (chol) then
                call ppcon(B, anorm, rcond, 'L')
                write(luout, '(4x,a,f15.8/)') 'Condition number of response matrix B: ', 1.0 / rcond
            else
                call spcon(B, ipiv, anorm, rcond, 'L')
                write(luout, '(4x,a,f15.8/)') 'Condition number of response matrix B: ', 1.0 / rcond
            end if
        end if

        call openfile('pelib_response_matrix.bin', lu, 'unknown', 'unformatted')
        rewind(lu)
        if (chol) then
            write(lu) B
        else
            write(lu) B, ipiv
        end if
        close(lu)

    end if

    if (present(Mkinds) .and. present(Qinds)) then
        F(1:3*npols,:) = Fks
        F(3*npols+1:3*npols+nsurp,:) = Vs
    else if (present(Mkinds) .and. .not. present(Qinds)) then
        F = Fks
    else if (.not. present(Mkinds) .and. present(Qinds)) then
        F = Vs
    end if

    if (chol) then
        call pptrs(B, F, 'L', info)
        if (info /= 0) then
            write(luout,*) 'INFO: Cholesky solver failed.'
            write(luout,*) 'INFO: Cholesky can be disabled with NOCHOL under .DIRECT.'
            error stop 'cannot solve for induced moments.'
        end if
        deallocate(B)
    else
        call sptrs(B, F, ipiv, 'L', info)
        if (info /= 0) then
            error stop 'cannot solve for induced moments.'
        end if
        deallocate(B, ipiv)
    end if

    if (pelib_sol) then
        if (response .and. pelib_noneq) then
            eps_fac = (epsinf - 1.0) / epsinf
        else if (fock .or. energy .or. (response .and. .not. pelib_noneq)) then
            eps_fac = (eps - 1.0) / eps
        end if
    end if

    if (present(Mkinds) .and. present(Qinds)) then
        Mkinds = F(1:3*npols,:)
        Qinds = eps_fac*F(3*npols+1:3*npols+nsurp,:)
    else if (present(Mkinds) .and. .not. present(Qinds)) then
        Mkinds = F
    else if (.not. present(Mkinds) .and. present(Qinds)) then
        Qinds = eps_fac*F
    end if

    if (fock) then
        if (present(Mkinds)) then
            call openfile('pelib_induced_moments.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Mkinds
            close(lu)
        end if
        if (present(Qinds)) then
            call openfile('pelib_induced_charges.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Qinds
            close(lu)
        end if

    end if

end subroutine direct_solver

!------------------------------------------------------------------------------

subroutine response_matrix(B)

    use pelib_options
    use pelib_constants
    use pelib_multipole_tools, only: damping_coefficient
    use pelib_blas_interfaces, only: nrm2
    use pelib_lapack_interfaces, only: pptrf, sptrf, pptri, sptri
    use pelib_cavity_generators, only: ntsatm

! TODO: Cutoff radius

    real(rp), dimension(:), intent(out) :: B

    logical(lp) :: exclude
    integer(ip) :: info
    integer(ip) :: i, j, k, l, m
    integer(ip), dimension(3) :: ipiv
    real(rp), parameter :: fourpi = 4.0 * pi
    real(rp), parameter :: i3 = 1.0 / 3.0
    real(rp), parameter :: i6 = 1.0 / 6.0
    real(rp) :: dism0
    real(rp) :: fE = 1.0
    real(rp) :: fT = 1.0
    real(rp) :: v, alpha_i, alpha_j
    real(rp) :: R, R3, R5, T
    real(rp), dimension(3) :: Rij
    real(rp), dimension(6) :: P11inv

    B = 0.0

    if (pelib_fixsol) then
        if (ntsatm ==  60) dism0 = 0.4 * aa2bohr
        if (ntsatm == 240) dism0 = 0.2 * aa2bohr
        if (ntsatm == 960) dism0 = 0.1 * aa2bohr
    end if

    m = 0
    do i = 1, nsites
        if (zeroalphas(i)) cycle
        P11inv = P11s(:,i)
        call pptrf(P11inv, 'L', info)
        if (info /= 0) then
            P11inv = P11s(:,i)
            call sptrf(P11inv, 'L', ipiv, info)
            if (info /= 0) then
                error stop 'could not factorize polarizability.'
            else if (chol) then
                chol = .false.
            end if
            call sptri(P11inv, ipiv, 'L')
        else
            call pptri(P11inv, 'L')
        end if
        if (pelib_ind_damp .or. pelib_amoeba_damp) then
            alpha_i = (P11s(1,i) + P11s(4,i) + P11s(6,i)) * i3
        end if
        do l = 3, 1, - 1
            do j = i, nsites
                if (zeroalphas(j)) cycle
                if (j == i) then
                    if (l == 3) then
                        do k = 1, l
                            B(m+k) = P11inv(k)
                        end do
                    else if (l == 2) then
                        do k = 1, l
                            B(m+k) = P11inv(3+k)
                        end do
                    else if (l == 1) then
                        do k = 1, l
                            B(m+k) = P11inv(5+k)
                        end do
                    end if
                    m = m + l
                else
                    if (pelib_nomb) then
                        m = m + 3
                        cycle
                    end if
                    exclude = .false.
                    do k = 1, pol_lexlst
                        if (pol_exclists(k,i) == pol_exclists(1,j)) then
                            exclude = .true.
                            exit
                        end if
                    end do
                    if (exclude) then
                        m = m + 3
                        cycle
                    end if
                    Rij = Rs(:,j) - Rs(:,i)
                    R = nrm2(Rij)
                    R3 = R**3
                    R5 = R**5
! TODO: cutoff radius
!                        if (R > cutoff) then
!                            m = m + 3
!                            cycle
!                        end if
                    if (pelib_amoeba_damp) then
                        alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                        call damping_coefficient(Rij, alpha_i, alpha_j, amoeba_damp, v)
                        fE = 1.0 - exp(-v**3)
                        fT = 1.0 - (1.0+v**3) * exp(-v**3)
                    else if (pelib_ind_damp) then
                        ! Thole damping
                        ! JPC A 102 (1998) 2399 and Mol. Sim. 32 (2006) 471
                        ! v = a * u, where a = 2.1304 (default) and u = R / (alpha_i * alpha_j)**(1/6)
                        ! fE = 1-(v^2/2+v+1)*exp(-v)
                        ! fT = 1-(v^3/6+v^2/2+v+1)*exp(-v)
                        alpha_j = (P11s(1,j) + P11s(4,j) + P11s(6,j)) * i3
                        call damping_coefficient(Rij, alpha_i, alpha_j, ind_damp, v)
                        fE = 1.0 - (0.5 * v**2 + v + 1.0) * exp(-v)
                        fT = fE - i6 * v**3 * exp(-v)
                    end if
                    if (l == 3) then
                        do k = 1, 3
                            T = 3.0 * Rij(1) * Rij(k) * fT / R5
                            if (k == 1) T = T - fE / R3
                            B(m+k) = - T
                        end do
                    else if (l == 2) then
                        do k = 1, 3
                            T = 3.0 * Rij(2) * Rij(k) * fT / R5
                            if (k == 2) T = T - fE / R3
                            B(m+k) = - T
                        end do
                    else if (l == 1) then
                        do k = 1, 3
                            T = 3.0 * Rij(3) * Rij(k) * fT / R5
                            if (k == 3) T = T - fE / R3
                            B(m+k) = - T
                        end do
                    end if
                    m = m + 3
                end if ! i /= j
            end do  ! do j = i, nsites
            if (pelib_sol) then
                do j = 1, nsurp
                    Rij = Rsp(:,j) - Rs(:,i)
                    R3 = nrm2(Rij)**3
                    if (l == 3) then
                        B(m+1) =  Rij(1) / R3
                        if (pelib_nomuq) B(m+1) = 0.0
                    else if (l == 2) then
                        B(m+1) =  Rij(2) / R3
                        if (pelib_nomuq) B(m+1) = 0.0
                    else if (l == 1) then
                        B(m+1) =  Rij(3) / R3
                        if (pelib_nomuq) B(m+1) = 0.0
                    end if
                    m = m + 1
                end do ! j = 1, nsurp
            end if ! pelib_sol
        end do  ! do l = 3, 1, - 1
    end do  ! do i = 1, nsites
    if (pelib_sol) then
        do i = 1, nsurp
            B(m+1) = 1.07 * sqrt(fourpi / Sa(i))
            m = m + 1
            do j = i + 1, nsurp
                Rij = Rsp(:,j) - Rsp(:,i)
                R = nrm2(Rij)
                if (pelib_fixsol) then
                    if (R > dism0) then
                        B(m+1) = ( 1.0 / R )
                    else
                        B(m+1) = (1.5 / dism0 - 0.5 * R**2 / dism0**3)
                    end if
                    m = m + 1
                else
                    B(m+1) = ( 1.0 / R )
                    m = m + 1
                end if
            end do ! j = i + 1, nsurp
        end do ! i = 1, nsurp
    end if

    if (pelib_debug) then
        do i = 1, (3 * npols + nsurp) * (3 * npols + nsurp + 1) / 2
            write (luout,*) 'Response matrix(i)',i, B(i)
        end do
    end if

end subroutine response_matrix

!------------------------------------------------------------------------------

subroutine d_local(D)

    use pelib_options
    use pelib_constants
    use pelib_blas_interfaces, only: nrm2, dot
    use pelib_lapack_interfaces, only: pptrf, sptrf, pptri, sptri
!    use pelib_cavity_generators, only: ntsatm

    real(rp), dimension(:,:), intent(out) :: D
    real(rp), dimension(:,:), allocatable :: Dtmp
    real(rp), parameter :: fourpi = 4.0 * pi
    real(rp), parameter :: twopi = 2.0 * pi
    real(rp), dimension(3) :: Rij
    real(rp) :: R, R3, Rtmp
    real(rp) :: x_atom, y_atom, z_atom
    real(rp) :: xi_ts, yi_ts, zi_ts
    real(rp) :: xj_ts, yj_ts, zj_ts
    real(rp) :: nrmj_x, nrmj_y, nrmj_z, deltaij
    integer(ip) :: i, j, atomi,atomj, id

    allocate(Dtmp(nsurp,nsurp))
    D = 0.0
    Dtmp = 0.0

    do i = 1, nsurp
        atomi = idatm(i)
        xi_ts   = Rsp(1,i)
        yi_ts   = Rsp(2,i)
        zi_ts   = Rsp(3,i)
        Dtmp(i,i) = - (1.07 * sqrt(Sa(i)/fourpi)*Sa(i))/(2.0*Radsp(atomi))
        do j = 1, nsurp
            if (i==j) cycle
            atomj = idatm(j)
            Rij = Rsp(:,i) - Rsp(:,j)
            R = nrm2(Rij)
            R3 = R**3

            x_atom = Rm(1,atomj)
            y_atom = Rm(2,atomj)
            z_atom = Rm(3,atomj)

            xj_ts   = Rsp(1,j)
            yj_ts   = Rsp(2,j)
            zj_ts   = Rsp(3,j)

            nrmj_x  = (xj_ts - x_atom)/Radsp(atomj)
            nrmj_y  = (yj_ts - y_atom)/Radsp(atomj)
            nrmj_z  = (zj_ts - z_atom)/Radsp(atomj)

            Rtmp = (xi_ts - xj_ts)*nrmj_x + (yi_ts - yj_ts)*nrmj_y + (zi_ts - zj_ts)*nrmj_z
            Dtmp(i,j) = (Sa(i)*Sa(j)*Rtmp)/(fourpi*R3)
        end do ! j = i + 1, nsurp
    end do ! i = 1, nsurp

    do i = 1, nsurp
        do j = 1, nsurp
            deltaij = 0.0
            if (i==j) then
                deltaij = Sa(i) * 0.5
            end if
            !if (i /= j) Dtmp(j,i) = 0.0
            D(i,j) = - fourpi*(deltaij - Dtmp(j,i))/Sa(i)
            !D(i,j) = - D(i,j)
        end do
    end do

    if (pelib_debug) then
        write(luout,*) 'nsurp * nsurp ', nsurp * nsurp
        id = 0
        do i = 1,  nsurp
            do j = 1, nsurp
                id = id +1
                write (luout,*) 'Response matrix(i) D',id, D(i,j)
            end do
        end do
    end if

end subroutine d_local

!------------------------------------------------------------------------------

subroutine diis_solver(Qinds, Vs, Mkinds, Fks)

    use pelib_mpi
    use pelib_options
    use pelib_constants
    use pelib_utils, only: openfile
    use pelib_cavity_generators, only: ntsatm
    use pelib_blas_interfaces, only: nrm2, dot

    real(rp), dimension(:,:), intent(out) :: Qinds
    real(rp), dimension(:,:), intent(in) :: Vs
    real(rp), dimension(:,:), intent(out), optional :: Mkinds
    real(rp), dimension(:,:), intent(in), optional :: Fks

    integer(ip) :: lu, itdiis
    integer(ip) :: i, j, k, n, m
    logical(lp) :: lexist, restart, converged
    integer(ip), parameter :: mxdiis = 100
    integer(ip), dimension(:), allocatable :: ipvt
    real(rp) :: factor, eps_fac, error
    real(rp) :: dism0, vtmp
    real(rp) :: R, R2, R3, oner3, Rji_x, Rji_y, Rji_z
    real(rp), dimension(3) :: Rji
    real(rp), dimension(nsurp) :: qfix, qnew, vfix2
    real(rp), dimension(:,:), allocatable :: tmpmat, dimat, field2
    real(rp), dimension(:), allocatable :: tmp
    real(rp), dimension(:,:,:), allocatable :: qrep

    if (any(lpol)) then
        if (.not. present(Mkinds) .or. .not. present(Fks)) then
            error stop 'in- and/or output pelib_options are missing in diis_solver'
        else
            allocate(field2(3*npols, ndens))
        end if
    end if
    if (myid == master) then
        allocate(tmpmat(mxdiis+1,mxdiis+1))
        allocate(dimat(mxdiis+1,mxdiis+1))
        allocate(tmp(mxdiis+1))
        allocate(qrep(nsurp,mxdiis,2))
        allocate(ipvt(mxdiis+1))
    end if

    if (myid == master) then
        inquire(file='pelib_induced_charges.bin', exist=lexist)
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(lexist, 1, lmpi, master, comm, ierr)
    end if
#endif

    if (lexist .and. (fock .or. energy) .and. (scf_cycle > 1 .or. pelib_restart)) then
        restart = .true.
    else
        restart = .false.
    end if

    if (restart) then
        if (myid == master) then
            call openfile('pelib_induced_charges.bin', lu, 'old', 'unformatted')
            rewind(lu)
            read(lu) qfix
            close(lu)
        end if
    else
        qfix = 0.0
    end if

    if (ntsatm ==  60) dism0 = 0.4 * aa2bohr
    if (ntsatm == 240) dism0 = 0.2 * aa2bohr
    if (ntsatm == 960) dism0 = 0.1 * aa2bohr
    factor =  1.0 / sqrt(4.0 * pi) / 1.07
    if (response .and. pelib_noneq) then
        eps_fac = (epsinf - 1.0) / epsinf
    else if (fock .or. energy .or. (response .and. .not. pelib_noneq)) then
        eps_fac = (eps - 1.0) / eps
    end if

    if (pelib_debug) then
        write(luout,*) 'diis_solver; eps_fac', eps_fac
        write(luout,*) 'diis_solver; factor', factor
    end if

    do n = 1, ndens

#if defined(VAR_MPI)
        if (nprocs > 1) then
            call mpi_bcast(Vs(1,n), nsurp, rmpi, master, comm, ierr)
            if (lpol(1,1)) call mpi_bcast(Fks(1,n), 3 * npols, rmpi, master, comm, ierr)
        end if
#endif

        if (.not. restart) then
            j = 1
            do i = surp_start, surp_finish
                qfix(j) = factor * sqrt(Sa(i)) * Vs(i,n)
                j = j + 1
            end do
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + surpdists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, dummy_int, rmpi, qfix(1), surpdists, displs, rmpi, master, comm, ierr)
            else if (myid /= master) then
                call mpi_gatherv(qfix(1), surpdists(myid), rmpi, dummy_real, surpdists, displs, rmpi, master, comm, ierr)
            end if
#endif
        end if

        do itdiis = 1, mxdiis
#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + surpdists(i-1)
                end do
                call mpi_scatterv(qfix(1), surpdists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
            else if (myid /= master) then
                call mpi_scatterv(dummy_real, surpdists, displs, rmpi, qfix(1), surpdists(myid), rmpi, master, comm, ierr)
            end if
#endif

!           Start Jacobi:
!           Mkinds_diis(:,itdiis+1) = M0*( F - M1* Mkinds_temp)
!     -- VFIX2: POTENTIAL AT NFFTS DUE TO SURFACE CHARGES --
            vfix2 = 0.0
            do i = 1, nsurp
                if (Sa(i) <=  zero) cycle
                k = 1
                do j = surp_start, surp_finish
                    if (i == j .or. Sa(j) <= zero) then
                        k = k + 1
                        cycle
                    end if
                    Rji = Rsp(:,i) - Rsp(:,j)
                    R = nrm2(Rji)
                    if (R > dism0) then
                        vfix2(i) = vfix2(i) + qfix(k) / R
                    else
                        vtmp = 1.5 / dism0 - 0.5 * R**2 / dism0**3
                        vfix2(i) = vfix2(i) + qfix(k) * vtmp
                    end if
                    k = k + 1
                end do
            end do

#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                call mpi_reduce(mpi_in_place, vfix2(1), nsurp, rmpi, mpi_sum, master, comm, ierr)
            else if (myid /= master) then
                call mpi_reduce(vfix2(1), dummy_real, nsurp, rmpi, mpi_sum, master, comm, ierr)
            end if
#endif

            if (lpol(1,1)) then
            ! calculate field at induced dipoles due to surface charges
                field2 = 0.0
                if (.not. pelib_nomuq) then
                    m = 1
                    do i = 1, nsites
                        if (zeroalphas(i)) cycle
                        k = 1
                        do j = surp_start, surp_finish
                            Rji = Rs(:,i) - Rsp(:,j)
                            R2 = Rji(1)**2 + Rji(2)**2 + Rji(3)**2
                            R3 = sqrt(R2) * R2
                            oner3 = 1.0 / R3
                            Rji_x =  Rji(1) * qfix(k) * eps_fac * oner3
                            Rji_y =  Rji(2) * qfix(k) * eps_fac * oner3
                            Rji_z =  Rji(3) * qfix(k) * eps_fac * oner3
                            field2(m,n) = field2(m,n) + Rji_x
                            field2(m+1,n) = field2(m+1,n) + Rji_y
                            field2(m+2,n) = field2(m+2,n) + Rji_z
                            k = k + 1
                        end do
                        m = m + 3
                    end do
                end if

#if defined(VAR_MPI)
                if (myid == master .and. nprocs > 1) then
                    call mpi_reduce(mpi_in_place, field2(1,n), 3 * npols, rmpi, mpi_sum, master, comm, ierr)
                else if (myid /= master) then
                    call mpi_reduce(field2(1,n), dummy_real, 3 * npols, rmpi, mpi_sum, master, comm, ierr)
                end if
#endif

! Get induced dipoles
                field2(:,n) = field2(:,n) + Fks(:,n)
                call iterative_solver(Mkinds(1:3*npols,n:n), field2(:,n:n))

! Potential at surface charge due to induced dipoles
#if defined(VAR_MPI)
                if (nprocs > 1) then
                    call mpi_bcast(Mkinds(1,n), 3 * npols, rmpi, master, comm, ierr)
                end if
                if (myid == master .and. nprocs > 1) then
                    displs(0) = 0
                    do i = 1, nprocs
                        displs(i) = displs(i-1) + surpdists(i-1)
                    end do
                    call mpi_scatterv(vfix2(1), surpdists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
                else if (myid /= master) then
                    call mpi_scatterv(dummy_real, surpdists, displs, rmpi, vfix2(1), surpdists(myid), rmpi, master, comm, ierr)
                end if
#endif

                if (.not. pelib_nomuq) then
                    k = 1
                    do i = surp_start, surp_finish
                        m = 1
                        do j = 1, nsites
                            if (zeroalphas(j)) cycle
                            Rji = Rsp(:,i) - Rs(:,j)
                            R2 = Rji(1)**2 + Rji(2)**2 + Rji(3)**2
                            R3 = sqrt(R2) * R2
                            oner3 = 1.0 / R3
                            vtmp =  dot(Mkinds(m:m+2,n), Rji) * oner3
                            vfix2(k) = vfix2(k) + vtmp
                            m = m + 3
                        end do
                        k = k + 1
                    end do
                end if

#if defined(VAR_MPI)
                if (myid == master .and. nprocs > 1) then
                    displs(0) = 0
                    do i = 1, nprocs
                        displs(i) = displs(i-1) + surpdists(i-1)
                    end do
                    call mpi_gatherv(mpi_in_place, dummy_int, rmpi, vfix2(1), surpdists, displs, rmpi, master, comm, ierr)
                else if (myid /= master) then
                    call mpi_gatherv(vfix2(1), surpdists(myid), rmpi, dummy_real, surpdists, displs, rmpi, master, comm, ierr)
                end if
#endif
            end if

#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + surpdists(i-1)
                end do
                call mpi_scatterv(vfix2(1), surpdists, displs, rmpi, mpi_in_place, dummy_int, rmpi, master, comm, ierr)
            else if (myid /= master) then
                call mpi_scatterv(dummy_real, surpdists, displs, rmpi, vfix2(1), surpdists(myid), rmpi, master, comm, ierr)
            end if
#endif

            j = 1
            do i = surp_start, surp_finish
                qnew(j) =  - factor * sqrt(Sa(i)) * (Vs(i,n) + vfix2(j))
                j = j + 1
            end do

#if defined(VAR_MPI)
            if (myid == master .and. nprocs > 1) then
                displs(0) = 0
                do i = 1, nprocs
                    displs(i) = displs(i-1) + surpdists(i-1)
                end do
                call mpi_gatherv(mpi_in_place, dummy_int, rmpi, qnew(1), surpdists, displs, rmpi, master, comm, ierr)
            else if (myid /= master) then
                call mpi_gatherv(qnew(1), surpdists(myid), rmpi, dummy_real, surpdists, displs, rmpi, master, comm, ierr)
            end if
#endif

!      Check convergance
            if (myid == master) then
                error = 0.0
                do i = 1, nsurp
                    error = error + abs(qfix(i) - qnew(i))
                end do
                if (error < nsurp * thrdiis) then
                    if (pelib_verbose) then
                        write (luout,'(4x,a,i2,a)') 'Induced surface charges converged in ', itdiis, ' iterations.'
                    end if
                    converged = .true.
                    Qinds(:,n) = eps_fac * qnew(:)
                else if (itdiis == mxdiis) then
                    write(luout, *) 'could not converge induced surface charges,'
                    write(luout, *) '       current threshold: ', itdiis
                    error stop 'could not converge induced surface charges.'
                else
                ! solve DIIS equations
                    converged = .false.
                    call fixdiis(nsurp, itdiis, mxdiis, qnew, qfix, dimat, qrep, tmp, tmpmat, ipvt, nsurp)
                    call dcopy(nsurp, qnew(1), 1, qfix(1), 1)
                end if
            end if
#if defined(VAR_MPI)
            if (nprocs > 1) then
                call mpi_bcast(converged, 1, lmpi, master, comm, ierr)
            end if
#endif
            if (converged) then
                exit
            end if

        end do ! itdiis
    end do ! n= 1, ndens

    if (fock) then
        if (myid == master) then
            call openfile('pelib_induced_charges.bin', lu, 'unknown', 'unformatted')
            rewind(lu)
            write(lu) Qinds
            close(lu)
        end if
    end if

end subroutine diis_solver

!-----------------------------------------------------------------------------

subroutine fixdiis(nffpar, nit, mxdiis, qout, qin, dimat, qrep, tmp, tmpmat, ipvt, nsize)

    use pelib_blas_interfaces, only: dot

    real(rp), dimension(:) :: qin, qout
    real(rp), dimension(:,:) :: dimat, tmpmat
    real(rp), dimension(:,:,:) :: qrep
    real(rp), dimension(:) :: tmp
    real(rp), dimension(:), allocatable :: bla, bla2
    integer(ip) :: nitmax, nit0, i0, nffme, info
! me used for parallelisation in games here we just set it to zero
    integer(ip) :: me, nit, i, j, nffpar, nsize, mxdiis, imax, imin
    integer(ip), dimension(:) :: ipvt

      nitmax = min(nit, mxdiis)
      nit0 = mod(nit-1,mxdiis) + 1

      me = 0
      imin=me*nffpar+1
      imax=min((me+1)*nffpar,nsize)
      nffme=imax-imin+1
!
!     -- store the charges
!        qrep(,,1)=qout
!        qrep(,,2)=qout-qin
!
      call dcopy(nffme,qout(imin),1,qrep(1,nit0,1),1)
      call dcopy(nffme,qout(imin),1,qrep(1,nit0,2),1)
      call daxpy(nffme,-1.0,qin(imin),1,qrep(1,nit0,2),1)
!
!     -- upgrade the interpolation matrix
!
      if(nit.gt.mxdiis) then
         do i=1,mxdiis-1
            do j=1,mxdiis-1
               dimat(i+1,j+1)=dimat(i+2,j+2)
            enddo
         enddo
      end if
!
      i0=nit0
      allocate(bla(nffme))
      allocate(bla2(nffme))
      bla = qrep(:,nit0,2)
      do i=nitmax,1,-1
         bla2 = qrep(:,i0,2)
         tmp(i)=dot(bla,bla2)
         i0=i0-1
         if(i0.eq.0) i0=i0+mxdiis
      enddo
      deallocate(bla,bla2)
!
      dimat(nitmax+1,1)=-1.0
      dimat(1,nitmax+1)=-1.0
      do i=nitmax,1,-1
         dimat(nitmax+1,i+1)=tmp(i)
         dimat(i+1,nitmax+1)=tmp(i)
      enddo
!
!     -- at the first iteration only matrix initialization
!
      if (nit.eq.1) then
         dimat(1,1)=0.0
         return
      end if
!
!     -- vector initialization
!
      tmp = 0.0
      tmp(1)=-1.0
!
!     -- copy the matrix (should be destroyed)
!
      do i=1,nitmax+1
         do j=1,nitmax+1
            tmpmat(i,j)=dimat(i,j)
         enddo
      enddo
!
!     -- solve the linear system
!
!      call dgetrf(mxdiis+1,mxdiis+1,tmpmat,mxdiis+1,ipvt,info)
      call dgefa(tmpmat(1,1),mxdiis+1,nitmax+1,ipvt(1),info)
      if (info.ne.0) then
         write(*,*) 'singular matrix in fixdiis.'
         error stop 'singular matrix in fixdiis.'
      end if
!      call dgetrs('n', mxdiis+1, 1, tmpmat, mxdiis+1, ipvt, tmp, mxdiis+1, info)
      call dgesl(tmpmat(1,1),mxdiis+1,nitmax+1,ipvt(1),tmp(1),0)
!
!     -- interpolate
!
      qout = 0.0
      i0=nit0
      do i=nitmax,1,-1
         call daxpy(nffme,tmp(i+1),qrep(1,i0,1),1,qout(imin),1)
         i0=i0-1
         if(i0.eq.0) i0=i0+mxdiis
      enddo
      return
!
end subroutine fixdiis

!------------------------------------------------------------------------------

end module pelib_induced_moments
