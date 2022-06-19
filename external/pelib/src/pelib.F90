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
module pelib

    use pelib_precision

    implicit none

    private

    logical(lp), public, save :: pelib_enabled = .false.

    ! public subroutines/functions
    public :: pelib_init, pelib_finalize, pelib_master
#if defined(VAR_MPI)
    public :: pelib_slave
#endif

contains

!------------------------------------------------------------------------------

subroutine pelib_init(lupri, coords, charges)

    ! Initialization routine for the PE library.

    use pelib_options
    use pelib_constants
    use pelib_integral_interfaces
    use pelib_blas_interfaces, only: nrm2
    use pelib_utils, only: charge2mass, charge2elem

    integer(ip) :: lupri
    real(rp), dimension(:), intent(in), optional :: charges
    real(rp), dimension(:,:), intent(in), optional :: coords

    integer(ip) :: i, j, k, l
    integer(ip) :: idx, jdx, kdx, nidx
    integer(ip), dimension(:), allocatable :: idxs
    logical(lp), dimension(:), allocatable :: redists
    logical(lp) :: lexist, skip
    real(rp) :: rclose, mass, totmas, redist_charge
    real(rp), dimension(3) :: Rsm, Rji
    character(len=80) :: fd

    if (allocated(Rm) .and. allocated(Zm)) then
        Rm(:,:) = coords
        synced = .false.
        scf_cycle = 0
        core_com = 0.0
        if (pelib_gauge) then
           !write(luout, '(a,i5)') 'INFO: The gauge for multipoles is changed  ', gauge_input
          core_com = gauge_input
        else
           totmas = 0.0
           do i = 1, nnucs
               mass = charge2mass(Zm(1,i))
               totmas = totmas + mass
               core_com = core_com + mass * Rm(:,i)
           end do
           core_com = core_com / totmas
        end if
        call init_pelib_integral_interfaces(core_com, Rm)
        return
    end if

    luout = lupri

    initialized = .true.

    if (present(coords) .and. present(charges)) then
        nnucs = size(charges)
        allocate(Rm(3,nnucs), Zm(1,nnucs), core_com(3))
        Rm(:,:) = coords
        Zm(1,:) = charges
        core_com = 0.0
        if (pelib_gauge) then
           !write(luout, '(a,i5)') 'INFO: The gauge for multipoles is changed  ', gauge_input
          core_com = gauge_input
        else
           totmas = 0.0
           do i = 1, nnucs
               mass = charge2mass(Zm(1,i))
               totmas = totmas + mass
               core_com = core_com + mass * Rm(:,i)
           end do
           core_com = core_com / totmas
        end if
        call init_pelib_integral_interfaces(core_com, Rm)
    else if (present(coords) .and. .not. present(charges)) then
        error stop 'pelib_init: coords present but charges missing'
    else if (.not. present(coords) .and. present(charges)) then
        error stop 'pelib_init: charges present but coords missing'
    end if

    ! setting up grid for MEP and CUBE calculation
    if ((pelib_mep .and. .not. mep_input) .or. pelib_cube) then
        origin(1) = minval(Rm(1,:)) - xsize
        origin(2) = minval(Rm(2,:)) - ysize
        origin(3) = minval(Rm(3,:)) - zsize
        step(1) = 1.0 / xgrid
        step(2) = 1.0 / ygrid
        step(3) = 1.0 / zgrid
        xsteps = int((maxval(Rm(1,:)) + xsize - origin(1)) / step(1)) + 1
        ysteps = int((maxval(Rm(2,:)) + ysize - origin(2)) / step(2)) + 1
        zsteps = int((maxval(Rm(3,:)) + zsize - origin(3)) / step(3)) + 1
        npoints = 0
        npoints = xsteps * ysteps * zsteps
        allocate(Rg(3,npoints))
        l = 1
        do i = 1, xsteps
            do j = 1, ysteps
                do k = 1, zsteps
                    Rg(1,l) = origin(1) + (i - 1) * step(1)
                    Rg(2,l) = origin(2) + (j - 1) * step(2)
                    Rg(3,l) = origin(3) + (k - 1) * step(3)
                    l = l + 1
                end do
            end do
        end do
    end if

    call read_potential(trim(potfile))

    if (pelib_sol) then
        call setup_solvent()
        call setup_cavity()
        if (pelib_diis) then
            thriter = thrdiis
        end if
        if (pelib_read_surf .or. .not. pelib_fixsol) call read_surface(trim(surfile))
    end if

    if (pelib_mul_damp .or. pelib_core_damp) then
        call setup_damping()
    end if

 10 format(4x,a)
 20 format(/6x,a)
 30 format(/8x,a)
 40 format(8x,a)
 50 format(6x,a,f8.4)
 60 format(6x,a,es8.1)

    write(luout, *)
    write(luout, 10) '.------------------------------------------------.'
    write(luout, 10)  '| Information from polarizable embedding library |'
    write(luout, 10) '`------------------------------------------------´'
    if (nsites > 0) then
        write(luout, '(/6x,a,i6)') 'Number of classical sites: ', nsites
    end if
    if (lmul(5)) then
        write(luout, 20) 'Multipole moments upto 5th order.'
    else if (lmul(4)) then
        write(luout, 20) 'Multipole moments upto 4th order.'
    else if (lmul(3)) then
        write(luout, 20) 'Multipole moments upto 3rd order.'
    else if (lmul(2)) then
        write(luout, 20) 'Multipole moments upto 2nd order.'
    else if (lmul(1)) then
        write(luout, 20) 'Multipole moments upto 1st order.'
    else if (lmul(0)) then
        write(luout, 20) 'Multipole moments upto 0th order.'
    end if
    if (pelib_polar) then
        if (lpol(1,1)) then
            write(luout, 20) 'Dipole-dipole polarizabilities.'
        end if
        if (pelib_gspol) then
            write(luout, 20) 'Dynamic response from environment will be&
                             & neglected during response calculation.'
        end if
        if (pelib_nomb) then
            write(luout, 20) 'Many-body interactions will be neglected.'
        end if
        if (pelib_iter) then
            if (pelib_diis) then
                write(luout, 20) 'Iterative DIIS solver for induced moments&
                                 & will be used'
                write(luout, 60) 'with convergence threshold:', thrdiis
            else
                write(luout, 20) 'Iterative solver for induced moments will&
                                 & be used'
                write(luout, 60) 'with convergence threshold:', thriter
            end if
            if (pelib_redthr) then
                write(luout, 20) 'Using reduced threshold in first few&
                                 & SCF iterations.'
            end if
        else
            write(luout, 20) 'Direct solver for induced moments will be used.'
        end if
        if (pelib_ind_damp) then
            write(luout, 20) 'Interactions between inducible moments will&
                             & be damped using Tholes scheme.'
            write(luout, 50) 'Damping coefficient:', ind_damp
        end if
        if (pelib_mul_damp) then
            write(luout, 20) 'Interactions between permanent and inducible&
                             & moments will be damped using Tholes scheme.'
            write(luout, 50) 'Damping coefficient:', mul_damp
        end if
        if (pelib_core_damp) then
            write(luout, 20) 'Interactions between electrons/nuclei and&
                                   & inducible moments will be damped using&
                                   & Tholes scheme.'
            write(luout, 50) 'Damping coefficient:', core_damp
            write(luout,'(6x,a)') 'Core-region polarizabilities:'
            do i = 1, nnucs
                write(luout,'(12x,a,1x,f8.4)') charge2elem(Zm(1,i)), core_alphas(1,i)
            enddo
        endif
        if (pelib_amoeba_damp) then
            write(luout, 20) 'Interactions will be damped using AMOEBAs Thole scheme.'
            write(luout, 50) 'Damping coefficient:', amoeba_damp
        end if
    end if
    if (pelib_fd) then
        if (pelib_fdes) then
            write(luout,'(/6x,a)') 'Electrostatic interactions with fragment&
                                   & densities are included'
        end if
        if (pelib_repuls) then
            write(luout,'(/6x,a,f6.3)') 'Repulsion operator will be used for&
                                        & fragment densities'
        end if
    end if
    if (pelib_sol) then
        write(luout,'(/6x,3a)') 'Continuum solvation in ', trim(solvent), &
                                & ' solvent modeled by FixSol.'
        write(luout,50) 'Static dielectric constant:', eps
        write(luout,50) 'Optical dielectric constant:', epsinf
        write(luout,'(6x,a,i8)') 'Number of surface points:', nsurp
        write(luout,'(6x,a,f12.3)') 'Surface area (in Aa**2):', sum(Sa) * bohr2aa**2
    end if
    if (pelib_restart) then
         write(luout, 20) 'Existing files will be used to restart if possible.'
    end if
    if (pelib_cube .and. .not. cube_field) then
        write(luout, 20) 'Cube files containing the potential from the&
                         & embedding potential will be written.'
    else if (pelib_cube .and. cube_field) then
        write(luout, 20) 'Cube files containing the potential and electric&
                         & field from the embedding potential will be written.'
    end if

    ! handling sites near border
    ! -----------------------------------------------
    if (pelib_border) then
        ! first locate all sites within given threshold of core nuclei
        allocate(idxs(nsites))
        idxs = 0; nidx = 0
        do i = 1, nnucs
            do j = 1, nsites
                lexist = .false.
                do k = 1, nidx
                    if (j == idxs(k)) then
                        lexist = .true.
                        exit
                    end if
                end do
                if (lexist) cycle
                Rsm = Rm(:,i) - Rs(:,j)
                if (nrm2(Rsm) <= Rmin) then
                    nidx = nidx + 1
                    idxs(nidx) = j
                end if
            end do
        end do

        if (border_type == 'REMOVE') then
            write(luout, *) ''
            do i = 1, nidx
                fd = '(8x,a,i6)'
                write(luout, fd) 'Removing all parameters on site:', idxs(i)
                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0
                endif
                if (lpol(1,1)) then
                    P11s(:,idxs(i)) = 0.0
                end if
            end do
        else if (border_type == 'REDIST') then
            allocate(redists(nsites))
            redists = .false.
            do i = 1, nidx
                rclose = 1.0e10
                do j = 1, nsites
                    lexist = .false.
                    do k = 1, nidx
                        if (j == idxs(k)) then
                            lexist = .true.
                            exit
                        end if
                    end do
                    if (lexist) cycle
                    Rji = Rs(:,idxs(i)) - Rs(:,j)
                    if (nrm2(Rji) <= rclose) then
                        rclose = nrm2(Rji)
                        idx = j
                    end if
                end do

                if (abs(redist_order) >= 1) then
                    if (lmul(0)) then
                        M0s(:,idx) = M0s(:,idx) + M0s(:,idxs(i)) / real(nredist, rp)
                    endif
                end if
                if (abs(redist_order) >= 2) then
                    if (lmul(1)) then
                        M1s(:,idx) = M1s(:,idx) + M1s(:,idxs(i)) / real(nredist, rp)
                    endif
                end if
                if (abs(redist_order) >= 3) then
                    if (lmul(2)) then
                        M2s(:,idx) = M2s(:,idx) + M2s(:,idxs(i)) / real(nredist, rp)
                    endif
                end if
                if (abs(redist_order) >= 4) then
                    if (lmul(3)) then
                        M3s(:,idx) = M3s(:,idx) + M3s(:,idxs(i)) / real(nredist, rp)
                    endif
                end if
                if (abs(redist_order) >= 5) then
                    if (lmul(4)) then
                        M4s(:,idx) = M4s(:,idx) + M4s(:,idxs(i)) / real(nredist, rp)
                    endif
                end if
                if (abs(redist_order) >= 6) then
                    if (lmul(5)) then
                        M5s(:,idx) = M5s(:,idx) + M5s(:,idxs(i)) / real(nredist, rp)
                    endif
                end if
                if (redist_order <= - 1) then
                    if (lpol(1,1)) then
                        P11s(:,idx) = P11s(:,idx) + P11s(:,idxs(i)) / real(nredist, rp)
                    end if
                end if

                if (nredist > 1) then
                    rclose = 1.0e10
                    do j = 1, nsites
                        if (j == idx) cycle
                        lexist = .false.
                        do k = 1, nidx
                            if (j == idxs(k)) then
                                lexist = .true.
                                exit
                            end if
                        end do
                        if (lexist) cycle
                        Rji = Rs(:,idxs(i)) - Rs(:,j)
                        if (nrm2(Rji) <= rclose) then
                            rclose = nrm2(Rji)
                            jdx = j
                        end if
                    end do
                    if (abs(redist_order) >= 1) then
                        if (lmul(0)) then
                            M0s(:,jdx) = M0s(:,jdx) + M0s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 2) then
                        if (lmul(1)) then
                            M1s(:,jdx) = M1s(:,jdx) + M1s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 3) then
                        if (lmul(2)) then
                            M2s(:,jdx) = M2s(:,jdx) + M2s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 4) then
                        if (lmul(3)) then
                            M3s(:,jdx) = M3s(:,jdx) + M3s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 5) then
                        if (lmul(4)) then
                            M4s(:,jdx) = M4s(:,jdx) + M4s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 6) then
                        if (lmul(5)) then
                            M5s(:,jdx) = M5s(:,jdx) + M5s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (redist_order <= - 1) then
                        if (lpol(1,1)) then
                            P11s(:,jdx) = P11s(:,jdx) + P11s(:,idxs(i)) / real(nredist, rp)
                        end if
                    end if
                end if

                if (nredist > 2) then
                    rclose = 1.0e10
                    do j = 1, nsites
                        if (j == idx .or. j == jdx) cycle
                        lexist = .false.
                        do k = 1, nidx
                            if (j == idxs(k)) then
                                lexist = .true.
                                exit
                            end if
                        end do
                        if (lexist) cycle
                        Rji = Rs(:,idxs(i)) - Rs(:,j)
                        if (nrm2(Rji) <= rclose) then
                            rclose = nrm2(Rji)
                            kdx = j
                        end if
                    end do
                    if (abs(redist_order) >= 1) then
                        if (lmul(0)) then
                            M0s(:,kdx) = M0s(:,kdx) + M0s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 2) then
                        if (lmul(1)) then
                            M1s(:,kdx) = M1s(:,kdx) + M1s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 3) then
                        if (lmul(2)) then
                            M2s(:,kdx) = M2s(:,kdx) + M2s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 4) then
                        if (lmul(3)) then
                            M3s(:,kdx) = M3s(:,kdx) + M3s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 5) then
                        if (lmul(4)) then
                            M4s(:,kdx) = M4s(:,kdx) + M4s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (abs(redist_order) >= 6) then
                        if (lmul(5)) then
                            M5s(:,kdx) = M5s(:,kdx) + M5s(:,idxs(i)) / real(nredist, rp)
                        endif
                    end if
                    if (redist_order <= - 1) then
                        if (lpol(1,1)) then
                            P11s(:,kdx) = P11s(:,kdx) + P11s(:,idxs(i)) / real(nredist, rp)
                        end if
                    end if
                end if

                if (lmul(0)) then
                    M0s(:,idxs(i)) = 0.0
                endif
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0
                endif
                if (lpol(1,1)) then
                    P11s(:,idxs(i)) = 0.0
                end if

                if (redist_order == 1) then
                    write(luout, 20) 'Redistributing multipoles upto 0th order'
                else if (redist_order == 2) then
                    write(luout, 20) 'Redistributing multipoles upto 1st order'
                else if (redist_order == 3) then
                    write(luout, 20) 'Redistributing multipoles upto 2nd order'
                else if (redist_order == 4) then
                    write(luout, 20) 'Redistributing multipoles upto 3rd order'
                else if (redist_order == 5) then
                    write(luout, 20) 'Redistributing multipoles upto 4th order'
                else if (redist_order == 6) then
                    write(luout, 20) 'Redistributing multipoles upto 5th order'
                else if (redist_order == - 1) then
                    write(luout, 20) 'Redistributing multipoles upto 0th&
                                     & order and polarizabilities'
                else if (redist_order == - 2) then
                    write(luout, 20) 'Redistributing multipoles upto 1st&
                                     & order and polarizabilities'
                else if (redist_order == - 3) then
                    write(luout, 20) 'Redistributing multipoles upto 2nd&
                                     & order and polarizabilities'
                else if (redist_order == - 4) then
                    write(luout, 20) 'Redistributing multipoles upto 3rd&
                                     & order and polarizabilities'
                else if (redist_order == - 5) then
                    write(luout, 20) 'Redistributing multipoles upto 4th&
                                     & order and polarizabilities'
                else if (redist_order == - 6) then
                    write(luout, 20) 'Redistributing multipoles upto 5th&
                                     & order and polarizabilities'
                end if
                write(luout, '(4x,a,i6)') 'from site:', idxs(i)
                fd = '(8x,a,3i6)'
                if (nredist == 3) then
                    write(luout, fd) 'to neighbouring sites:', idx, jdx, kdx
                else if (nredist == 2) then
                    write(luout, fd) 'to neighbouring sites:', idx, jdx
                else if (nredist == 1) then
                    write(luout, fd) 'to neighbouring site:', idx
                end if
                redists(idx) = .true.
                if (nredist > 1) redists(jdx) = .true.
                if (nredist > 2) redists(kdx) = .true.
            end do
            if (abs(redist_order) >= 1) then
                if (lmul(0)) then
                    write(luout, 30) ' Resulting monopoles: '
                    write(luout, 40) '----------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,f9.4)'
                            write(luout, fd) elems(i), i, M0s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 2) then
                if (lmul(1)) then
                    write(luout, 30) ' Resulting dipoles: '
                    write(luout, 40) '--------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,3f9.4)'
                            write(luout, fd) elems(i), i, M1s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 3) then
                if (lmul(2)) then
                    write(luout, 30) ' Resulting quadrupoles: '
                    write(luout, 40) '------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,6f9.4)'
                            write(luout, fd) elems(i), i, M2s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 4) then
                if (lmul(3)) then
                    write(luout, 30) ' Resulting octopoles: '
                    write(luout, 40) '----------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,10f9.4)'
                            write(luout, fd) elems(i), i, M3s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 5) then
                if (lmul(4)) then
                    write(luout, 30) ' Resulting hexadecapoles: '
                    write(luout, 40) '--------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,15f9.4)'
                            write(luout, fd) elems(i), i, M4s(:,i)
                        end if
                    end do
                end if
            end if
            if (abs(redist_order) >= 6) then
                if (lmul(5)) then
                    write(luout, 30) ' Resulting ditriacontapoles: '
                    write(luout, 40) '-----------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,21f9.4)'
                            write(luout, fd) elems(i), i, M5s(:,i)
                        end if
                    end do
                end if
            end if
            if (redist_order <= - 1) then
                if (lpol(1,1)) then
                    write(luout, 30) ' Resulting polarizabilities: '
                    write(luout, 40) '-----------------------------'
                    do i = 1, nsites
                        if (redists(i)) then
                            fd = '(8x,a,1x,i6,2x,6f9.4)'
                            write(luout, fd) elems(i), i, P11s(:,i)
                        end if
                    end do
                end if
            end if
            deallocate(redists)
        else if (border_type == 'CHGRED') then
            write(luout, *) ''
            redist_charge = 0.0
            do i = 1, nidx
                fd = '(6x,a,i6)'
                if (lmul(0)) then
                    write(luout, fd) 'Redistributing charges from site:', idxs(i)
                    redist_charge = redist_charge + M0s(1,idxs(i))
                    M0s(:,idxs(i)) = 0.0
                end if
                if (lmul(1)) then
                    M1s(:,idxs(i)) = 0.0
                endif
                if (lmul(2)) then
                    M2s(:,idxs(i)) = 0.0
                endif
                if (lmul(3)) then
                    M3s(:,idxs(i)) = 0.0
                endif
                if (lmul(4)) then
                    M4s(:,idxs(i)) = 0.0
                endif
                if (lmul(5)) then
                    M5s(:,idxs(i)) = 0.0
                endif
                if (lpol(1,1)) then
                    P11s(:,idxs(i)) = 0.0
                end if
            end do
            redist_charge = redist_charge / (nsites - nidx)
            skip = .false.
            do i = 1, nsites
                do j = 1, nidx
                    if (i == idxs(j)) then
                        skip = .true.
                        exit
                    end if
                end do
                if (skip) then
                    skip = .false.
                    cycle
                end if
                M0s(1,i) = M0s(1,i) + redist_charge
            end do
        end if
        deallocate(idxs)
        write(luout, *) ''
    end if

    if (lmul(0)) then
        redist_charge = 0.0
        do i = 1, nsites
            redist_charge = redist_charge + M0s(1,i)
        end do
        fd = '(6x,a,f9.4)'
        write(luout, fd) 'Total charge of classical region: ', redist_charge
    end if

    write(luout, *) ''

    ! number of polarizabilities different from zero
    if (any(lpol)) then
        allocate(zeroalphas(nsites))
        do i = 1, nsites
            if (maxval(abs(P11s(:,i))) <= zero) then
                zeroalphas(i) = .true.
            else
                zeroalphas(i) = .false.
                npols = npols + 1
            end if
        end do
    end if

end subroutine pelib_init

!------------------------------------------------------------------------------

subroutine pelib_finalize()

    use pelib_options
    use pelib_integral_interfaces

    initialized = .false.

    if (allocated(Rm)) deallocate(Rm)
    if (allocated(Zm)) deallocate(Zm)
    if (allocated(Rg)) deallocate(Rg)
    if (allocated(elems)) deallocate(elems)
    if (allocated(Rs)) deallocate(Rs)
    if (allocated(Zs)) deallocate(Zs)
    if (allocated(M0s)) deallocate(M0s)
    if (allocated(M1s)) deallocate(M1s)
    if (allocated(M2s)) deallocate(M2s)
    if (allocated(M3s)) deallocate(M3s)
    if (allocated(M4s)) deallocate(M4s)
    if (allocated(M5s)) deallocate(M5s)
    if (allocated(P11s)) deallocate(P11s)
    if (allocated(mul_exclists)) deallocate(mul_exclists)
    if (allocated(pol_exclists)) deallocate(pol_exclists)
    if (allocated(zeroalphas)) deallocate(zeroalphas)
    if (allocated(Epe)) deallocate(Epe)
    if (allocated(Ees)) deallocate(Ees)
    if (allocated(Epol)) deallocate(Epol)
    if (allocated(displs)) deallocate(displs)
    if (allocated(sitedists)) deallocate(sitedists)
    if (allocated(siteloops)) deallocate(siteloops)
    if (allocated(poldists)) deallocate(poldists)
    if (allocated(cubeloops)) deallocate(cubeloops)
    if (allocated(cubedists)) deallocate(cubedists)

    call finalize_pelib_integral_interfaces

end subroutine pelib_finalize

!------------------------------------------------------------------------------

! This routine serves as an example of how input reading can be performed in a
! host program.
!subroutine pelib_input_reader(word, luinp)
!
!    use pelib_options
!    use pelib_constants
!    use pelib_utils, only: chcase
!    use pelib_cavity_generators, only: ntsatm
!
!    character(len=*), intent(inout) :: word
!    integer(ip), intent(in) :: luinp
!
!    character(len=80) :: option
!    character(len=2) :: auoraa
!    integer(ip) :: i, j
!    real(rp), dimension(2) :: temp
!
!    do
!        read(luinp, '(a80)') option
!        call chcase(option)
!
!        ! Read potential (optionally from potfile)
!        if (trim(option(2:7)) == 'POTENT') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, '(a240)') potfile
!            end if
!        ! direct solver for induced moments
!        else if (trim(option(2:7)) == 'DIRECT') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, '(a80)') option
!                call chcase(option)
!                if (trim(option(1:7)) == 'NOCHOL') then
!                    chol = .false.
!                end if
!            end if
!            pelib_iter = .false.
!        ! iterative solver for induced moments (default)
!        else if (trim(option(2:7)) == 'ITERAT') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) thriter
!            end if
!            pelib_iter = .true.
!        else if (trim(option(2:7)) == 'DIIS T') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) thrdiis
!            end if
!            thriter = thrdiis
!        ! use reduced threshold in iterative induced moments solver
!        else if (trim(option(2:7)) == 'REDTHR') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) redlvl
!            end if
!            pelib_redthr = .true.
!        ! handling sites near quantum-classical border
!         else if (trim(option(2:7)) == 'BORDER') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, '(a)', advance='no') border_type
!                backspace(luinp)
!                call chcase(border_type)
!                if ((trim(border_type) /= 'REMOVE') .and.&
!                  & (trim(border_type) /= 'REDIST') .and.&
!                  & (trim(border_type) /= 'CHGRED')) then
!                    error stop 'unknown handling of border sites'
!                else if (trim(border_type) == 'REMOVE') then
!                    read(luinp, *) border_type, Rmin, auoraa
!                else if (trim(border_type) == 'REDIST') then
!                    read(luinp, *) border_type, redist_order, Rmin, auoraa, nredist
!                    if ((nredist > 3) .or. (nredist < 1)) then
!                        error stop 'parameters can only be distributed to&
!                             & minimum one site and maximum three sites'
!                    end if
!                else if (trim(border_type) == 'CHGRED') then
!                    read(luinp, *) border_type, Rmin, auoraa
!                else
!                    error stop 'unrecognized input in .BORDER option'
!                end if
!                call chcase(auoraa)
!                if (trim(auoraa) == 'AA') Rmin = Rmin * aa2bohr
!            end if
!            pelib_border = .true.
!        ! damp electric field from induced multipoles
!        else if (trim(option(2:7)) == 'DAMP I') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) ind_damp
!            end if
!            pelib_ind_damp = .true.
!        ! damp electric field from permanent multipoles
!        else if (trim(option(2:7)) == 'DAMP M') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) mul_damp
!            end if
!            pelib_mul_damp = .true.
!        ! damp electric field from core region
!        else if (trim(option(2:7)) == 'DAMP C') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) core_damp
!                ! attempt to optionally read a custom specification of the polarizabilities
!                read(luinp, '(a80)') option
!                backspace(luinp)
!                if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!                   & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                    read(luinp, *) j
!                    allocate(core_alphas(6,j))
!                    core_alphas = 0.0
!                    do i = 1, j
!                        read(luinp, *) core_alphas(1,i)
!                        core_alphas(4,i) = core_alphas(1,i)
!                        core_alphas(6,i) = core_alphas(1,i)
!                    enddo
!                end if
!            end if
!            pelib_core_damp = .true.
!        ! damp electric fields using AMOEABA-style Thole damping
!        else if (trim(option(2:7)) == 'DAMP A') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) amoeba_damp
!            end if
!            pelib_amoeba_damp = .true.
!        ! the old deprecated DAMP option
!        else if (trim(option(2:7)) == 'DAMP') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) ind_damp
!            end if
!            pelib_ind_damp = .true.
!            write(luout, *) 'INFO: the .DAMP option is deprecated, please use .DAMP INDUCED'
!        ! neglect dynamic response from environment
!        else if (trim(option(2:7)) == 'GSPOL') then
!            pelib_gspol = .true.
!        ! neglect many-body interactions
!        else if (trim(option(2:7)) == 'NOMB') then
!            pelib_nomb = .true.
!        ! Use existing files for restart
!        else if (trim(option(2:7)) == 'RESTAR') then
!            pelib_restart = .true.
!        ! calculate intermolecular two-electron integrals
!        else if (trim(option(2:7)) == 'TWOINT') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, '(a240)') h5pdefile
!            end if
!            pelib_twoints = .true.
!        ! save density matrix
!        else if (trim(option(2:7)) == 'SAVE D') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, '(a240)') h5pdefile
!            end if
!            pelib_savden = .true.
!        ! disable exchange repulsion
!        else if (trim(option(2:7)) == 'NO REP') then
!            pelib_repuls = .false.
!        ! electrostatics and exchange repulsion from fragment densities
!        else if (trim(option(2:7)) == 'PDE') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, '(a240)') h5pdefile
!            end if
!            pelib_fd = .true.
!        ! skip electrostatics from fragment densities
!        else if (trim(option(2:7)) == 'NO FD') then
!            pelib_fdes = .false.
!        ! request calculation of effective dipole integrals
!        else if (trim(option(2:7)) == 'EEF') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) ncrds
!                allocate(crds(3,ncrds))
!                do i = 1, ncrds
!                    read(luinp, *) (crds(j,i), j = 1, 3)
!                end do
!            end if
!            pelib_lf = .true.
!        ! provide LJ parameters for the QM region
!        else if (trim(option(2:7)) == 'LJ') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) qmLJsites
!                allocate(qmLJs(2,qmLJsites))
!                do i = 1, qmLJsites
!                    read(luinp, *) (temp(j), j = 1, 2)
!                    qmLJs(1,i) = temp(1) * 2.0
!                    qmLJs(2,i) = temp(2)
!                end do
!                lvdw = .true.
!             end if
!        ! skip QM calculations, i.e. go directly into PE library
!        else if (trim(option(2:7)) == 'SKIPQM') then
!            pelib_skipqm = .true.
!        ! Write cube files
!        else if (trim(option(2:7)) == 'CUBE') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!              & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                do
!                    read(luinp, '(a80)') option
!                    call chcase(option)
!                    if (trim(option(1:7)) == 'COARSE') then
!                        xgrid = 3
!                        ygrid = 3
!                        zgrid = 3
!                    else if (trim(option(1:7)) == 'MEDIUM') then
!                        xgrid = 6
!                        ygrid = 6
!                        zgrid = 6
!                    else if (trim(option(1:7)) == 'FINE') then
!                        xgrid = 12
!                        ygrid = 12
!                        zgrid = 12
!                    else if (trim(option(1:7)) == 'GRID') then
!                        read(luinp, *) xsize, xgrid, ysize, ygrid, zsize, zgrid
!                    else if (trim(option(1:7)) == 'FIELD') then
!                        cube_field = .true.
!                    else if (option(1:1) == '.' .or. option(1:1) == '*') then
!                        backspace(luinp)
!                        exit
!                    else if (option(1:1) == '!' .or. option(1:1) == '#') then
!                        cycle
!                    else
!                        error stop 'unknown input under .CUBE option'
!                    end if
!                end do
!            end if
!            pelib_cube = .true.
!        ! evaluate molecular electrostatic potential
!        else if (trim(option(2:7)) == 'MEP') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                do
!                    read(luinp, '(a80)') option
!                    call chcase(option)
!                    if (trim(option(1:7)) == 'COARSE') then
!                        xgrid = 3
!                        ygrid = 3
!                        zgrid = 3
!                    else if (trim(option(1:7)) == 'MEDIUM') then
!                        xgrid = 6
!                        ygrid = 6
!                        zgrid = 6
!                    else if (trim(option(1:7)) == 'FINE') then
!                        xgrid = 12
!                        ygrid = 12
!                        zgrid = 12
!                    else if (trim(option(1:7)) == 'GRID') then
!                        read(luinp, *) xsize, xgrid, ysize, ygrid, zsize, zgrid
!                    else if (trim(option(1:7)) == 'FIELD') then
!                        mep_field = .true.
!                    else if (trim(option(1:7)) == 'EXTFLD') then
!                        read(luinp, *) (extfld(i), i = 1, 3)
!                        mep_extfld = .true.
!                    else if (trim(option(1:7)) == 'LOCFLD') then
!                        read(luinp, *) lf_component
!                        mep_lf = .true.
!                    else if (trim(option(1:7)) == 'SKIPQM') then
!                        mep_qmcube = .false.
!                    else if (trim(option(1:7)) == 'SKIPMUL') then
!                        mep_mulcube = .false.
!                    else if (trim(option(1:7)) == 'INPUT') then
!                        mep_input = .true.
!                        read(luinp, '(a240)') h5gridfile
!                    else if (option(1:1) == '.' .or. option(1:1) == '*') then
!                        backspace(luinp)
!                        exit
!                    else if (option(1:1) == '!' .or. option(1:1) == '#') then
!                        cycle
!                    else
!                        error stop 'unknown option present in .MEP section.'
!                    end if
!                end do
!            end if
!            pelib_mep = .true.
!        ! continuum solvation calculation
!        else if (trim(option(2:7)) == 'SOLVAT') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) solvent
!            else
!                solvent = 'H2O'
!            end if
!            pelib_sol = .true.
!            pelib_diis = .true.
!            pelib_fixsol = .true.
!            pelib_polar = .true.
!            chol = .false.
!        ! equilibrium solvation (non-equilibrium is default)
!        else if (trim(option(2:7)) == 'EQSOL') then
!            pelib_noneq = .false.
!        ! Turn off interaction between induced dipoles and surface charges
!        else if (trim(option(2:7)) == 'NOMUQ') then
!            pelib_nomuq = .true.
!        else if (trim(option(2:7)) == 'NOVMU') then
!            pelib_novmu = .true.
!        ! specify surface file
!        else if (trim(option(2:7)) == 'SURFAC') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, '(a240)') surfile
!            end if
!            pelib_read_surf = .true.
!        ! FixSol solvation using FIXPVA2 tesselation (optional: number of tessera per atom)
!        else if (trim(option(2:7)) == 'NTESS ') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) ntsatm
!                if ((ntsatm .ne. 60) .and. (ntsatm .ne. 240) .and. (ntsatm .ne. 960)) then
!                   error stop 'number of tessera per atom must be equal to 60, 240 or 960'
!                end if
!            end if
!        ! apply external electric field
!        else if (trim(option(2:7)) == 'FIELD') then
!            pelib_field = .true.
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) (efields(i), i = 1, 3)
!            end if
!        ! verbose output
!        else if (trim(option(2:7)) == 'VERBOS') then
!            pelib_verbose = .true.
!        ! debug output
!        else if (trim(option(2:7)) == 'DEBUG') then
!            pelib_debug = .true.
!            pelib_verbose = .true.
!        ! isotropic polarizabilities
!        else if (trim(option(2:7)) == 'ISOPOL') then
!            pelib_isopol = .true.
!        ! zero out the polarizabilities
!        else if (trim(option(2:7)) == 'ZEROPO') then
!            pelib_zeropol = .true.
!        ! zero out higher-order multipoles
!        else if (trim(option(2:7)) == 'ZEROMU') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                read(luinp, *) zeromul_order
!                if (zeromul_order < 0) then
!                    error stop 'ZEROMUL order cannot be negative'
!                end if
!            end if
!            pelib_zeromul = .true.
!        ! skip calculation of multipole-multipole interaction energy
!        else if (trim(option(2:7)) == 'SKIPMU') then
!            pelib_skipmul = .true.
!        else if (trim(option(2:7)) == 'GAUGE') then
!            read(luinp, '(a80)') option
!            backspace(luinp)
!            if ((option(1:1) /= '.') .and. (option(1:1) /= '*') .and.&
!               & (option(1:1) /= '!') .and. (option(1:1) /= '#')) then
!                allocate(gauge_input(3))
!                gauge_input = 0.0
!                do i = 1,3
!                    read(luinp, *) gauge_input(i)
!                enddo
!            end if
!            pelib_gauge = .true.
!        else if (option(1:1) == '*') then
!            word = option
!            exit
!        else if (option(1:1) == '!' .or. option(1:1) == '#') then
!            cycle
!        else
!            write(luout, *) 'unknown option:', option
!            error stop 'unknown option in *PEQM section'
!        end if
!    end do
!
!! check options
!    if (pelib_diis .and. .not. pelib_iter) then ! Assume that you want to use a direct solver for FixSol
!        pelib_diis = .false.
!    end if
!
!    if (pelib_sol .and. ndens > 1) then
!        error stop 'Continuum solvation only implemented for ndens = 1'
!    end if
!
!    if (pelib_sol .and. pelib_iter .and. .not. pelib_fixsol) then
!        error stop '.SOLV without .FIXSOL requires .DIRECT'
!    end if
!
!    if (pelib_mep .and. pelib_cube) then
!        error stop '.MEP and .CUBE are not compatible.'
!    end if
!
!end subroutine pelib_input_reader

!------------------------------------------------------------------------------

subroutine read_potential(filename)

    use pelib_options
    use pelib_constants
    use pelib_utils, only: chcase, openfile, elem2charge

    character(len=*) :: filename

    integer(ip) :: i, j, s
    integer(ip) :: nlines
    integer(ip) :: lupot
    integer(ip) :: multype
    integer(ip), dimension(2) :: poltype
    integer(ip), dimension(:), allocatable :: itemp
    real(rp) :: trace
    real(rp), dimension(21) :: temp
    character(len=2) :: bohroraa
    character(len=80) :: word
    logical(lp) :: lexist

    inquire(file=filename, exist=lexist)
    if (lexist) then
        call openfile(filename, lupot, 'old', 'formatted')
    else if (pelib_sol) then
        return
    else
        write(luout, *) 'potential input file not found'
        error stop 'potential input file not found'
    end if

    write(luout, '(/2x,a)') 'Reading potential input file for polarizable embedding calculation'
    do
        read(lupot, *, end=112) word
        call chcase(word)
        if (trim(word) == '@COORDINATES') then
            read(lupot, *) nsites
            read(lupot, *) bohroraa
            allocate(elems(nsites), Zs(1,nsites), Rs(3,nsites))
            do i = 1, nsites
                read(lupot, *) elems(i), (Rs(j,i), j = 1, 3)
                Zs(1,i) = elem2charge(elems(i))
            end do
            call chcase(bohroraa)
            if (bohroraa == 'AA') then
                Rs = Rs * aa2bohr
            end if
            exit
        else
            cycle
        end if
112     error stop 'no coordinates found in potential input file'
    end do

    rewind(lupot)

    do
        read(lupot, *, end=100) word
        call chcase(word)
        if (word(1:1) /= '@') then
            cycle
        else if (trim(word) == '@COORDINATES') then
            cycle
        else if (trim(word) == '@MULTIPOLES') then
            do
                read(lupot, *, end=100) word
                call chcase(word)
                if (word(1:1) == '@') then
                    backspace(lupot)
                    exit
                else if (trim(word) /= 'ORDER' .and. trim(word) /= 'EXCLISTS') then
                    cycle
                else if (trim(word) == 'ORDER') then
                    backspace(lupot)
                    read(lupot, *) word, multype
                    if (multype == 0) then
                        if (pelib_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(0) = .true.
                        allocate(M0s(1,nsites))
                        M0s = 0.0
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, temp(1)
                            M0s(1,s) = temp(1)
                        end do
                    else if (multype == 1) then
                        if (pelib_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(1) = .true.
                        allocate(M1s(3,nsites))
                        M1s = 0.0
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 3)
                            M1s(:,s) = temp(1:3)
                        end do
                    else if (multype == 2) then
                        if (pelib_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(2) = .true.
                        allocate(M2s(6,nsites))
                        M2s = 0.0
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 6)
                            ! remove trace
                            trace = (temp(1) + temp(4) + temp(6)) / 3.0
                            temp(1) = temp(1) - trace
                            temp(4) = temp(4) - trace
                            temp(6) = temp(6) - trace
                            M2s(:,s) = temp(1:6)
                        end do
                    else if (multype == 3) then
                        if (pelib_zeromul .and. (zeromul_order <= multype)) cycle
                        lmul(3) = .true.
                        allocate(M3s(10,nsites))
                        M3s = 0.0
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 10)
                            trace = (temp(1) + temp(4) + temp(6)) / 5.0
                            temp(1) = temp(1) - 3.0 * trace
                            temp(4) = temp(4) - trace
                            temp(6) = temp(6) - trace
                            trace = (temp(2) + temp(7) + temp(9)) / 5.0
                            temp(2) = temp(2) - trace
                            temp(7) = temp(7) - 3.0 * trace
                            temp(9) = temp(9) - trace
                            trace = (temp(3) + temp(8) + temp(10)) / 5.0
                            temp(3) = temp(3) - trace
                            temp(8) = temp(8) - trace
                            temp(10) = temp(10) - 3.0 * trace
                            M3s(:,s) = temp(1:10)
                        end do
                    else if (multype == 4) then
                        if (pelib_zeromul .and. (zeromul_order <= multype)) cycle
                        write(luout, *) 'WARNING: results will be wrong if&
                                        & non-traceless hexadecapoles&
                                        & (16-poles) are used'
                        lmul(4) = .true.
                        allocate(M4s(15,nsites))
                        M4s = 0.0
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 15)
                            M4s(:,s) = temp(1:15)
                        end do
                    else if (multype == 5) then
                        if (pelib_zeromul .and. (zeromul_order <= multype)) cycle
                        write(luout, *) 'WARNING: results will be wrong if&
                                        & non-traceless ditriacontapoles&
                                        & (32-poles) are used'
                        lmul(5) = .true.
                        allocate(M5s(21,nsites))
                        M5s = 0.0
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 21)
                            M5s(:,s) = temp(1:21)
                        end do
                    else
                        write(luout, *) 'unsupported or unrecognized multipole order'
                        write(luout, *) multype
                        error stop 'unsupported or unrecognized multipole order'
                    end if
                else if (trim(word) == 'EXCLISTS') then
                    if (pelib_zeropol) cycle
                    read(lupot, *) nlines, mul_lexlst
                    allocate(mul_exclists(mul_lexlst,nsites))
                    allocate(itemp(mul_lexlst))
                    mul_exclists = 0
                    do i = 1, nsites
                        mul_exclists(1,i) = i
                    end do
                    if (mul_lexlst > 1) then
                        do i = 1, nlines
                            itemp = 0
                            read(lupot, *) s, (itemp(j), j = 2, mul_lexlst)
                            mul_exclists(2:mul_lexlst,s) = itemp(2:mul_lexlst)
                        end do
                    end if
                    deallocate(itemp)
                else
                    write(luout, *) 'unknown keyword in @MULTIPOLES section'
                    write(luout, *) trim(word)
                    error stop 'unknown keyword in @MULTIPOLES section'
                end if
            end do
        else if (trim(word) == '@POLARIZABILITIES') then
            do
                read(lupot, *, end=100) word
                call chcase(word)
                if (word(1:1) == '@') then
                    backspace(lupot)
                    exit
                else if (trim(word) /= 'ORDER' .and. trim(word) /= 'EXCLISTS') then
                    cycle
                else if (trim(word) == 'ORDER') then
                    backspace(lupot)
                    read(lupot, *) word, poltype
                    if (poltype(1) == 1 .and. poltype(2) == 1) then
                        if (pelib_zeropol) cycle
                        lpol(1,1) = .true.
                        pelib_polar = .true.
                        if (.not. allocated(P11s)) then
                            allocate(P11s(6,nsites))
                            P11s = 0.0
                        end if
                        read(lupot, *) nlines
                        do i = 1, nlines
                            read(lupot, *) s, (temp(j), j = 1, 6)
                            if (pelib_isopol) then
                                trace = temp(1) + temp(4) + temp(6)
                                temp = 0.0
                                temp(1) = trace / 3.0
                                temp(4) = temp(1)
                                temp(6) = temp(1)
                            end if
                            P11s(:,s) = temp(1:6)
                        end do
                    else if (word(1:1) == '!' .or. word(1:1) == '#') then
                        cycle
                    else
                        write(luout, *) 'unsupported or unrecognized polarizability order'
                        write(luout, *) poltype
                        error stop 'unsupported or unrecognized polarizability'
                    end if
                else if (trim(word) == 'EXCLISTS') then
                    if (pelib_zeropol) cycle
                    read(lupot, *) nlines, pol_lexlst
                    allocate(pol_exclists(pol_lexlst,nsites))
                    allocate(itemp(pol_lexlst))
                    pol_exclists = 0
                    do i = 1, nsites
                        pol_exclists(1,i) = i
                    end do
                    if (pol_lexlst > 1) then
                        do i = 1, nlines
                            itemp = 0
                            read(lupot, *) s, (itemp(j), j = 2, pol_lexlst)
                            pol_exclists(2:pol_lexlst,s) = itemp(2:pol_lexlst)
                        end do
                    end if
                    deallocate(itemp)
                else
                    write(luout, *) 'unknown keyword in @POLARIZABILITIES section'
                    write(luout, *) trim(word)
                    error stop 'unknown keyword in @POLARIZABILITIES section'
                end if
            end do
        else if (trim(word) == '@LJ') then
            !lvdw = .true.
            allocate(LJs(2,nsites))
            LJs = 0.0
            read(lupot,*) nlines
            do i = 1, nlines
                read(lupot,*) s, (temp(j), j = 1, 2)
                LJs(1,s) = temp(1) * 2.0
                LJs(2,s) = temp(2)
            end do
        else
            write(luout, *) 'unknown section in potential input file'
            write(luout, *) trim(word)
            error stop 'unknown section in potential input file'
        end if
    end do

100 continue

    close(lupot)

    ! default exclusion list (everything polarizes everything)
    if (.not. allocated(pol_exclists)) then
        pol_lexlst = 1
        allocate(pol_exclists(pol_lexlst,nsites))
        do i = 1, nsites
            pol_exclists(1,i) = i
        end do
        if (any(lpol)) then
            write(luout, *) 'WARNING: no exclusion lists provided which means all polarizable'
            write(luout, *) '         sites are polarized by all other induced multipoles'
        end if
    end if
    ! default multipole exclusion list is the same as the polarizability exclusion list (backwards compatibility)
    if (.not. allocated(mul_exclists)) then
        mul_lexlst = pol_lexlst
        allocate(mul_exclists(mul_lexlst,nsites))
        mul_exclists = pol_exclists
    end if

end subroutine read_potential

!------------------------------------------------------------------------------

subroutine read_surface(filename)

    use pelib_options
    use pelib_constants
    use pelib_utils, only: chcase, openfile
    use pelib_blas_interfaces, only: nrm2

    character(len=*) :: filename

    logical(lp) :: lexist
    integer(ip) :: i, j, lusurf
    character(len=2) :: bohroraa
    real(rp), dimension(3) :: Rji

    inquire(file=filename, exist=lexist)
    if (lexist) then
        call openfile(filename, lusurf, 'old', 'formatted')
    else
        write(luout, *) 'surface file not found: ', filename
        error stop 'surface file not found'
    end if

    read(lusurf, *) nsurp
    read(lusurf, *) bohroraa

    allocate(Rsp(3,nsurp))
    allocate(Sa(nsurp))

    do i = 1, nsurp
        read(lusurf, *) (Rsp(j,i), j = 1, 3), Sa(i)
    end do

    close(lusurf)

    call chcase(bohroraa)
    if (bohroraa == 'AA') then
        Rsp = Rsp * aa2bohr
        Sa = Sa * aa2bohr**2
    end if

    do i = 1, nsurp
       do j = 1, nsites
          Rji = Rs(:,j) - Rsp(:,i)
          if (nrm2(Rji) < 1.2 ) then
              write(luout, '(a,f12.8)') 'WARNING: Cavity to close to classical site:', nrm2(Rji)
              write(luout, '(a,f12.8)') 'Surface point:', Rsp(:,i)
              write(luout, '(a,f12.8)') 'Classical site:', Rs(:,j)
          end if
       end do
    end do

!    if (pelib_debug) then
!       write(luout, *) 'Rsp in read_surface, number of surface points:',nsurp
!       do i = 1, nsurp
!              write(luout, *) i, Rsp(:,i)
!       end do
!       write(luout, *) 'Sa in read_surface'
!       do i = 1, nsurp
!          write (luout, *) i, Sa(i)
!       end do
!       write(luout, *) 'Rsp in read_surface in AU'
!       do i = 1, nsurp
!              write (luout, *) Rsp(:,i)
!       end do
!       write(luout, *) 'Sa in read_surface in AU'
!       do i = 1, nsurp
!          write (luout, *) Sa(i)
!       end do
!    end if

end subroutine read_surface

!------------------------------------------------------------------------------

subroutine pelib_master(runtype, triang, ndim, nmats, denmats, idx, fckmats, expvals)

    use pelib_mpi
    use pelib_options
    use pelib_analysis_tools
    use pelib_operators

    character(*), intent(in) :: runtype
    logical(ip), intent(in), optional :: triang
    integer(ip), intent(in), optional :: ndim
    integer(ip), intent(in), optional :: nmats
    integer(ip), intent(in), optional :: idx
    real(rp), dimension(:), intent(in), optional :: denmats
    real(rp), dimension(:), intent(out), optional :: fckmats
    real(rp), dimension(:), intent(out), optional :: expvals

    real(rp), dimension(:), allocatable :: dummy

    ! get root id and number of MPI processes
#if defined(VAR_MPI)
    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)
    master = myid
    if (master /= 0) then
        error stop 'a slave cannot be master'
    end if
#else
    myid = 0
    nprocs = 1
    master = myid
#endif

    ! in-/output consistency checks
    if ((present(denmats) .or. present(fckmats)) .and. .not. (present(ndim) .and. present(triang))) then
        error stop 'triang and ndim are required input if either denmats or fckmats are provided'
    end if

    if (present(denmats) .and. .not. present(nmats)) then
        error stop 'nmats is required input if denmats is provided'
    end if

    if (present(expvals) .and. .not. present(denmats)) then
        error stop 'expvals requires denmats'
    end if

    ! set variables according to input
    if (present(triang)) then
        trimat = logical(triang, lp)
    else
        trimat = .true.
    end if

    if (.not. trimat) then
        error stop 'only triangular matrices in packed storage available'
    end if

    if (present(ndim)) then
        nbas = ndim
        nnbas = nbas * (nbas + 1) / 2
        n2bas = nbas * nbas
    else
        nbas = 0
        nnbas = 0
        n2bas = 0
    end if

    if (present(nmats)) then
        ndens = nmats
    else
        ndens = 1
    end if

    if (present(idx)) then
        if (ndens > 1) then
            error stop 'state_idx only valid for one density'
        else
            state_idx = idx
        end if
    else
        if (ndens > 1) then
            state_idx = -1
        else
            state_idx = 0
        end if
    end if

    if (present(denmats)) then
        lden = .true.
    else
        lden = .false.
    end if

    if (present(fckmats)) then
        fckmats = 0.0
        lfck = .true.
    else
        lfck = .false.
    end if

    if (present(expvals)) then
        expvals = 0.0
        lexp = .true.
    else
        lexp = .false.
    end if

    ! checks according to runtype
    if (present(denmats)) then
        if (trimat) then
            if (size(denmats) < ndens * nnbas) then
                error stop 'denmats too small'
            end if
        else if (.not. trimat) then
            if (size(denmats) < ndens * n2bas) then
                error stop 'denmats too small'
            end if
        end if
    end if

    fock = .false.
    energy = .false.
    response = .false.
    london = .false.
    molgrad = .false.
    mep = .false.
    effdipole = .false.
    cube = .false.

    if (runtype == 'full_fock') then
        fock = .true.
        scf_cycle = scf_cycle + 1
        if (.not. present(fckmats)) then
            error stop 'fckmats required for runtype=full_fock'
        else if (.not. present(expvals)) then
            error stop 'expvals required for runtype=full_fock'
        end if
        if (size(fckmats) < ndens * nnbas) then
            error stop 'fckmats too small'
        else if (size(expvals) < ndens) then
            error stop 'expvals too small'
        end if
    else if (runtype == 'print_energy') then
        energy = .true.
        if (.not. present(denmats)) then
            error stop 'denmats required for runtype=print_energy'
        end if
    else if (runtype == 'dynamic_response') then
        if (pelib_gspol) return
        response = .true.
        if (.not. present(denmats)) then
            error stop 'denmats required for runtype=dynamic_response'
        else if (.not. present(fckmats)) then
            error stop 'fckmats requires for runtype=dynamic_response'
        end if
        if (size(fckmats) < ndens * nnbas) then
            error stop 'fckmats too small'
        end if
    else if (runtype == 'magnetic_gradient') then
        if (ndens > 1) then
            error stop 'magnetic gradient not implemented for more than one density matrix'
        end if
        london = .true.
        if (.not. present(fckmats)) then
            error stop 'fckmats required for runtype=magnetic_gradient'
        end if
        if (size(fckmats) < 3 * nnbas) then
            error stop 'fckmats too small'
        end if
    else if (runtype == 'molecular_gradient') then
        if (ndens > 1) then
            error stop 'molecular gradient not implemented for more than one density matrix'
        end if
        molgrad = .true.
        if (present(expvals)) then
            if (size(expvals) < 3 * nnucs) then
                error stop 'expvals too small'
            end if
        end if
        if (present(fckmats)) then
            if (size(fckmats) < 3 * nnucs * nnbas) then
                error stop 'fckmats too small'
            end if
        end if
    else if (runtype == 'mep') then
        if (ndens > 1) then
            error stop 'MEP not implemented for more than one density matrix'
        end if
        mep = .true.
        if (mep_qmcube .and. .not. present(denmats)) then
            error stop 'denmats required for runtype=mep'
        end if
        if (mep_input .and. .not. present(denmats)) then
            error stop 'denmats required for runtype=mep'
        end if
    else if (runtype == 'cube') then
        if (ndens > 1) then
            error stop 'CUBE not implemented for more than one density&
   &              matrix (ndens>1)'
        end if
        if (state_idx < 0) then
            error stop 'CUBE not implemented for more than one density&
   &              matrix (state_idx<0)'
        end if
        cube = .true.
    else if (runtype == 'effdipole') then
        if (present(fckmats)) then
            if (size(fckmats) < 3 * nnbas) then
                error stop 'integrals matrix too small'
            end if
            effdipole = .true.
        end if
    else
        error stop 'unknown runtype in pelib_master'
    end if

    if (nprocs == 1) then
        site_start = 1
        site_finish = nsites
        surp_start = 1
        surp_finish = nsurp
        cube_start = 1
        cube_finish = npoints
    end if

#if defined(VAR_MPI)
    if (nprocs > 1) then
        call mpi_bcast(trimat, 1, lmpi, master, comm, ierr)
        call mpi_bcast(ndens, 1, impi, master, comm, ierr)
        call mpi_bcast(nbas, 1, impi, master, comm, ierr)
        call mpi_bcast(nnbas, 1, impi, master, comm, ierr)
        call mpi_bcast(n2bas, 1, impi, master, comm, ierr)
        call mpi_bcast(lden, 1, lmpi, master, comm, ierr)
        call mpi_bcast(lfck, 1, lmpi, master, comm, ierr)
        call mpi_bcast(lexp, 1, lmpi, master, comm, ierr)
        call mpi_bcast(state_idx, 1, lmpi, master, comm, ierr)

        if (lden) then
            call mpi_bcast(denmats(1), nnbas * ndens, rmpi, master, comm, ierr)
        end if

        if (fock) then
            call mpi_bcast(scf_cycle, 1, impi, master, comm, ierr)
        end if

        call mpi_bcast(synced, 1, lmpi, master, comm, ierr)
        if (.not. synced) then
            call mpi_sync()
        end if
    end if
#endif

    if (fock) then
        call get_fock(denmats=denmats, fckmats=fckmats, energies=expvals)
    else if (energy) then
        if (present(expvals)) then
            call get_fock(denmats=denmats, energies=expvals)
        else
            allocate(dummy(ndens))
            call get_fock(denmats=denmats, energies=dummy)
            deallocate(dummy)
        end if
        if (pelib_cube) then
            if (ndens > 1) then
                error stop 'CUBE not implemented for more than one density matrix (ndens>1)'
            end if
            if (state_idx < 0) then
                error stop 'CUBE not implemented for more than one&
   &                                    density matrix (state_idx<0)'
            end if
            call compute_potential_cube()
            if (cube_field) call compute_field_cube()
        end if
    else if (response) then
        call get_fock(denmats=denmats, fckmats=fckmats)
    else if (london) then
        call get_magnetic_gradient(fckmats)
    else if (molgrad) then
        if (present(expvals) .and. present(fckmats)) then
            call get_molecular_gradient(denmats=denmats, fckmats=fckmats, molgrads=expvals)
        else if (present(expvals) .and. .not. present(fckmats)) then
            call get_molecular_gradient(denmats=denmats, molgrads=expvals)
        else if (present(fckmats) .and. .not. present(expvals)) then
            call get_molecular_gradient(fckmats=fckmats)
        end if
    else if (mep) then
        if (mep_input) then
            call compute_mep_input(denmats)
        else if (mep_qmcube) then
            call compute_mep_cube(denmats)
            if (mep_field) call compute_mef_cube(denmats)
        else
            call compute_mep_cube()
            if (mep_field) call compute_mef_cube()
        end if
    else if (effdipole) then
        call get_effdipole(fckmats)
    else if (cube) then
        call compute_potential_cube(denmats)
    end if

end subroutine pelib_master

!------------------------------------------------------------------------------

#if defined(VAR_MPI)
subroutine pelib_slave(runtype)

    use pelib_mpi
    use pelib_options
    use pelib_analysis_tools
    use pelib_operators

    character(*), intent(in) :: runtype

    real(rp), dimension(:), allocatable :: denmats
    real(rp), dimension(:), allocatable :: fckmats
    real(rp), dimension(:), allocatable :: expvals

    call mpi_comm_rank(comm, myid, ierr)
    call mpi_comm_size(comm, nprocs, ierr)
    master = 0

    fock = .false.
    energy = .false.
    response = .false.
    london = .false.
    molgrad = .false.
    mep = .false.
    effdipole = .false.
    cube = .false.

    if (runtype == 'full_fock') then
        fock = .true.
    else if (runtype == 'print_energy') then
        energy = .true.
    else if (runtype == 'dynamic_response') then
        response = .true.
    else if (runtype == 'magnetic_gradient') then
        london = .true.
    else if (runtype == 'molecular_gradient') then
        molgrad = .true.
    else if (runtype == 'mep') then
        mep = .true.
    else if (runtype == 'effdipole') then
        effdipole = .true.
    else if (runtype == 'cube') then
        cube = .true.
    else
        error stop 'unknown runtype in pelib_slave'
    end if

    call mpi_bcast(trimat, 1, lmpi, master, comm, ierr)
    call mpi_bcast(ndens, 1, impi, master, comm, ierr)
    call mpi_bcast(nbas, 1, impi, master, comm, ierr)
    call mpi_bcast(nnbas, 1, impi, master, comm, ierr)
    call mpi_bcast(n2bas, 1, impi, master, comm, ierr)
    call mpi_bcast(lden, 1, lmpi, master, comm, ierr)
    call mpi_bcast(lfck, 1, lmpi, master, comm, ierr)
    call mpi_bcast(lexp, 1, lmpi, master, comm, ierr)
    call mpi_bcast(state_idx, 1, lmpi, master, comm, ierr)

    if (lden) then
        allocate(denmats(ndens*nnbas))
        call mpi_bcast(denmats(1), nnbas * ndens, rmpi, master, comm, ierr)
    end if

    if (fock .or. response) then
        allocate(fckmats(ndens*nnbas))
    else if (london) then
        allocate(fckmats(3*n2bas))
    else if (molgrad .and. lfck) then
        allocate(fckmats(3*nnucs*nnbas))
    else if (effdipole) then
        allocate(fckmats(3*nnbas))
    end if

    if (fock .or. energy) then
        allocate(expvals(ndens))
    else if (molgrad .and. lexp) then
        allocate(expvals(3*nnucs))
    end if

    if (lfck) fckmats = 0.0
    if (lexp) expvals = 0.0

    if (fock) then
        call mpi_bcast(scf_cycle, 1, impi, master, comm, ierr)
    end if

    call mpi_bcast(synced, 1, lmpi, master, comm, ierr)
    if (.not. synced) then
        call mpi_sync()
    end if

    if (fock) then
        call get_fock(denmats=denmats, fckmats=fckmats, energies=expvals)
        deallocate(denmats)
        deallocate(fckmats)
        deallocate(expvals)
    else if (energy) then
        call get_fock(denmats=denmats, energies=expvals)
        if (pelib_cube) then
            if (ndens > 1) then
                error stop 'CUBE not implemented for more than 1 density matrix'
            end if
            call compute_potential_cube()
            if (cube_field) call compute_field_cube()
        end if
        deallocate(denmats)
        deallocate(expvals)
    else if (response) then
        call get_fock(denmats=denmats, fckmats=fckmats)
        deallocate(denmats)
        deallocate(fckmats)
    else if (london) then
        call get_magnetic_gradient(fckmats)
        deallocate(fckmats)
    else if (molgrad) then
        if (lexp .and. lfck) then
            call get_molecular_gradient(denmats=denmats, fckmats=fckmats, molgrads=expvals)
        else if (lexp .and. .not. lfck) then
            call get_molecular_gradient(denmats=denmats, molgrads=expvals)
        else if (lfck .and. .not. lexp) then
            call get_molecular_gradient(fckmats=fckmats)
        end if
        deallocate(denmats)
        deallocate(expvals)
    else if (mep) then
        if (mep_input) then
            call compute_mep_input(denmats)
            deallocate(denmats)
        else if (mep_qmcube) then
            call compute_mep_cube(denmats)
            if (mep_field) call compute_mef_cube(denmats)
            deallocate(denmats)
        else
            call compute_mep_cube()
            if (mep_field) call compute_mef_cube()
        end if
    else if (effdipole) then
        call get_effdipole(fckmats)
        deallocate(fckmats)
    else if (cube) then
        if (pelib_cube) then
            call compute_potential_cube(denmats)
            deallocate(denmats)
        end if
    end if

end subroutine pelib_slave
#endif

!------------------------------------------------------------------------------

subroutine setup_solvent()

    use pelib_options

    integer(ip) :: i
    logical(lp) :: not_found = .true.
    character(len=8), dimension(18) :: solvents
    character(len=8), dimension(18) :: solvent_names
    real(rp), dimension(18) :: epslist, epsinflist

    solvents = ['H2O     ', 'CH3OH   ', 'C2H5OH  ', 'CHCL3   ', &
                'CH2CL2  ', 'C2H4CL2 ', 'CCL4    ', 'C6H6    ', &
                'C6H5CH3 ', 'C6H5CL  ', 'CH3NO2  ', 'C7H16   ', &
                'C6H12   ', 'C6H5NH2 ', 'CH3COCH3', 'THF     ', &
                'DMSO    ', 'CH3CN   ']
    solvent_names = ['WATER   ', 'METHANOL', 'ETHANOL ', 'CLFORM  ', &
                     'METHYLCL', '12DICLET', 'TETRACLC', 'BENZENE ', &
                     'TOLUENE ', 'CLBENZ  ', 'NITROMET', 'N-HEPTAN', &
                     'CYCLOHEX', 'ANILINE ', 'ACETONE ', 'TETHYDFU', &
                     'DIMETSOX', 'ACETONIT']
    ! water, methanol, ethanol, chloroform, methylenechloride,
    ! 1,2-dichloroethane, carbon tetrachloride, benzene, toluene,
    ! chlorobenzene, nitromethane, n-heptane, cyclohexane, aniline,
    ! acetone, tetrahydrofurane, dimethylsulfoxide, acetonitrile
    epslist = [78.39, 32.63, 24.55, 4.9, 8.93, 10.36, 2.228, 2.247, &
               2.379, 5.621, 38.20, 1.92, 2.023, 6.89, 20.7, 7.58, &
               46.7, 36.64]
    epsinflist = [1.776, 1.758, 1.847, 2.085, 2.02, 2.085, 2.129, 2.244, &
                  2.232, 2.32, 1.904, 1.918, 2.028, 2.506, 1.841, 1.971, &
                  2.179, 1.806]

    do i = 1, size(solvents)
        if (trim(solvent) == trim(solvents(i)) .or. trim(solvent) == trim(solvent_names(i))) then
            if (eps < 0.0) then
                eps = epslist(i)
            end if
            if (epsinf < 0.0) then
                epsinf = epsinflist(i)
            end if
            not_found = .false.
        end if
    end do

    if (not_found) then
        write(luout,*) 'Unknown solvent ', trim(solvent)
        write(luout,*) 'Currently available solvents are'
        do i = 1, size(solvents)
            write(luout,*) trim(solvents(i)), ' or the equivalent name ', trim(solvent_names(i))
        end do
        error stop 'ERROR: unknown solvent'
    end if

end subroutine setup_solvent

!------------------------------------------------------------------------------

subroutine setup_cavity()

    use pelib_options
    use pelib_cavity_generators, only: fixtes

    real(rp), dimension(:), allocatable :: all_charges

    allocate(all_coords(3,nnucs+nsites))
    allocate(all_charges(nnucs+nsites))

    all_coords(1:3,1:nnucs) = Rm(1:3,1:nnucs)
    all_charges(1:nnucs) = Zm(1,1:nnucs)
    if (nsites > 0) then
        all_coords(1:3,nnucs+1:nnucs+nsites) = Rs(1:3,1:nsites)
        all_charges(nnucs+1:nnucs+nsites) = Zs(1,1:nsites)
    end if

    call fixtes(all_coords, all_charges)

end subroutine setup_cavity

!------------------------------------------------------------------------------

subroutine setup_damping()

    use pelib_options

    integer(ip) :: i

    ! check if there are polarizabilities on all sites needed for damping field
    ! from permanent multipoles
    if (pelib_mul_damp) then
        do i = 1, nsites
            if (maxval(abs(P11s(:,i))) <= zero) then
                error stop 'damping field from permanent multipoles requires polarizabilities on all sites'
            end if
        end do
    end if

    ! if core polarizabilities were not read in, then attempt to use values from
    ! table 7 in Piet Th. van Duijnen and Marcel Swart, JPC A 102, 2399, 1998
    ! DOI: 10.1021/jp980221f
    if (pelib_core_damp .and. .not. allocated(core_alphas)) then
        allocate(core_alphas(6,nnucs))
        core_alphas = 0.0
        do i = 1, nnucs
            if (nint(Zm(1,i)) == 1) then
                core_alphas(1,i) = 2.7927
            end if
            if (nint(Zm(1,i)) == 6) then
                core_alphas(1,i) = 8.6959
            end if
            if (nint(Zm(1,i)) == 7) then
                core_alphas(1,i) = 6.5565
            end if
            if (nint(Zm(1,i)) == 8) then
                core_alphas(1,i) = 5.7494
            end if
            if (nint(Zm(1,i)) == 9) then
                core_alphas(1,i) = 3.0013
            end if
            if (nint(Zm(1,i)) == 16) then
                core_alphas(1,i) = 16.6984
            end if
            if (nint(Zm(1,i)) == 17) then
                core_alphas(1,i) = 16.1979
            end if
            if (nint(Zm(1,i)) == 35) then
                core_alphas(1,i) = 23.5714
            end if
            if (nint(Zm(1,i)) == 53) then
                core_alphas(1,i) = 36.9880
            end if
            if (core_alphas(1,i) <= zero) then
                write(luout, '(a,i5)') 'isotropic polarizability not available for core atom ', i
                error stop 'damping electric field from core failed'
            end if
            core_alphas(4,i) = core_alphas(1,i)
            core_alphas(6,i) = core_alphas(1,i)
        enddo
    else if (pelib_core_damp .and. allocated(core_alphas)) then
        if (size(core_alphas, 2) /= nnucs) then
            error stop 'number of core polarizabilities given in input does not match number of core atoms'
        end if
    end if

end subroutine setup_damping

!------------------------------------------------------------------------------

end module pelib
