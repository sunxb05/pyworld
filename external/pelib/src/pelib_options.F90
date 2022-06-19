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
module pelib_options

    use pelib_precision

    implicit none

    ! options
    logical(lp), public, save :: pelib_gspol = .false.
    logical(lp), public, save :: pelib_iter = .true.
    logical(lp), public, save :: pelib_diis = .false.
    logical(lp), public, save :: pelib_redthr = .false.
    logical(lp), public, save :: pelib_border = .false.
    logical(lp), public, save :: pelib_ind_damp = .false.
    logical(lp), public, save :: pelib_mul_damp = .false.
    logical(lp), public, save :: pelib_core_damp = .false.
    logical(lp), public, save :: pelib_amoeba_damp = .false.
    logical(lp), public, save :: pelib_nomb = .false.
    logical(lp), public, save :: pelib_polar = .false.
    logical(lp), public, save :: pelib_cube = .false.
    logical(lp), public, save :: pelib_mep = .false.
    logical(lp), public, save :: pelib_skipqm = .false.
    logical(lp), public, save :: pelib_twoints = .false.
    logical(lp), public, save :: pelib_repuls = .true.
    logical(lp), public, save :: pelib_savden = .false.
    logical(lp), public, save :: pelib_fd = .false.
    logical(lp), public, save :: pelib_fdes = .true.
    logical(lp), public, save :: pelib_sol = .false.
    logical(lp), public, save :: pelib_noneq = .true.
    logical(lp), public, save :: pelib_restart = .false.
    logical(lp), public, save :: pelib_verbose = .false.
    logical(lp), public, save :: pelib_debug = .false.
    logical(lp), public, save :: pelib_fixsol = .false.
    logical(lp), public, save :: pelib_nomuq = .false.
    logical(lp), public, save :: pelib_novmu = .false.
    logical(lp), public, save :: pelib_read_surf = .false.
    logical(lp), public, save :: pelib_zeromul = .false.
    logical(lp), public, save :: pelib_zeropol = .false.
    logical(lp), public, save :: pelib_isopol = .false.
    logical(lp), public, save :: pelib_gauge = .false.
    logical(lp), public, save :: pelib_lf = .false.
    logical(lp), public, save :: pelib_field = .false.
    logical(lp), public, save :: pelib_skipmul = .true.
    logical(lp), public, save :: pelib_fmm = .false.
    logical(lp), public, save :: pelib_old_field = .false.

    ! runtypes
    logical(lp) :: fock = .false.
    logical(lp) :: mixed = .false.
    logical(lp) :: energy = .false.
    logical(lp) :: response = .false.
    logical(lp) :: molgrad = .false.
    logical(lp) :: mep = .false.
    logical(lp) :: london = .false.
    logical(lp) :: effdipole = .false.
    logical(lp) :: cube = .false.

    ! in- and output logicals, i.e. denmats, fckmats and expvals
    logical(lp) :: lden = .false.
    logical(lp) :: lfck = .false.
    logical(lp) :: lexp = .false.

    ! matrix type
    logical(lp) :: trimat = .true.

    ! filenames
    character(len=240) :: potfile = 'POTENTIAL.INP'
    character(len=240) :: surfile = 'SURFACE.INP'
    character(len=240) :: h5pdefile = 'standard.h5'
    character(len=240) :: h5gridfile = 'grid.h5'

    logical(lp), save :: synced = .false.
    logical(lp), save :: initialized = .false.
    integer(ip), save :: site_start, site_finish
    integer(ip), save :: surp_start, surp_finish
    integer(ip), save :: cube_start, cube_finish
    integer(ip), dimension(:), save, allocatable :: siteloops
    integer(ip), dimension(:), save, allocatable :: surploops
    integer(ip), dimension(:), save, allocatable :: cubeloops
    integer(ip), dimension(:), save, allocatable :: poldists
    integer(ip), dimension(:), save, allocatable :: sitedists
    integer(ip), dimension(:), save, allocatable :: surpdists
    integer(ip), dimension(:), save, allocatable :: cubedists
    integer(ip), dimension(:), save, allocatable :: displs

    ! logical unit for output file (default is stdout)
    integer(ip), save :: luout = 6

    ! integer used to check memory allocation status,
    ! i.e. allocate(array(len), stat=mem_stat)
    integer(ip) :: mem_stat

    ! dummy variables
    integer(ip) :: dummy_int
    real(rp) :: dummy_real
    integer(ip), dimension(1) :: dummy_int_array
    real(rp), dimension(1) :: dummy_real_array

    ! constants, thresholds and stuff
    real(rp), parameter :: zero = 1.0e-8
    integer(ip), save :: scf_cycle = 0
    integer(ip), save :: redlvl = 0
    real(rp), save :: thriter = 1.0e-8
    real(rp), save :: thrdiis = 1.0e-8
    real(rp), save :: ind_damp = 2.1304
    real(rp), save :: mul_damp = 2.1304
    real(rp), save :: core_damp = 2.1304
    real(rp), save :: amoeba_damp = 0.39
    real(rp), save :: gauss_factor = 1.0
    real(rp), dimension(3), save :: repfacs = 1.0
    real(rp), save :: Rmin = 2.2
    integer(ip), save :: nredist = 1
    integer(ip), save :: redist_order = 1
    character(len=6), save :: border_type = 'REDIST'
    ! use Cholesky factorization of classical response matrix
    logical(lp), save :: chol = .true.
    ! solvent and dielectric constant (defaults to water)
    character(len=80) :: solvent
    real(rp), save :: eps = -huge(0.0_rp)
    real(rp), save :: epsinf = -huge(0.0_rp)
    ! the order from which multipoles are zeroed
    integer(ip), save :: zeromul_order = 1
    ! fmm constants
    real(rp), save :: fmm_theta = 0.4
    integer(ip), save :: fmm_expansion_order = 6
    integer(ip), save :: fmm_ncrit = 64

    ! polarizable embedding potential info
    ! ------------------------------------
    ! total number of classical sites
    integer(ip), save :: nsites = 0
    ! number of polarizable sites
    integer(ip), save :: npols = 0
    ! number of surface points
    integer(ip), save :: nsurp = 0
    ! number fragment densities
    integer(ip), save :: nfds = 0
    ! number of nuclei in fragment density
    integer(ip), save :: fdnucs = 0
    ! multipole exclusion list length
    integer(ip), save :: mul_lexlst = 0
    ! polarizability exclusion list length
    integer(ip), save :: pol_lexlst = 0
    ! number of density matrices
    integer(ip) :: ndens = 0
    ! number of basis functions in core fragment
    integer(ip), save :: nbas
    ! size of packed matrices
    integer(ip), save :: nnbas
    ! size of full matrices
    integer(ip), save :: n2bas
    ! number of nuclei in core region
    integer(ip), save :: nnucs = 0
    ! number of LJ sites in core region
    integer(ip), save :: qmLJsites = 0
    ! index of electronic state
    integer(ip), save :: state_idx = -1

    ! Thole qm damping stuff
    ! ----------------------
    integer(ip), save :: nqmpols = 0
    real(rp), dimension(:,:), allocatable, save :: qmpols
    real(rp), save :: qmdamp = 2.1304

    ! specifies what type of parameters are present
    ! lmul(0): monopoles, lmul(1): dipoles etc.
    logical(lp), dimension(0:5), save :: lmul = .false.
    ! lpol(1): dipole-dipole polarizabilities
    logical(lp), dimension(0:2,0:2), save :: lpol = .false.
    ! lhypol(1): dipole-dipole-dipole polarizabilities/1st hyperpolarizability
!    logical(lp), dimension(1), save :: lhypol
!    ! lvdw: LJ parameters
    logical(lp), save :: lvdw = .false.

    ! charges, areas, coordinates, elements and exclusion lists
    ! site nuclear charges
    real(rp), dimension(:,:), allocatable, save :: Zs
    ! fragment density nuclear charges
    real(rp), dimension(:,:), allocatable, save :: Zfd
    ! core fragment nuclear charges
    real(rp), dimension(:,:), allocatable, save :: Zm
    ! surface point areas
    real(rp), dimension(:), allocatable, save :: Sa
    ! Atom ID
    integer(ip), dimension(:), save, allocatable :: idatm
    ! Radii of the sphere from which each tessera stems from
    real(rp), dimension(:), allocatable, save :: Radsp
    ! site coordinates
    real(rp), dimension(:,:), allocatable, save :: Rs
    ! fragment density nuclear coordinates
    real(rp), dimension(:,:), allocatable, save :: Rfd
    ! core fragment nuclear coordinates
    real(rp), dimension(:,:), allocatable, save :: Rm
    ! all nuclear coordinates
    !real(rp), dimension(:,:), allocatable, save :: Rm_full
    real(rp), dimension(:,:), allocatable, save :: all_coords
    ! surface point coordinates
    real(rp), dimension(:,:), allocatable, save :: Rsp
    ! The normal vector to each surface point
    real(rp), dimension(:,:), allocatable, save :: nrmsp
    ! core center-of-mass
    real(rp), dimension(:), allocatable, save :: core_com
    ! core when changing gauge
    real(rp), dimension(:), allocatable, save :: gauge_input
    ! core polarizabilities
    real(rp), dimension(:,:), allocatable, save :: core_alphas
    ! site elements
    character(len=2), dimension(:), allocatable, save :: elems
    ! multipole exclusion list
    integer(ip), dimension(:,:), allocatable, save :: mul_exclists
    ! polarizability exclusion list
    integer(ip), dimension(:,:), allocatable, save :: pol_exclists

    ! energy contributions
    ! total
    real(rp), dimension(:), allocatable, save :: Epe
    ! electrostatic
    real(rp), dimension(:,:), allocatable, save :: Ees
    ! polarization
    real(rp), dimension(:,:), allocatable, save :: Epol
    ! continuum solvation
    real(rp), dimension(:,:), allocatable, save :: Esol
    ! fragment density
    real(rp), dimension(:,:), allocatable, save :: Efd
    ! LJ
    real(rp), dimension(:), allocatable, save :: Elj

    ! multipole moments
    ! monopoles, dipoles, quadrupoles, octopoles, etc.
    real(rp), dimension(:,:), allocatable, save :: M0s, M1s, M2s, M3s, M4s, M5s
    ! polarizabilities
    ! symmetric dipole-dipole polarizabilities
    real(rp), dimension(:,:), allocatable, save :: P11s
    ! .true. if P11 > 0 else .false.
    logical(lp), dimension(:), allocatable, save :: zeroalphas

    ! LJ parameters for MM region - from pot file
    real(rp), dimension(:,:), allocatable, save :: LJs
    ! LJ paramters for QM region - from dal file
    real(rp), dimension(:,:), allocatable, save :: qmLJs

    ! CUBE stuff
    ! ---------

    ! options for MEP
    ! create QM cubes
    logical(lp), save :: mep_qmcube = .true.
    ! create multipole cubes
    logical(lp), save :: mep_mulcube = .true.
    ! external electric field
    logical(lp), save :: mep_extfld = .false.
    real(rp), dimension(3), save :: extfld = huge(0.0)
    ! calculate electric field
    logical(lp), save :: mep_field = .false.
    ! lf component
    logical(lp), save :: mep_lf = .false.
    integer(ip), save :: lf_component = huge(0)
    ! grid file
    logical(lp), save :: mep_input = .false.

    ! options for CUBE
    ! calculate electric field
    logical(lp), save :: cube_field = .false.

    ! general cube information
    ! number of grid points
    integer(ip), save :: npoints
    ! grid points
    real(rp), dimension(:,:), allocatable, save :: Rg
    ! CUBE file origin and step sizes
    real(rp), dimension(3), save :: origin, step
    ! grid density in x, y and z direction
    integer(ip), save :: xgrid = 6
    integer(ip), save :: ygrid = 6
    integer(ip), save :: zgrid = 6
    ! number of steps in x, y and z direction
    integer(ip), save :: xsteps
    integer(ip), save :: ysteps
    integer(ip), save :: zsteps
    ! box size relative to molecule size
    real(rp), save :: xsize = 8.0
    real(rp), save :: ysize = 8.0
    real(rp), save :: zsize = 8.0

    ! Internal field stuff and locfld stuff
    ! Coordinates on which potential and field are calculated
    real(rp), dimension(:,:), allocatable, save :: crds
    ! Number of coordinates (length of crds) / 3
    integer(ip), save :: ncrds = 0

    ! External field stuff
    ! Field strengths
    real(rp), dimension(3), save :: efields = 0.0

    ! QM/MM interaction scheme options
    ! Scheme to use
    character(len=16), public, save :: qmmm_interaction_scheme = 'EXACT'
                                                              ! 'FMM_FOCK'
                                                              ! 'SINGLE_CENTER'
                                                              ! 'ESPF'
                                                              ! 'CESPF'
    ! For ESPF/CESPF, which solver to use for the least-squares fitting
    character(len=16), public, save :: ESPF_solver = 'SPTRS'  ! or SVD
    ! tolerance / threshold under which singular values are discarded for the SVD solver
    real(rp), public, save :: svdtol = 1e-9

    ! For ESPF/CESPF, which grid to use for building the least-squares matrix
    ! SOLVENT uses as grid positions the coordinates of the MM atoms (from the .pot file)
    ! VDW creates a VDW grid and uses the surface points as grid positions
    character(len=16), public, save :: ESPF_grid_type = 'SOLVENT' 
    ! For SOLVENT grid, cutoff within which MM positions will be used as grid points
    real(rp), public, save :: qmmm_interaction_rfit = 20.0 ! bohr

    ! For ESPF/CESPF, add bond midpoints as expansion points
    logical(lp), public, save :: ESPF_bond_midpoints = .false.

    ! Save near-field (within Rcut or theta) in memory
    logical(lp), public, save :: incore_nearfield = .false.

    ! Cutoff beyond which to use the approximate interaction (ESPF/CESPF/SINGLE_CENTER)
    real(rp), public, save :: qmmm_interaction_rcut = 20.0 ! bohr

    ! theta (same meaning as in standard FMM) for the FMM-like FMM_FOCK interaction type 
    real(rp), save :: qmmm_interaction_theta = 0.4

    ! expansion order for the approximate QM/MM interaction,
    ! max order of multipoles for ESPF/CESPF/SINGLE_CENTER
    ! max order of expansion/multipoles for FMM_FOCK
    integer(ip), public, save :: qmmm_interaction_order = 0   

    ! allocatables for approximate QM/MM coupling
    ! integral storage for FMM_FOCK or near-field t1 integrals
    real(rp), dimension(:,:,:,:), allocatable, save :: Fel_ints_fmm, near_field_integrals
    ! coordinates of near/far field MM positions
    real(rp), allocatable, public, save, dimension(:,:) :: mm_coordinates, mm_coordinates_far
    real(rp), allocatable, target, public, save, dimension(:,:) :: mm_coordinates_near
    ! multipole interaction operators (traced with density to give QM multipoles) for ESPF/CESPF/SINGLE_CENTER
    real(rp), allocatable, public, save, dimension(:,:,:) :: multipole_interaction_operators  ! nnbas x num_operators x multipole_length
    ! positions of multipole interaction operators
    real(rp), allocatable, public, save, dimension(:,:) :: multipole_operator_centers  ! 3 x num_operators 
    ! index array for near/far field MM positions
    integer(ip), allocatable, public, save, dimension(:) :: idx_near_field, idx_far_field
end module pelib_options
