module pelib_interaction_tools
    use pelib_precision
    implicit none
    private

    public build_multipole_interaction_operators
    public get_near_field_idx
    public get_multipole_field
contains

pure function xyz2idx(x, y, z) result(idx)
    implicit none
    integer, intent(in) :: x, y, z
    integer :: k
    integer :: idx
    k = x + y + z
    ! number of components before order k is k*(k+1)*(k+2)/6
    ! (y**2 + 2*y*z + y + z**2 + 3*z)/2 + 1 is the symmetry-packed index of the current slice
    idx = k * (k + 1) * (k + 2) / 6 + (y**2 + 2 * y * z + y + z**2 + 3 * z) / 2 + 1
end function xyz2idx

pure function factorial(n)
    integer, intent(in) :: n
    integer :: i
    real(8) :: factorial
    factorial = 1
    do i = n, 1, -1
        factorial = factorial * i
    end do
end function

pure function binom(n, k)
    integer, intent(in) :: n, k
    real(8) :: binom
    binom = factorial(n) / (factorial(k) * factorial(n - k))
end function binom

pure function trinom(i, j, k)
    integer, intent(in) :: i, j, k
    real(8) :: trinom
    trinom = factorial(i + j + k) / (factorial(i) * factorial(j) * factorial(k))
end function trinom

subroutine solve_svd(A, b, svdtol)
    external dgesvd
    real(rp), intent(in), dimension(:,:) :: A
    real(rp), intent(inout), dimension(:,:) :: b
    real(rp), intent(in), optional :: svdtol
    integer(ip) :: N, M, lda, ldb, ldu, ldvt, lwork, info, incx, nrhs, k, i
    real(rp), allocatable, dimension(:,:) :: U, VT
    real(rp), allocatable, dimension(:) :: S, work
    real(rp), allocatable, dimension(:,:) :: x
    real(rp) :: alpha, beta
    real(rp) :: tol

    if (present(svdtol)) then
        tol = svdtol
    else
        tol = 1e-9
    end if

    if (size(A,1) /= size(A,2)) error stop "Not implemented"

    N = size(A, 1)
    M = N
    K = N
    lda = N
    ldb = N
    ldu = N
    ldvt = N
    nrhs = size(b, 2)

    allocate(U(N,N))
    allocate(VT(N,N))
    allocate(S(N))
    allocate(work(1))
    lwork = -1

    ! query size
    call dgesvd('S', 'S', M, N, A, lda, S, U, ldu, VT, ldvt, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    ! factorize
    call dgesvd('S', 'S', M, N, A, lda, S, U, ldu, VT, ldvt, work, lwork, info)

    ! A = U @ diag(S) @ VT
    ! Ainv = VT.T @ diag(1/S) @ U.T
    do i = 1, N
        if (S(i) > tol) then
            U(:, i) = U(:, i) / S(i)
        else
            U(:, i) = 0.0_rp
        end if
    end do

    alpha = 1.0_rp
    beta = 0.0_rp
    call dgemm("T", "T", M, N, K, alpha, VT, ldvt, U, ldu, beta, A, lda)

    ! apply A^{+} b => x
    allocate(x(size(b,1), size(b,2)))
    alpha = 1.0_rp
    beta = 0.0_rp
    call dgemm("N", "N", M, nrhs, K, alpha, A, lda, b, ldb, beta, x, ldb)
    b = x
    deallocate(U, VT, S, work, x)
end subroutine solve_svd

subroutine vdw_grid(charges, coordinates, density, r_min, delta_r, num_surfaces, grid)
    use pelib_blas_interfaces, only: nrm2
    integer(ip), dimension(:), intent(in) :: charges
    real(rp), dimension(:,:), intent(in) :: coordinates
    real(rp), intent(in) :: density, r_min, delta_r
    integer(ip), intent(in) :: num_surfaces
    real(rp), dimension(:,:), allocatable, intent(out) :: grid
    real(rp), dimension(:,:), allocatable :: tmpgrid
    integer(ip) :: N, i, j, k, l, m, num_points
    real(rp) :: r, h, theta, phi, r_scale, distance, x, y, z, grid_spacing
    logical :: near_other_atom, near_other_point
    real(rp), parameter :: pi = 4.0_rp*atan2(1.0_rp,1.0_rp)
    real(rp), parameter, dimension(118) :: vdw_radii =             &
      [2.26767119, 2.64561638, 3.4393013 , 2.89128076, 3.6282739 , &
       3.21253418, 2.92907528, 2.8723835 , 2.7778972 , 2.91017802, &
       4.28967799, 3.26922596, 3.47709582, 3.96842458, 3.40150678, &
       3.40150678, 3.30702048, 3.55268486, 5.19674647, 4.36526703, &
       3.98732184, 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 3.08025336, 2.64561638, 2.62671912, &
       3.5337876 , 3.98732184, 3.49599308, 3.59047938, 3.49599308, &
       3.8172465 , 5.72586975, 4.70541771, 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       3.08025336, 3.2503287 , 2.98576706, 3.64717116, 4.1007054 , &
       3.89283554, 3.89283554, 3.74165746, 4.08180814, 6.48176014, &
       5.06446565, 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 3.30702048, 3.13694514, 2.92907528, &
       3.70386294, 3.8172465 , 3.9117328 , 3.7227602 , 3.8172465 , &
       4.15739717, 6.57624644, 5.34792455, 4.        , 4.        , &
       4.        , 3.51489034, 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        , 4.        , 4.        , &
       4.        , 4.        , 4.        ]  ! missing radii set to 4 bohr 


    N = size(coordinates, 2)

    r_scale = r_min
    do j = 1, num_surfaces
        do i = 1, N
            num_points = num_points + int(density * 4 * pi * (r_scale * vdw_radii(i))**2)
        end do
        r_scale = r_scale + delta_r
    end do
    r_scale = r_min

    allocate(tmpgrid(3, num_points))
    l = 1
    do j = 1, num_surfaces
        do i = 1, N
            num_points = nint(density * 4 * pi * (r_scale * vdw_radii(charges(i)))**2)
            grid_spacing = 0.5_rp * 3.809 / sqrt(real(num_points)) * r_scale * vdw_radii(charges(i))
            do k = 0, num_points - 1
                h = 2.0_rp * k / (num_points) - 1.0_rp
                theta = acos(h)
                if (k == 0 .or. k == num_points - 1) then 
                    phi = 0.0_rp
                else
                    x = num_points * (1-h**2)
                    phi = phi + 3.6_rp / sqrt(x)
                    phi = mod(phi, 2*pi)
                end if
                r = r_scale * vdw_radii(charges(i))
                x = r * cos(phi) * sin(theta)
                y = r * sin(phi) * sin(theta)
                z = r * cos(theta)
                tmpgrid(1, l) = coordinates(1, i) + x
                tmpgrid(2, l) = coordinates(2, i) + y
                tmpgrid(3, l) = coordinates(3, i) + z
                ! keep point only if not near other atom
                near_other_atom = .false.
                do m = 1, N
                    distance = nrm2(tmpgrid(:, l) - coordinates(:, m))
                    if (distance < r_scale * 0.99_rp * vdw_radii(charges(m))) then
                        near_other_atom = .true.
                        exit
                    end if
                end do
                if (.not. near_other_atom) then
                    near_other_point = .false.
                    do m = 1, l - 1
                        distance = nrm2(tmpgrid(:, l) - tmpgrid(:, m))
                        if (distance < grid_spacing) then
                            near_other_point = .true.
                            exit
                        end if
                    end do
                    if (.not. near_other_point) then
                        l = l + 1
                    end if
                end if
            end do
        end do
        r_scale = r_scale + delta_r
    end do
    l = l - 1
    allocate(grid(3, l))
    grid(:, 1:l) = tmpgrid(:, 1:l)
    deallocate(tmpgrid)
end subroutine vdw_grid


subroutine build_multipole_interaction_operators
    use pelib_blas_interfaces, only: nrm2
    use pelib_mpi
    use pelib_options
    use pelib_lapack_interfaces, only: sptrf, sptrs, sptri, spcon
    use T_tensor, only: Tn
    use pelib_integral_interfaces, only: tk_integrals, multipole_integrals
    use pelib_multipole_tools, only: prefactors
    integer(ip) :: l, k, ncomp, multipole_size, matrix_dimension
    integer(ip) :: N, M
    integer(ip) :: i, j, idx, info, p, mu, nu
    real(rp), allocatable, dimension(:, :, :) :: T
    real(rp), allocatable, dimension(:, :) :: A, T_slice, weights
    real(rp), allocatable, dimension(:) :: Ap, dx, dy, dz
    real(rp), dimension(:,:,:), allocatable :: Fel_ints, integrals
    real(rp), dimension(:,:), allocatable :: untransformed_operators
    real(rp), dimension(3) :: multipole_origin
    integer(ip), dimension(:), allocatable :: ipiv
    integer(ip), dimension(:), allocatable :: near_field_indices, far_field_indices
    real(rp), dimension(:), allocatable :: factors
    real(rp) :: factor
    real(rp) :: x, y, z, monomial, distance, bond_threshold
    integer(ip) :: kx, ky, kz, lx, ly, lz, lidx, kidx
    integer(ip) :: work_start, work_stop, work_size
    real(rp), dimension(:,:), allocatable, target :: grid
    real(rp), dimension(:,:), pointer :: espf_grid
    integer(ip), dimension(:), allocatable :: Zi


    if (.not. allocated(mm_coordinates)) then 
        allocate(mm_coordinates(3, npols))
        k = 1
        do i = 1, nsites
            if (zeroalphas(i)) cycle
            mm_coordinates(:, k) = Rs(:, i)
            k = k + 1
        end do
    end if
    if (.not. allocated(idx_near_field)) then
        call get_near_field_idx(Rm, mm_coordinates, qmmm_interaction_rcut, idx_near_field, idx_far_field)
        allocate(mm_coordinates_near(3, size(idx_near_field)))
        allocate(mm_coordinates_far(3, size(idx_far_field)))
        do i = 1, size(idx_near_field)
            mm_coordinates_near(:, i) = mm_coordinates(:, idx_near_field(i))
        end do
        do i = 1, size(idx_far_field)
            mm_coordinates_far(:, i) = mm_coordinates(:, idx_far_field(i))
        end do
    end if

    multipole_size = (qmmm_interaction_order+1)*(qmmm_interaction_order+2)*(qmmm_interaction_order+3)/6

    if (trim(qmmm_interaction_scheme) == 'SINGLE_CENTER') then
        allocate(multipole_interaction_operators(nnbas, multipole_size, 1))
        allocate(multipole_operator_centers(3, 1))

        multipole_origin = core_com
        multipole_operator_centers(:, 1) = multipole_origin
        l = 1
        do k = 0, qmmm_interaction_order
            ncomp = (k+1)*(k+2) / 2
            call multipole_integrals(multipole_origin, multipole_interaction_operators(:, l:l+ncomp-1, 1:1))
            l = l + ncomp
        end do
    else if (trim(qmmm_interaction_scheme) == 'ESPF' .or. trim(qmmm_interaction_scheme) == 'CESPF') then
        ! build Aij
        !  1   N N+1    
        ! ( Aqq  Aqmu  ... )
        ! ( Amuq Amumu ... )
        ! ( ...  ...   ... )
        k = 0
        bond_threshold = 3.2 ! bohr ~ 1.7 Å
        if (ESPF_bond_midpoints) then
            do i = 1, size(Rm, 2)
                do j = i + 1, size(Rm, 2)
                    distance = nrm2(Rm(:, i) - Rm(:,j))
                    if (distance < bond_threshold) then
                        k = k + 1
                    end if
                end do
            end do
        end if

        N = size(Rm, 2) + k
        allocate(multipole_operator_centers(3, N))
        do i = 1, size(Rm, 2)
            multipole_operator_centers(:, i) = Rm(:, i)
        end do
        if (ESPF_bond_midpoints) then
            k = size(Rm, 2) + 1
            do i = 1, size(Rm, 2)
                do j = i + 1, size(Rm, 2)
                    distance = nrm2(Rm(:, i) - Rm(:,j))
                    if (distance < bond_threshold) then
                        multipole_operator_centers(:, k) = 0.5_rp * (Rm(:, i) + Rm(:, j))
                        k = k + 1
                    end if
                end do
            end do
        end if

        if (trim(ESPF_grid_type) == 'SOLVENT') then
            if (qmmm_interaction_rcut == qmmm_interaction_rfit) then
                espf_grid => mm_coordinates_near
            else
                call get_near_field_idx(Rm, mm_coordinates, qmmm_interaction_rfit, near_field_indices, far_field_indices)
                allocate(grid(3, size(near_field_indices)))
                do i = 1, size(near_field_indices)
                    grid(:, i) = mm_coordinates(:, near_field_indices(i))
                end do
                espf_grid => grid
            end if
        else if (trim(ESPF_grid_type) == 'VDW') then
            allocate(Zi(size(Zm, 2)))
            do i = 1, size(Zm, 2)
                Zi(i) = nint(Zm(1, i))
            end do
            !                    density  rmin    delta_r num_surfaces
            call vdw_grid(Zi, Rm, 1.0_rp, 2.0_rp, 1.0_rp, 5, grid)
            espf_grid => grid
        else
            error stop "Invalid ESPF grid type"
        end if
        M = size(espf_grid, 2)

        if (trim(qmmm_interaction_scheme) == 'CESPF') then
            matrix_dimension = (N+1)*multipole_size
        else
            matrix_dimension = (N)*multipole_size
        end if

        allocate(multipole_interaction_operators(nnbas, multipole_size, N))
        allocate(A(matrix_dimension, matrix_dimension))
        allocate(T(N, M, multipole_size)) 
        allocate(T_slice(N, multipole_size)) 
        allocate(dx(N))
        allocate(dy(N))
        allocate(dz(N))

        do i = 1, M
            dx = espf_grid(1, i) - multipole_operator_centers(1, :) 
            dy = espf_grid(2, i) - multipole_operator_centers(2, :) 
            dz = espf_grid(3, i) - multipole_operator_centers(3, :) 
            call Tn(qmmm_interaction_order, dx, dy, dz, T_slice)
            T(:, i, :) = T_slice
        end do

        allocate(factors(multipole_size))
        l = 1
        do k = 0, qmmm_interaction_order
            ncomp = (k+1)*(k+2)/2
            call prefactors(factors(l:l+ncomp-1))
            l = l + ncomp
        end do

        do k = 1, multipole_size
            T(:,:,k) = T(:,:,k) * factors(k)
        end do

        ! esp blocks
        A = 0.0_rp
        do k = 1, multipole_size
            do l = 1, multipole_size
                A(N*(k-1)+1:N*k, N*(l-1)+1:N*l) = matmul(T(:,:,k), transpose(T(:,:,l)))
            end do
        end do

        !constraint blocks
        if (trim(qmmm_interaction_scheme) == 'CESPF') then
            do k = 0, qmmm_interaction_order
                do l = 0, k
                    do kx = k, 0, -1
                        do ky = k, 0, -1
                            do kz = k, 0, -1
                                if (kx + ky + kz /= k) cycle
                                kidx = xyz2idx(kx, ky, kz)
                                do lx = kx, 0, -1
                                    do ly = ky, 0, -1
                                        do lz = kz, 0, -1
                                            if (lx + ly + lz /= l) cycle
                                            lidx = xyz2idx(lx, ly, lz)
                                            do i = 1, N
                                                x = multipole_operator_centers(1, i) - core_com(1)
                                                y = multipole_operator_centers(2, i) - core_com(2)
                                                z = multipole_operator_centers(3, i) - core_com(3)
                                                monomial = x**(kx - lx) * y**(ky - ly) * z**(kz - lz)
                                                A(N*multipole_size + kidx, (lidx-1)*N + i) = monomial
                                                A((lidx-1)*N + i, N*multipole_size + kidx) = monomial
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end if

        ! pack A
        ! and factorize
        if (trim(ESPF_solver) == 'SPTRS') then
            allocate(Ap((matrix_dimension)*(matrix_dimension+1)/2))
            k = 1
            do i = 1, matrix_dimension
                do j = i, matrix_dimension
                    Ap(k) = A(i, j)
                    k = k + 1
                end do
            end do
            allocate(ipiv(matrix_dimension))
            call sptrf(Ap, 'L', ipiv, info)
        end if

        ! get untransformed opererators
        ! esp part
        allocate(Fel_ints(nnbas, multipole_size, 1))
        allocate(untransformed_operators(matrix_dimension, nnbas))
        untransformed_operators = 0.0_rp
        ! parallelized over M grid points (p)
        work_size = M / nprocs
        work_start = 1 + myid * work_size
        work_stop = (myid + 1) * work_size
        if (myid == nprocs - 1) work_stop = M
        do p = work_start, work_stop
            call Tk_integrals('potential_derivative', espf_grid(:, p), Fel_ints(:, 1:1, 1:1))
            k = 1
            do j = 1, multipole_size
                do i = 1, N
                    untransformed_operators(k, :) = untransformed_operators(k, :) - T(i,p,j) * Fel_ints(:, 1, 1)
                    k = k + 1
                end do
            end do
        end do
        deallocate(Fel_ints)
#ifdef VAR_MPI
        call mpi_allreduce(MPI_IN_PLACE, untransformed_operators(1,1), size(untransformed_operators), rmpi, MPI_SUM, comm, ierr)
#endif
        ! constraint part - multipole integrals
        if (trim(qmmm_interaction_scheme) == 'CESPF') then
            multipole_origin = core_com
            l = 1
            allocate(integrals(nnbas, multipole_size, 1))
            do k = 0, qmmm_interaction_order
                ncomp = (k+1)*(k+2) / 2
                call multipole_integrals(multipole_origin, integrals(:, l:l+ncomp-1, 1:1))
                l = l + ncomp
            end do
            do j = 1, multipole_size
                untransformed_operators(N*multipole_size + j, :) = integrals(:, j, 1)
            end do
            deallocate(integrals)
        end if
        ! first N*multipole size operators are the multipole operators
        ! the rest is lagrange multiplier operators

        ! parallelized over nnbas
        work_size = nnbas / nprocs
        work_start = 1 + myid * work_size
        work_stop = (myid + 1) * work_size
        if (myid == nprocs - 1) work_stop = nnbas

        if (trim(ESPF_solver) == 'SVD') then
            call solve_svd(A, untransformed_operators(:, work_start:work_stop))
        else if (trim(ESPF_solver) == 'SPTRS') then
            call sptrs(Ap, untransformed_operators(:, work_start:work_stop), ipiv, 'L', info)
            if (info /= 0) error stop "Failed to solve ESP linear system"
            deallocate(ipiv)
        else
            error stop "Invalid ESPF solver."
        end if

        multipole_interaction_operators = 0.0_rp
        do j = 1, multipole_size
            do i = 1, N
                multipole_interaction_operators(work_start:work_stop, j, i) = untransformed_operators((j-1)*N+i, work_start:work_stop)
            end do
        end do
#ifdef VAR_MPI
        call mpi_allreduce(MPI_IN_PLACE, multipole_interaction_operators(1, 1, 1), size(multipole_interaction_operators), rmpi, MPI_SUM, comm, ierr)
#endif

        deallocate(untransformed_operators)
    else
        print *, qmmm_interaction_scheme, "is not available."
        error stop "Not Implemented"
    end if
end subroutine build_multipole_interaction_operators

subroutine get_near_field_idx(target_coordinates, source_coordinates, Rcut, near_field_indices, far_field_indices)
    use pelib_blas_interfaces, only: nrm2
    use pelib_options

    real(rp), dimension(:,:), intent(in) :: target_coordinates, source_coordinates
    real(rp), intent(in) :: Rcut
    integer(ip), allocatable, dimension(:), intent(out) :: near_field_indices, far_field_indices
    logical, allocatable, dimension(:) :: is_near
    integer(ip) :: num_source, num_target, num_far_idx, num_near_idx
    integer(ip) :: i, j, idx_near, idx_far
    real(rp), dimension(3) :: r

    num_source = size(source_coordinates, 2)
    num_target = size(target_coordinates, 2)
    allocate(is_near(num_source))
    is_near(:) = .false.
    num_near_idx = 0
    do i = 1, num_source
        do j = 1, num_target
            r = source_coordinates(:,i) - target_coordinates(:,j)
            if (nrm2(r) < Rcut) then
                is_near(i) = .true.
                num_near_idx = num_near_idx + 1
                exit
            end if
        end do
    end do
    num_far_idx = num_source - num_near_idx
    allocate(near_field_indices(num_near_idx))
    allocate(far_field_indices(num_far_idx))
    idx_near = 1
    idx_far = 1
    do i = 1, num_source
        if (is_near(i)) then
            near_field_indices(idx_near) = i
            idx_near = idx_near + 1
        else
            far_field_indices(idx_far) = i
            idx_far = idx_far + 1
        end if
    end do
    deallocate(is_near)
end subroutine get_near_field_idx

subroutine get_multipole_field(field, source_coordinates, target_coordinates, multipoles, max_field_order)
    use pelib_mpi
    use T_tensor, only: Tn
    real(rp), intent(out) :: field(:, :, :)
    real(rp), dimension(:,:), intent(in) :: source_coordinates, target_coordinates
    real(rp), dimension(:,:,:), intent(in) :: multipoles
    integer(ip), intent(in) :: max_field_order

    integer(ip) :: N, M, ndens
    integer(ip) :: i, j, k
    integer(ip) :: max_multipole_order, max_tensor_order, field_length, tensor_length
    real(rp), dimension(:), allocatable :: dx, dy, dz
    real(rp), dimension(:,:), allocatable :: T
    integer(ip) :: work_start, work_stop, work_size, iwork
    integer(ip) :: field_order, multipole_order, tx, ty, tz, sx, sy, sz, field_index, tensor_index, multipole_index
    real(rp) :: taylor, symfac
    N = size(source_coordinates, 2)
    M = size(target_coordinates, 2)
    ndens = size(multipoles, 3)

    ! parallelization over N*M
    ! this hurts vectorization
    ! maybe make two variants of get_multipole_field to go only over either N or M


    work_size = N*M / nprocs
    work_start = 1 + myid * work_size
    work_stop = (myid + 1) * work_size
    if (myid == nprocs - 1) work_stop = N*M

    allocate(dx(1))
    allocate(dy(1))
    allocate(dz(1))

    max_multipole_order = NINT((size(multipoles, 2) * 6)**(1./3.)) - 2
    max_tensor_order = max_field_order + max_multipole_order
    field_length = (max_field_order + 1) * (max_field_order + 2) * (max_field_order + 3) / 6
    tensor_length = (max_tensor_order + 1) * (max_tensor_order + 2) * (max_tensor_order + 3) / 6

    allocate(T(1, tensor_length))

    field = 0.0_rp
    do iwork = work_start, work_stop
        i = 1 + ((iwork-1) / M)
        j = 1 + mod(iwork - 1, M)
        dx = target_coordinates(1, j) - source_coordinates(1, i)
        dy = target_coordinates(2, j) - source_coordinates(2, i)
        dz = target_coordinates(3, j) - source_coordinates(3, i)
        call Tn(max_tensor_order, dx, dy, dz, T)
        do field_order = 0, max_field_order
            do multipole_order = 0, max_multipole_order
                if (field_order + multipole_order > max_tensor_order) cycle
                taylor = (-1.0)**(multipole_order + 1) / factorial(multipole_order)
                do sx = multipole_order, 0, -1
                    do sy = multipole_order, 0, -1
                        do sz = multipole_order, 0, -1
                            if (sx + sy + sz /= multipole_order) cycle
                            multipole_index = xyz2idx(sx, sy, sz)
                            do tx = field_order, 0, -1
                                do ty = field_order, 0, -1
                                    do tz = field_order, 0, -1
                                        if (tx + ty + tz /= field_order) cycle
                                        field_index = xyz2idx(tx, ty, tz)
                                        tensor_index = xyz2idx(sx + tx, sy + ty, sz + tz)
                                        symfac = trinom(sx, sy, sz)
                                        do k = 1, ndens
                                            field(j, field_index, k) = field(j, field_index, k) + taylor * symfac * multipoles(i, multipole_index, k) * T(1, tensor_index)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

#ifdef VAR_MPI
    if (myid == 0) then
        call mpi_reduce(MPI_IN_PLACE, field(1, 1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
    else
        call mpi_reduce(field(1, 1, 1), field(1, 1, 1), size(field), MPI_REAL8, MPI_SUM, 0, comm, ierr)
    end if
#endif
end subroutine get_multipole_field

end module pelib_interaction_tools
