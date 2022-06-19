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
module pelib_multipole_tools

    use pelib_precision
    use pelib_options, only: pelib_amoeba_damp

    implicit none

    private

    public :: multipole_interaction_energy, multipole_derivative, multipole_derivative_damped
    public :: T, Td, Tk_tensor, Tk_damped_tensor, prefactors, damping_coefficient

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476
    ! C^(n)_ij coefficients for calculating T(k) tensor elements
    integer(ip), dimension(:,:,:), allocatable, save :: Cnij

contains

!------------------------------------------------------------------------------

subroutine multipole_interaction_energy(E, Mi, Rij, Mj)

    real(rp), intent(out) :: E
    real(rp), dimension(:), intent(in) :: Mi
    real(rp), dimension(:), intent(in) :: Mj
    real(rp), dimension(3), intent(in) :: Rij

    integer(ip) :: a_order
    integer(ip) :: b_order
    integer(ip) :: a_x, a_y, a_z
    integer(ip) :: b_x, b_y, b_z
    integer(ip) :: a_norm, a_fact, a_idx
    integer(ip) :: b_norm, b_fact, b_idx
    real(rp) :: taylor

    E = 0.0

    ! the order of the first multipole
    a_order = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Mi), rp)) - 1.0)) - 1

    ! the order of the second multipole
    b_order = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Mj), rp)) - 1.0)) - 1

    do a_x = a_order, 0, -1
        do a_y = a_order, 0, -1
            do a_z = a_order, 0, -1
                if (a_x + a_y + a_z /= a_order) cycle
                a_norm = a_x + a_y + a_z
                a_fact = factorial(a_x) * factorial(a_y) * factorial(a_z)
                a_idx = xyz2idx(a_x, a_y, a_z)
                do b_x = b_order, 0, -1
                    do b_y = b_order, 0, -1
                        do b_z = b_order, 0, -1
                            if (b_x + b_y + b_z /= b_order) cycle
                            b_norm = b_x + b_y + b_z
                            b_fact = factorial(b_x) * factorial(b_y) * factorial(b_z)
                            b_idx = xyz2idx(b_x, b_y, b_z)
                            taylor = (-1)**(a_norm+b_norm) / (a_fact * b_fact)
                            E = E + taylor * Mi(a_idx) * T(Rij, a_x+b_x, a_y+b_y, a_z+b_z) * Mj(b_idx)
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine multipole_interaction_energy

subroutine multipole_derivative(Fi, Rji, Mkj)

    real(rp), dimension(:), intent(inout) :: Fi
    real(rp), dimension(3), intent(in) :: Rji
    real(rp), dimension(:), intent(in) :: Mkj

    integer(ip) :: i, j, k, l, m
    integer(ip) :: a, b, c, x, y, z
    real(rp) :: taylor, symfac
    real(rp), dimension(:), allocatable :: Tk

    ! the order of the incoming multipole moment
    k = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Mkj), rp)) - 1.0)) - 1

    ! the order of the derivative
    l = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Fi), rp)) - 1.0)) - 1

    if (mod(k + l, 2) == 0) then
        taylor = 1.0 / real(factorial(k), rp)
    else if (mod(k + l, 2) /= 0) then
        taylor = - 1.0 / real(factorial(k), rp)
    end if

    allocate(Tk((k + l + 1) * (k + l + 2) / 2))
    call Tk_tensor(Tk, Rji)

    do x = k + l, 0, - 1
        do y = k + l, 0, - 1
            do z = k + l, 0, - 1
                if (x + y + z /= k + l) cycle
                i = xyz2idx(x, y, z)
                do a = x, 0, -1
                    do b = y, 0, -1
                        do c = z, 0, -1
                            if (a + b + c /= k) cycle
                            j = xyz2idx(a, b, c)
                            m = xyz2idx(x-a, y-b, z-c)
                            symfac = real(trinom(a, b, c), rp)
                            Fi(m) = Fi(m) + taylor * symfac * Tk(i) * Mkj(j)
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine multipole_derivative

!------------------------------------------------------------------------------

subroutine multipole_derivative_damped(Fi, Rji, Mkj, alpha_i, alpha_j, factor)

    real(rp), dimension(3), intent(inout) :: Fi
    real(rp), dimension(3), intent(in) :: Rji
    real(rp), dimension(:), intent(in) :: Mkj
    real(rp), intent(in) :: alpha_i, alpha_j, factor

    integer(ip) :: i, j, k, l, m
    integer(ip) :: a, b, c, x, y, z
    real(rp) :: taylor, symfac
    real(rp), dimension(:), allocatable :: Tk

    ! the order of the incoming multipole moment
    k = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Mkj), rp)) - 1.0)) - 1

    ! the order of the derivative
    l = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Fi), rp)) - 1.0)) - 1

    if (mod(k + l, 2) == 0) then
        taylor = 1.0 / real(factorial(k), rp)
    else if (mod(k + l, 2) /= 0) then
        taylor = - 1.0 / real(factorial(k), rp)
    end if

    allocate(Tk((k + l + 1) * (k + l + 2) / 2))
    call Tk_damped_tensor(Tk, Rji, alpha_i, alpha_j, factor)

    do x = k + l, 0, - 1
        do y = k + l, 0, - 1
            do z = k + l, 0, - 1
                if (x + y + z /= k + l) cycle
                i = xyz2idx(x, y, z)
                do a = x, 0, -1
                    do b = y, 0, -1
                        do c = z, 0, -1
                            if (a + b + c /= k) cycle
                            j = xyz2idx(a, b, c)
                            m = xyz2idx(x-a, y-b, z-c)
                            symfac = real(trinom(a, b, c), rp)
                            Fi(m) = Fi(m) + taylor * symfac * Tk(i) * Mkj(j)
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine multipole_derivative_damped

!------------------------------------------------------------------------------

subroutine Tk_tensor(Tk, Rij)

    real(rp), dimension(:), intent(out) :: Tk
    real(rp), dimension(3), intent(in) :: Rij

    integer(ip) :: k, i
    integer(ip) :: x, y, z

    k = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Tk), rp)) - 1.0)) - 1

    i = 1
    do x = k, 0, - 1
        do y = k, 0, - 1
            do z = k, 0, - 1
                if (x + y + z /= k) cycle
                Tk(i) = T(Rij, x, y, z)
                i = i + 1
            end do
        end do
    end do

end subroutine Tk_tensor

!------------------------------------------------------------------------------

subroutine Tk_damped_tensor(Tk, Rij, alpha_i, alpha_j, factor)

    real(rp), dimension(:), intent(out) :: Tk
    real(rp), dimension(3), intent(in) :: Rij
    real(rp), intent(in) :: alpha_i, alpha_j, factor

    integer(ip) :: k, i
    integer(ip) :: x, y, z

    k = int(0.5 * (sqrt(1.0 + 8.0 * real(size(Tk), rp)) - 1.0)) - 1

    i = 1
    do x = k, 0, - 1
        do y = k, 0, - 1
            do z = k, 0, - 1
                if (x + y + z /= k) cycle
                Tk(i) = Td(Rij, x, y, z, alpha_i, alpha_j, factor)
                i = i + 1
            end do
        end do
    end do

end subroutine Tk_damped_tensor

!------------------------------------------------------------------------------

function T(Rij, x, y, z)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    use pelib_blas_interfaces, only: nrm2

    integer(ip), intent(in) :: x, y, z
    real(rp), dimension(3), intent(in) :: Rij

    integer(ip) :: l, m, n
    real(rp) :: T
    real(rp) :: R, Cx, Cy, Cz

    if (.not. allocated(Cnij)) call Tk_coefficients

    R = nrm2(Rij)

    T = 0.0

    do l = 0, x
        Cx = Cnij(1,x,l) * (Rij(1) / R)**l
        do m = 0, y
            Cy = Cx * Cnij(l+x+1,y,m) * (Rij(2) / R)**m
            do n = 0, z
                Cz = Cy * Cnij(l+x+m+y+1,z,n) * (Rij(3) / R)**n
                T = T + Cz
            end do
        end do
    end do

    T = T / R**(x + y + z + 1)

end function T

!------------------------------------------------------------------------------

function Td(Rij, x, y, z, alpha_i, alpha_j, factor)

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    use pelib_blas_interfaces, only: nrm2

    integer(ip), intent(in) :: x, y, z
    real(rp), intent(in) :: alpha_i, alpha_j, factor
    real(rp), dimension(3), intent(in) :: Rij

    integer(ip) :: l, m, n, k
    real(rp) :: Td
    real(rp) :: R, Cx, Cy, Cz
    real(rp) :: v, fV, fE, fT, fD
    real(rp) :: exp_factor

    k = x + y + z

    R = nrm2(Rij)

    if (.not. allocated(Cnij)) call Tk_coefficients

    if (pelib_amoeba_damp) then
        call damping_coefficient(Rij, alpha_i, alpha_j, factor, v)
        if (k == 0) then
            ! Missing the upper incomplete gamma function
            error stop 'damping does not work for zeroth-order interaction tensors'
        end if
        exp_factor = exp(-v**3)
        if (k >= 1) then
            fE = 1.0 - exp_factor
        end if
        if (k >= 2) then
            fT = 1.0 - (1.0 + v**3) * exp_factor
        end if
        if (k >= 3) then
            fD = 1.0 - (1.0 + v**3 + (3.0/5.0)*v**6) * exp_factor
        end if
        if (k > 3) then
            error stop 'damping only works upto third-order interaction tensors'
        end if
    else
        ! Thole damping
        ! JPC A 102 (1998) 2399 and Mol. Sim. 32 (2006) 471
        ! v = factor * u, where factor = 2.1304 (default) and u = R / (alpha_i * alpha_j)**(1/6)
        call damping_coefficient(Rij, alpha_i, alpha_j, factor, v)
        exp_factor = exp(-v)
        if (k >= 0) then
            fV = 1.0 - (0.5 * v + 1.0) * exp_factor
        end if
        if (k >= 1) then
            fE = fV  - (0.5 * v**2 + 0.5 * v) * exp_factor
        end if
        if (k >= 2) then
            fT = fE - (1.0 / 6.0) * v**3 * exp_factor
        end if
        if (k >= 3) then
            fD = fT - (1.0 / 30.0) * v**4 * exp_factor
        end if
        if (k > 3) then
            error stop 'damping only works upto third-order interaction tensors'
        end if
    end if

    Td = 0.0

    do l = 0, x
        Cx = Cnij(1,x,l) * (Rij(1) / R)**l
        do m = 0, y
            Cy = Cx * Cnij(l+x+1,y,m) * (Rij(2) / R)**m
            do n = 0, z
                Cz = Cy * Cnij(l+x+m+y+1,z,n) * (Rij(3) / R)**n
                if (l + m + n == 0) then
                    if (k == 0) then
                        Cz = Cz * fV
                    else if (k == 2) then
                        Cz = Cz * fE
                    end if
                else if (l + m + n == 1) then
                    if (k == 1) then
                        Cz = Cz * fE
                    else if (k == 3) then
                        Cz = Cz * fT
                    end if
                else if (l + m + n == 2) then
                    if (k == 2) then
                        Cz = Cz * fT
                    end if
                else if (l + m + n == 3) then
                    if (k == 3) then
                        Cz = Cz * fD
                    end if
                end if
                Td = Td + Cz
            end do
        end do
    end do

    Td = Td / R**(x + y + z + 1)

end function Td

!------------------------------------------------------------------------------

subroutine Tk_coefficients

    ! C. E. Dykstra, J. Comp. Chem., 9 (1988), 476

    integer(ip) :: i, j, k, n

! TODO
!    i = max(mulorder, polorder)
!    allocate(Cnij(2*i+3,0:i+1,0:i+1))
    allocate(Cnij(2*5+3,0:5+1,0:5+1))

    Cnij = 0
    Cnij(:,0,0) = 1
    do n = 1, 2*5+3
        if (mod(n,2) == 0) cycle
        do i = 1, 5+1
            if (mod(i,2) /= 0) then
                k = i - 1
            else if (mod(i,2) == 0) then
                k = i
            end if
            do j = 0, i
                if (mod(i+j,2) /= 0) cycle
                if (j == 0) then
                    Cnij(n,i,j) = Cnij(n,i-1,j+1)
                else if (j /= i) then
                    Cnij(n,i,j) = (j + 1) * Cnij(n,i-1,j+1)
                    Cnij(n,i,j) = Cnij(n,i,j) - (n + k) * Cnij(n,i-1,j-1)
                    k = k + 2
                else if (j == i) then
                    Cnij(n,i,j) = - (n + k) * Cnij(n,i-1,j-1)
                end if
            end do
        end do
    end do

end subroutine Tk_coefficients

!------------------------------------------------------------------------------

function xyz2idx(x, y, z) result(idx)

    integer(ip), intent(in) :: x, y, z

    integer(ip) :: k, a, b, c
    integer(ip) :: idx

    k = x + y + z
    idx = 1
    do a = k, 0, -1
        do b = k, 0, -1
            do c = k, 0, -1
                if (a + b + c /= k) cycle
                if (a /= x .or. b /= y .or. c /= z) then
                    idx = idx + 1
                else
                    return
                end if
            end do
        end do
    end do

end function xyz2idx

!------------------------------------------------------------------------------

subroutine prefactors(factors)

    real(rp), dimension(:), intent(out) :: factors

    integer(ip) :: k
    real(rp) :: taylor

    k = int(0.5 * (sqrt(1.0 + 8.0 * real(size(factors), rp)) - 1.0)) - 1

    if (mod(k,2) == 0) then
        taylor = 1.0 / real(factorial(k), rp)
    else if (mod(k,2) /= 0) then
        taylor = - 1.0 / real(factorial(k), rp)
    end if

    call symmetry_factors(factors)

    factors = taylor * factors

end subroutine prefactors

!------------------------------------------------------------------------------

subroutine symmetry_factors(factors)

    real(rp), dimension(:), intent(out) :: factors

    integer(ip) :: idx, x, y, z, k

    k = int(0.5 * (sqrt(1.0 + 8.0 * real(size(factors), rp)) - 1.0)) - 1

    idx = 1
    do x = k, 0, - 1
        do y = k, 0, - 1
            do z = k, 0, - 1
                if (x + y + z /= k) cycle
                factors(idx) = real(trinom(x, y, z), rp)
                idx = idx + 1
            end do
        end do
     end do

end subroutine symmetry_factors

!------------------------------------------------------------------------------

function trinom(i, j, k)

    ! trinomial coefficient

    integer(ip), intent(in) :: i, j, k

    integer(ip) :: trinom

    trinom = factorial(i+j+k) / (factorial(i) * factorial(j) * factorial(k))

end function trinom

!------------------------------------------------------------------------------

subroutine damping_coefficient(Rij, alpha_i, alpha_j, factor, coeff)

    use pelib_options
    use pelib_options
    use pelib_blas_interfaces, only: nrm2

    real(rp), dimension(3), intent(in) :: Rij
    real(rp), intent(in) :: alpha_i, alpha_j
    real(rp), intent(in) :: factor
    real(rp), intent(out) :: coeff

    real(rp), parameter :: i3 = 1.0 / 3.0
    real(rp), parameter :: i6 = 1.0 / 6.0

    ! Thole damping
    ! JPC A 102 (1998) 2399 and Mol. Sim. 32 (2006) 471
    ! coeff = a * u , where a = 2.1304 (default) and u = R / (alpha_i * alpha_j)**(1/6)
    !
    ! Thole-style AMOEBA damping added
    ! Ref?
    ! coeff = a^(1/3) * u, where a = 0.39 (default) and u same as above

    if ((alpha_i < zero) .or. (alpha_j < zero)) then
        error stop 'damping electric fields requires non-zero polarizabilities'
    end if
    if (pelib_amoeba_damp) then
        coeff = factor**(i3) * nrm2(Rij) * (alpha_i * alpha_j)**(-i6)
    else
        coeff = factor * nrm2(Rij) * (alpha_i * alpha_j)**(-i6)
    end if

end subroutine damping_coefficient

!------------------------------------------------------------------------------

recursive function factorial(n) result(nfact)

    ! Clive Page, http://www.star.le.ac.uk/~cgp/f90course/f90.html

    integer(ip), intent(in) :: n

    integer(ip) :: nfact

    if (n > 0) then
        nfact = n * factorial(n-1)
    else
        nfact = 1
    end if

end function factorial

!------------------------------------------------------------------------------

end module pelib_multipole_tools
