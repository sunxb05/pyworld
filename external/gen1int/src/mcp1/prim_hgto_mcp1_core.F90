!!  gen1int: compute one-electron integrals using rotational London atomic-orbitals
!!  Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen
!! 
!!  gen1int is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!! 
!!  gen1int is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!! 
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file contains the recurrence relations of primitive Hermite Gaussians.
!!
!!  2010-10-10, Bin Gao:
!!  * first version

! (a) prim_hgto_carmom: Cartesian multipole moments
! (b) prim_hgto_delta:  Delta function
! (c) prim_hgto_nucpot: nuclear attraction potential
! (d) prim_hgto_isdpot: inverse square distance potential
! (e) prim_hgto_gaupot: Gaussian charge potential
! (f) prim_hgto_dso:    diamagnetic spin-orbit coupling
! (g) prim_hgto_ecp:    effective core potential
! (h) prim_hgto_mcp1:   model core potential (Version 1)

!  subroutine prim_hgto_nucpot(min_order_bra, max_order_bra,   &
!                              min_order_ket, max_order_ket,   & 
!                              order_c, order_mom, order_elec, &
!                              gint_matrix)
!    implicit none
!    integer, intent(in) :: min_order_bra
!    integer, intent(in) :: max_order_bra
!    integer, intent(in) :: min_order_ket
!    integer, intent(in) :: max_order_ket
!    integer, intent(in) :: order_c
!    integer, intent(in) :: order_mom
!    integer, intent(in) :: order_elec
!    real(8), intent(out) :: &
!      gint_matrix(min_order_bra:max_order_bra, &
!                  min_order_ket:max_order_ket, &
!                  order_c)
!!f2py intent(in) :: min_order_bra
!!f2py intent(in) :: max_order_bra
!!f2py intent(in) :: min_order_ket
!!f2py intent(in) :: max_order_ket
!!f2py intent(in) :: order_c
!!f2py intent(in) :: order_mom
!!f2py intent(in) :: order_elec
!!f2py intent(out) :: gint_matrix
!    call dboys_vector(order_boys, arg_boys, val_boys)
!#ifdef ROW_MAJOR
!    ! row-major order, ket increases rapidly
!#else
!    ! column-major order, bra increases rapidly
!#endif
!  end subroutine prim_hgto_nucpot
