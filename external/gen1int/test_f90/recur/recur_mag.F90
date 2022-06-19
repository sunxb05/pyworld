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
!!  Recursive functions related to magnetic derivatives.
!!
!!  2012-03-19, Bin Gao
!!  * first version

#include "stdout.h"

!> \brief recursive functions related to magnetic derivatives
!> \author Bin Gao
!> \date 2012-03-19
module recur_mag
  use xkind
  implicit none
  real(REALK), private, save :: half_recip_b       !\f$1/2b\f$, where \f$b\f$ is the exponent on ket center
  real(REALK), private, save :: ket_wrt_london(3)  !relative coordinates of ket center w.r.t. ...
                                                   !London phase factor origin

  public :: recur_mag_create
  public :: recur_london_mom_hgto

  contains

  !> \brief prepares the prerequisite quantities used in the following recursive functions
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param coord_ket is the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param london_origin contains the coordinates of London phase factor origin
  subroutine recur_mag_create(coord_ket, exponent_ket, london_origin)
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: london_origin(3)
    integer ixyz  !incremental recorder over xyz components
    ! computes the half of the reciprocal of exponent on ket center
    half_recip_b = 0.5_REALK/exponent_ket
    ! computes the relative coordinates of ket center w.r.t. London phase factor origin
    do ixyz = 1, 3
      ket_wrt_london(ixyz) = coord_ket(ixyz)-london_origin(ixyz)
    end do
  end subroutine recur_mag_create

  !> \brief recursive function for \fn(london_mom_hgto)
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param order_mom contains the orders of xyz components of Cartesian multipole moment
  !> \param order_hket contains the orders of xyz components of HGTO
  !> \return lmom_pint is the primitive integral with required Cartesian multipole moment
  recursive function recur_london_mom_hgto(order_mom, order_hket) result(lmom_pint)
    real(REALK) lmom_pint
    integer, intent(in) :: order_mom(3)
    integer, intent(in) :: order_hket(3)
    real(REALK) xyz_function  !function returning artificial value for given orders of xyz components
    ! recurrence relation along x-direction
    if (order_mom(1)>0) then
      lmom_pint = ket_wrt_london(1)                                          &
                * recur_london_mom_hgto((/order_mom(1)-1,order_mom(2:3)/),   &
                                        order_hket)                          &
                + half_recip_b*real(order_hket(1),REALK)                     &
                * recur_london_mom_hgto((/order_mom(1)-1,order_mom(2:3)/),   &
                                        (/order_hket(1)-1,order_hket(2:3)/)) &
                + recur_london_mom_hgto((/order_mom(1)-1,order_mom(2:3)/),   &
                                        (/order_hket(1)+1,order_hket(2:3)/))
    ! recurrence relation along y-direction
    else if (order_mom(2)>0) then
      lmom_pint = ket_wrt_london(2)                                                      &
                * recur_london_mom_hgto((/order_mom(1),order_mom(2)-1,order_mom(3)/),    &
                                        order_hket)                                      &
                + half_recip_b*real(order_hket(2),REALK)                                 &
                * recur_london_mom_hgto((/order_mom(1),order_mom(2)-1,order_mom(3)/),    &
                                        (/order_hket(1),order_hket(2)-1,order_hket(3)/)) &
                + recur_london_mom_hgto((/order_mom(1),order_mom(2)-1,order_mom(3)/),    &
                                        (/order_hket(1),order_hket(2)+1,order_hket(3)/)) 
    ! recurrence relation along z-direction
    else if (order_mom(3)>0) then
      lmom_pint = ket_wrt_london(3)                                          &
                * recur_london_mom_hgto((/order_mom(1:2),order_mom(3)-1/),   &
                                        order_hket)                          &
                + half_recip_b*real(order_hket(3),REALK)                     &
                * recur_london_mom_hgto((/order_mom(1:2),order_mom(3)-1/),   &
                                        (/order_hket(1:2),order_hket(3)-1/)) &
                + recur_london_mom_hgto((/order_mom(1:2),order_mom(3)-1/),   &
                                        (/order_hket(1:2),order_hket(3)+1/)) 
    else if (any(order_mom<0)) then
      lmom_pint = 0.0_REALK
    ! gets the artificial value for given orders of xyz components
    else
      lmom_pint = xyz_function(order_hket)
    end if
  end function recur_london_mom_hgto

end module recur_mag
