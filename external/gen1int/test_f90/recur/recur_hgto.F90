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
!!  Recursive functions related to HGTOs.
!!
!!  2012-03-19, Bin Gao
!!  * first version

#include "stdout.h"

!> \brief recursive functions related to HGTOs
!> \author Bin Gao
!> \date 2012-03-19
module recur_hgto
  use xkind
  implicit none

  public :: recur_hgto_to_cgto
  public :: recur_hgto_value

  contains

  !> \brief recursive function for \fn(hgto_to_cgto)
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param half_recip_expnt is the half of the reciprocal of exponent
  !> \param order_cgto contains the orders of xyz components of CGTO
  !> \param order_hgto contains the orders of xyz components of HGTO
  !> \return cgto_pint is the primitive CGTO integral
  recursive function recur_hgto_to_cgto(half_recip_expnt, order_cgto, order_hgto) &
                     result(cgto_pint)
    real(REALK) cgto_pint
    real(REALK), intent(in) :: half_recip_expnt
    integer, intent(in) :: order_cgto(3)
    integer, intent(in) :: order_hgto(3)
    real(REALK) xyz_function  !function returning artificial value for given orders of xyz components
    ! recurrence relation along x-direction
    if (order_cgto(1)>0) then
      cgto_pint = half_recip_expnt                                         &
                * (recur_hgto_to_cgto(half_recip_expnt,                    &
                                      (/order_cgto(1)-1,order_cgto(2:3)/), &
                                      (/order_hgto(1)+1,order_hgto(2:3)/)) &
                + real(order_cgto(1)-1,REALK)                              &
                * recur_hgto_to_cgto(half_recip_expnt,                     &
                                     (/order_cgto(1)-2,order_cgto(2:3)/),  &
                                     order_hgto))
    ! recurrence relation along y-direction
    else if (order_cgto(2)>0) then
      cgto_pint = half_recip_expnt                                                     &
                * (recur_hgto_to_cgto(half_recip_expnt,                                &
                                      (/order_cgto(1),order_cgto(2)-1,order_cgto(3)/), &
                                      (/order_hgto(1),order_hgto(2)+1,order_hgto(3)/)) &
                + real(order_cgto(2)-1,REALK)                                          &
                * recur_hgto_to_cgto(half_recip_expnt,                                 &
                                     (/order_cgto(1),order_cgto(2)-2,order_cgto(3)/),  &
                                     order_hgto))
    ! recurrence relation along z-direction
    else if (order_cgto(3)>0) then
      cgto_pint = half_recip_expnt                                         &
                * (recur_hgto_to_cgto(half_recip_expnt,                    &
                                      (/order_cgto(1:2),order_cgto(3)-1/), &
                                      (/order_hgto(1:2),order_hgto(3)+1/)) &
                + real(order_cgto(3)-1,REALK)                              &
                * recur_hgto_to_cgto(half_recip_expnt,                     &
                                     (/order_cgto(1:2),order_cgto(3)-2/),  &
                                     order_hgto))
    else if (any(order_cgto<0)) then
      cgto_pint = 0.0_REALK
    ! gets the artificial value for given orders of xyz components
    else
      cgto_pint = (1.0_REALK/half_recip_expnt)**sum(order_hgto) &
                * xyz_function(order_hgto)
    end if
  end function recur_hgto_to_cgto

  !> \brief recursive function for a given order HGTO
  !> \author Bin Gao
  !> \date 2012-02-10
  !> \param order_hgto contains the orders of xyz components of HGTO
  !> \param half_nr_expnt is the half of the negative reciprocal of exponent of HGTO
  !> \param elec_coord contains the coordinates of electron
  !> \return val_hgto is the value of the given order HGTO
  recursive function recur_hgto_value(order_hgto, half_nr_expnt, elec_coord) &
                     result(val_hgto)
    real(REALK) val_hgto
    integer, intent(in) :: order_hgto(3)
    real(REALK), intent(in) :: half_nr_expnt
    real(REALK), intent(in) :: elec_coord(3)
    ! recurrence relation along x-direction
    if (order_hgto(1)>0) then
      val_hgto = elec_coord(1)                                         &
               * recur_hgto_value((/order_hgto(1)-1,order_hgto(2:3)/), &
                                  half_nr_expnt, elec_coord)           &
               + real(order_hgto(1)-1,REALK)*half_nr_expnt             &
               * recur_hgto_value((/order_hgto(1)-2,order_hgto(2:3)/), &
                                  half_nr_expnt, elec_coord)
    ! recurrence relation along y-direction
    else if (order_hgto(2)>0) then
      val_hgto = elec_coord(2)                                                     &
               * recur_hgto_value((/order_hgto(1),order_hgto(2)-1,order_hgto(3)/), &
                                  half_nr_expnt, elec_coord)                       &
               + real(order_hgto(2)-1,REALK)*half_nr_expnt                         &
               * recur_hgto_value((/order_hgto(1),order_hgto(2)-2,order_hgto(3)/), &
                                  half_nr_expnt, elec_coord)
    ! recurrence relation along z-direction
    else if (order_hgto(3)>0) then
      val_hgto = elec_coord(3)                                         &
               * recur_hgto_value((/order_hgto(1:2),order_hgto(3)-1/), &
                                  half_nr_expnt, elec_coord)           &
               + real(order_hgto(3)-1,REALK)*half_nr_expnt             &
               * recur_hgto_value((/order_hgto(1:2),order_hgto(3)-2/), &
                                  half_nr_expnt, elec_coord)
    ! zero integral
    else if (any(order_hgto<0)) then
      val_hgto = 0.0_REALK
    ! zeroth order HGTO
    else
      val_hgto = elec_coord(1)**2+elec_coord(2)**2+elec_coord(3)**2
      val_hgto = exp(0.5_REALK/half_nr_expnt*val_hgto)
    end if
  end function recur_hgto_value

end module recur_hgto
