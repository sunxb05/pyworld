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
!!  Recursive functions of Dirac delta function integrals.
!!
!!  2012-03-18, Bin Gao
!!  * first version

#include "stdout.h"

!> \brief recursive functions of Dirac delta function integrals
!> \author Bin Gao
!> \date 2012-03-18
module recur_delta
  use xkind
  implicit none
  real(REALK), private, save :: zero_pint            !Gauss-Weierstrass transform of Dirac delta function
  real(REALK), private, save :: half_nr_a            !half of the negative reciprocal of exponent on bra center
  real(REALK), private, save :: half_nr_b            !half of the negative reciprocal of exponent on ket center
  real(REALK), private, save :: neg_twice_p          !negative of the twice total exponent
  real(REALK), private, save :: delta_wrt_bra(3)     !relative coordinates of Dirac delta function origin ...
                                                     !w.r.t. bra center
  real(REALK), private, save :: delta_wrt_ket(3)     !relative coordinates of Dirac delta function origin ...
                                                     !w.r.t. ket center
  real(REALK), private, save :: delta_wrt_cc(3)      !relative coordinates of Dirac delta function origin ...
                                                     !w.r.t. center-of-charge
  real(REALK), private, save :: delta_wrt_diporg(3)  !relative coordinates of Dirac delta function origin ...
                                                     !w.r.t. dipole origin

  public :: recur_delta_create
  public :: recur_delta_geom
  public :: recur_delta_hbra
  public :: recur_delta_hket
  public :: recur_delta_moment

  contains

  !> \brief prepares the prerequisite quantities used in the following recursive functions
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param coord_bra is the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param coord_ket is the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param delta_origin contains the coordinates of Dirac delta function origin
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant
  !> \param order_elec is the order of electronic derivatives
  subroutine recur_delta_create(coord_bra, exponent_bra, coord_ket, exponent_ket, &
                                delta_origin, dipole_origin, scal_const, order_elec)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: delta_origin(3)
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_elec
    real(REALK) neg_total_expnt    !negative total exponent
    real(REALK) recip_total_expnt  !reciprocal of total exponent
    real(REALK) neg_reduced_expnt  !negative reduced exponent
    real(REALK) sd_bra_ket         !square of the relative distance between bra and ket centers
    real(REALK) sd_delta_cc        !square of the relative distance between delta function ...
                                   !and center-of-charge
    integer ixyz                   !incremental recorder over xyz components
    ! sets the half of the negative reciprocal of exponents
    half_nr_a = -0.5_REALK/exponent_bra
    half_nr_b = -0.5_REALK/exponent_ket
    ! computes the negative total exponent, reciprocal of total exponent and negative reduced exponent
    neg_total_expnt = -(exponent_bra+exponent_ket)
    recip_total_expnt = 1.0_REALK/neg_total_expnt
    neg_reduced_expnt = exponent_bra*exponent_ket*recip_total_expnt
    recip_total_expnt = -recip_total_expnt
    ! sets the negative of twice the total exponents
    neg_twice_p = 2.0_REALK*neg_total_expnt
    ! computes the relative coordinates, and squares of the relative distance
    ! between bra and ket centers, delta function and center-of-charge
    sd_bra_ket = 0.0_REALK
    sd_delta_cc = 0.0_REALK
    do ixyz = 1, 3
      delta_wrt_bra(ixyz) = delta_origin(ixyz)-coord_bra(ixyz)
      delta_wrt_ket(ixyz) = delta_origin(ixyz)-coord_ket(ixyz)
      delta_wrt_cc(ixyz) = delta_origin(ixyz)-(exponent_bra*coord_bra(ixyz) &
                         + exponent_ket*coord_ket(ixyz))*recip_total_expnt
      delta_wrt_diporg(ixyz) = delta_origin(ixyz)-dipole_origin(ixyz)
      sd_bra_ket = sd_bra_ket+(coord_ket(ixyz)-coord_bra(ixyz))**2
      sd_delta_cc = sd_delta_cc+delta_wrt_cc(ixyz)**2
    end do
    ! computes the Gauss-Weierstrass transform of Dirac delta function
    zero_pint = scal_const*(-exponent_ket-exponent_ket)**order_elec &
              * exp(neg_reduced_expnt*sd_bra_ket+neg_total_expnt*sd_delta_cc)
  end subroutine recur_delta_create

  !> \brief recursive function for \fn(delta_geom)
  !> \author Bin Gao
  !> \date 2012-02-12
  !> \param order_geo contains the orders of xyz components of geometric derivatives
  !> \return geo_pint is the geometric derivative on Dirac delta function
  recursive function recur_delta_geom(order_geo) result(geo_pint)
    real(REALK) geo_pint
    integer, intent(in) :: order_geo(3)
    ! recurrence relation along x-direction
    if (order_geo(1)>0) then
      geo_pint = neg_twice_p*(delta_wrt_cc(1)                        &
               * recur_delta_geom((/order_geo(1)-1,order_geo(2:3)/)) &
               + real(order_geo(1)-1,REALK)                          &
               * recur_delta_geom((/order_geo(1)-2,order_geo(2:3)/)))
    ! recurrence relation along y-direction
    else if (order_geo(2)>0) then
      geo_pint = neg_twice_p*(delta_wrt_cc(2)                                   &
               * recur_delta_geom((/order_geo(1),order_geo(2)-1,order_geo(3)/)) &
               + real(order_geo(2)-1,REALK)                                     &
               * recur_delta_geom((/order_geo(1),order_geo(2)-2,order_geo(3)/)))
    ! recurrence relation along z-direction
    else if (order_geo(3)>0) then
      geo_pint = neg_twice_p*(delta_wrt_cc(3)                        &
               * recur_delta_geom((/order_geo(1:2),order_geo(3)-1/)) &
               + real(order_geo(3)-1,REALK)                          &
               * recur_delta_geom((/order_geo(1:2),order_geo(3)-2/)))
    ! zero integral
    else if (any(order_geo<0)) then
      geo_pint = 0.0_REALK
    ! zeroth order geometric derivative
    else
      geo_pint = zero_pint
    end if
  end function recur_delta_geom

  !> \brief recursive function for \fn(delta_hket) but on bra center
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \param order_geo contains the orders of xyz components of geometric derivatives
  !> \return hbra_pint is the primitive integral with requied order of HGTOs on bra center
  recursive function recur_delta_hbra(order_hbra, order_geo) result(hbra_pint)
    real(REALK) hbra_pint
    integer, intent(in) :: order_hbra(3)
    integer, intent(in) :: order_geo(3)
    ! recurrence relation along x-direction
    if (order_hbra(1)>0) then
      hbra_pint = delta_wrt_bra(1)                                                 &
                * recur_delta_hbra((/order_hbra(1)-1,order_hbra(2:3)/), order_geo) &
                + real(order_geo(1),REALK)                                         &
                * recur_delta_hbra((/order_hbra(1)-1,order_hbra(2:3)/),            &
                                   (/order_geo(1)-1,order_geo(2:3)/))              &
                + real(order_hbra(1)-1,REALK)*half_nr_a                            &
                * recur_delta_hbra((/order_hbra(1)-2,order_hbra(2:3)/), order_geo)
    ! recurrence relation along y-direction
    else if (order_hbra(2)>0) then
      hbra_pint = delta_wrt_bra(2)                                                  &
                * recur_delta_hbra((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/), &
                                   order_geo)                                       &
                + real(order_geo(2),REALK)                                          &
                * recur_delta_hbra((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/), &
                                   (/order_geo(1),order_geo(2)-1,order_geo(3)/))    &
                + real(order_hbra(2)-1,REALK)*half_nr_a                             &
                * recur_delta_hbra((/order_hbra(1),order_hbra(2)-2,order_hbra(3)/), &
                                   order_geo)
    ! recurrence relation along z-direction
    else if (order_hbra(3)>0) then
      hbra_pint = delta_wrt_bra(3)                                                 &
                * recur_delta_hbra((/order_hbra(1:2),order_hbra(3)-1/), order_geo) &
                + real(order_geo(3),REALK)                                         &
                * recur_delta_hbra((/order_hbra(1:2),order_hbra(3)-1/),            &
                                   (/order_geo(1:2),order_geo(3)-1/))              &
                + real(order_hbra(3)-1,REALK)*half_nr_a                            &
                * recur_delta_hbra((/order_hbra(1:2),order_hbra(3)-2/), order_geo)
    ! zero integral
    else if (any(order_hbra<0) .or. any(order_geo<0)) then
      hbra_pint = 0.0_REALK
    ! recurrence relations on the order of geometric derivatives
    else
      hbra_pint = recur_delta_geom(order_geo)
    end if
  end function recur_delta_hbra

  !> \brief recursive function for \fn(delta_hket)
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \param order_hket contains the orders of xyz components of HGTO on ket center
  !> \param order_geo contains the orders of xyz components of geometric derivatives
  !> \return hket_pint is the primitive integral with requied order of HGTOs on ket center
  recursive function recur_delta_hket(order_hbra, order_hket, order_geo) &
                     result(hket_pint)
    real(REALK) hket_pint
    integer, intent(in) :: order_hbra(3)
    integer, intent(in) :: order_hket(3)
    integer, intent(in) :: order_geo(3)
    ! recurrence relation along x-direction
    if (order_hket(1)>0) then
      hket_pint = delta_wrt_ket(1)                                                  &
                * recur_delta_hket(order_hbra, (/order_hket(1)-1,order_hket(2:3)/), &
                                   order_geo)                                       &
                + real(order_geo(1),REALK)                                          &
                * recur_delta_hket(order_hbra, (/order_hket(1)-1,order_hket(2:3)/), &
                                   (/order_geo(1)-1,order_geo(2:3)/))               &
                + real(order_hket(1)-1,REALK)*half_nr_b                             &
                * recur_delta_hket(order_hbra, (/order_hket(1)-2,order_hket(2:3)/), &
                                   order_geo)
    ! recurrence relation along y-direction
    else if (order_hket(2)>0) then
      hket_pint = delta_wrt_ket(2)                                                  &
                * recur_delta_hket(order_hbra,                                      &
                                   (/order_hket(1),order_hket(2)-1,order_hket(3)/), &
                                   order_geo)                                       &
                + real(order_geo(2),REALK)                                          &
                * recur_delta_hket(order_hbra,                                      &
                                   (/order_hket(1),order_hket(2)-1,order_hket(3)/), &
                                   (/order_geo(1),order_geo(2)-1,order_geo(3)/))    &
                + real(order_hket(2)-1,REALK)*half_nr_b                             &
                * recur_delta_hket(order_hbra,                                      &
                                   (/order_hket(1),order_hket(2)-2,order_hket(3)/), &
                                   order_geo)
    ! recurrence relation along z-direction
    else if (order_hket(3)>0) then
      hket_pint = delta_wrt_ket(3)                                                  &
                * recur_delta_hket(order_hbra, (/order_hket(1:2),order_hket(3)-1/), &
                                   order_geo)                                       &
                + real(order_geo(3),REALK)                                          &
                * recur_delta_hket(order_hbra, (/order_hket(1:2),order_hket(3)-1/), &
                                   (/order_geo(1:2),order_geo(3)-1/))               &
                + real(order_hket(3)-1,REALK)*half_nr_b                             &
                * recur_delta_hket(order_hbra, (/order_hket(1:2),order_hket(3)-2/), &
                                   order_geo)
    ! zero integral
    else if (any(order_hket<0) .or. any(order_geo<0)) then
      hket_pint = 0.0_REALK
    ! recurrence relations on the order of HGTOs on bra center
    else
      hket_pint = recur_delta_hbra(order_hbra, order_geo)
    end if
  end function recur_delta_hket

  !> \brief recursive function for \fn(delta_moment)
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \param order_hket contains the orders of xyz components of HGTO on ket center
  !> \param order_mom contains the orders of xyz components of Cartesian multipole moment
  !> \param order_geo contains the orders of xyz components of geometric derivatives
  !> \return hmom_pint is the primitive Cartesian multipole moment integral with requied
  !>         order of geometric derivatives on Dirac delta function
  recursive function recur_delta_moment(order_hbra, order_hket, order_mom, order_geo) &
                     result(hmom_pint)
    real(REALK) hmom_pint
    integer, intent(in) :: order_hbra(3)
    integer, intent(in) :: order_hket(3)
    integer, intent(in) :: order_mom(3)
    integer, intent(in) :: order_geo(3)
    ! recurrence relation along x-direction
    if (order_mom(1)>0) then
      hmom_pint = delta_wrt_diporg(1)                                   &
                * recur_delta_moment(order_hbra, order_hket,            &
                                     (/order_mom(1)-1,order_mom(2:3)/), &
                                     order_geo)                         &
                + real(order_geo(1),REALK)                              &
                * recur_delta_moment(order_hbra, order_hket,            &
                                     (/order_mom(1)-1,order_mom(2:3)/), &
                                     (/order_geo(1)-1,order_geo(2:3)/))
    ! recurrence relation along y-direction
    else if (order_mom(2)>0) then
      hmom_pint = delta_wrt_diporg(2)                                              &
                * recur_delta_moment(order_hbra, order_hket,                       &
                                     (/order_mom(1),order_mom(2)-1,order_mom(3)/), &
                                     order_geo)                                    &
                + real(order_geo(2),REALK)                                         &
                * recur_delta_moment(order_hbra, order_hket,                       &
                                     (/order_mom(1),order_mom(2)-1,order_mom(3)/), &
                                     (/order_geo(1),order_geo(2)-1,order_geo(3)/))
    ! recurrence relation along z-direction
    else if (order_mom(3)>0) then
      hmom_pint = delta_wrt_diporg(3)                                   &
                * recur_delta_moment(order_hbra, order_hket,            &
                                     (/order_mom(1:2),order_mom(3)-1/), &
                                     order_geo)                         &
                + real(order_geo(3),REALK)                              &
                * recur_delta_moment(order_hbra, order_hket,            &
                                     (/order_mom(1:2),order_mom(3)-1/), &
                                     (/order_geo(1:2),order_geo(3)-1/))
    ! zero integral
    else if (any(order_mom<0) .or. any(order_geo<0)) then
      hmom_pint = 0.0_REALK
    ! recurrence relations on the order of HGTOs on ket center
    else
      hmom_pint = recur_delta_hket(order_hbra, order_hket, order_geo)
    end if
  end function recur_delta_moment

end module recur_delta
