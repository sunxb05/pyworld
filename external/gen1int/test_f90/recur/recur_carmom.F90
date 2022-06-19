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
!!  Recursive functions of Cartesian multipole moment integrals.
!!
!!  2012-03-19, Bin Gao
!!  * first version

#include "stdout.h"

!> \brief recursive functions of Cartesian multipole moment integrals
!> \author Bin Gao
!> \date 2012-03-19
module recur_carmom
  use xkind
  implicit none
! constant Pi
#include "private/pi.h"
  real(REALK), private, save :: zero_pint         !Gauss-Weierstrass transform of \f$1\f$
  real(REALK), private, save :: half_recip_nup    !\f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
  real(REALK), private, save :: neg_ratio_braket  !negative ratio between exponents of bra and ket centers
  real(REALK), private, save :: hrp_total_expnt   !half of the reciprocal of total exponent
  real(REALK), private, save :: cc_wrt_bra(3)     !relative coordinates of center-of-charge w.r.t. bra center
  real(REALK), private, save :: cc_wrt_diporg(3)  !relative coordinates of center-of-charge w.r.t. dipole origin

  public :: recur_carmom_create
  public :: recur_carmom_hbra
  public :: recur_carmom_hrr_ket
  public :: recur_carmom_moment

  contains

  !> \brief prepares the prerequisite quantities used in the following recursive functions
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param coord_bra is the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param coord_ket is the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for operator
  !> \param order_elec is the order of electronic derivatives
  subroutine recur_carmom_create(coord_bra, exponent_bra, coord_ket, exponent_ket, &
                                 dipole_origin, scal_const, order_elec)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_elec
    real(REALK) recip_total_expnt  !reciprocal of total exponent
    real(REALK) reduced_expnt      !reduced exponent
    real(REALK) coord_cc           !xyz coordinates of center-of-charge
    real(REALK) sd_bra_ket         !square of the relative distance between bra and ket centers
    integer ixyz                   !incremental recorder over xyz components
    ! computes the reciprocal of total exponent
    recip_total_expnt = 1.0_REALK/(exponent_bra+exponent_ket)
    ! computes the half of the reciprocal of total exponent
    hrp_total_expnt = 0.5_REALK*recip_total_expnt
    ! computes the negative ratio between exponents of bra and ket centers
    neg_ratio_braket = -exponent_bra/exponent_ket
    ! computes \f$-\frac{1}{2p_{ij}}\frac{b_{j\lambda}}{a_{i\kappa}}\f$
    half_recip_nup = -hrp_total_expnt*exponent_ket/exponent_bra
    ! computes the relative coordinates and square of the relative distance between bra and ket centers
    sd_bra_ket = 0.0_REALK
    do ixyz = 1, 3
      coord_cc = recip_total_expnt*(exponent_bra*coord_bra(ixyz) &
               + exponent_ket*coord_ket(ixyz))
      cc_wrt_bra(ixyz) = coord_cc-coord_bra(ixyz)
      cc_wrt_diporg(ixyz) = coord_cc-dipole_origin(ixyz)
      sd_bra_ket = sd_bra_ket+(coord_ket(ixyz)-coord_bra(ixyz))**2
    end do
    ! computes the reduced exponent
    reduced_expnt = exponent_bra*exponent_ket*recip_total_expnt
    ! computes the Gauss-Weierstrass transform of \f$1\f$
    zero_pint = scal_const*(-exponent_ket-exponent_ket)**order_elec &
              * exp(-reduced_expnt*sd_bra_ket)*sqrt((PI*recip_total_expnt)**3)
  end subroutine recur_carmom_create

  !> \brief recursive function for \fn(carmom_hbra)
  !> \author Bin Gao
  !> \date 2012-02-12
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \return hbra_pint is the primitive zeroth order Cartesian multipole moment
  !>         integral with requied order of HGTO on bra center
  recursive function recur_carmom_hbra(order_hbra) result(hbra_pint)
    real(REALK) hbra_pint
    integer, intent(in) :: order_hbra(3)
    ! recurrence relation along x-direction
    if (order_hbra(1)>0) then
      hbra_pint = cc_wrt_bra(1)                                          &
                * recur_carmom_hbra((/order_hbra(1)-1,order_hbra(2:3)/)) &
                + half_recip_nup*real(order_hbra(1)-1,REALK)             &
                * recur_carmom_hbra((/order_hbra(1)-2,order_hbra(2:3)/))
    ! recurrence relation along y-direction
    else if (order_hbra(2)>0) then
      hbra_pint = cc_wrt_bra(2)                                                      &
                * recur_carmom_hbra((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/)) &
                + half_recip_nup*real(order_hbra(2)-1,REALK)                         &
                * recur_carmom_hbra((/order_hbra(1),order_hbra(2)-2,order_hbra(3)/))
    ! recurrence relation along z-direction
    else if (order_hbra(3)>0) then
      hbra_pint = cc_wrt_bra(3)                                          &
                * recur_carmom_hbra((/order_hbra(1:2),order_hbra(3)-1/)) &
                + half_recip_nup*real(order_hbra(3)-1,REALK)             &
                * recur_carmom_hbra((/order_hbra(1:2),order_hbra(3)-2/))
    ! zero integral
    else if (any(order_hbra<0)) then
      hbra_pint = 0.0_REALK
    ! zeroth order of Cartesian multipole moment
    else
      hbra_pint = zero_pint
    end if
  end function recur_carmom_hbra

  !> \brief recursive function for \fn(carmom_hrr_ket)
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \param order_hket contains the orders of xyz components of HGTO on ket center
  !> \return hket_pint is the zeroth order Cartesian multipole moment integral with requied
  !>         order of HGTOs on bra and ket centers
  recursive function recur_carmom_hrr_ket(order_hbra, order_hket) result(hket_pint)
    real(REALK) hket_pint
    integer, intent(in) :: order_hbra(3)
    integer, intent(in) :: order_hket(3)
    ! recurrence relation along x-direction
    if (order_hket(1)>0) then
      hket_pint = neg_ratio_braket                                          &
                * recur_carmom_hrr_ket((/order_hbra(1)+1,order_hbra(2:3)/), &
                                       (/order_hket(1)-1,order_hket(2:3)/))
    ! recurrence relation along y-direction
    else if (order_hket(2)>0) then
      hket_pint = neg_ratio_braket                                                      &
                * recur_carmom_hrr_ket((/order_hbra(1),order_hbra(2)+1,order_hbra(3)/), &
                                       (/order_hket(1),order_hket(2)-1,order_hket(3)/))
    ! recurrence relation along z-direction
    else if (order_hket(3)>0) then
      hket_pint = neg_ratio_braket                                          &
                * recur_carmom_hrr_ket((/order_hbra(1:2),order_hbra(3)+1/), &
                                       (/order_hket(1:2),order_hket(3)-1/))
    ! zero integral
    else if (any(order_hket<0)) then
      hket_pint = 0.0_REALK
    ! recurrence relations on the order of HGTO on bra center
    else
      hket_pint = recur_carmom_hbra(order_hbra)
    end if
  end function recur_carmom_hrr_ket

  !> \brief recursive function for \fn(carmom_moment)
  !> \author Bin Gao
  !> \date 2012-02-16
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \param order_hket contains the orders of xyz components of HGTO on ket center
  !> \param order_mom contains the orders of xyz components of Cartesian multipole moment
  !> \return hmom_pint is the Cartesian multipole moment integral with requied
  !>         order of HGTOs on bra and ket centers
  recursive function recur_carmom_moment(order_hbra, order_hket, order_mom) &
                     result(hmom_pint)
    real(REALK) hmom_pint
    integer, intent(in) :: order_hbra(3)
    integer, intent(in) :: order_hket(3)
    integer, intent(in) :: order_mom(3)
    ! recurrence relation along x-direction
    if (order_mom(1)>0) then
      hmom_pint = cc_wrt_diporg(1)                                         &
                * recur_carmom_moment(order_hbra, order_hket,              &
                                      (/order_mom(1)-1,order_mom(2:3)/))   &
                + hrp_total_expnt*(real(order_hbra(1),REALK)               &
                * recur_carmom_moment((/order_hbra(1)-1,order_hbra(2:3)/), &
                                      order_hket,                          &
                                      (/order_mom(1)-1,order_mom(2:3)/))   &
                + real(order_hket(1),REALK)                                &
                * recur_carmom_moment(order_hbra,                          &
                                      (/order_hket(1)-1,order_hket(2:3)/), &
                                      (/order_mom(1)-1,order_mom(2:3)/))   &
                + real(order_mom(1)-1,REALK)                               &
                * recur_carmom_moment(order_hbra, order_hket,              &
                                      (/order_mom(1)-2,order_mom(2:3)/)))
    ! recurrence relation along y-direction
    else if (order_mom(2)>0) then
      hmom_pint = cc_wrt_diporg(2)                                                     &
                * recur_carmom_moment(order_hbra, order_hket,                          &
                                      (/order_mom(1),order_mom(2)-1,order_mom(3)/))    &
                + hrp_total_expnt*(real(order_hbra(2),REALK)                           &
                * recur_carmom_moment((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/), &
                                      order_hket,                                      &
                                      (/order_mom(1),order_mom(2)-1,order_mom(3)/))    &
                + real(order_hket(2),REALK)                                            &
                * recur_carmom_moment(order_hbra,                                      &
                                      (/order_hket(1),order_hket(2)-1,order_hket(3)/), &                         
                                      (/order_mom(1),order_mom(2)-1,order_mom(3)/))    &
                + real(order_mom(2)-1,REALK)                                           &
                * recur_carmom_moment(order_hbra, order_hket,                          &
                                      (/order_mom(1),order_mom(2)-2,order_mom(3)/)))
    ! recurrence relation along z-direction
    else if (order_mom(3)>0) then
      hmom_pint = cc_wrt_diporg(3)                                         &
                * recur_carmom_moment(order_hbra, order_hket,              &
                                      (/order_mom(1:2),order_mom(3)-1/))   &
                + hrp_total_expnt*(real(order_hbra(3),REALK)               &
                * recur_carmom_moment((/order_hbra(1:2),order_hbra(3)-1/), &
                                      order_hket,                          &
                                      (/order_mom(1:2),order_mom(3)-1/))   &
                + real(order_hket(3),REALK)                                &
                * recur_carmom_moment(order_hbra,                          &
                                      (/order_hket(1:2),order_hket(3)-1/), &                         
                                      (/order_mom(1:2),order_mom(3)-1/))   &
                + real(order_mom(3)-1,REALK)                               &
                * recur_carmom_moment(order_hbra, order_hket,              &
                                      (/order_mom(1:2),order_mom(3)-2/)))
    ! zero integral
    else if (any(order_mom<0) .or. any(order_hket<0) .or. any(order_hbra<0)) then
      hmom_pint = 0.0_REALK
    ! recurrence relations on the order of HGTO on ket center
    else
      hmom_pint = recur_carmom_hrr_ket(order_hbra, order_hket)
    end if
  end function recur_carmom_moment

end module recur_carmom
