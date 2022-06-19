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
!!  Recursive functions of nuclear attraction potential integrals.
!!
!!  2012-03-19, Bin Gao
!!  * first version

#include "stdout.h"

!> \brief recursive functions of nuclear attraction potential integrals
!> \author Bin Gao
!> \date 2012-03-19
module recur_nucpot
  use xkind
  implicit none
! constant Pi
#include "private/pi.h"
  real(REALK), private, save :: half_recip_p      !half of the reciprocal of total exponent \f$p_{ij}\f$
  real(REALK), private, save :: half_neg_rp       !half of the negative reciprocal of total exponent \f$p_{ij}\f$
  real(REALK), private, save :: bra_to_ket        !ratio of exponent on bra center to that on ket center
  real(REALK), private, save :: ket_to_bra        !ratio of exponent on ket center to that on bra center
  real(REALK), private, save :: cc_wrt_bra(3)     !relative coordinates of center-of-charge w.r.t. bra center
  real(REALK), private, save :: cc_wrt_ket(3)     !relative coordinates of center-of-charge w.r.t. ket center
  real(REALK), private, save :: nucorg_wrt_cc(3)  !relative coordinates of nuclear potential origin ...
                                                  !w.r.t. center-of-charge
  real(REALK), private, save :: cc_wrt_diporg(3)  !relative coordinates of center-of-charge w.r.t. dipole origin
  integer, private, save :: max_order_aux         !maximum order of auxiliary integrals
  real(REALK), private, allocatable, save :: zero_pint(:)
                                                  !Gauss-Weierstrass transform of nuclear attraction potential

  public :: recur_nucpot_create
  public :: recur_nucpot_destroy
  public :: recur_nucpot_geom
  public :: recur_nucpot_hket
  public :: recur_nucpot_hbra
  public :: recur_nucpot_moment

  contains

  !> \brief prepares the prerequisite quantities used in the following recursive functions
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param coord_bra is the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param coord_ket is the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param nucpot_origin contains the coordinates of nuclear potential origin
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for operator
  !> \param order_elec is the order of electronic derivatives
  !> \param order_aux is the maximum order of auxiliary integrals
  !> \param gaupot_expt is the exponent used in the Gaussian broadening function of the charge
  subroutine recur_nucpot_create(coord_bra, exponent_bra, coord_ket, exponent_ket, &
                                 nucpot_origin, dipole_origin, scal_const,         &
                                 order_elec, order_aux, gaupot_expt)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: nucpot_origin(3)
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: order_elec
    integer, intent(in) :: order_aux
    real(REALK), optional, intent(in) :: gaupot_expt
    real(REALK) total_expnt        !total exponent
    real(REALK) recip_total_expnt  !reciprocal of total exponent
    real(REALK) coord_cc           !xyz coordinates of center-of-charge
    real(REALK) sd_bra_ket         !square of the relative distance between bra and ket centers
    real(REALK) reduced_expnt      !reduced exponent
    real(REALK) prefact_zero       !prefactor of integrals of zero order GTOs and operators
    real(REALK) arg_aux_fun        !argument of auxiliary functions
    real(REALK) scale_factor       !scale factor
    integer ixyz                   !incremental recorder over xyz components
    integer iorder                 !incremental recorder over auxiliary functions
    integer ierr                   !error information
    ! computes the total exponent and its reciprocal
    total_expnt = exponent_bra+exponent_ket
    recip_total_expnt = 1.0_REALK/total_expnt
    ! computes the half of the (negative) reciprocal of total exponent \f$p_{ij}\f$
    half_recip_p = 0.5_REALK*recip_total_expnt
    half_neg_rp = -half_recip_p
    ! computes the ratios of exponents
    bra_to_ket = exponent_bra/exponent_ket
    ket_to_bra = exponent_ket/exponent_bra
    ! computes the relative coordinates and square of the relative distance between bra and ket centers
    sd_bra_ket = 0.0_REALK
    do ixyz = 1, 3
      coord_cc = recip_total_expnt*(exponent_bra*coord_bra(ixyz) &
               + exponent_ket*coord_ket(ixyz))
      cc_wrt_bra(ixyz) = coord_cc-coord_bra(ixyz)
      cc_wrt_ket(ixyz) = coord_cc-coord_ket(ixyz)
      nucorg_wrt_cc(ixyz) = nucpot_origin(ixyz)-coord_cc
      cc_wrt_diporg(ixyz) = coord_cc-dipole_origin(ixyz)
      sd_bra_ket = sd_bra_ket+(coord_ket(ixyz)-coord_bra(ixyz))**2
    end do
    ! computes the reduced exponent
    reduced_expnt = exponent_bra*exponent_ket*recip_total_expnt
    ! computes the prefactor of auxiliary integrals
    prefact_zero = scal_const*(-exponent_ket-exponent_ket)**order_elec &
                 * exp(-reduced_expnt*sd_bra_ket)
    ! computes the auxiliary integrals for Gaussian charge potential
    if (present(gaupot_expt)) then
      ! computes the argument of Boys functions
      arg_aux_fun = 0.0_REALK
      do ixyz = 1, 3
        arg_aux_fun = arg_aux_fun+nucorg_wrt_cc(ixyz)**2
      end do
      scale_factor = gaupot_expt/(gaupot_expt+total_expnt)
      total_expnt = scale_factor*total_expnt
      arg_aux_fun = total_expnt*arg_aux_fun
      ! sets the maximum order of Boys functions
      max_order_aux = order_aux
      allocate(zero_pint(0:max_order_aux), stat=ierr)
      if (ierr/=0)                           &
      call error_stop("recur_nucpot_create", &
                      "failed to allocate zero_pint", max_order_aux+1)
      ! computes the Boys functions
      call aux_boys_vec(0, max_order_aux, arg_aux_fun, zero_pint)
      arg_aux_fun = prefact_zero*(PI+PI)*recip_total_expnt*sqrt(scale_factor)
      total_expnt = -2.0_REALK*total_expnt
      ! loops over the order of Boys functions
      do iorder = 0, max_order_aux
        zero_pint(iorder) = arg_aux_fun*zero_pint(iorder)
        arg_aux_fun = arg_aux_fun*total_expnt
      end do
    ! computes the auxiliary integrals for nuclear attraction potential
    else
      ! computes the argument of Boys functions
      arg_aux_fun = 0.0_REALK
      do ixyz = 1, 3
        arg_aux_fun = arg_aux_fun+nucorg_wrt_cc(ixyz)**2
      end do
      arg_aux_fun = total_expnt*arg_aux_fun
      ! sets the maximum order of Boys functions
      max_order_aux = order_aux
      allocate(zero_pint(0:max_order_aux), stat=ierr)
      if (ierr/=0)                           &
      call error_stop("recur_nucpot_create", &
                      "failed to allocate zero_pint", max_order_aux+1)
      ! computes the Boys functions
      call aux_boys_vec(0, max_order_aux, arg_aux_fun, zero_pint)
      arg_aux_fun = prefact_zero*(PI+PI)*recip_total_expnt
      total_expnt = -2.0_REALK*total_expnt
      ! loops over the order of Boys functions
      do iorder = 0, max_order_aux
        zero_pint(iorder) = arg_aux_fun*zero_pint(iorder)
        arg_aux_fun = arg_aux_fun*total_expnt
      end do
    end if
  end subroutine recur_nucpot_create

  !> \brief cleans quantities used in the recursive functions
  !> \author Bin Gao
  !> \date 2011-07-27
  subroutine recur_nucpot_destroy()
    deallocate(zero_pint)
  end subroutine recur_nucpot_destroy

  !> \brief recursive function for \fn(nucpot_geom)
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param order_geo contains the orders of xyz components of geometric derivative on
  !>        nuclear potential center
  !> \param order_boys is the order of Boys function
  !> \return geo_pint is the geometric derivative on nuclear potential center
  recursive function recur_nucpot_geom(order_geo, order_boys) result(geo_pint)
    real(REALK) geo_pint
    integer, intent(in) :: order_geo(3)
    integer, intent(in) :: order_boys
    ! recurrence relation along x-direction
    if (order_geo(1)>0) then
      geo_pint = nucorg_wrt_cc(1)                                     &
               * recur_nucpot_geom((/order_geo(1)-1,order_geo(2:3)/), &
                                   order_boys+1)                      &
               + real(order_geo(1)-1,REALK)                           &
               * recur_nucpot_geom((/order_geo(1)-2,order_geo(2:3)/), &
                                   order_boys+1)
    ! recurrence relation along y-direction
    else if (order_geo(2)>0) then
      geo_pint = nucorg_wrt_cc(2)                                                &
               * recur_nucpot_geom((/order_geo(1),order_geo(2)-1,order_geo(3)/), &
                                   order_boys+1)                                 &
               + real(order_geo(2)-1,REALK)                                      &
               * recur_nucpot_geom((/order_geo(1),order_geo(2)-2,order_geo(3)/), &
                                   order_boys+1)
    ! recurrence relation along z-direction
    else if (order_geo(3)>0) then
      geo_pint = nucorg_wrt_cc(3)                                     &
               * recur_nucpot_geom((/order_geo(1:2),order_geo(3)-1/), &
                                   order_boys+1)                      &
               + real(order_geo(3)-1,REALK)                           &
               * recur_nucpot_geom((/order_geo(1:2),order_geo(3)-2/), &
                                   order_boys+1)
    ! zero integral
    else if (any(order_geo<0)) then
      geo_pint = 0.0_REALK
    ! auxiliary integral
    else
      if (order_boys<=max_order_aux) then
        geo_pint = zero_pint(order_boys)
      else
        geo_pint = 0.0_REALK  !gets rid of compiler's warning
        call error_stop("recur_nucpot_geom", &
                        "no required order of auxiliary integrals", order_boys)
      end if
    end if
  end function recur_nucpot_geom

  !> \brief recursive function for \fn(nucpot_hket)
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param order_hket contains the orders of xyz components of HGTO on ket center
  !> \param order_geo contains the orders of xyz components of geometric derivatives
  !> \return hket_pint is the primitive integral with requied order of HGTOs on ket center
  recursive function recur_nucpot_hket(order_hket, order_geo) result(hket_pint)
    real(REALK) hket_pint
    integer, intent(in) :: order_hket(3)
    integer, intent(in) :: order_geo(3)
    ! recurrence relation along x-direction
    if (order_hket(1)>0) then
      hket_pint = cc_wrt_ket(1)                                          &
                * recur_nucpot_hket((/order_hket(1)-1,order_hket(2:3)/), &
                                    order_geo)                           &
                + half_neg_rp*(bra_to_ket*real(order_hket(1)-1,REALK)    &
                * recur_nucpot_hket((/order_hket(1)-2,order_hket(2:3)/), &
                                    order_geo)                           &
                + recur_nucpot_hket((/order_hket(1)-1,order_hket(2:3)/), &
                                    (/order_geo(1)+1,order_geo(2:3)/)))
    ! recurrence relation along y-direction
    else if (order_hket(2)>0) then
      hket_pint = cc_wrt_ket(2)                                                      &
                * recur_nucpot_hket((/order_hket(1),order_hket(2)-1,order_hket(3)/), &
                                    order_geo)                                       &
                + half_neg_rp*(bra_to_ket*real(order_hket(2)-1,REALK)                &
                * recur_nucpot_hket((/order_hket(1),order_hket(2)-2,order_hket(3)/), &
                                    order_geo)                                       &
                + recur_nucpot_hket((/order_hket(1),order_hket(2)-1,order_hket(3)/), &
                                    (/order_geo(1),order_geo(2)+1,order_geo(3)/)))
    ! recurrence relation along z-direction
    else if (order_hket(3)>0) then
      hket_pint = cc_wrt_ket(3)                                          &
                * recur_nucpot_hket((/order_hket(1:2),order_hket(3)-1/), &
                                    order_geo)                           &
                + half_neg_rp*(bra_to_ket*real(order_hket(3)-1,REALK)    &
                * recur_nucpot_hket((/order_hket(1:2),order_hket(3)-2/), &
                                    order_geo)                           &
                + recur_nucpot_hket((/order_hket(1:2),order_hket(3)-1/), &
                                    (/order_geo(1:2),order_geo(3)+1/)))
    ! zero integral
    else if (any(order_hket<0)) then
      hket_pint = 0.0_REALK
    ! recurrence relations on the order of geometric derivatives
    else
      hket_pint = recur_nucpot_geom(order_geo, 0)
    end if
  end function recur_nucpot_hket

  !> \brief recursive function for \fn(nucpot_hbra)
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \param order_hket contains the orders of xyz components of HGTO on ket center
  !> \param order_geo contains the orders of xyz components of geometric derivatives
  !> \return hbra_pint is the primitive integral with requied order of HGTOs on bra center
  recursive function recur_nucpot_hbra(order_hbra, order_hket, order_geo) &
                     result(hbra_pint)
    real(REALK) hbra_pint
    integer, intent(in) :: order_hbra(3)
    integer, intent(in) :: order_hket(3)
    integer, intent(in) :: order_geo(3)
    ! recurrence relation along x-direction
    if (order_hbra(1)>0) then
      hbra_pint = cc_wrt_bra(1)                                          &
                * recur_nucpot_hbra((/order_hbra(1)-1,order_hbra(2:3)/), &
                                    order_hket, order_geo)               &
                + half_neg_rp*(ket_to_bra*real(order_hbra(1)-1,REALK)    &
                * recur_nucpot_hbra((/order_hbra(1)-2,order_hbra(2:3)/), &
                                    order_hket, order_geo)               &
                + recur_nucpot_hbra((/order_hbra(1)-1,order_hbra(2:3)/), &
                                    order_hket,                          &
                                    (/order_geo(1)+1,order_geo(2:3)/))   &
                - real(order_hket(1),REALK)                              &
                * recur_nucpot_hbra((/order_hbra(1)-1,order_hbra(2:3)/), &
                                    (/order_hket(1)-1,order_hket(2:3)/), &
                                    order_geo))
    ! recurrence relation along y-direction
    else if (order_hbra(2)>0) then
      hbra_pint = cc_wrt_bra(2)                                                      &
                * recur_nucpot_hbra((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/), &
                                    order_hket, order_geo)                           &
                + half_neg_rp*(ket_to_bra*real(order_hbra(2)-1,REALK)                &
                * recur_nucpot_hbra((/order_hbra(1),order_hbra(2)-2,order_hbra(3)/), &
                                    order_hket, order_geo)                           &
                + recur_nucpot_hbra((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/), &
                                    order_hket,                                      &
                                    (/order_geo(1),order_geo(2)+1,order_geo(3)/))    &
                - real(order_hket(2),REALK)                                          &
                * recur_nucpot_hbra((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/), &
                                    (/order_hket(1),order_hket(2)-1,order_hket(3)/), &
                                    order_geo))
    ! recurrence relation along z-direction
    else if (order_hbra(3)>0) then
      hbra_pint = cc_wrt_bra(3)                                          &
                * recur_nucpot_hbra((/order_hbra(1:2),order_hbra(3)-1/), &
                                    order_hket, order_geo)               &
                + half_neg_rp*(ket_to_bra*real(order_hbra(3)-1,REALK)    &
                * recur_nucpot_hbra((/order_hbra(1:2),order_hbra(3)-2/), &
                                    order_hket, order_geo)               &
                + recur_nucpot_hbra((/order_hbra(1:2),order_hbra(3)-1/), &
                                    order_hket,                          &
                                    (/order_geo(1:2),order_geo(3)+1/))   &
                - real(order_hket(3),REALK)                              &
                * recur_nucpot_hbra((/order_hbra(1:2),order_hbra(3)-1/), &
                                    (/order_hket(1:2),order_hket(3)-1/), &
                                    order_geo))
    ! zero integral
    else if (any(order_hbra<0) .or. any(order_hket<0)) then
      hbra_pint = 0.0_REALK
    ! recurrence relations on the order of HGTOs on ket center
    else
      hbra_pint = recur_nucpot_hket(order_hket, order_geo)
    end if
  end function recur_nucpot_hbra

  !> \brief recursive function for \fn(nucpot_moment)
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param order_hbra contains the orders of xyz components of HGTO on bra center
  !> \param order_hket contains the orders of xyz components of HGTO on ket center
  !> \param order_geo contains the orders of xyz components of geometric derivatives
  !> \param order_mom contains the orders of xyz components of Cartesian multipole moment
  !> \return hmom_pint is the nuclear attraction integral with required orders of Cartesian
  !>         multipole moment integral and HGTOs on bra and ket centers
  recursive function recur_nucpot_moment(order_hbra, order_hket, order_mom, order_geo) &
                     result(hmom_pint)
    real(REALK) hmom_pint
    integer, intent(in) :: order_hbra(3)
    integer, intent(in) :: order_hket(3)
    integer, intent(in) :: order_mom(3)
    integer, intent(in) :: order_geo(3)
    ! recurrence relation along x-direction
    if (order_mom(1)>0) then
      hmom_pint = cc_wrt_diporg(1)                                         &
                * recur_nucpot_moment(order_hbra, order_hket,              &
                                      (/order_mom(1)-1,order_mom(2:3)/),   &
                                      order_geo)                           &
                + half_recip_p*(real(order_hbra(1),REALK)                  &
                * recur_nucpot_moment((/order_hbra(1)-1,order_hbra(2:3)/), &
                                      order_hket,                          &
                                      (/order_mom(1)-1,order_mom(2:3)/),   &
                                       order_geo)                          &
                + real(order_hket(1),REALK)                                &
                * recur_nucpot_moment(order_hbra,                          &
                                      (/order_hket(1)-1,order_hket(2:3)/), &
                                      (/order_mom(1)-1,order_mom(2:3)/),   &
                                      order_geo)                           &
                + real(order_mom(1)-1,REALK)                               &
                * recur_nucpot_moment(order_hbra, order_hket,              &
                                      (/order_mom(1)-2,order_mom(2:3)/),   &
                                      order_geo)                           &
                - recur_nucpot_moment(order_hbra, order_hket,              &
                                      (/order_mom(1)-1,order_mom(2:3)/),   &
                                      (/order_geo(1)+1,order_geo(2:3)/)))
    ! recurrence relation along y-direction
    else if (order_mom(2)>0) then
      hmom_pint = cc_wrt_diporg(2)                                                     &
                * recur_nucpot_moment(order_hbra, order_hket,                          &
                                      (/order_mom(1),order_mom(2)-1,order_mom(3)/),    &
                                      order_geo)                                       &
                + half_recip_p*(real(order_hbra(2),REALK)                              &
                * recur_nucpot_moment((/order_hbra(1),order_hbra(2)-1,order_hbra(3)/), &
                                      order_hket,                                      & 
                                      (/order_mom(1),order_mom(2)-1,order_mom(3)/),    &
                                       order_geo)                                      &
                + real(order_hket(2),REALK)                                            &
                * recur_nucpot_moment(order_hbra,                                      &
                                      (/order_hket(1),order_hket(2)-1,order_hket(3)/), &
                                      (/order_mom(1),order_mom(2)-1,order_mom(3)/),    &
                                      order_geo)                                       &
                + real(order_mom(2)-1,REALK)                                           &
                * recur_nucpot_moment(order_hbra, order_hket,                          &
                                      (/order_mom(1),order_mom(2)-2,order_mom(3)/),    &
                                      order_geo)                                       &
                - recur_nucpot_moment(order_hbra, order_hket,                          &
                                      (/order_mom(1),order_mom(2)-1,order_mom(3)/),    &
                                      (/order_geo(1),order_geo(2)+1,order_geo(3)/)))
    ! recurrence relation along z-direction
    else if (order_mom(3)>0) then
      hmom_pint = cc_wrt_diporg(3)                                         &
                * recur_nucpot_moment(order_hbra, order_hket,              &
                                      (/order_mom(1:2),order_mom(3)-1/),   &
                                      order_geo)                           &
                + half_recip_p*(real(order_hbra(3),REALK)                  &
                * recur_nucpot_moment((/order_hbra(1:2),order_hbra(3)-1/), &
                                      order_hket,                          & 
                                      (/order_mom(1:2),order_mom(3)-1/),   & 
                                       order_geo)                          &
                + real(order_hket(3),REALK)                                &
                * recur_nucpot_moment(order_hbra,                          &
                                      (/order_hket(1:2),order_hket(3)-1/), &
                                      (/order_mom(1:2),order_mom(3)-1/),   &
                                      order_geo)                           &
                + real(order_mom(3)-1,REALK)                               &
                * recur_nucpot_moment(order_hbra, order_hket,              &
                                      (/order_mom(1:2),order_mom(3)-2/),   &
                                      order_geo)                           &
                - recur_nucpot_moment(order_hbra, order_hket,              &
                                      (/order_mom(1:2),order_mom(3)-1/),   &
                                      (/order_geo(1:2),order_geo(3)+1/)))
    ! zero integral
    else if (any(order_hbra<0) .or. any(order_hket<0) .or. any(order_mom<0)) then
      hmom_pint = 0.0_REALK
    ! recurrence relations on the order of HGTOs on bra center
    else
      hmom_pint = recur_nucpot_hbra(order_hbra, order_hket, order_geo)
    end if
  end function recur_nucpot_moment

end module recur_nucpot
