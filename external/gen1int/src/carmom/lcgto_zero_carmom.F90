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
!!  This file calculates the Cartesian multipole moment integrals
!!  using contracted London Cartesian Gaussians at zero fields.
!!
!!  2011-03-24, Bin Gao:
!!  * first version

  !> \brief calculates the Cartesian multipole moment integrals using
  !>        contracted London Cartesian Gaussians at zero fields
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param idx_bra is the atomic index of bra center
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param idx_ket is the atomic index of ket center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_elec is the order of electronic derivatives
  !> \param order_mag is the order of magnetic derivatives
  !> \param num_cents is the number of differentiated centers
  !> \param idx_cent contains the indices of atomic centers
  !> \param order_cent contains the order of derivatives of the corresponding atomic centers
  !> \param num_cart_bra is the number of Cartesian GTOs on bra,
  !>        equals to \f$(\var(angular_bra)+1)(\var(angular_bra)+2)/2\f$
  !> \param num_cart_ket is the number of Cartesian GTOs on ket,
  !>        equals to \f$(\var(angular_ket)+1)(\var(angular_ket)+2)/2\f$
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the contracted Cartesian multipole moment integrals
  subroutine lcgto_zero_carmom(idx_bra, coord_bra, angular_bra, num_prim_bra, &
                               exponent_bra, num_contr_bra, contr_coef_bra,   &
                               idx_ket, coord_ket, angular_ket, num_prim_ket, &
                               exponent_ket, num_contr_ket, contr_coef_ket,   &
                               order_mag_bra, order_mag_ket, order_mag_total, &
                               order_ram_bra, order_ram_ket, order_ram_total, &
                               order_geo_bra, order_geo_ket,                  &
                               num_cents, idx_cent, order_cent,               &
                               idx_diporg, dipole_origin, scal_const,         &
                               order_geo_mom, order_mom, order_elec,          &
                               num_cart_bra, num_cart_ket, num_opt, contr_ints)
    implicit none
    integer, intent(in) :: idx_bra
    real(8), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(8), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(8), intent(in) :: contr_coef_bra(num_prim_bra,num_contr_bra)
    integer, intent(in) :: idx_ket
    real(8), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(8), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(8), intent(in) :: contr_coef_ket(num_prim_ket,num_contr_ket)
    real(8), intent(in) :: dipole_origin(3)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_elec
    integer, intent(in) :: order_mag
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(2)
    integer, intent(in) :: order_cent(2)
    integer, intent(in) :: num_cints
    real(8), intent(out) :: contr_ints(num_cints)
!f2py intent(in) :: idx_bra
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: idx_ket
!f2py intent(in) :: coord_ket
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(in) :: contr_coef_ket
!f2py depend(num_prim_ket) :: contr_coef_ket
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_elec
!f2py intent(in) :: order_mag
!f2py intent(in) :: num_cents
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py intent(in) :: num_cints
!f2py intent(out) :: contr_ints
!f2py depend(num_cints) :: contr_ints
#if defined(DEBUG)
    ! dumps the contracted GTOs, operators and derivatives to check
    call dump_int_info("LCGTO", idx_bra, coord_bra, angular_bra, &
                       num_prim_bra, exponent_bra,               &
                       num_contr_bra, contr_coef_bra,            &
                       "LCGTO", idx_ket, coord_ket, angular_ket, &
                       num_prim_ket, exponent_ket,               &
                       num_contr_ket, contr_coef_ket,            &
                       "carmom", 1, (/0/), dipole_origin,        &
                       (/order_mom/), order_elec, order_mag,     &
                       2, num_cents, idx_cent, order_cent)
#endif
    ! too many differentiated centers, or the centers of bra and ket coincide
    if (num_cents>2 .or. (num_cents>0 .and. idx_bra==idx_ket)) then
      contr_ints = 0.0D+00
    else
      select case(num_cents)
      ! no geometric derivatives
      case(0)
        call lcgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                      num_prim_bra, exponent_bra,    &
                                      num_contr_bra, contr_coef_bra, &
                                      coord_ket, angular_ket,        &
                                      num_prim_ket, exponent_ket,    &
                                      num_contr_ket, contr_coef_ket, &
                                      dipole_origin, order_mom,      &
                                      order_elec, order_mag,         &
                                      0, 0, num_cints, contr_ints)
      ! one-center geometric derivatives
      case(1)
        if (idx_cent(1)==idx_bra) then
          call lcgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                        num_prim_bra, exponent_bra,    &
                                        num_contr_bra, contr_coef_bra, &
                                        coord_ket, angular_ket,        &
                                        num_prim_ket, exponent_ket,    &
                                        num_contr_ket, contr_coef_ket, &
                                        dipole_origin, order_mom,      &
                                        order_elec, order_mag,         &
                                        order_cent(1), 0, num_cints, contr_ints)
        else if (idx_cent(1)==idx_ket) then
          call lcgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                        num_prim_bra, exponent_bra,    &
                                        num_contr_bra, contr_coef_bra, &
                                        coord_ket, angular_ket,        &
                                        num_prim_ket, exponent_ket,    &
                                        num_contr_ket, contr_coef_ket, &
                                        dipole_origin, order_mom,      &
                                        order_elec, order_mag,         &
                                        0, order_cent(1), num_cints, contr_ints)
        else
          contr_ints = 0.0D+00
        end if
      ! two-center geometric derivatives
      case(2)
        if (idx_cent(1)==idx_bra .and. idx_cent(2)==idx_ket) then
          call lcgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                        num_prim_bra, exponent_bra,    &
                                        num_contr_bra, contr_coef_bra, &
                                        coord_ket, angular_ket,        &
                                        num_prim_ket, exponent_ket,    &
                                        num_contr_ket, contr_coef_ket, &
                                        dipole_origin, order_mom,      &
                                        order_elec, order_mag,         &
                                        order_cent(1), order_cent(2),  &
                                        num_cints, contr_ints)
        else if (idx_cent(2)==idx_bra .and. idx_cent(1)==idx_ket) then
          call lcgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                        num_prim_bra, exponent_bra,    &
                                        num_contr_bra, contr_coef_bra, &
                                        coord_ket, angular_ket,        &
                                        num_prim_ket, exponent_ket,    &
                                        num_contr_ket, contr_coef_ket, &
                                        dipole_origin, order_mom,      &
                                        order_elec, order_mag,         &
                                        order_cent(2), order_cent(1),  &
                                        num_cints, contr_ints)
        else
          contr_ints = 0.0D+00
        end if
      end select
    end if
    return
  end subroutine lcgto_zero_carmom

  !> \brief recurrence relations of Cartesian multipole moment integrals using
  !>        contracted London Cartesian Gaussians at zero field
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra is the angular number of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket is the angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_elec is the order of electronic derivatives
  !> \param order_mag is the order of magnetic derivatives
  !> \param order_geo_bra is the order of geometric derivatives with respect to the bra center
  !> \param order_geo_ket is the order of geometric derivatives with respect to the ket center
  !> \param num_cints is the number of contracted Cartesian multipole moment integrals
  !> \return contr_ints contains the Cartesian multipole moment integrals
  subroutine lcgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                      num_prim_bra, exponent_bra,    &
                                      num_contr_bra, contr_coef_bra, &
                                      coord_ket, angular_ket,        &
                                      num_prim_ket, exponent_ket,    &
                                      num_contr_ket, contr_coef_ket, &
                                      dipole_origin, order_mom,      &
                                      order_elec, &
                                      order_geo_bra, order_geo_ket,  &
                                      order_mag_bra, order_mag_ket, order_ram_bra, order_ram_ket,
                                      num_cints, contr_ints)
    implicit none
    real(8), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra
    integer, intent(in) :: num_prim_bra
    real(8), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(8), intent(in) :: contr_coef_bra(num_prim_bra,num_contr_bra)
    real(8), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket
    integer, intent(in) :: num_prim_ket
    real(8), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(8), intent(in) :: contr_coef_ket(num_prim_ket,num_contr_ket)
    real(8), intent(in) :: dipole_origin(3)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_elec
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_mag
    integer, intent(in) :: order_mag_bra
    integer, intent(in) :: order_mag_ket
    integer, intent(in) :: order_ram
    integer, intent(in) :: order_ram_bra
    integer, intent(in) :: order_ram_ket
    integer, intent(in) :: num_cints
    real(8), intent(out) :: contr_ints(num_cints)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(in) :: contr_coef_ket
!f2py depend(num_prim_ket) :: contr_coef_ket
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_elec
!f2py intent(in) :: order_mag
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(in) :: num_cints
!f2py intent(out) :: contr_ints
!f2py depend(num_cints) :: contr_ints
    integer iprim, jprim                            !incremental recorders over primitives
    integer icontr, jcontr                          !incremental recorders over contractions
    integer ierr                                    !error information
    if (order_mom>0) then
      if (order_ram>0) then
      else
      end if
    else
      if (order_ram>0) then
        do iram = 0, order_ram
        end do
      else
        call contr_cgto_carmom_recurr(coord_bra, angular_bra,        &
                                      num_prim_bra, exponent_bra,    &
                                      num_contr_bra, contr_coef_bra, &
                                      coord_ket, angular_ket,        &
                                      num_prim_ket, exponent_ket,    &
                                      num_contr_ket, contr_coef_ket, &
                                      order_elec, dipole_origin,     &
                                      order_mom, order_geo_mom,      &
                                      order_geo_bra, order_geo_ket,  &
                                      num_cart_bra, num_cart_ket,    &
                                      num_opt, contr_ints)
      end if
    end if
    return
  end subroutine lcgto_zero_carmom_recurr

  !> \brief recurrence relations of Cartesian multipole moment integrals using
  !>        contracted Cartesian Gaussians
  !> \author Bin Gao
  !> \date 2011-03-24
  !> \param coord_bra contains the coordinates of bra center
  !> \param angular_bra contains the minimum and maximum angular numbers of bra center (s=0, p=1, d=2, ...)
  !> \param num_prim_bra is the number of primitive Gaussians of bra center
  !> \param exponent_bra contains the exponents of primitive Gaussians of bra center
  !> \param num_contr_bra is the number of contractions of bra center
  !> \param contr_coef_bra contains the contraction coefficients of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param angular_ket contains the minimum and maximum angular number of ket center
  !> \param num_prim_ket is the number of primitive Gaussians of ket center
  !> \param exponent_ket contains the exponents of primitive Gaussians of ket center
  !> \param num_contr_ket is the number of contractions of ket center
  !> \param contr_coef_ket contains the contraction coefficients of ket center
  !> \param order_elec is the order of electronic derivatives
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_mom is the order of geometric derivatives on Cartesian multipole moments
  !> \param order_geo_bra contains the minimum and maximum orders of geometric derivatives
  !>        with respect to the bra center
  !> \param order_geo_ket contains the minimum and maximum orders of geometric derivatives
  !>        with respect to the ket center
  !> \param dim_cart_bra is the number of Cartesian GTOs on bra
  !> \param dim_cart_ket is the number of Cartesian GTOs on ket
  !> \param num_opt is the number of operators including derivatives
  !> \return contr_ints contains the contracted Cartesian multipole moment integrals
  subroutine contr_cgto_carmom_recurr(coord_bra, angular_bra,        &
                                      num_prim_bra, exponent_bra,    &
                                      num_contr_bra, contr_coef_bra, &
                                      coord_ket, angular_ket,        &
                                      num_prim_ket, exponent_ket,    &
                                      num_contr_ket, contr_coef_ket, &
                                      order_elec, dipole_origin,     &
                                      order_mom, order_geo_mom,      &
                                      order_geo_bra, order_geo_ket,  &
                                      dim_cart_bra, dim_cart_ket,    &
                                      num_opt, contr_ints)
    implicit none
    real(8), intent(in) :: coord_bra(3)
    integer, intent(in) :: angular_bra(2)
    integer, intent(in) :: num_prim_bra
    real(8), intent(in) :: exponent_bra(num_prim_bra)
    integer, intent(in) :: num_contr_bra
    real(8), intent(in) :: contr_coef_bra(num_prim_bra,num_contr_bra)
    real(8), intent(in) :: coord_ket(3)
    integer, intent(in) :: angular_ket(2)
    integer, intent(in) :: num_prim_ket
    real(8), intent(in) :: exponent_ket(num_prim_ket)
    integer, intent(in) :: num_contr_ket
    real(8), intent(in) :: contr_coef_ket(num_prim_ket,num_contr_ket)
    integer, intent(in) :: order_elec
    real(8), intent(in) :: dipole_origin(3)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_mom
    integer, intent(in) :: order_geo_bra(2)
    integer, intent(in) :: order_geo_ket(2)
    integer, intent(in) :: dim_cart_bra
    integer, intent(in) :: dim_cart_ket
    integer, intent(in) :: num_opt
    real(8), intent(out) :: contr_ints(dim_cart_bra,num_contr_bra, &
                                       dim_cart_ket,num_contr_ket,num_opt)
!f2py intent(in) :: coord_bra
!f2py intent(in) :: angular_bra
!f2py intent(hide) :: num_prim_bra
!f2py intent(in) :: exponent_bra
!f2py intent(hide) :: num_contr_bra
!f2py intent(in) :: contr_coef_bra
!f2py depend(num_prim_bra) :: contr_coef_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: angular_ket
!f2py intent(hide) :: num_prim_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: num_contr_ket
!f2py intent(in) :: contr_coef_ket
!f2py depend(num_prim_ket) :: contr_coef_ket
!f2py intent(in) :: order_elec
!f2py intent(in) :: dipole_origin
!f2py intent(in) :: order_mom
!f2py intent(in) :: order_geo_mom
!f2py intent(in) :: order_geo_bra
!f2py intent(in) :: order_geo_ket
!f2py intent(in) :: dim_cart_bra
!f2py intent(in) :: dim_cart_ket
!f2py intent(in) :: num_opt
!f2py intent(out) :: contr_ints
!f2py depend(dim_cart_bra) :: contr_ints
!f2py depend(num_contr_bra) :: contr_ints
!f2py depend(dim_cart_ket) :: contr_ints
!f2py depend(num_contr_ket) :: contr_ints
!f2py depend(num_opt) :: contr_ints
    integer, parameter :: STDOUT = 6  !IO of standard output
    integer order_herm_bra(2)         !range of orders of Hermite Gaussians on bra
    integer order_herm_ket(2)         !range of orders of Hermite Gaussians on ket
    integer low_order_mom             !order of Cartesian multipole moment by considering its derivative
    integer dim_herm_bra              !number of Hermite Gaussians on bra
    integer dim_herm_ket              !number of Hermite Gaussians on ket
    integer num_elec                  !number of xyz components of electronic derivatives
    integer num_low_mom               !number of xyz components of lower order Cartesian multipole moment
    integer num_uniq_ops              !number of unique operators after considering the geometric
                                      !derivative of Cartesian multipole moment
    integer num_geo_bra               !number of geometric derivatives on bra
    integer num_geo_ket               !number of geometric derivatives on ket
    real(8), allocatable :: herm_pints(:,:,:)        !primitive Hermite integrals
    real(8), allocatable :: herm_ket_pints(:,:,:,:)  !primitive Hermite (ket) integrals
    integer prod_op_geo                              !product of number of derivatives and operators
    real(8), allocatable :: cart_pints(:,:,:,:,:)    !primitive Cartesian integrals
    integer dim_contr_cgto            !dimension of contracted CGTOs
    integer start_opt                 !start address of the first operator
    integer iprim, jprim              !incremental recorders over primitives
    integer icontr, jcontr            !incremental recorders over contractions
    integer icart, jcart              !incremental recorders over Cartesian Gaussians
    integer iop                       !incremental recorder over operators
    integer ierr                      !error information
    ! computes the minimum and maximum orders of Hermite Gaussians on bra and ket
    order_herm_bra(1) = order_geo_bra(1) !FIXME +mod(angular_bra,2)
    order_herm_bra(2) = order_geo_bra(2)+angular_bra(2)
    order_herm_ket(1) = order_geo_ket(1) !FIXME +mod(angular_ket,2)
    order_herm_ket(2) = order_geo_ket(2)+angular_ket(2)
    ! computes the order of Cartesian multipole moment after considering its geometric derivative
    low_order_mom = order_mom-order_geo_mom
    ! computes the number of primitive Hermite Gaussians used in recurrence relations
    if (order_herm_bra(1)==0) then
      dim_herm_bra = (order_herm_bra(2)+1)*(order_herm_bra(2)+2)*(order_herm_bra(2)+3)/6
    else
      dim_herm_bra = ((order_herm_bra(2)+1)*(order_herm_bra(2)+2)*(order_herm_bra(2)+3) &
                   -  order_herm_bra(1)*(order_herm_bra(1)+1)*(order_herm_bra(1)+2))/6
    end if
    if (order_herm_ket(1)==0) then
      dim_herm_ket = (order_herm_ket(2)+1)*(order_herm_ket(2)+2)*(order_herm_ket(2)+3)/6
    else
      dim_herm_ket = ((order_herm_ket(2)+1)*(order_herm_ket(2)+2)*(order_herm_ket(2)+3) &
                   -  order_herm_ket(1)*(order_herm_ket(1)+1)*(order_herm_ket(1)+2))/6
    end if
    ! computes the number of operators, which is more efficient than using \var(TRI_SIZE)
    num_elec = (order_elec+1)*(order_elec+2)/2
    num_low_mom = (low_order_mom+1)*(low_order_mom+2)/2
    num_uniq_ops = num_elec*num_low_mom
    ! computes the number of geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
    ! allocates the memory for primitive integrals
    allocate(herm_pints(dim_herm_bra,dim_herm_ket,num_uniq_ops), stat=ierr)
    if (ierr/=0) then
      write(STDOUT,100) "contr_cgto_carmom_recurr: failed to allocate", &
                        dim_herm_bra*dim_herm_ket*num_uniq_ops,         &
                        "primitive Hermite integrals!"
      stop
    end if
    allocate(herm_ket_pints(dim_cart_bra,dim_herm_ket,num_uniq_ops,num_geo_bra), &
             stat=ierr)
    if (ierr/=0) then
      write(STDOUT,100) "contr_cgto_carmom_recurr: failed to allocate",     &
                        dim_cart_bra*dim_herm_ket*num_uniq_ops*num_geo_bra, &
                        "primitive Cartesian-Hermite integrals!"
      stop
    end if
    prod_op_geo = num_uniq_ops*num_geo_bra*num_geo_ket
    ! computes the start address of the first operator
    if (order_geo_mom==0) then
      start_opt = 1
    else
      dim_contr_cgto = dim_cart_bra*num_contr_bra*dim_cart_ket*num_contr_ket
      start_opt = num_opt-prod_op_geo+1
    end if
    allocate(cart_pints(dim_cart_bra,dim_cart_ket,start_opt:num_opt, &
                        num_prim_bra,num_prim_ket), stat=ierr)
    if (ierr/=0) then
      write(STDOUT,100) "contr_cgto_carmom_recurr: failed to allocate",  &
        dim_cart_bra*dim_cart_ket*prod_op_geo*num_prim_bra*num_prim_ket, &
        "primitive Cartesian integrals!"
      stop
    end if
#if defined(DEBUG)
    write(STDOUT,100) "Number of primitive Hermite integrals:", &
                      dim_herm_bra*dim_herm_ket*num_uniq_ops
    write(STDOUT,100) "Number of primitive Cartesian-Hermite integrals:", &
                      dim_cart_bra*dim_herm_ket*num_uniq_ops*num_geo_bra
    write(STDOUT,100) "Number of primitive Cartesian integrals:", &
                      dim_cart_bra*dim_cart_ket*prod_op_geo*num_prim_bra*num_prim_ket
#endif
    ! computes the primitive integrals
    do jprim = 1, num_prim_ket
      do iprim = 1, num_prim_bra
        ! computes the primitive Hermite integrals
        call prim_hgto_carmom(order_herm_bra, coord_bra, exponent_bra(iprim),    &
                              order_herm_ket, coord_ket, exponent_ket(jprim),    &
                              order_elec, dipole_origin, low_order_mom,          &
                              dim_herm_bra, dim_herm_ket, num_elec, num_low_mom, &
                              herm_pints)
        ! transforms Hermite Gaussians on bra to Cartesian ones first
        call hgto_to_lcgto_ket(angular_bra,            &
                              (/order_geo_bra,order_geo_bra/),        &
                              exponent_bra(iprim), 1, dim_herm_bra,   &
                              num_uniq_ops, dim_herm_ket, herm_pints, &
                              dim_cart_bra, num_geo_bra, herm_ket_pints)
        ! then Hermite Gaussians on ket
        call hgto_to_lcgto_ket(angular_ket,               &
                              (/order_geo_ket,order_geo_ket/),           &
                              exponent_ket(jprim), dim_cart_bra,         &
                              dim_herm_ket, num_uniq_ops, num_geo_bra,   &
                              herm_ket_pints, dim_cart_ket, num_geo_ket, &
                              cart_pints(:,:,:,iprim,jprim))
      end do
    end do
    ! cleans
    deallocate(herm_pints)
    deallocate(herm_ket_pints)
    ! constructs contracted integrals (without geometric derivative on Cartesian multipole moment),
    ! in the order of contr_ints(dim_cart_bra,num_contr_bra,dim_cart_ket,num_contr_ket, &
    !                            num_elec,num_low_mom,num_geo_bra,num_geo_ket)
    contr_ints(:,:,:,:,start_opt:num_opt) = 0.0
    do jprim = 1, num_prim_ket
      do iprim = 1, num_prim_bra
        do iop = start_opt, num_opt
          do jcontr = 1, num_contr_ket
            do jcart = 1, dim_cart_ket
              do icontr = 1, num_contr_bra
                do icart = 1, dim_cart_bra
                  !FIXME: uses intermediate results to reduce the number of loops
                  contr_ints(icart,icontr,jcart,jcontr,iop)     &
                    = contr_ints(icart,icontr,jcart,jcontr,iop) &
                    + contr_coef_bra(iprim,icontr)              & !FIXME
                    * contr_coef_ket(jprim,jcontr)              & !FIXME
                    * cart_pints(icart,jcart,iop,iprim,jprim)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    ! cleans
    deallocate(cart_pints)
    ! constructs the geometric derivatives of Cartesian multipole moments
    ! in the order of contr_ints(dim_cart_bra,num_contr_bra,dim_cart_ket,num_contr_ket, &
    !                            num_elec,num_mom,num_geo_mom,num_geo_bra,num_geo_ket)
    if (order_geo_mom>0) then
      call carmom_deriv(order_mom, low_order_mom, order_geo_mom, &
                        .false., dim_contr_cgto*num_elec,        &
                        num_low_mom, num_geo_bra*num_geo_ket,    &
                        contr_ints(:,:,:,:,start_opt:num_opt),   &
                        (order_mom+1)*(order_mom+2)/2,           &
                        (order_geo_mom+1)*(order_geo_mom+2)/2,   &
                        contr_ints(:,:,:,:,1:start_opt-1))
    end if
    return
100 format("contr_cgto_carmom_recurr>> ",A,I8,1X,A)
  end subroutine contr_cgto_carmom_recurr
