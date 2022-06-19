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
!!  using contracted London spherical Gaussians at zero field.
!!
!!  2011-03-24, Bin Gao:
!!  * first version

  !> \brief calculates the Cartesian multipole moment integrals using
  !>        contracted London spherical Gaussians at zero field
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
  !> \param num_cints is the number of contracted Cartesian multipole moment integrals
  !> \return contr_ints contains the contracted Cartesian multipole moment integrals
  subroutine lsgto_zero_carmom(idx_bra, coord_bra, angular_bra, &
                               num_prim_bra, exponent_bra,      &
                               num_contr_bra, contr_coef_bra,   &
                               idx_ket, coord_ket, angular_ket, &
                               num_prim_ket, exponent_ket,      &
                               num_contr_ket, contr_coef_ket,   &
                               dipole_origin, order_mom,        &
                               order_elec, order_mag,           &
                               num_cents, idx_cent, order_cent, &
                               num_cints, contr_ints)
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
    call dump_int_info("LSGTO", idx_bra, coord_bra, angular_bra, &
                       num_prim_bra, exponent_bra,               &
                       num_contr_bra, contr_coef_bra,            &
                       "LSGTO", idx_ket, coord_ket, angular_ket, &
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
        call lsgto_zero_carmom_recurr(coord_bra, angular_bra,        &
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
          call lsgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                        num_prim_bra, exponent_bra,    &
                                        num_contr_bra, contr_coef_bra, &
                                        coord_ket, angular_ket,        &
                                        num_prim_ket, exponent_ket,    &
                                        num_contr_ket, contr_coef_ket, &
                                        dipole_origin, order_mom,      &
                                        order_elec, order_mag,         &
                                        order_cent(1), 0, num_cints, contr_ints)
        else if (idx_cent(1)==idx_ket) then
          call lsgto_zero_carmom_recurr(coord_bra, angular_bra,        &
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
          call lsgto_zero_carmom_recurr(coord_bra, angular_bra,        &
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
          call lsgto_zero_carmom_recurr(coord_bra, angular_bra,        &
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
  end subroutine lsgto_zero_carmom

  !> \brief recurrence relations of Cartesian multipole moment integrals using
  !>        contracted London spherical Gaussians at zero field
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
  !> \return contr_ints contains the Cartesian multipole moment integrals
  subroutine lsgto_zero_carmom_recurr(coord_bra, angular_bra,        &
                                      num_prim_bra, exponent_bra,    &
                                      num_contr_bra, contr_coef_bra, &
                                      coord_ket, angular_ket,        &
                                      num_prim_ket, exponent_ket,    &
                                      num_contr_ket, contr_coef_ket, &
                                      dipole_origin, order_mom,      &
                                      order_elec, order_mag,         &
                                      order_geo_bra, order_geo_ket,  &
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
    integer, intent(in) :: order_mag
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
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
    return
  end subroutine lsgto_zero_carmom_recurr
