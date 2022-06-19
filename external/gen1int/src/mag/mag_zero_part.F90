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
!!  This file contains the recurrence relations of constructing magnetic
!!  and total rotational angular momentum derivatives (at zero filed) with
!!  rotational London atomic-orbitals.
!!
!!  2010-11-01, Bin Gao
!!  * first version

!  !> \brief constructs total magnetic or total rotational angular momentum
!  !>        derivatives (at zero field, with rotational London atomic-orbitals)
!  !>        from their partial couterparts, using the recurrence relation of
!  !>        \f$\{i+1,jk\}=\frac{i+1}{i+1-K_{0}}[\{ij,k+1\}+\{i,j+1,k\}]\f$
!  !> \author Bin Gao
!  !> \date 2010-11-01
!  !> \param order_total is the order (\f$i\f$) of total derivative
!  !> \param order_operator is the order \f$K_{0}\f$, derivative on the operator
!  !> \param part_deriv contains the partial derivatives \f$\{0j_{1}k_{1}\}\f$,
!  !>        where \f$j_{1}+k_{1}=i+j+k\f$, \f$j\le j_{1}\le j+i\f$, and
!  !>        \f$k\le k_{1}\le k+i\f$; on exit, it is destroyed
!  !> \return total_deriv is \f$\{ijk\}\f$
!  subroutine mag_deriv_total(order_total, order_operator, part_deriv, total_deriv)
!    implicit none
!    integer, intent(in) :: order_total
!    integer, intent(in) :: order_operator
!    real(8), intent(inout) :: part_deriv(order_total+1)
!    real(8), intent(out) :: total_deriv(dim_cints,num_bra_mag,num_ket_mag,num_tot_mag)
!!f2py intent(in) :: order_total
!!f2py intent(in) :: order_operator
!!f2py intent(inout) :: part_deriv
!!f2py depend(order_total) :: part_deriv
!!f2py intent(out) :: total_deriv
!    integer iorder, jorder    !incremental recorders
!    real(8) curr_tot_order    !current total order
!    real(8) curr_denominator  !\f$i+1-K_{0}\f$
!    real(8) scal_factor       !scale factor \f$\frac{i+1}{i+1-K_{0}}\f$
!    curr_tot_order = dfloat(order_operator)+1.0D+00
!    curr_denominator = 1.0D+00
!    do iorder = 0, order_total-2
!      scal_factor = curr_tot_order/curr_denominator
!      ! \f$\{i+1,jk\}=\frac{i+1}{i+1-K_{0}}[\{ij,k+1\}+\{i,j+1,k\}]\f$
!      do jorder = 1, order_total-iorder
!        part_deriv(jorder) = scal_factor*(part_deriv(jorder)+part_deriv(jorder+1))
!      end do
!      curr_tot_order = curr_tot_order+1.0D+00
!      curr_denominator = curr_denominator+1.0D+00
!    end do
!    scal_factor = curr_tot_order/curr_denominator
!    total_deriv = scal_factor*(part_deriv(1)+part_deriv(2))
!    return
!  end subroutine mag_deriv_total

  !> \brief constructs partial magnetic derivatives (at zero field, with London
  !>        atomic-orbitals)
  !> \author Bin Gao
  !> \date 2011-02-07
  !> \param order_aux is the order of auxiliary integrals
  !> \param low_geo is the lower order of geometric derivatives
  !> \param up_geo is the upper order of geometric derivatives
  !> \param order_mag is the required order of magnetic derivatives
  !> \param val_aux(order_aux,low_geo:up_geo)
  !> \return val_mag(low_geo+:up_geo,order_mag)
  subroutine drecurr_lgto_zero_geo_mag(dparams, order_aux, low_geo, up_geo, order_mag)
    implicit none
    real(8), intent(in) :: dparams
    integer, intent(in) :: order_aux
    integer, intent(in) :: low_geo
    integer, intent(in) :: up_geo
    integer, intent(in) :: order_mag
!>>>>!f2py intent(in) :: dparams
!>>>>!f2py intent(in) :: order_aux
!>>>>!f2py intent(in) :: low_geo
!>>>>!f2py intent(in) :: up_geo
!>>>>!f2py intent(in) :: order_mag
!>>>>! number of xyz components for given orders
!>>>>#include "tri_size.h"
!>>>>    integer curr_order_aux  !current order of auxiliary integrals
!>>>>    integer curr_size_aux   !current size of xyz triangular auxiliary integrals
!>>>>    integer next_low_geo    !order of next lower order geometric derivatives
!>>>>    integer curr_size_glow  !current size of xyz triangular lower order geometric derivatives
!>>>>    integer curr_size_gup   !current size of xyz triangular upper order geometric derivatives
!>>>>    integer iorder, jorder  !incremental recorders of order of derivatives
!>>>>    integer iterm, jterm    !incremental recorders of xyz components of derivatives
!>>>>    ! magnetic derivatives only
!>>>>    if (up_geo==0) then
!>>>>      curr_order_aux = order_aux
!>>>>      ! recurrence relations of getting required order of magnetic derivatives
!>>>>      ! from auxiliary integrals (order decreasing), order by order
!>>>>      do iorder = 1, order_mag
!>>>>        ! recurrence relations along x, y and z directions, respectively
!>>>>        do iterm = iorder, 1, -1
!>>>>          do jterm = 1, iterm
!>>>>            call drecurr_lgto_zero_aux_x(curr_order_aux)
!>>>>          end do
!>>>>          call drecurr_lgto_zero_aux_y(curr_order_aux)
!>>>>        end do
!>>>>        call drecurr_lgto_zero_aux_z(curr_order_aux)
!>>>>        ! decreases the order of auxiliary integrals
!>>>>        curr_order_aux = curr_order_aux-1
!>>>>      end do
!>>>>    ! mixed magnetic and geometric derivatives
!>>>>    else
!>>>>      ! computes the sizes of xyz triangles
!>>>>      iorder = min(order_aux, low_geo)
!>>>>      jorder = max(order_aux, up_geo)
!>>>>      curr_order_aux = order_aux
!>>>>      ! recurrence relations of getting required order of magnetic derivatives
!>>>>      ! from auxiliary integrals (order decreasing), order by order
!>>>>      do iorder = 1, order_mag
!>>>>        next_low_geo = low_geo+iorder
!>>>>        curr_size_aux = TRI_SIZE(curr_order_aux)
!>>>>        ! recurrence relations along x, y and z directions, respectively
!>>>>        do iterm = iorder, 1, -1
!>>>>          do jterm = 1, iterm
!>>>>            ! recurrence relations by considering geometric derivatives
!>>>>            do jorder = next_low_geo, up_geo
!>>>>              curr_size_glow = TRI_SIZE(next_low_geo-1)
!>>>>              curr_size_gup = TRI_SIZE(next_low_geo)
!>>>>
!>>>>              ! recurrence relations along x direction
!>>>>              call drecurr_lgto_zero_aux_x(curr_order_aux)
!>>>>            end do
!>>>>          end do
!>>>>          ! recurrence relations by considering geometric derivatives
!>>>>          do jorder = next_low_geo, up_geo
!>>>>            curr_size_glow = TRI_SIZE(next_low_geo-1)
!>>>>            curr_size_gup = TRI_SIZE(next_low_geo)
!>>>>            ! recurrence relations along y direction
!>>>>            
!>>>>          end do
!>>>>        end do
!>>>>        ! recurrence relations by considering geometric derivatives
!>>>>        do jorder = next_low_geo, up_geo
!>>>>          curr_size_glow = TRI_SIZE(next_low_geo-1)
!>>>>          curr_size_gup = TRI_SIZE(next_low_geo)
!>>>>          ! recurrence relations along z direction
!>>>>          
!>>>>        end do
!>>>>        ! decreases the order of auxiliary integrals
!>>>>        curr_order_aux = curr_order_aux-1
!>>>>      end do
!>>>>    end if
    return
  end subroutine drecurr_lgto_zero_geo_mag

  subroutine drecurr_lgto_zero_mag_x()
    implicit none
    return
  end subroutine drecurr_lgto_zero_mag_x

  subroutine drecurr_lgto_zero_mag_y()
    implicit none
    return
  end subroutine drecurr_lgto_zero_mag_y

  subroutine drecurr_lgto_zero_mag_z()
    implicit none
    return
  end subroutine drecurr_lgto_zero_mag_z

  !> \brief two-term recurrence relation by transferring the xyz triangluar auxiliary integrals
  !>        to magnetic derivatives along x direction
  !> \author Bin Gao
  !> \date 2011-02-07
  !> \param order_aux is the order of auxiliary integrals
  subroutine drecurr_lgto_zero_aux_x(order_aux)
    implicit none
    integer, intent(in) :: order_aux
!f2py intent(in) :: order_aux
    integer iorder          !incremental recorder of order of derivatives
    integer iterm, jterm    !incremental recorders of xyz components of derivatives
    integer up_order_aux    !order of auxiliary integrals plus 1
    write(6,*) 'x direction'
    up_order_aux = order_aux+1
    iterm = 1
    do iorder = 1, order_aux
      do jterm = iterm, iterm+order_aux-iorder
        write(6,*) jterm, jterm+iorder, jterm+up_order_aux
      end do
      iterm = iterm+up_order_aux-iorder
    end do
    return
  end subroutine drecurr_lgto_zero_aux_x

  !> \brief two-term recurrence relation by transferring the xyz triangluar auxiliary integrals
  !>        to magnetic derivatives along y direction
  !> \author Bin Gao
  !> \date 2011-02-07
  !> \param order_aux is the order of auxiliary integrals
  subroutine drecurr_lgto_zero_aux_y(order_aux)
    implicit none
    integer, intent(in) :: order_aux
!f2py intent(in) :: order_aux
    integer iorder, jorder  !incremental recorders of order of derivatives
    integer iterm, jterm    !incremental recorders of xyz components of derivatives
    integer up_order_aux    !order of auxiliary integrals plus 1
    write(6,*) 'y direction'
    up_order_aux = order_aux+1
    iterm = 1
    do iorder = 1, order_aux
      jorder = iorder-1
      do jterm = iterm, iterm+order_aux-iorder
        write(6,*) jterm, jterm+jorder, jterm+up_order_aux
      end do
      iterm = iterm+up_order_aux-iorder
    end do
    return
  end subroutine drecurr_lgto_zero_aux_y

  !> \brief two-term recurrence relation by transferring the xyz triangluar auxiliary integrals
  !>        to magnetic derivatives along z direction
  !> \author Bin Gao
  !> \date 2011-02-07
  !> \param order_aux is the order of auxiliary integrals
  subroutine drecurr_lgto_zero_aux_z(order_aux)
    implicit none
    integer, intent(in) :: order_aux
!f2py intent(in) :: order_aux
    integer iorder, jorder  !incremental recorders of order of derivatives
    integer iterm, jterm    !incremental recorders of xyz components of derivatives
    integer up_order_aux    !order of auxiliary integrals plus 1
    write(6,*) 'z direction'
    up_order_aux = order_aux+1
    iterm = 1
    do iorder = 1, order_aux
      jorder = iorder-1
      do jterm = iterm, iterm+order_aux-iorder
        write(6,*) jterm, jterm+jorder, jterm+iorder
      end do
      iterm = iterm+up_order_aux-iorder
    end do
    return
  end subroutine drecurr_lgto_zero_aux_z

! 3 terms
  subroutine drecurr_lgto_zero_aux_ygeo_x()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_ygeo_x

  subroutine drecurr_lgto_zero_aux_zgeo_x()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_zgeo_x

  subroutine drecurr_lgto_zero_aux_zgeo_y()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_zgeo_y

  subroutine drecurr_lgto_zero_aux_xgeo_y()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_xgeo_y

  subroutine drecurr_lgto_zero_aux_xgeo_z()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_xgeo_z

  subroutine drecurr_lgto_zero_aux_ygeo_z()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_ygeo_z

! 4 terms
  subroutine drecurr_lgto_zero_aux_geo_x()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_geo_x

  subroutine drecurr_lgto_zero_aux_geo_y()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_geo_y

  subroutine drecurr_lgto_zero_aux_geo_z()
    implicit none
    return
  end subroutine drecurr_lgto_zero_aux_geo_z

!  !> \brief constructs partial derivatives with respect to the total rotational
!  !>        angular momentum (at zero field, with rotational London atomic-orbitals)
!  !> \author Bin Gao
!  !> \date 2010-11-01
!  !> \param
!  !> \return
!  subroutine drecurr_lgto_zero_geo_rot()
!    implicit none
!    return
!  end subroutine drecurr_lgto_zero_geo_rot
