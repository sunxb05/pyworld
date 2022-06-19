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
!!  This file recovers Cartesian multipole moments (using CGTOs) at the center of
!!  London phase factor.
!!
!!  2011-10-24, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief recovers the HGTOs on bra center in nuclear attraction potential integrals
  !> \author Bin Gao
  !> \date 2011-10-24
  !> \param order_herm_bra contains the minimum and maximum orders of HGTOs on bra center to return
  !> \param order_herm_ket contains the minimum and maximum orders of HGTOs on ket center to return
  !> \param orders_derv contains the minimum and maximum orders of geometric derivatives to return
  !> \param coord_ket contains the coordinates of bra center
  !> \param exponent_ket is the exponent of HGTOs of bra center
  !> \param coord_ket contains the coordinates of bra center
  !> \param exponent_ket is the exponent of HGTOs of bra center
  !> \param dim_hket is the dimension of HGTOs on ket center
  !> \param dim_hket_derv is the dimension of geometric derivatives
  !> \param hket_pints contains the nuclear attraction integrals with zeroth order
  !>        HGTOs on bra center
  !> \param dim_herm_bra is the dimension of HGTOs of bra center afterwards
  !> \param dim_herm_ket is the dimension of HGTOs of ket center afterwards
  !> \param dim_hbra_derv is the dimension of geometric derivatives afterwards
  !> \return hbra_pints contains the integrals with required HGTOs on bra center
  subroutine london_mom_cgto(order_herm_bra, order_herm_ket, orders_derv,       &
                         coord_bra, exponent_bra, coord_ket, exponent_ket,  &
                         dim_hket, dim_hket_derv, hket_pints, dim_herm_bra, &
                         dim_herm_ket, dim_hbra_derv, hbra_pints)
    use xkind
    implicit none
    integer, intent(in) :: order_herm_bra(2)
    integer, intent(in) :: order_herm_ket(2)
    integer, intent(in) :: orders_derv(2)
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    integer, intent(in) :: dim_hket
    integer, intent(in) :: dim_hket_derv
    real(REALK), intent(in) :: hket_pints(dim_hket,dim_hket_derv)
    integer, intent(in) :: dim_herm_bra
    integer, intent(in) :: dim_herm_ket
    integer, intent(in) :: dim_hbra_derv
    real(REALK), intent(out) :: hbra_pints(dim_herm_bra,dim_herm_ket,dim_hbra_derv)
!f2py intent(in) :: order_herm_bra
!f2py intent(in) :: order_herm_ket
!f2py intent(in) :: orders_derv
!f2py intent(in) :: coord_bra
!f2py intent(in) :: exponent_bra
!f2py intent(in) :: coord_ket
!f2py intent(in) :: exponent_ket
!f2py intent(hide) :: dim_hket
!f2py intent(hide) :: dim_hket_derv
!f2py intent(in) :: hket_pints
!f2py intent(in) :: dim_herm_bra
!f2py intent(in) :: dim_herm_ket
!f2py intent(in) :: dim_hbra_derv
!f2py intent(out) :: hbra_pints
!f2py depend(dim_herm_bra) :: hbra_pints
!f2py depend(dim_herm_ket) :: hbra_pints
!f2py depend(dim_hbra_derv) :: hbra_pints
    integer max_low_derv       !maximum order of lower order geometric derivatives
    real(REALK) neg_hrp_expt   !negative half of the reciprocal of exponent
    real(REALK) ket_to_bra     !ratio of exponent on ket center to that on bra center
    real(REALK) cc_wrt_bra(3)  !relative coordinates of center-of-charge w.r.t. bra center
    integer max_order_hbra     !current maximum order of HGTOs on bra center
    integer min_order_hket     !current minimum order of HGTOs on ket center
    integer dim_low_hbra       !dimensions of HGTOs on bra center of intermediate primitive integrals
    integer dim_up_hbra
    integer dim_cur_hket       !dimensions of HGTOs on ket cetner of intermediate primitive integrals
    integer dim_low_hket
    integer dim_up_hket
    integer num_low_derv       !number of lower order geometric derivatives
    integer num_up_derv        !number of higher order geometric derivatives
    integer strt_up_derv       !start address of upper order geometric derivatives
    integer end_up_derv        !end address of upper order geometric derivatives
    integer strt_low_derv      !start address of lower order geometric derivatives
    integer end_low_derv       !end address of lower order geometric derivatives
    integer mod2_recurr        !module of number of recurrence relations w.r.t. 2
    integer order_derv         !order of geometric derivatives
    real(REALK), allocatable :: int1_pints(:,:,:)  !intermediate primitive integrals
    real(REALK), allocatable :: int2_pints(:,:,:)
    integer strt_hbra_derv     !start address of geometric derivatives in returned integrals
    integer end_hbra_derv      !end address of geometric derivatives in returned integrals
    integer strt_hket_zero     !start address of HGTOs on ket center in integrals with zeroth
                               !order HGTO on bra center
    integer strt_cur_hket      !start address of HGTOs on ket center in current intermediate integrals
    integer strt_up_hket       !start address of HGTOs on ket center in upper order geometric derivatives
    integer strt_s_uket        !start addresses of HGTOs on ket center in integrals with zeroth
    integer strt_s_lket        !order HGTO on bra center, and upper and lower order geometric derivatives
    integer strt_addr_hbra     !start address of HGTOs on bra center in returned integrals
    integer strt_up_hbra       !start address of HGTOs on bra center in intermediate integrals
    integer ixyz               !incremental recorder over xyz components of coordinates
    integer ierr               !error information
#if defined(XTIME)
    real(REALK) curr_time      !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    select case(order_herm_bra(2))
    ! returns s-shell HGTO
    case(0)
      hbra_pints(1,:,:) = hket_pints
    ! maximum returned HGTOs are other shells, at least p-shell
    case default
      ! reciprocal of total exponent
      neg_hrp_expt = 1.0_REALK/(exponent_bra+exponent_ket)
      ! relative coordinates of center-of-charge w.r.t. ket center
      do ixyz = 1, 3
        cc_wrt_bra(ixyz) = neg_hrp_expt*(exponent_bra*coord_bra(ixyz) &
                         + exponent_ket*coord_ket(ixyz))-coord_bra(ixyz)
      end do
      neg_hrp_expt = -0.5_REALK*neg_hrp_expt
      ! ratio of exponent on ket center to that on bra center
      ket_to_bra = exponent_ket/exponent_bra
      ! \var(max_low_derv)-1 order geometric derivatives
      max_low_derv = orders_derv(2)+order_herm_bra(2)     !maximum order of geometric derivatives
      end_up_derv = dim_hket_derv                         !end address of upper order geometric derivatives
      num_up_derv = (max_low_derv+1)*(max_low_derv+2)/2   !number of upper order geometric derivatives
      end_low_derv = end_up_derv-num_up_derv              !end address of lower order geometric derivatives
      strt_up_derv = end_low_derv+1                       !start address of upper order geometric derivatives
      max_low_derv = max_low_derv-1                       !updates the maximum order of geometric derivatives
      num_low_derv = (max_low_derv+1)*(max_low_derv+2)/2  !number of lower order geometric derivatives
      strt_low_derv = end_low_derv-num_low_derv+1         !end address of \var(orders_derv(2)) order derivatives
#if defined(DEBUG)
      write(STDOUT,100) "DERV-s-HGTO/upper/start/end:", &
                        max_low_derv+1, strt_up_derv, end_up_derv
      write(STDOUT,100) "DERV-s-HGTO/upper/start/end:", &
                        max_low_derv, strt_low_derv, end_low_derv
#endif
      ! initializes current maximum order of HGTOs on bra center
      max_order_hbra = 0
      ! initializes the dimensions of HGTOs on bra center
      dim_up_hbra = 3
      ! sets current minimum order of HGTOs on ket center
      min_order_hket = order_herm_ket(1)
      ! sets the dimensions of HGTOs on ket center
      if (min_order_hket==0) then
        dim_cur_hket = (order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                     * (order_herm_ket(2)+3)/6
        dim_low_hket = dim_cur_hket
      else
        dim_cur_hket = ((order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                     *  (order_herm_ket(2)+3)                       &
                     - min_order_hket*(min_order_hket+1)*(min_order_hket+2))/6
        dim_low_hket = dim_cur_hket+min_order_hket*(min_order_hket+1)/2
      end if
      ! allocates faked integrals
      dim_low_hbra = 1
      dim_up_hket = dim_cur_hket
      allocate(int1_pints(dim_low_hbra,dim_up_hket,num_up_derv), stat=ierr)
      if (ierr/=0) stop "london_mom_cgto>> failed to allocate int1_pints/0!"
      ! allocates memory for the intermediate results
      allocate(int2_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
      if (ierr/=0) stop "london_mom_cgto>> failed to allocate int2_pints/0!"
      ! gets the intermediate integrals for order of HGTOs on bra center up to 
      ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(max_low_derv)
      call sub_london_mom_cgto(max_low_derv, order_herm_ket, min_order_hket, &
                           max_order_hbra, cc_wrt_bra, neg_hrp_expt,     &
                           ket_to_bra, dim_cur_hket, num_up_derv,        &
                           hket_pints(dim_hket-dim_cur_hket+1:dim_hket,  &
                                      strt_up_derv:end_up_derv),         &
                           dim_low_hket, num_low_derv,                   &
                           hket_pints(dim_hket-dim_low_hket+1:dim_hket,  &
                                      strt_low_derv:end_low_derv),       &
                           dim_low_hbra, dim_up_hket, int1_pints,        &
                           dim_up_hbra, int2_pints)
      deallocate(int1_pints)
      ! initializes the dimension of lower order HGTOs on bra center
      dim_low_hbra = 0
      ! loops over the orders of geometric derivatives not returned and
      ! the maximum returned order HGTOs on bra center, and current maximum
      ! order of HGTOs \var(max_order_hbra) needs to update each iteration
      mod2_recurr = mod(max_low_derv-orders_derv(2),2)
      do order_derv = max_low_derv-1, orders_derv(2)+mod2_recurr, -2
        ! (1) \var(order_derv) order geometric derivatives
        ! updates the numbers of lower and upper order geometric derivatives
        num_up_derv = num_low_derv
        num_low_derv = (order_derv+1)*(order_derv+2)/2
        ! updates the start and end addresses of lower and upper order geometric derivatives
        end_up_derv = end_low_derv
        strt_up_derv = strt_low_derv
        end_low_derv = strt_low_derv-1
        strt_low_derv = end_low_derv-num_low_derv+1
#if defined(DEBUG)
        write(STDOUT,100) "DERV-HGTO/loop/1/1/upper/start/end:", &
                          order_derv+1, strt_up_derv, end_up_derv
        write(STDOUT,100) "DERV-HGTO/loop/1/1/lower/start/end:", &
                          order_derv, strt_low_derv, end_low_derv
#endif
        ! updates current maximum order of HGTOs on bra center
        max_order_hbra = max_order_hbra+1
        ! updates the dimensions of HGTOs on bra center
        dim_low_hbra = dim_up_hbra
        dim_up_hbra = dim_low_hbra+(max_order_hbra+2)*(max_order_hbra+3)/2
        ! sets current minimum order of HGTOs on ket center
        min_order_hket = max(order_herm_ket(1)-max_order_hbra,0)
        ! sets the dimensions of HGTOs on ket center
        dim_up_hket = dim_cur_hket
        if (min_order_hket==0) then
          dim_cur_hket = (order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                       * (order_herm_ket(2)+3)/6
          dim_low_hket = dim_cur_hket
        else
          dim_cur_hket = ((order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                       *  (order_herm_ket(2)+3)                       &
                       - min_order_hket*(min_order_hket+1)*(min_order_hket+2))/6
          dim_low_hket = dim_cur_hket+min_order_hket*(min_order_hket+1)/2
        end if
        ! allocates memory for the intermediate results
        allocate(int1_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
        if (ierr/=0) stop "london_mom_cgto>> failed to allocate int1_pints/loop/1!"
        ! gets the intermediate integrals for order of HGTOs on bra center up to 
        ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_derv)
        call sub_london_mom_cgto(order_derv, order_herm_ket, min_order_hket,  &
                             max_order_hbra, cc_wrt_bra, neg_hrp_expt,    &
                             ket_to_bra, dim_cur_hket, num_up_derv,       &
                             hket_pints(dim_hket-dim_cur_hket+1:dim_hket, &
                                        strt_up_derv:end_up_derv),        &
                             dim_low_hket, num_low_derv,                  &
                             hket_pints(dim_hket-dim_low_hket+1:dim_hket, &
                                        strt_low_derv:end_low_derv),      &
                             dim_low_hbra, dim_up_hket, int2_pints,       &
                             dim_up_hbra, int1_pints)
        deallocate(int2_pints)
        ! (2) \var(order_derv)-1 order geometric derivatives
        ! updates the numbers of lower and upper order geometric derivatives
        num_up_derv = num_low_derv
        num_low_derv = order_derv*(order_derv+1)/2
        ! updates the start and end addresses of lower and upper order geometric derivatives
        end_up_derv = end_low_derv
        strt_up_derv = strt_low_derv
        end_low_derv = strt_low_derv-1
        strt_low_derv = end_low_derv-num_low_derv+1
#if defined(DEBUG)
        write(STDOUT,100) "DERV-HGTO/loop/1/2/upper/start/end:", &
                          order_derv, strt_up_derv, end_up_derv
        write(STDOUT,100) "DERV-HGTO/loop/1/2/lower/start/end:", &
                          order_derv-1, strt_low_derv, end_low_derv
#endif
        ! updates current maximum order of HGTOs on bra center
        max_order_hbra = max_order_hbra+1
        ! updates the dimensions of HGTOs on bra center
        dim_low_hbra = dim_up_hbra
        dim_up_hbra = dim_low_hbra+(max_order_hbra+2)*(max_order_hbra+3)/2
        ! sets current minimum order of HGTOs on ket center
        min_order_hket = max(order_herm_ket(1)-max_order_hbra,0)
        ! sets the dimensions of HGTOs on ket center
        dim_up_hket = dim_cur_hket
        if (min_order_hket==0) then
          dim_cur_hket = (order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                       * (order_herm_ket(2)+3)/6
          dim_low_hket = dim_cur_hket
        else
          dim_cur_hket = ((order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                       *  (order_herm_ket(2)+3)                       &
                       - min_order_hket*(min_order_hket+1)*(min_order_hket+2))/6
          dim_low_hket = dim_cur_hket+min_order_hket*(min_order_hket+1)/2
        end if
        ! allocates memory for the intermediate results
        allocate(int2_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
        if (ierr/=0) stop "london_mom_cgto>> failed to allocate int2_pints/loop/1!"
        ! gets the intermediate integrals for order of HGTOs on bra center up to 
        ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_derv)-1
        call sub_london_mom_cgto(order_derv-1, order_herm_ket, min_order_hket, &
                             max_order_hbra, cc_wrt_bra, neg_hrp_expt,     &
                             ket_to_bra, dim_cur_hket, num_up_derv,        &
                             hket_pints(dim_hket-dim_cur_hket+1:dim_hket,  &
                                        strt_up_derv:end_up_derv),         &
                             dim_low_hket, num_low_derv,                   &
                             hket_pints(dim_hket-dim_low_hket+1:dim_hket,  &
                                        strt_low_derv:end_low_derv),       &
                             dim_low_hbra, dim_up_hket, int1_pints,        &
                             dim_up_hbra, int2_pints)
        deallocate(int1_pints)
      end do
      ! we first need to perform recurrence relations for \var(orders_derv(2)) order
      ! geometric derivatives, current maximum order of HGTOs on bra center \var(max_order_hbra)
      ! still needs to update
      if (mod2_recurr==1) then
        ! updates the numbers of lower and upper order geometric derivatives
        num_up_derv = num_low_derv
        num_low_derv = (orders_derv(2)+1)*(orders_derv(2)+2)/2
        ! updates the start and end addresses of lower and upper order geometric derivatives
        end_up_derv = end_low_derv
        strt_up_derv = strt_low_derv
        end_low_derv = strt_low_derv-1
        strt_low_derv = end_low_derv-num_low_derv+1
#if defined(DEBUG)
        write(STDOUT,100) "DERV-HGTO/1/upper/start/end:", &
                          orders_derv(2)+1, strt_up_derv, end_up_derv
        write(STDOUT,100) "DERV-HGTO/1/lower/start/end:", &
                          orders_derv(2), strt_low_derv, end_low_derv
#endif
        ! updates current maximum order of HGTOs on bra center
        max_order_hbra = max_order_hbra+1
        ! updates the dimensions of HGTOs on bra center
        dim_low_hbra = dim_up_hbra
        dim_up_hbra = dim_low_hbra+(max_order_hbra+2)*(max_order_hbra+3)/2
        ! sets current minimum order of HGTOs on ket center
        min_order_hket = max(order_herm_ket(1)-max_order_hbra,0)
        ! sets the dimensions of HGTOs on ket center
        dim_up_hket = dim_cur_hket
        if (min_order_hket==0) then
          dim_cur_hket = (order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                       * (order_herm_ket(2)+3)/6
          dim_low_hket = dim_cur_hket
        else
          dim_cur_hket = ((order_herm_ket(2)+1)*(order_herm_ket(2)+2) &
                       *  (order_herm_ket(2)+3)                       &
                       - min_order_hket*(min_order_hket+1)*(min_order_hket+2))/6
          dim_low_hket = dim_cur_hket+min_order_hket*(min_order_hket+1)/2
        end if
        ! allocates memory for the intermediate results
        allocate(int1_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
        if (ierr/=0) stop "london_mom_cgto>> failed to allocate int1_pints/1!"
        ! gets the intermediate integrals for order of HGTOs on bra center up to 
        ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(orders_derv(2))
        call sub_london_mom_cgto(orders_derv(2), order_herm_ket, min_order_hket, &
                             max_order_hbra, cc_wrt_bra, neg_hrp_expt,       &
                             ket_to_bra, dim_cur_hket, num_up_derv,          &
                             hket_pints(dim_hket-dim_cur_hket+1:dim_hket,    &
                                        strt_up_derv:end_up_derv),           &
                             dim_low_hket, num_low_derv,                     &
                             hket_pints(dim_hket-dim_low_hket+1:dim_hket,    &
                                        strt_low_derv:end_low_derv),         &
                             dim_low_hbra, dim_up_hket, int2_pints,          &
                             dim_up_hbra, int1_pints)
        deallocate(int2_pints)
        ! initializes the end address of geometric derivatives in returned integrals
        end_hbra_derv = dim_hbra_derv
        ! sets the start address of geometric derivatives in returned integrals
        strt_hbra_derv = end_hbra_derv-num_low_derv+1
        ! checks if s-shell HGTO returned
        if (order_herm_bra(1)==0) then
          strt_hket_zero = dim_hket-dim_herm_ket+1
          hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
            = hket_pints(strt_hket_zero:dim_hket,strt_low_derv:end_low_derv)
          ! sets the start address of HGTOs in returned integrals
          strt_addr_hbra = 2
          ! sets the start address of HGTOs in intermediate integrals
          strt_up_hbra = dim_up_hbra-dim_herm_bra+2
        else
          ! sets the start address of HGTOs in returned integrals
          strt_addr_hbra = 1
          ! sets the start address of HGTOs in intermediate integrals
          strt_up_hbra = dim_up_hbra-dim_herm_bra+1
        end if
        ! assigns the returned integrals for other HGTOs and \var(orders_derv(2)) order
        ! geometric derivatives
        strt_cur_hket = dim_cur_hket-dim_herm_ket+1
        hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
          = int1_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
        ! updates the end address of geometric derivatives in returned integrals
        end_hbra_derv = strt_hbra_derv-1
        ! sets the start address of HGTOs on ket center in upper order geometric derivatives
        strt_up_hket = dim_cur_hket-dim_up_hket+1
        ! sets the start addresses of HGTOs on ket center in integrals with zeroth
        ! order HGTO on bra center, and upper and lower order geometric derivatives
        strt_s_uket = dim_hket-dim_cur_hket+1
        strt_s_lket = dim_hket-dim_low_hket+1
        ! loops over other returned orders of geometric derivatives, current maximum
        ! order of HGTOs on bra center \var(max_order_hbra) does not need to update
        mod2_recurr = mod(orders_derv(2)-orders_derv(1),2)
        do order_derv = orders_derv(2)-1, orders_derv(1)+mod2_recurr, -2
          ! (1) \var(order_derv) order geometric derivatives
          ! updates the numbers of lower and upper order geometric derivatives
          num_up_derv = num_low_derv
          num_low_derv = (order_derv+1)*(order_derv+2)/2
          ! updates the start and end addresses of lower and upper order geometric derivatives
          end_up_derv = end_low_derv
          strt_up_derv = strt_low_derv
          end_low_derv = strt_low_derv-1
          strt_low_derv = end_low_derv-num_low_derv+1
#if defined(DEBUG)
          write(STDOUT,100) "DERV-HGTO/loop/2/1/upper/start/end:", &
                            order_derv+1, strt_up_derv, end_up_derv
          write(STDOUT,100) "DERV-HGTO/loop/2/1/lower/start/end:", &
                            order_derv, strt_low_derv, end_low_derv
#endif
          ! allocates memory for the intermediate results
          allocate(int2_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
          if (ierr/=0) stop "london_mom_cgto>> failed to allocate int2_pints/loop/2!"
          ! gets the intermediate integrals for order of HGTOs on bra center up to 
          ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_derv)
          call sub_london_mom_cgto(order_derv, order_herm_ket, min_order_hket,  &
                               max_order_hbra, cc_wrt_bra, neg_hrp_expt,    &
                               ket_to_bra, dim_cur_hket, num_up_derv,       &
                               hket_pints(strt_s_uket:dim_hket,             &
                                          strt_up_derv:end_up_derv),        &
                               dim_low_hket, num_low_derv,                  &
                               hket_pints(strt_s_lket:dim_hket,             &
                                          strt_low_derv:end_low_derv),      &
                               dim_low_hbra, dim_up_hket,                   &
                               int1_pints(1:dim_low_hbra,                   &
                                          strt_up_hket:dim_cur_hket,:),     &
                               dim_up_hbra, int2_pints)
          deallocate(int1_pints)
          ! sets the start address of geometric derivatives in returned integrals
          strt_hbra_derv = end_hbra_derv-num_low_derv+1
          ! checks if s-shell HGTO returned
          if (order_herm_bra(1)==0) then
            hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
              = hket_pints(strt_hket_zero:dim_hket,strt_low_derv:end_low_derv)
          end if
          ! assigns the returned integrals for other HGTOs
          hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
            = int2_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
          ! updates the end address of geometric derivatives in returned integrals
          end_hbra_derv = strt_hbra_derv-1
          ! (2) \var(order_derv)-1 order geometric derivatives
          ! updates the numbers of lower and upper order geometric derivatives
          num_up_derv = num_low_derv
          num_low_derv = order_derv*(order_derv+1)/2
          ! updates the start and end addresses of lower and upper order geometric derivatives
          end_up_derv = end_low_derv
          strt_up_derv = strt_low_derv
          end_low_derv = strt_low_derv-1
          strt_low_derv = end_low_derv-num_low_derv+1
#if defined(DEBUG)
          write(STDOUT,100) "DERV-HGTO/loop/2/2/upper/start/end:", &
                            order_derv, strt_up_derv, end_up_derv
          write(STDOUT,100) "DERV-HGTO/loop/2/2/lower/start/end:", &
                            order_derv-1, strt_low_derv, end_low_derv
#endif
          ! allocates memory for the intermediate results
          allocate(int1_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
          if (ierr/=0) stop "london_mom_cgto>> failed to allocate int1_pints/loop/2!"
          ! gets the intermediate integrals for order of HGTOs on bra center up to 
          ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_derv)-1
          call sub_london_mom_cgto(order_derv-1, order_herm_ket, min_order_hket, &
                               max_order_hbra, cc_wrt_bra, neg_hrp_expt,     &
                               ket_to_bra, dim_cur_hket, num_up_derv,        &
                               hket_pints(strt_s_uket:dim_hket,              &
                                          strt_up_derv:end_up_derv),         &
                               dim_low_hket, num_low_derv,                   &
                               hket_pints(strt_s_lket:dim_hket,              &
                                          strt_low_derv:end_low_derv),       &
                               dim_low_hbra, dim_up_hket,                    &
                               int2_pints(1:dim_low_hbra,                    &
                                          strt_up_hket:dim_cur_hket,:),      &
                               dim_up_hbra, int1_pints)
          deallocate(int2_pints)
          ! sets the start address of geometric derivatives in returned integrals
          strt_hbra_derv = end_hbra_derv-num_low_derv+1
          ! checks if s-shell HGTO returned
          if (order_herm_bra(1)==0) then
            hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
              = hket_pints(strt_hket_zero:dim_hket,strt_low_derv:end_low_derv)
          end if
          ! assigns the returned integrals for other HGTOs
          hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
            = int1_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
          ! updates the end address of geometric derivatives in returned integrals
          end_hbra_derv = strt_hbra_derv-1
        end do
        ! performs recurrence relations for \var(orders_derv(1)) order geometric derivatives
        if (mod2_recurr==1) then
          ! updates the numbers of lower and upper order geometric derivatives
          num_up_derv = num_low_derv
          num_low_derv = (orders_derv(1)+1)*(orders_derv(1)+2)/2
          ! updates the start and end addresses of upper order geometric derivatives
          end_up_derv = end_low_derv
          strt_up_derv = strt_low_derv
#if defined(DEBUG)
          write(STDOUT,100) "DERV-HGTO/last/upper/start/end:", &
                            orders_derv(1)+1, strt_up_derv, end_up_derv
          write(STDOUT,100) "DERV-HGTO/last/lower/start/end:", &
                            orders_derv(1), 1, num_low_derv
#endif
          ! allocates memory for the intermediate results
          allocate(int2_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
          if (ierr/=0) stop "london_mom_cgto>> failed to allocate int2_pints/last!"
          ! gets the intermediate integrals for order of HGTOs on bra center up to 
          ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(orders_derv(1))
          call sub_london_mom_cgto(orders_derv(1), order_herm_ket, min_order_hket, &
                               max_order_hbra, cc_wrt_bra, neg_hrp_expt,       &
                               ket_to_bra, dim_cur_hket, num_up_derv,          &
                               hket_pints(strt_s_uket:dim_hket,                &
                                          strt_up_derv:end_up_derv),           &
                               dim_low_hket, num_low_derv,                     &
                               hket_pints(strt_s_lket:dim_hket,                &
                                          1:num_low_derv),                     &
                               dim_low_hbra, dim_up_hket,                      &
                               int1_pints(1:dim_low_hbra,                      &
                                          strt_up_hket:dim_cur_hket,:),        &
                               dim_up_hbra, int2_pints)
          deallocate(int1_pints)
          ! sets the start address of geometric derivatives in returned integrals
          strt_hbra_derv = end_hbra_derv-num_low_derv+1
          ! checks if s-shell HGTO returned
          if (order_herm_bra(1)==0) then
            hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
              = hket_pints(strt_hket_zero:dim_hket,1:num_low_derv)
          end if
          ! assigns the returned integrals for other HGTOs
          hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
            = int2_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
          deallocate(int2_pints)
        else
          deallocate(int1_pints)
        end if
      else
        ! initializes the end address of geometric derivatives in returned integrals
        end_hbra_derv = dim_hbra_derv
        ! sets the start address of geometric derivatives in returned integrals
        strt_hbra_derv = end_hbra_derv-num_low_derv+1
        ! checks if s-shell HGTO returned
        if (order_herm_bra(1)==0) then
          strt_hket_zero = dim_hket-dim_herm_ket+1
          hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
            = hket_pints(strt_hket_zero:dim_hket,strt_low_derv:end_low_derv)
          ! sets the start address of HGTOs in returned integrals
          strt_addr_hbra = 2
          ! sets the start address of HGTOs in intermediate integrals
          strt_up_hbra = dim_up_hbra-dim_herm_bra+2
        else
          ! sets the start address of HGTOs in returned integrals
          strt_addr_hbra = 1
          ! sets the start address of HGTOs in intermediate integrals
          strt_up_hbra = dim_up_hbra-dim_herm_bra+1
        end if
        ! assigns the returned integrals for other HGTOs and \var(orders_derv(2))
        ! order geometric derivatives
        strt_cur_hket = dim_cur_hket-dim_herm_ket+1
        hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
          = int2_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
        ! updates the end address of geometric derivatives in returned integrals
        end_hbra_derv = strt_hbra_derv-1
        ! sets the start address of HGTOs on ket center in upper order geometric derivatives
        strt_up_hket = dim_cur_hket-dim_up_hket+1
        ! sets the start addresses of HGTOs on ket center in integrals with zeroth
        ! order HGTO on bra center, and upper and lower order geometric derivatives
        strt_s_uket = dim_hket-dim_cur_hket+1
        strt_s_lket = dim_hket-dim_low_hket+1
        ! loops over other returned orders of geometric derivatives, current maximum
        ! order of HGTOs on bra center \var(max_order_hbra) does not need to update
        mod2_recurr = mod(orders_derv(2)-orders_derv(1),2)
        do order_derv = orders_derv(2)-1, orders_derv(1)+mod2_recurr, -2
          ! (1) \var(order_derv) order geometric derivatives
          ! updates the numbers of lower and upper order geometric derivatives
          num_up_derv = num_low_derv
          num_low_derv = (order_derv+1)*(order_derv+2)/2
          ! updates the start and end addresses of lower and upper order geometric derivatives
          end_up_derv = end_low_derv
          strt_up_derv = strt_low_derv
          end_low_derv = strt_low_derv-1
          strt_low_derv = end_low_derv-num_low_derv+1
#if defined(DEBUG)
          write(STDOUT,100) "DERV-HGTO/loop/3/1/upper/start/end:", &
                            order_derv+1, strt_up_derv, end_up_derv
          write(STDOUT,100) "DERV-HGTO/loop/3/1/lower/start/end:", &
                            order_derv, strt_low_derv, end_low_derv
#endif
          ! allocates memory for the intermediate results
          allocate(int1_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
          if (ierr/=0) stop "london_mom_cgto>> failed to allocate int1_pints/loop/3!"
          ! gets the intermediate integrals for order of HGTOs on bra center up to 
          ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_derv)
          call sub_london_mom_cgto(order_derv, order_herm_ket, min_order_hket,  &
                               max_order_hbra, cc_wrt_bra, neg_hrp_expt,    &
                               ket_to_bra, dim_cur_hket, num_up_derv,       &
                               hket_pints(strt_s_uket:dim_hket,             &
                                          strt_up_derv:end_up_derv),        &
                               dim_low_hket, num_low_derv,                  &
                               hket_pints(strt_s_lket:dim_hket,             &
                                          strt_low_derv:end_low_derv),      &
                               dim_low_hbra, dim_up_hket,                   &
                               int2_pints(1:dim_low_hbra,                   &
                                          strt_up_hket:dim_cur_hket,:),     &
                               dim_up_hbra, int1_pints)
          deallocate(int2_pints)
          ! sets the start address of geometric derivatives in returned integrals
          strt_hbra_derv = end_hbra_derv-num_low_derv+1
          ! checks if s-shell HGTO returned
          if (order_herm_bra(1)==0) then
            hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
              = hket_pints(strt_hket_zero:dim_hket,strt_low_derv:end_low_derv)
          end if
          ! assigns the returned integrals for other HGTOs
          hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
            = int1_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
          ! updates the end address of geometric derivatives in returned integrals
          end_hbra_derv = strt_hbra_derv-1
          ! (2) \var(order_derv)-1 order geometric derivatives
          ! updates the numbers of lower and upper order geometric derivatives
          num_up_derv = num_low_derv
          num_low_derv = order_derv*(order_derv+1)/2
          ! updates the start and end addresses of lower and upper order geometric derivatives
          end_up_derv = end_low_derv
          strt_up_derv = strt_low_derv
          end_low_derv = strt_low_derv-1
          strt_low_derv = end_low_derv-num_low_derv+1
#if defined(DEBUG)
          write(STDOUT,100) "DERV-HGTO/loop/3/2/upper/start/end:", &
                            order_derv, strt_up_derv, end_up_derv
          write(STDOUT,100) "DERV-HGTO/loop/3/2/lower/start/end:", &
                            order_derv-1, strt_low_derv, end_low_derv
#endif
          ! allocates memory for the intermediate results
          allocate(int2_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
          if (ierr/=0) stop "london_mom_cgto>> failed to allocate int2_pints/loop/3!"
          ! gets the intermediate integrals for order of HGTOs on bra center up to 
          ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(order_derv)-1
          call sub_london_mom_cgto(order_derv-1, order_herm_ket, min_order_hket, &
                               max_order_hbra, cc_wrt_bra, neg_hrp_expt,     &
                               ket_to_bra, dim_cur_hket, num_up_derv,        &
                               hket_pints(strt_s_uket:dim_hket,              &
                                          strt_up_derv:end_up_derv),         &
                               dim_low_hket, num_low_derv,                   &
                               hket_pints(strt_s_lket:dim_hket,              &
                                          strt_low_derv:end_low_derv),       &
                               dim_low_hbra, dim_up_hket,                    &
                               int1_pints(1:dim_low_hbra,                    &
                                          strt_up_hket:dim_cur_hket,:),      &
                               dim_up_hbra, int2_pints)
          deallocate(int1_pints)
          ! sets the start address of geometric derivatives in returned integrals
          strt_hbra_derv = end_hbra_derv-num_low_derv+1
          ! checks if s-shell HGTO returned
          if (order_herm_bra(1)==0) then
            hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
              = hket_pints(strt_hket_zero:dim_hket,strt_low_derv:end_low_derv)
          end if
          ! assigns the returned integrals for other HGTOs
          hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
            = int2_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
          ! updates the end address of geometric derivatives in returned integrals
          end_hbra_derv = strt_hbra_derv-1
        end do
        ! performs recurrence relations for \var(orders_derv(1)) order geometric derivatives
        if (mod2_recurr==1) then
          ! updates the numbers of lower and upper order geometric derivatives
          num_up_derv = num_low_derv
          num_low_derv = (orders_derv(1)+1)*(orders_derv(1)+2)/2
          ! updates the start and end addresses of upper order geometric derivatives
          end_up_derv = end_low_derv
          strt_up_derv = strt_low_derv
#if defined(DEBUG)
          write(STDOUT,100) "DERV-HGTO/last/upper/start/end:", &
                            orders_derv(1)+1, strt_up_derv, end_up_derv
          write(STDOUT,100) "DERV-HGTO/last/lower/start/end:", &
                            orders_derv(1), 1, num_low_derv
#endif
          ! allocates memory for the intermediate results
          allocate(int1_pints(dim_up_hbra,dim_cur_hket,num_low_derv), stat=ierr)
          if (ierr/=0) stop "london_mom_cgto>> failed to allocate int1_pints/last!"
          ! gets the intermediate integrals for order of HGTOs on bra center up to 
          ! \var(max_order_hbra)+1 and order of geometric derivatives as \var(orders_derv(1))
          call sub_london_mom_cgto(orders_derv(1), order_herm_ket, min_order_hket, &
                               max_order_hbra, cc_wrt_bra, neg_hrp_expt,       &
                               ket_to_bra, dim_cur_hket, num_up_derv,          &
                               hket_pints(strt_s_uket:dim_hket,                &
                                          strt_up_derv:end_up_derv),           &
                               dim_low_hket, num_low_derv,                     &
                               hket_pints(strt_s_lket:dim_hket,                &
                                          1:num_low_derv),                     &
                               dim_low_hbra, dim_up_hket,                      &
                               int2_pints(1:dim_low_hbra,                      &
                                          strt_up_hket:dim_cur_hket,:),        &
                               dim_up_hbra, int1_pints)
          deallocate(int2_pints)
          ! sets the start address of geometric derivatives in returned integrals
          strt_hbra_derv = end_hbra_derv-num_low_derv+1
          ! checks if s-shell HGTO returned
          if (order_herm_bra(1)==0) then
            hbra_pints(1,:,strt_hbra_derv:end_hbra_derv) &
              = hket_pints(strt_hket_zero:dim_hket,1:num_low_derv)
          end if
          ! assigns the returned integrals for other HGTOs
          hbra_pints(strt_addr_hbra:dim_herm_bra,:,strt_hbra_derv:end_hbra_derv) &
            = int1_pints(strt_up_hbra:dim_up_hbra,strt_cur_hket:dim_cur_hket,:)
          deallocate(int1_pints)
        else
          deallocate(int2_pints)
        end if
      end if
    end select
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "london_mom_cgto", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("london_mom_cgto>> ",A,I6,2I8)
110 format("london_mom_cgto>> ",I8,2X,A,I8)
#endif
  end subroutine london_mom_cgto

  !> \brief sub-recurrence relations by recovering upper order HGTOs on bra center
  !> \author Bin Gao
  !> \date 2011-10-24
  !> \param order_cur_derv is current order of geometric derivatives
  !> \param order_herm_ket contains the minimum and maximum orders HGTOs on ket center
  !> \param min_order_hket is current minimum order of HGTOs on ket center
  !> \param max_order_hbra is current maximum order of HGTOs on bra center
  !> \param cc_wrt_bra contains the relative coordinates of center-of-charge w.r.t. bra center
  !> \param neg_hrp_expt is the negative half reciprocal of the total exponent
  !> \param ket_to_bra is the ratio of exponent on ket center to that on bra center
  !> \param dim_cur_hket is the dimension of HGTOs on ket center for integrals with
  !>        upper order geometric derivatives and zeroth order HGTOs on bra center
  !> \param num_up_derv is the number of upper order geometric derivatives
  !> \param up_derv_sint contains the integrals with s-shell HGTO on bra center and
  !>        upper order geometric derivatives
  !> \param dim_low_hket is the dimension of HGTOs on ket center for integrals with
  !>        lower order geometric derivatives and zeroth order HGTOs on bra center
  !> \param num_low_derv is the number lower order geometric derivatives
  !> \param low_derv_sint contains the integrals with s-shell HGTO on bra center and
  !>        lower order geometric derivatives
  !> \param dim_low_hbra is the dimension of lower order HGTOs on bra center
  !> \param dim_up_hket is the dimension of HGTOs on ket center with upper order
  !>        geometric derivatives and lower order HGTOs on bra center
  !> \param up_derv_bint contains the integrals with lower order HGTOs on bra center and
  !>        upper order geometric derivatives
  !> \param dim_up_hbra is the dimension of upper order HGTOs on bra center
  !> \return low_derv_bint contains the integrals with upper order HGTOs on bra center and
  !>         lower order geometric derivatives
  subroutine sub_london_mom_cgto(order_cur_derv, order_herm_ket, min_order_hket, &
                             max_order_hbra, cc_wrt_bra, neg_hrp_expt,       &
                             ket_to_bra, dim_cur_hket, num_up_derv,          &
                             up_derv_sint, dim_low_hket, num_low_derv,       &
                             low_derv_sint, dim_low_hbra, dim_up_hket,       &
                             up_derv_bint, dim_up_hbra, low_derv_bint)
    use xkind
    implicit none
    integer, intent(in) :: order_cur_derv
    integer, intent(in) :: order_herm_ket(2)
    integer, intent(in) :: min_order_hket
    integer, intent(in) :: max_order_hbra
    real(REALK), intent(in) :: cc_wrt_bra(3)
    real(REALK), intent(in) :: neg_hrp_expt
    real(REALK), intent(in) :: ket_to_bra
    integer, intent(in) :: dim_cur_hket
    integer, intent(in) :: num_up_derv
    real(REALK), intent(in) :: up_derv_sint(dim_cur_hket,num_up_derv)
    integer, intent(in) :: dim_low_hket
    integer, intent(in) :: num_low_derv
    real(REALK), intent(in) :: low_derv_sint(dim_low_hket,num_low_derv)
    integer, intent(in) :: dim_low_hbra
    integer, intent(in) :: dim_up_hket
    real(REALK), intent(in) :: up_derv_bint(dim_low_hbra,dim_up_hket,num_up_derv)
    integer, intent(in) :: dim_up_hbra
    real(REALK), intent(inout) :: low_derv_bint(dim_up_hbra,dim_cur_hket,num_low_derv)
!f2py intent(in) :: order_cur_derv
!f2py intent(in) :: order_herm_ket
!f2py intent(in) :: min_order_hket
!f2py intent(in) :: max_order_hbra
!f2py intent(in) :: cc_wrt_bra
!f2py intent(in) :: neg_hrp_expt
!f2py intent(in) :: ket_to_bra
!f2py intent(hide) :: dim_cur_hket
!f2py intent(hide) :: num_up_derv
!f2py intent(in) :: up_derv_sint
!f2py intent(hide) :: dim_low_hket
!f2py intent(hide) :: num_low_derv
!f2py intent(in) :: low_derv_sint
!f2py intent(hide) :: dim_low_hbra
!f2py intent(hide) :: dim_up_hket
!f2py intent(in) :: up_derv_bint
!f2py depend(num_up_derv) :: up_derv_bint
!f2py intent(hide) :: dim_up_hbra
!f2py intent(inout) :: low_derv_bint
!f2py depend(dim_cur_hket) :: low_derv_bint
!f2py depend(num_low_derv) :: low_derv_bint
    integer addr_up_derv    !addresses of upper order geometric derivatives
    integer addr_up_derv_y
    integer addr_up_derv_z
    integer addr_cur_derv   !address of current order geometric derivatives
    integer iderv, jderv    !incremental recorders over geometric derivatives
    integer min_cur_hket    !current minimum order of HGTOs on ket center
    logical zero_cur_hket   !if the minimum order of HGTOs on ket center is zeroth order
    integer max_cur_hbra    !current maximum order of HGTOs on bra center for given order HGTOs on ket center
    integer base_low_hket   !base address of lower order HGTOs on ket center
    integer base_s_hket     !base address of lower order HGTOs on ket center in integrals with
                            !s-shell HGTO on bra center and lower order geometric derivatives
    integer addr_cur_hket   !address of current order HGTOs on ket center
    integer addr_s_hket     !address of current order HGTOs on ket center in integrals with
                            !s-shell HGTO on bra center and lower order geometric derivatives
    integer addr_b_hket     !address of current order HGTOs on ket center in integrals with
                            !lower order HGTOs on bra center and lower order geometric derivatives
    integer addr_low_xket   !addresses of lower order HGTOs on ket center
    integer addr_low_yket
    integer addr_low_zket
    integer addr_s_xket     !addresses of lower order HGTOs on ket center in integrals with
    integer addr_s_yket     !s-shell HGTO on bra center and lower order geometric derivatives
    integer addr_s_zket
    integer order_hket      !order of HGTOs on ket center
    integer order_xket      !orders of xyz components of HGTOs on ket center
    integer order_yket
    integer order_zket
    integer addr_up_hbra    !address of upper order HGTOs on bra center
    integer addr_cur_hbra   !address of current order HGTOs on bra center
    integer addr_low_hbra   !address of lower order HGTOs on bra center
    integer order_hbra      !order of HGTOs on bra center
    integer ibra, jbra      !incremental recorder over HGTOs on bra center
#if defined(XTIME)
    real(REALK) curr_time   !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    if (min_order_hket==0) then
      zero_cur_hket = .true.
      min_cur_hket = min_order_hket+1
    else
      zero_cur_hket = .false.
      min_cur_hket = min_order_hket
    end if
    addr_up_derv = 0
    addr_cur_derv = 0
    ! loops over xyz components of geometric derivatives
    do iderv = order_cur_derv, 0, -1
      do jderv = 0, iderv
        addr_up_derv = addr_up_derv+1
        addr_cur_derv = addr_cur_derv+1
        addr_cur_hket = 0
        addr_b_hket = 0
        ! zeroth order HGTOs on ket center
        if (zero_cur_hket) then
          addr_cur_hket = addr_cur_hket+1
          addr_s_hket = 1
          ! sets current maximum order of HGTOs on bra center for zeroth order HGTOs on ket center
          max_cur_hbra = max_order_hbra-order_herm_ket(1)
          ! px HGTO on bra center
          low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv)
          ! py HGTO on bra center
          addr_up_derv_y = addr_up_derv+1
          low_derv_bint(2,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_y)
          ! pz HGTO on bra center
          addr_up_derv_z = addr_up_derv+iderv+2
          low_derv_bint(3,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_z)
          ! d-shell on bra center
          if (max_cur_hbra>0) then
            addr_b_hket = addr_b_hket+1
            ! dxx
            low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(1,addr_b_hket,addr_up_derv))
            ! dxy
            low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_y)
            ! dyy
            low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(2,addr_b_hket,addr_up_derv_y))
            ! dxz
            low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_z)
            ! dyz
            low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(2,addr_b_hket,addr_up_derv_z)
            ! dzz
            low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(3,addr_b_hket,addr_up_derv_z))
          end if
          if (max_cur_hbra>1) then
            addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
            addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
            addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
            ! loops over other current order of HGTOs on bra center, starting from d-shell
            do order_hbra = 2, max_cur_hbra
              addr_up_hbra = addr_up_hbra+1
              addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
              write(STDOUT,100) "orders:", order_hbra, order_hket, order_cur_derv
              write(STDOUT,100) "x-direction"
              write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
              write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv
              write(STDOUT,100) "------------------------------------"
#endif
              ! x-direction
              low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                = cc_wrt_bra(1)                                              &
                * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv))
              ! y-direction
              addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
              write(STDOUT,100) "y-direction"
              write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
              write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
              write(STDOUT,100) "------------------------------------"
#endif
              low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                = cc_wrt_bra(2)                                            &
                * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt                                             &
                * up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1)
              do ibra = 1, order_hbra
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                write(STDOUT,100) "------------------------------------"
#endif
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                  = cc_wrt_bra(2)                                                 &
                  * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                  * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                  + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1))
              end do
              ! z-direction
#if defined(DEBUG)
              write(STDOUT,100) "z-direction"
#endif
              addr_cur_hbra = addr_cur_hbra-1
              do jbra = 0, order_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                write(STDOUT,100) "------------------------------------"
#endif
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                  = cc_wrt_bra(3)                                            &
                  * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt                                             &
                  * up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2)
              end do
              do ibra = 1, order_hbra
                do jbra = 0, order_hbra-ibra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
                  addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                    = cc_wrt_bra(3)                                            &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                &
                    * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv) &
                    + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2))
                end do
              end do
            end do
          end if
        else
          addr_s_hket = min_order_hket*(min_order_hket+1)/2
        end if
        base_low_hket = 0
        base_s_hket = 0
        ! loops over other current order HGTOs on ket center
        do order_hket = min_cur_hket, order_herm_ket(2)
          addr_cur_hket = addr_cur_hket+1
          addr_s_hket = addr_s_hket+1
          addr_s_xket = base_s_hket+1
          ! sets current maximum order of HGTOs on bra center for order \var(order_hket)
          ! HGTOs on ket center
          max_cur_hbra = min(max_order_hbra,max_order_hbra-order_herm_ket(1)+order_hket)
          ! (1) x...x component of upper order HGTOs on ket center
          !
          ! px HGTO on bra center
          low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv) &
            - real(order_hket,REALK)*low_derv_sint(addr_s_xket,addr_cur_derv))
          ! py HGTO on bra center
          addr_up_derv_y = addr_up_derv+1
          low_derv_bint(2,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_y)
          ! pz HGTO on bra center
          addr_up_derv_z = addr_up_derv+iderv+2
          low_derv_bint(3,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_z)
          ! d-shell on bra center
          if (max_cur_hbra>0) then
            addr_low_xket = base_low_hket+1
            addr_b_hket = addr_b_hket+1
            ! dxx
            low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(1,addr_b_hket,addr_up_derv)                          &
              - real(order_hket,REALK)*low_derv_bint(1,addr_low_xket,addr_cur_derv))
            ! dxy
            low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_y)
            ! dyy
            low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(2,addr_b_hket,addr_up_derv_y))
            ! dxz
            low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_z)
            ! dyz
            low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(2,addr_b_hket,addr_up_derv_z)
            ! dzz
            low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(3,addr_b_hket,addr_up_derv_z))
            if (max_cur_hbra>1) then
              addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
              addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
              addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
              ! loops over other current order of HGTOs on bra center, starting from d-shell
              do order_hbra = 2, max_cur_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "orders:", order_hbra, order_hket, order_cur_derv
                write(STDOUT,100) "x-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv
                write(STDOUT,110) addr_up_hbra, addr_low_xket, addr_cur_derv
                write(STDOUT,100) "------------------------------------"
#endif
                ! x-direction
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                  = cc_wrt_bra(1)                                              &
                  * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                  + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                  * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                  + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv)       &
                  - real(order_hket,REALK)                                     &
                  * low_derv_bint(addr_cur_hbra,addr_low_xket,addr_cur_derv))
                ! y-direction
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "y-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
                write(STDOUT,100) "------------------------------------"
#endif
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                  = cc_wrt_bra(2)                                            &
                  * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt                                             &
                  * up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1)
                do ibra = 1, order_hbra
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                    = cc_wrt_bra(2)                                                 &
                    * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                    * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                    + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1))
                end do
                ! z-direction
#if defined(DEBUG)
                write(STDOUT,100) "z-direction"
#endif
                addr_cur_hbra = addr_cur_hbra-1
                do jbra = 0, order_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                    = cc_wrt_bra(3)                                            &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt                                             &
                    * up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2)
                end do
                do ibra = 1, order_hbra
                  do jbra = 0, order_hbra-ibra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
                    addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                      = cc_wrt_bra(3)                                            &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                      + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                &
                      * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv) &
                      + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2))
                  end do
                end do
              end do
            end if
          end if
          ! (2) x...xy to xy...y components of upper order HGTOs on ket center
          addr_low_yket = base_low_hket
          addr_s_yket = base_s_hket
          do order_yket = 1, order_hket-1
            addr_cur_hket = addr_cur_hket+1
            addr_s_hket = addr_s_hket+1
            addr_s_xket = addr_s_xket+1
            addr_s_yket = addr_s_yket+1
            ! sets the order along x-direction
            order_xket = order_hket-order_yket
            ! px HGTO on bra center
            low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
              = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv) &
              - real(order_xket,REALK)*low_derv_sint(addr_s_xket,addr_cur_derv))
            ! py HGTO on bra center
            addr_up_derv_y = addr_up_derv+1
            low_derv_bint(2,addr_cur_hket,addr_cur_derv)                 &
              = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
              + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_y) &
              - real(order_yket,REALK)*low_derv_sint(addr_s_yket,addr_cur_derv))
            ! pz HGTO on bra center
            addr_up_derv_z = addr_up_derv+iderv+2
            low_derv_bint(3,addr_cur_hket,addr_cur_derv)               &
              = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_z)
            ! d-shell on bra center
            if (max_cur_hbra>0) then
              addr_low_xket = addr_low_xket+1
              addr_low_yket = addr_low_yket+1
              addr_b_hket = addr_b_hket+1
              ! dxx
              low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(1,addr_b_hket,addr_up_derv)                          &
                - real(order_xket,REALK)*low_derv_bint(1,addr_low_xket,addr_cur_derv))
              ! dxy
              low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_y)   &
                - real(order_yket,REALK)*low_derv_bint(1,addr_low_yket,addr_cur_derv))
              ! dyy
              low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(2,addr_b_hket,addr_up_derv_y)                        &
                - real(order_yket,REALK)*low_derv_bint(2,addr_low_yket,addr_cur_derv))
              ! dxz
              low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_z)
              ! dyz
              low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*up_derv_bint(2,addr_b_hket,addr_up_derv_z)
              ! dzz
              low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(3,addr_b_hket,addr_up_derv_z))
              if (max_cur_hbra>1) then
                addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                ! loops over other current order of HGTOs on bra center, starting from d-shell
                do order_hbra = 2, max_cur_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "orders:", order_hbra, order_xket, order_cur_derv
                  write(STDOUT,100) "x-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv
                  write(STDOUT,110) addr_cur_hbra, addr_low_xket, addr_cur_derv
                  write(STDOUT,100) "------------------------------------"
#endif
                  ! x-direction
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                    = cc_wrt_bra(1)                                              &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                    + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                    * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                    + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv)       &
                    - real(order_xket,REALK)                                     &
                    * low_derv_bint(addr_cur_hbra,addr_low_xket,addr_cur_derv))
                  ! y-direction
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "y-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
                  write(STDOUT,110) addr_up_hbra, addr_low_yket, addr_cur_derv
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                    = cc_wrt_bra(2)                                            &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt                                             &
                    * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1)  &
                    - real(order_yket,REALK)                                   &
                    * low_derv_bint(addr_cur_hbra,addr_low_yket,addr_cur_derv))
                  do ibra = 1, order_hbra
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_derv
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                      = cc_wrt_bra(2)                                                 &
                      * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                      + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                      * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                      + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1)   &
                      - real(order_yket,REALK)                                        &
                      * low_derv_bint(addr_cur_hbra+ibra,addr_low_yket,addr_cur_derv))
                  end do
                  ! z-direction
#if defined(DEBUG)
                  write(STDOUT,100) "z-direction"
#endif
                  addr_cur_hbra = addr_cur_hbra-1
                  do jbra = 0, order_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                      = cc_wrt_bra(3)                                            &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                      + neg_hrp_expt                                             &
                      * up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2)
                  end do
                  do ibra = 1, order_hbra
                    do jbra = 0, order_hbra-ibra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
                      addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                      write(STDOUT,100) "------------------------------------"
#endif
                      low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                        = cc_wrt_bra(3)                                            &
                        * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                        + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                &
                        * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv) &
                        + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2))
                    end do
                  end do
                end do
              end if
            end if
          end do
          ! (3) y...y component of upper order HGTOs on ket center
          addr_cur_hket = addr_cur_hket+1
          addr_s_hket = addr_s_hket+1
          addr_s_yket = addr_s_yket+1
          ! px HGTO on bra center
          low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv)
          ! py HGTO on bra center
          addr_up_derv_y = addr_up_derv+1
          low_derv_bint(2,addr_cur_hket,addr_cur_derv)                 &
            = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
            + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_y) &
            - real(order_hket,REALK)*low_derv_sint(addr_s_yket,addr_cur_derv))
          ! pz HGTO on bra center
          addr_up_derv_z = addr_up_derv+iderv+2
          low_derv_bint(3,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_z)
          ! d-shell on bra center
          if (max_cur_hbra>0) then
            addr_low_yket = addr_low_yket+1
            addr_b_hket = addr_b_hket+1
            ! dxx
            low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(1,addr_b_hket,addr_up_derv))
            ! dxy
            low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_y)   &
              - real(order_hket,REALK)*low_derv_bint(1,addr_low_yket,addr_cur_derv))
            ! dyy
            low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(2,addr_b_hket,addr_up_derv_y)                        &
              - real(order_hket,REALK)*low_derv_bint(2,addr_low_yket,addr_cur_derv))
            ! dxz
            low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_z)
            ! dyz
            low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(2,addr_b_hket,addr_up_derv_z)
            ! dzz
            low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(3,addr_b_hket,addr_up_derv_z))
            if (max_cur_hbra>1) then
              addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
              addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
              addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
              ! loops over other current order of HGTOs on bra center, starting from d-shell
              do order_hbra = 2, max_cur_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "orders:", order_hbra, order_xket, order_cur_derv
                write(STDOUT,100) "x-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv
                write(STDOUT,100) "------------------------------------"
#endif
                ! x-direction
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                  = cc_wrt_bra(1)                                              &
                  * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                  + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                  * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                  + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv))
                ! y-direction
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "y-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
                write(STDOUT,110) addr_cur_hbra, addr_low_yket, addr_cur_derv
                write(STDOUT,100) "------------------------------------"
#endif
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                  = cc_wrt_bra(2)                                            &
                  * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt                                             &
                  * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1)  &
                  - real(order_hket,REALK)                                   &
                  * low_derv_bint(addr_cur_hbra,addr_low_yket,addr_cur_derv))
                do ibra = 1, order_hbra
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_derv
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                    = cc_wrt_bra(2)                                                 &
                    * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                    * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                    + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1)   &
                    - real(order_hket,REALK)                                        &
                    * low_derv_bint(addr_cur_hbra+ibra,addr_low_yket,addr_cur_derv))
                end do
                ! z-direction
#if defined(DEBUG)
                write(STDOUT,100) "z-direction"
#endif
                addr_cur_hbra = addr_cur_hbra-1
                do jbra = 0, order_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                    = cc_wrt_bra(3)                                            &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt                                             &
                    * up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2)
                end do
                do ibra = 1, order_hbra
                  do jbra = 0, order_hbra-ibra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
                    addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                      = cc_wrt_bra(3)                                            &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                      + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                &
                      * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv) &
                      + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2))
                  end do
                end do
              end do
            end if
          end if
          ! (4) x...xz to yz...z components of upper order HGTOs on ket center
          addr_low_xket = base_low_hket+order_hket
          addr_low_zket = base_low_hket
          addr_s_xket = base_s_hket+order_hket
          addr_s_zket = base_s_hket
          do order_zket = 1, order_hket-1
            addr_cur_hket = addr_cur_hket+1
            addr_s_hket = addr_s_hket+1
            addr_s_xket = addr_s_xket+1
            addr_s_zket = addr_s_zket+1
            ! (4.1) x...xz...z component
            order_xket = order_hket-order_zket
            ! px HGTO on bra center
            low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
              = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv) &
              - real(order_xket,REALK)*low_derv_sint(addr_s_xket,addr_cur_derv))
            ! py HGTO on bra center
            addr_up_derv_y = addr_up_derv+1
            low_derv_bint(2,addr_cur_hket,addr_cur_derv)               &
              = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_y)
            ! pz HGTO on bra center
            addr_up_derv_z = addr_up_derv+iderv+2
            low_derv_bint(3,addr_cur_hket,addr_cur_derv)                 &
              = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
              + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_z) &
              - real(order_zket,REALK)*low_derv_sint(addr_s_zket,addr_cur_derv))
            ! d-shell on bra center
            if (max_cur_hbra>0) then
              addr_low_xket = addr_low_xket+1
              addr_low_zket = addr_low_zket+1
              addr_b_hket = addr_b_hket+1
              ! dxx
              low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(1,addr_b_hket,addr_up_derv)                          &
                - real(order_xket,REALK)*low_derv_bint(1,addr_low_xket,addr_cur_derv))
              ! dxy
              low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_y)
              ! dyy
              low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(2,addr_b_hket,addr_up_derv_y))
              ! dxz
              low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_z)   &
                - real(order_zket,REALK)*low_derv_bint(1,addr_low_zket,addr_cur_derv))
              ! dyz
              low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*(up_derv_bint(2,addr_b_hket,addr_up_derv_z)   &
                - real(order_zket,REALK)*low_derv_bint(2,addr_low_zket,addr_cur_derv))
              ! dzz
              low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(3,addr_b_hket,addr_up_derv_z)                        &
                - real(order_zket,REALK)*low_derv_bint(3,addr_low_zket,addr_cur_derv))
              if (max_cur_hbra>1) then
                addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                ! loops over other current order of HGTOs on bra center, starting from d-shell
                do order_hbra = 2, max_cur_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "orders:", order_hbra, order_xket, order_cur_derv
                  write(STDOUT,100) "x-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_low_xket, addr_up_derv
                  write(STDOUT,100) "------------------------------------"
#endif
                  ! x-direction
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                    = cc_wrt_bra(1)                                              &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                    + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                    * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                    + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv)       &
                    - real(order_xket,REALK)                                     &
                    * low_derv_bint(addr_cur_hbra,addr_low_xket,addr_cur_derv))
                  ! y-direction
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "y-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                    = cc_wrt_bra(2)                                            &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt                                             &
                    * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1))
                  do ibra = 1, order_hbra
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                      = cc_wrt_bra(2)                                                 &
                      * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                      + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                      * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                      + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1))
                  end do
                  ! z-direction
#if defined(DEBUG)
                  write(STDOUT,100) "z-direction"
#endif
                  addr_cur_hbra = addr_cur_hbra-1
                  do jbra = 0, order_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                    write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                      = cc_wrt_bra(3)                                                 &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)      &
                      + neg_hrp_expt                                                  &
                      * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                      - real(order_zket,REALK)                                        &
                      * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                  end do
                  do ibra = 1, order_hbra
                    do jbra = 0, order_hbra-ibra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
                      addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                      write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                      write(STDOUT,100) "------------------------------------"
#endif
                      low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)          &
                        = cc_wrt_bra(3)                                                &
                        * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)     &
                        + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                    &
                        * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv)     &
                        + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                        - real(order_zket,REALK)                                       &
                        * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                    end do
                  end do
                end do
              end if
            end if
            ! (4.2) x...xyz...z to xy...yz...z components of upper order HGTOs on ket center
            do order_yket = 1, order_hket-(order_zket+1)
              addr_cur_hket = addr_cur_hket+1
              addr_s_hket = addr_s_hket+1
              addr_s_xket = addr_s_xket+1
              addr_s_yket = addr_s_yket+1
              addr_s_zket = addr_s_zket+1
              order_xket = order_xket-1
              ! px HGTO on bra center
              low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
                = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv) &
                - real(order_xket,REALK)*low_derv_sint(addr_s_xket,addr_cur_derv))
              ! py HGTO on bra center
              addr_up_derv_y = addr_up_derv+1
              low_derv_bint(2,addr_cur_hket,addr_cur_derv)                 &
                = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
                + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_y) &
                - real(order_yket,REALK)*low_derv_sint(addr_s_yket,addr_cur_derv))
              ! pz HGTO on bra center
              addr_up_derv_z = addr_up_derv+iderv+2
              low_derv_bint(3,addr_cur_hket,addr_cur_derv)                 &
                = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
                + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_z) &
                - real(order_zket,REALK)*low_derv_sint(addr_s_zket,addr_cur_derv))
              ! d-shell on bra center
              if (max_cur_hbra>0) then
                addr_low_xket = addr_low_xket+1
                addr_low_yket = addr_low_yket+1
                addr_low_zket = addr_low_zket+1
                addr_b_hket = addr_b_hket+1
                ! dxx
                low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
                  = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
                  + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                  + up_derv_bint(1,addr_b_hket,addr_up_derv)                          &
                  - real(order_xket,REALK)*low_derv_bint(1,addr_low_xket,addr_cur_derv))
                ! dxy
                low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
                  = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_y)   &
                  - real(order_yket,REALK)*low_derv_bint(1,addr_low_yket,addr_cur_derv))
                ! dyy
                low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
                  = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
                  + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                  + up_derv_bint(2,addr_b_hket,addr_up_derv_y)                        &
                  - real(order_yket,REALK)*low_derv_bint(2,addr_low_yket,addr_cur_derv))
                ! dxz
                low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
                  = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_z)   &
                  - real(order_zket,REALK)*low_derv_bint(1,addr_low_zket,addr_cur_derv))
                ! dyz
                low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
                  = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt*(up_derv_bint(2,addr_b_hket,addr_up_derv_z)   &
                  - real(order_zket,REALK)*low_derv_bint(2,addr_low_zket,addr_cur_derv))
                ! dzz
                low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
                  = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
                  + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                  + up_derv_bint(3,addr_b_hket,addr_up_derv_z)                        &
                  - real(order_zket,REALK)*low_derv_bint(3,addr_low_zket,addr_cur_derv))
                if (max_cur_hbra>1) then
                  addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                  addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                  addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                  ! loops over other current order of HGTOs on bra center, starting from d-shell
                  do order_hbra = 2, max_cur_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,100) "orders:", order_hbra, order_xket, order_cur_derv
                    write(STDOUT,100) "x-direction"
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_low_xket, addr_up_derv
                    write(STDOUT,100) "------------------------------------"
#endif
                    ! x-direction
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                      = cc_wrt_bra(1)                                              &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                      + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                      * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                      + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv)       &
                      - real(order_xket,REALK)                                     &
                      * low_derv_bint(addr_cur_hbra,addr_low_xket,addr_cur_derv))
                    ! y-direction
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,100) "y-direction"
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
                    write(STDOUT,110) addr_cur_hbra, addr_low_yket, addr_cur_derv
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                      = cc_wrt_bra(2)                                            &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                      + neg_hrp_expt                                             &
                      * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1)  &
                      - real(order_yket,REALK)                                   &
                      * low_derv_bint(addr_cur_hbra,addr_low_yket,addr_cur_derv))
                    do ibra = 1, order_hbra
                      addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                      write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_derv
                      write(STDOUT,100) "------------------------------------"
#endif
                      low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                        = cc_wrt_bra(2)                                                 &
                        * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                        + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                        * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                        + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1)   &
                        - real(order_yket,REALK)                                        &
                        * low_derv_bint(addr_cur_hbra+ibra,addr_low_yket,addr_cur_derv))
                    end do
                    ! z-direction
#if defined(DEBUG)
                    write(STDOUT,100) "z-direction"
#endif
                    addr_cur_hbra = addr_cur_hbra-1
                    do jbra = 0, order_hbra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                      write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                      write(STDOUT,100) "------------------------------------"
#endif
                      low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                        = cc_wrt_bra(3)                                                 &
                        * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)      &
                        + neg_hrp_expt                                                  &
                        * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                        - real(order_zket,REALK)                                        &
                        * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                    end do
                    do ibra = 1, order_hbra
                      do jbra = 0, order_hbra-ibra
                        addr_up_hbra = addr_up_hbra+1
                        addr_cur_hbra = addr_cur_hbra+1
                        addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                        write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                        write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                        write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                        write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                        write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                        write(STDOUT,100) "------------------------------------"
#endif
                        low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)          &
                          = cc_wrt_bra(3)                                                &
                          * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)     &
                          + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                    &
                          * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv)     &
                          + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                          - real(order_zket,REALK)                                       &
                          * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                      end do
                    end do
                  end do
                end if
              end if
            end do
            ! (4.3) y...yz...z component of upper order HGTOs on ket center
            addr_cur_hket = addr_cur_hket+1
            addr_s_hket = addr_s_hket+1
            addr_s_yket = addr_s_yket+1
            addr_s_zket = addr_s_zket+1
            order_yket = order_hket-order_zket
            ! px HGTO on bra center
            low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
              = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv)
            ! py HGTO on bra center
            addr_up_derv_y = addr_up_derv+1
            low_derv_bint(2,addr_cur_hket,addr_cur_derv)                 &
              = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
              + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_y) &
              - real(order_yket,REALK)*low_derv_sint(addr_s_yket,addr_cur_derv))
            ! pz HGTO on bra center
            addr_up_derv_z = addr_up_derv+iderv+2
            low_derv_bint(3,addr_cur_hket,addr_cur_derv)                 &
              = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
              + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_z) &
              - real(order_zket,REALK)*low_derv_sint(addr_s_zket,addr_cur_derv))
            ! d-shell on bra center
            if (max_cur_hbra>0) then
              addr_low_yket = addr_low_yket+1
              addr_low_zket = addr_low_zket+1
              addr_b_hket = addr_b_hket+1
              ! dxx
              low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(1,addr_b_hket,addr_up_derv))
              ! dxy
              low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_y)   &
                - real(order_yket,REALK)*low_derv_bint(1,addr_low_yket,addr_cur_derv))
              ! dyy
              low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(2,addr_b_hket,addr_up_derv_y)                        &
                - real(order_yket,REALK)*low_derv_bint(2,addr_low_yket,addr_cur_derv))
              ! dxz
              low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_z)   &
                - real(order_zket,REALK)*low_derv_bint(1,addr_low_zket,addr_cur_derv))
              ! dyz
              low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
                = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
                + neg_hrp_expt*(up_derv_bint(2,addr_b_hket,addr_up_derv_z)   &
                - real(order_zket,REALK)*low_derv_bint(2,addr_low_zket,addr_cur_derv))
              ! dzz
              low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
                = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
                + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
                + up_derv_bint(3,addr_b_hket,addr_up_derv_z)                        &
                - real(order_zket,REALK)*low_derv_bint(3,addr_low_zket,addr_cur_derv))
              if (max_cur_hbra>1) then
                addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
                addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
                addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
                ! loops over other current order of HGTOs on bra center, starting from d-shell
                do order_hbra = 2, max_cur_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "orders:", order_hbra, order_xket, order_cur_derv
                  write(STDOUT,100) "x-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,100) "------------------------------------"
#endif
                  ! x-direction
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                    = cc_wrt_bra(1)                                              &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                    + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                    * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                    + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv))
                  ! y-direction
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,100) "y-direction"
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
                  write(STDOUT,110) addr_cur_hbra, addr_low_yket, addr_cur_derv
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                    = cc_wrt_bra(2)                                            &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt                                             &
                    * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1)  &
                    - real(order_yket,REALK)                                   &
                    * low_derv_bint(addr_cur_hbra,addr_low_yket,addr_cur_derv))
                  do ibra = 1, order_hbra
                    addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                    write(STDOUT,110) addr_cur_hbra+ibra, addr_low_yket, addr_cur_derv
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                      = cc_wrt_bra(2)                                                 &
                      * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                      + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                      * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                      + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1)   &
                      - real(order_yket,REALK)                                        &
                      * low_derv_bint(addr_cur_hbra+ibra,addr_low_yket,addr_cur_derv))
                  end do
                  ! z-direction
#if defined(DEBUG)
                  write(STDOUT,100) "z-direction"
#endif
                  addr_cur_hbra = addr_cur_hbra-1
                  do jbra = 0, order_hbra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                    write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                      = cc_wrt_bra(3)                                                 &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)      &
                      + neg_hrp_expt                                                  &
                      * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                      - real(order_zket,REALK)                                        &
                      * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                  end do
                  do ibra = 1, order_hbra
                    do jbra = 0, order_hbra-ibra
                      addr_up_hbra = addr_up_hbra+1
                      addr_cur_hbra = addr_cur_hbra+1
                      addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                      write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                      write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                      write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                      write(STDOUT,100) "------------------------------------"
#endif
                      low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)          &
                        = cc_wrt_bra(3)                                                &
                        * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)     &
                        + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                    &
                        * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv)     &
                        + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                        - real(order_zket,REALK)                                       &
                        * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                    end do
                  end do
                end do
              end if
            end if
          end do
          ! (5) z...z component of upper order HGTOs on ket center
          addr_cur_hket = addr_cur_hket+1
          addr_s_hket = addr_s_hket+1
          addr_s_zket = addr_s_zket+1
          ! px HGTO on bra center
          low_derv_bint(1,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(1)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv)
          ! py HGTO on bra center
          addr_up_derv_y = addr_up_derv+1
          low_derv_bint(2,addr_cur_hket,addr_cur_derv)               &
            = cc_wrt_bra(2)*low_derv_sint(addr_s_hket,addr_cur_derv) &
            + neg_hrp_expt*up_derv_sint(addr_cur_hket,addr_up_derv_y)
          ! pz HGTO on bra center
          addr_up_derv_z = addr_up_derv+iderv+2
          low_derv_bint(3,addr_cur_hket,addr_cur_derv)                 &
            = cc_wrt_bra(3)*low_derv_sint(addr_s_hket,addr_cur_derv)   &
            + neg_hrp_expt*(up_derv_sint(addr_cur_hket,addr_up_derv_z) &
            - real(order_hket,REALK)*low_derv_sint(addr_s_zket,addr_cur_derv))
          ! d-shell on bra center
          if (max_cur_hbra>0) then
            addr_low_zket = addr_low_zket+1
            addr_b_hket = addr_b_hket+1
            ! dxx
            low_derv_bint(4,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(1)*low_derv_bint(1,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(1,addr_b_hket,addr_up_derv))
            ! dxy
            low_derv_bint(5,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(2)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*up_derv_bint(1,addr_b_hket,addr_up_derv_y)
            ! dyy
            low_derv_bint(6,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(2)*low_derv_bint(2,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(2,addr_b_hket,addr_up_derv_y))
            ! dxz
            low_derv_bint(7,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(1,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*(up_derv_bint(1,addr_b_hket,addr_up_derv_z)   &
              - real(order_hket,REALK)*low_derv_bint(1,addr_low_zket,addr_cur_derv))
            ! dyz
            low_derv_bint(8,addr_cur_hket,addr_cur_derv)                   &
              = cc_wrt_bra(3)*low_derv_bint(2,addr_cur_hket,addr_cur_derv) &
              + neg_hrp_expt*(up_derv_bint(2,addr_b_hket,addr_up_derv_z)   &
              - real(order_hket,REALK)*low_derv_bint(2,addr_low_zket,addr_cur_derv))
            ! dzz
            low_derv_bint(9,addr_cur_hket,addr_cur_derv)                          &
              = cc_wrt_bra(3)*low_derv_bint(3,addr_cur_hket,addr_cur_derv)        &
              + neg_hrp_expt*(ket_to_bra*low_derv_sint(addr_s_hket,addr_cur_derv) &
              + up_derv_bint(3,addr_b_hket,addr_up_derv_z)                        &
              - real(order_hket,REALK)*low_derv_bint(3,addr_low_zket,addr_cur_derv))
            if (max_cur_hbra>1) then
              addr_up_hbra = 9   !base address of the f-shell HGTOs on bra center
              addr_cur_hbra = 3  !base address of the d-shell HGTOs on bra center
              addr_low_hbra = 0  !base address of the p-shell HGTOs on bra center
              ! loops over other current order of HGTOs on bra center, starting from d-shell
              do order_hbra = 2, max_cur_hbra
                addr_up_hbra = addr_up_hbra+1
                addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "orders:", order_hbra, order_xket, order_cur_derv
                write(STDOUT,100) "x-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_low_hbra+1, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,100) "------------------------------------"
#endif
                ! x-direction
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)        &
                  = cc_wrt_bra(1)                                              &
                  * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)   &
                  + neg_hrp_expt*(real(order_hbra,REALK)*ket_to_bra            &
                  * low_derv_bint(addr_low_hbra+1,addr_cur_hket,addr_cur_derv) &
                  + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv))
                ! y-direction
                addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                write(STDOUT,100) "y-direction"
                write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+1
                write(STDOUT,100) "------------------------------------"
#endif
                low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)      &
                  = cc_wrt_bra(2)                                            &
                  * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv) &
                  + neg_hrp_expt                                             &
                  * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+1))
                do ibra = 1, order_hbra
                  addr_up_hbra = addr_up_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_low_hbra+ibra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra+ibra, addr_cur_hket, addr_up_derv+1
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                    = cc_wrt_bra(2)                                                 &
                    * low_derv_bint(addr_cur_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                    + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                     &
                    * low_derv_bint(addr_low_hbra+ibra,addr_cur_hket,addr_cur_derv) &
                    + up_derv_bint(addr_cur_hbra+ibra,addr_b_hket,addr_up_derv+1))
                end do
                ! z-direction
#if defined(DEBUG)
                write(STDOUT,100) "z-direction"
#endif
                addr_cur_hbra = addr_cur_hbra-1
                do jbra = 0, order_hbra
                  addr_up_hbra = addr_up_hbra+1
                  addr_cur_hbra = addr_cur_hbra+1
#if defined(DEBUG)
                  write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                  write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                  write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                  write(STDOUT,100) "------------------------------------"
#endif
                  low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)           &
                    = cc_wrt_bra(3)                                                 &
                    * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)      &
                    + neg_hrp_expt                                                  &
                    * (up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                    - real(order_hket,REALK)                                        &
                    * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                end do
                do ibra = 1, order_hbra
                  do jbra = 0, order_hbra-ibra
                    addr_up_hbra = addr_up_hbra+1
                    addr_cur_hbra = addr_cur_hbra+1
                    addr_low_hbra = addr_low_hbra+1
#if defined(DEBUG)
                    write(STDOUT,110) addr_up_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_low_hbra, addr_cur_hket, addr_cur_derv
                    write(STDOUT,110) addr_cur_hbra, addr_cur_hket, addr_up_derv+iderv+2
                    write(STDOUT,110) addr_cur_hbra, addr_low_zket, addr_cur_derv
                    write(STDOUT,100) "------------------------------------"
#endif
                    low_derv_bint(addr_up_hbra,addr_cur_hket,addr_cur_derv)          &
                      = cc_wrt_bra(3)                                                &
                      * low_derv_bint(addr_cur_hbra,addr_cur_hket,addr_cur_derv)     &
                      + neg_hrp_expt*(real(ibra,REALK)*ket_to_bra                    &
                      * low_derv_bint(addr_low_hbra,addr_cur_hket,addr_cur_derv)     &
                      + up_derv_bint(addr_cur_hbra,addr_b_hket,addr_up_derv+iderv+2) &
                      - real(order_hket,REALK)                                       &
                      * low_derv_bint(addr_cur_hbra,addr_low_zket,addr_cur_derv))
                  end do
                end do
              end do
            end if
          end if
          ! updates the base addresses of lower order HGTOs on ket center
          base_low_hket = addr_low_zket
          base_s_hket = addr_s_zket
        end do
      end do
      addr_up_derv = addr_up_derv+1
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "sub_london_mom_cgto", STDOUT)
#endif
    return
#if defined(DEBUG)
100 format("sub_london_mom_cgto>> ",A,3I6)
110 format("sub_london_mom_cgto>> ","HBRA",I8,4X,"HKET",I8,4X,"DERV",I8)
#endif
  end subroutine sub_london_mom_cgto
