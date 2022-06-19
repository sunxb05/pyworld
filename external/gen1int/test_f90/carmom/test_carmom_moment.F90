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
!!  Tests subroutines in carmom_moment.F90.
!!
!!  2012-02-21, Bin Gao
!!  * first version

  !> \brief tests subroutines carmom_moment.F90
  !> \author Bin Gao
  !> \date 2012-02-21
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_carmom_moment(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer, parameter :: ORDER_ELEC = 2  !order of electronic derivatives
    integer orders_hgto_bra(2)            !orders of HGTOs on bra center
    integer orders_hgto_ket(2)            !orders of HGTOs on ket center
    integer order_mom                     !order of Cartesian multipole moments
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_real_array(log_text="Coordinates of bra center:", &
                             log_real=COORD_BRA, fmt_real="F16.8", io_log=io_log)
    call html_log_real_number("Orbital exponent of bra center:", &
                              EXPONENT_BRA, "F16.8", io_log)
    call html_log_real_array(log_text="Coordinates of ket center:", &
                             log_real=COORD_KET, fmt_real="F16.8", io_log=io_log)
    call html_log_real_number("Orbital exponent of ket center:", &
                              EXPONENT_KET, "F16.8", io_log)
    call html_log_int_number("Order of electronic derivatives:", &
                             ORDER_ELEC, "I3", io_log)
    call html_log_real_number("Scaling constant on Cartesian multipole moments:", &
                              SCAL_CONST, "F16.8", io_log)
    call html_log_real_array(log_text="Coordinates of dipole origin:", &
                             log_real=DIPORG, fmt_real="F16.8", io_log=io_log)
    ! the following tests are for individual cases in \fn(carmom_moment)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,0/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,0/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,0/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,0/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,1/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,1/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,1/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,1/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/0,2/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/0,2/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/0,2/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/0,2/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,3/)
    orders_hgto_ket = (/0,3/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,3/)
    orders_hgto_ket = (/0,3/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,3/)
    orders_hgto_ket = (/0,3/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,3/)
    orders_hgto_ket = (/0,3/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,6/)
    orders_hgto_ket = (/0,6/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,6/)
    orders_hgto_ket = (/0,6/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,6/)
    orders_hgto_ket = (/0,6/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,6/)
    orders_hgto_ket = (/0,6/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/1,1/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/1,1/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/1,1/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/1,1/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/1,2/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/1,2/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/1,2/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/1,2/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,3/)
    orders_hgto_ket = (/1,3/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,3/)
    orders_hgto_ket = (/1,3/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,3/)
    orders_hgto_ket = (/1,3/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,3/)
    orders_hgto_ket = (/1,3/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,6/)
    orders_hgto_ket = (/1,6/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,6/)
    orders_hgto_ket = (/1,6/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,6/)
    orders_hgto_ket = (/1,6/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,6/)
    orders_hgto_ket = (/1,6/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    orders_hgto_ket = (/2,2/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    orders_hgto_ket = (/2,2/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    orders_hgto_ket = (/2,2/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    orders_hgto_ket = (/2,2/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,3/)
    orders_hgto_ket = (/2,3/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,3/)
    orders_hgto_ket = (/2,3/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,3/)
    orders_hgto_ket = (/2,3/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,3/)
    orders_hgto_ket = (/2,3/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/2,6/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/2,6/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/2,6/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/2,6/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/3,7/)
    order_mom = 0
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/3,7/)
    order_mom = 1
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/3,7/)
    order_mom = 2
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/3,7/)
    order_mom = 6
    call sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                ORDER_ELEC, io_log, test_failed)
    return
  end subroutine test_carmom_moment

  !> \brief tests subroutines in carmom_moment.F90 with specific orders of HGTOs
  !>        and order of Cartesian multipole moments
  !> \author Bin Gao
  !> \date 2012-02-21
  !> \param orders_hgto_bra contains the orders of HGTOs on bra center
  !> \param orders_hgto_ket contains the orders of HGTOs on ket center
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_elec is the order of electronic derivatives
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_carmom_moment(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                    order_elec, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions of Cartesian multipole moment integrals
    use recur_carmom
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_elec
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer orders_hbra(2)                         !orders of HGTOs on bra center for \fn(carmom_hbra)
    integer dim_hbra                               !dimension of integrals from \fn(carmom_hbra)
    real(REALK), allocatable :: hbra_pints(:)      !integrals from \fn(carmom_hbra)
    integer orders_hbra_zero(2)                    !orders of HGTOs on bra center for \fn(carmom_hrr_ket)
    integer orders_hket_zero(2)                    !orders of HGTOs on ket center for \fn(carmom_hrr_ket)
    integer dim_hbra_zero                          !dimensions of integrals from \fn(carmom_hrr_ket)
    integer dim_hket_zero
    real(REALK), allocatable :: hket_pints(:,:)    !integrals from \fn(carmom_hrr_ket)
    integer dim_hgto_bra                           !dimensions of integrals from \fn(carmom_moment)
    integer dim_hgto_ket
    integer num_mom                                !number of xyz components of Cartesian multipole moments
    real(REALK), allocatable :: hmom_pints(:,:,:)  !integrals from \fn(carmom_moment)
    integer order_hbra                             !order of HGTO on bra center
    integer x_hbra, y_hbra, z_hbra                 !orders of xyz components of HGTO on bra center
    integer order_hket                             !order of HGTO on ket center
    integer x_hket, y_hket, z_hket                 !orders of xyz components of HGTO on ket center
    integer x_mom, y_mom, z_mom                    !orders of xyz components of Cartesian multipole moments
    real(REALK) recur_pint                         !results from recursive function \fn(recur_carmom_moment)
    integer addr_hbra, addr_hket, addr_mom         !addresses of \var(hmom_pints)
    integer ierr                                   !error information
    real(REALK) begin_time                         !begin of CPU time
    real(REALK) end_time                           !end of CPU time
    logical different                              !if result from \fn(carmom_moment) is different from reference
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on bra center:",            &
                            log_int=orders_hgto_bra, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(log_text="Orders of HGTOs on ket center:",            &
                            log_int=orders_hgto_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_number("Order of Cartesian multipole moments:", &
                             order_mom, "I3", io_log)
    ! gets the begin time
    call xtimer_set(begin_time)
    ! sets the orders of HGTOs for \fn(carmom_hrr_ket)
    orders_hbra_zero(1) = max(0,orders_hgto_bra(1)-order_mom)
    orders_hbra_zero(2) = orders_hgto_bra(2)
    orders_hket_zero(1) = max(0,orders_hgto_ket(1)-order_mom)
    orders_hket_zero(2) = orders_hgto_ket(2)
    ! sets the orders of HGTOs on bra center for \fn(carmom_hbra)
    orders_hbra = orders_hbra_zero+orders_hket_zero
    ! allocates memory for the integrals from \fn(carmom_hbra)
    if (orders_hbra(1)>0) then
      dim_hbra = ((orders_hbra(2)+1)*(orders_hbra(2)+2)*(orders_hbra(2)+3) &
               -  orders_hbra(1)*(orders_hbra(1)+1)*(orders_hbra(1)+2))/6
    else
      dim_hbra = (orders_hbra(2)+1)*(orders_hbra(2)+2)*(orders_hbra(2)+3)/6
    end if
    allocate(hbra_pints(dim_hbra), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hbra_pints:", &
                               log_int=dim_hbra, fmt_int="I12",           &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! calculates the integrals using \fn(carmom_hbra)
    call carmom_hbra(orders_hbra, COORD_BRA, EXPONENT_BRA, &
                     COORD_KET, EXPONENT_KET, order_elec,  &
                     SCAL_CONST, dim_hbra, hbra_pints)
    ! allocates memory for the integrals from \fn(carmom_hrr_ket)
    if (orders_hbra_zero(1)>0) then
      dim_hbra_zero = ((orders_hbra_zero(2)+1)*(orders_hbra_zero(2)+2) &
                       *(orders_hbra_zero(2)+3)                        &
                    -  orders_hbra_zero(1)*(orders_hbra_zero(1)+1)     &
                       *(orders_hbra_zero(1)+2))/6
    else
      dim_hbra_zero = (orders_hbra_zero(2)+1)*(orders_hbra_zero(2)+2) &
                    * (orders_hbra_zero(2)+3)/6
    end if
    if (orders_hket_zero(1)>0) then
      dim_hket_zero = ((orders_hket_zero(2)+1)*(orders_hket_zero(2)+2) &
                       *(orders_hket_zero(2)+3)                        &
                    -  orders_hket_zero(1)*(orders_hket_zero(1)+1)     &
                       *(orders_hket_zero(1)+2))/6
    else
      dim_hket_zero = (orders_hket_zero(2)+1)*(orders_hket_zero(2)+2) &
                    * (orders_hket_zero(2)+3)/6
    end if
    allocate(hket_pints(dim_hbra_zero,dim_hket_zero), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hket_pints:",          &
                               log_int=dim_hbra_zero*dim_hket_zero, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      deallocate(hbra_pints)
      test_failed = .true.
      return
    end if
    ! calculates the integrals using \fn(carmom_hrr_ket)
    call carmom_hrr_ket(orders_hbra_zero, EXPONENT_BRA, &
                        orders_hket_zero, EXPONENT_KET, &
                        dim_hbra, hbra_pints,           &
                        dim_hbra_zero, dim_hket_zero, hket_pints)
    deallocate(hbra_pints)
    ! allocates memory for the integrals from \fn(carmom_moment)
    if (orders_hgto_bra(1)>0) then
      dim_hgto_bra = ((orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                      *(orders_hgto_bra(2)+3)                       &
                   -  orders_hgto_bra(1)*(orders_hgto_bra(1)+1)     & 
                      *(orders_hgto_bra(1)+2))/6
    else
      dim_hgto_bra = (orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                   * (orders_hgto_bra(2)+3)/6
    end if
    if (orders_hgto_ket(1)>0) then
      dim_hgto_ket = ((orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                      *(orders_hgto_ket(2)+3)                       &
                   -  orders_hgto_ket(1)*(orders_hgto_ket(1)+1)     & 
                      *(orders_hgto_ket(1)+2))/6
    else
      dim_hgto_ket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                   * (orders_hgto_ket(2)+3)/6
    end if
    num_mom = (order_mom+1)*(order_mom+2)/2
    allocate(hmom_pints(dim_hgto_bra,dim_hgto_ket,num_mom), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hmom_pints:", &
                               log_int=dim_hgto_bra*dim_hgto_ket*num_mom, &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(hket_pints)
      test_failed = .true.
      return
    end if
    ! calculates the integrals using \fn(carmom_moment)
    call carmom_moment(orders_hgto_bra, COORD_BRA, EXPONENT_BRA,        &
                       orders_hgto_ket, COORD_KET, EXPONENT_KET,        &
                       DIPORG, order_mom, dim_hbra_zero, dim_hket_zero, &
                       hket_pints, dim_hgto_bra, dim_hgto_ket, hmom_pints)
    deallocate(hket_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(                                             &
            log_text="Time (s) used for carmom_hbra, carmom_hrr_ket "//    &
                     "and carmom_moment:",                                 &
            log_real=end_time-begin_time, fmt_real="F10.4", io_log=io_log, &
            font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the data used for the module \fn(recur_carmom)
    call recur_carmom_create(COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET, &
                             DIPORG, SCAL_CONST, order_elec)
    ! initializes the address of integrals from \fn(carmom_moment)
    addr_mom = 0
    do z_mom = 0, order_mom
      do y_mom = 0, order_mom-z_mom
        x_mom = order_mom-(z_mom+y_mom)
        addr_mom = addr_mom+1
        ! initializes the address of integrals from \fn(carmom_moment)
        addr_hket = 0
        ! loops over the orders of HGTOs on ket center
        do order_hket = orders_hgto_ket(1), orders_hgto_ket(2)
          do z_hket = 0, order_hket
            do y_hket = 0, order_hket-z_hket
              x_hket = order_hket-(z_hket+y_hket)
              addr_hket = addr_hket+1
              ! initializes the address of integrals from \fn(carmom_moment)
              addr_hbra = 0
              ! loops over the orders of HGTOs on bra center
              do order_hbra = orders_hgto_bra(1), orders_hgto_bra(2)
                do z_hbra = 0, order_hbra
                  do y_hbra = 0, order_hbra-z_hbra
                    x_hbra = order_hbra-(z_hbra+y_hbra)
                    ! integral from \fn(recur_carmom_moment)
                    recur_pint = recur_carmom_moment((/x_hbra,y_hbra,z_hbra/), &
                                                     (/x_hket,y_hket,z_hket/), &
                                                     (/x_mom,y_mom,z_mom/))
                    addr_hbra = addr_hbra+1
                    call check_difference(recur_pint,                               &
                                          hmom_pints(addr_hbra,addr_hket,addr_mom), &
                                          different)
                    if (different) then
                      call html_log_int_array(log_text="Orders of HGTO on bra:", &
                                              log_int=(/x_hbra,y_hbra,z_hbra/),  &
                                              fmt_int="I3", io_log=io_log)
                      call html_log_int_array(log_text="Orders of HGTO on ket:", &
                                              log_int=(/x_hket,y_hket,z_hket/),  &
                                              fmt_int="I3", io_log=io_log)
                      call html_log_int_array(                                  &
                             log_text="Orders of Cartesian multipole moments:", &
                             log_int=(/x_mom,y_mom,z_mom/), fmt_int="I3",       &
                             io_log=io_log)
                      call html_log_real_number(                             &
                             log_text="Reference from recur_carmom_moment:", &
                             log_real=recur_pint, fmt_real="Es20.12",        &
                             io_log=io_log, font_color="blue")
                      call html_log_real_number(                                &
                             log_text="Result from carmom_moment:",             &
                             log_real=hmom_pints(addr_hbra,addr_hket,addr_mom), &
                             fmt_real="Es20.12", io_log=io_log, font_color="red")
                      test_failed = .true.
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_carmom_moment:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",    &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(hmom_pints)
    return
  end subroutine sub_test_carmom_moment
