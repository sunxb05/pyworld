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
!!  Tests subroutines in carmom_hrr_ket.F90.
!!
!!  2012-02-15, Bin Gao
!!  * first version

  !> \brief tests subroutines carmom_hrr_ket.F90
  !> \author Bin Gao
  !> \date 2012-02-15
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_carmom_hrr_ket(io_log, test_failed)
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
    ! the following tests are for individual cases in \fn(carmom_hrr_ket)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,0/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,1/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/0,2/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/0,3/)
    orders_hgto_ket = (/0,3/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/0,6/)
    orders_hgto_ket = (/0,6/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/1,1/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/1,2/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/1,3/)
    orders_hgto_ket = (/1,3/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/1,6/)
    orders_hgto_ket = (/1,6/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    orders_hgto_ket = (/2,2/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/2,3/)
    orders_hgto_ket = (/2,3/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/2,6/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    orders_hgto_ket = (/3,7/)
    call sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, ORDER_ELEC, &
                                 io_log, test_failed)
    return
  end subroutine test_carmom_hrr_ket

  !> \brief tests subroutine carmom_hrr_ket with specific orders of HGTOs
  !> \author Bin Gao
  !> \date 2012-02-15
  !> \param orders_hgto_bra contains the orders of HGTOs on bra center
  !> \param orders_hgto_ket contains the orders of HGTOs on ket center
  !> \param order_elec is the order of electronic derivatives
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_carmom_hrr_ket(orders_hgto_bra, orders_hgto_ket, order_elec, &
                                     io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions of Cartesian multipole moment integrals
    use recur_carmom
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer orders_hbra(2)                       !orders of HGTOs on bra center for \fn(carmom_hbra)
    integer dim_hbra                             !dimension of integrals from \fn(carmom_hbra)
    real(REALK), allocatable :: hbra_pints(:)    !integrals from \fn(carmom_hbra)
    integer dim_hgto_bra                         !dimensions of integrals from \fn(carmom_hrr_ket)
    integer dim_hgto_ket
    real(REALK), allocatable :: hket_pints(:,:)  !integrals from \fn(carmom_hrr_ket)
    integer order_hbra                           !order of HGTO on bra center
    integer x_hbra, y_hbra, z_hbra               !orders of xyz components of HGTO on bra center
    integer order_hket                           !order of HGTO on ket center
    integer x_hket, y_hket, z_hket               !orders of xyz components of HGTO on ket center
    real(REALK) recur_pint                       !results from recursive function \fn(recur_carmom_hrr_ket)
    integer addr_hbra, addr_hket                 !addresses of \var(hket_pints)
    integer ierr                                 !error information
    real(REALK) begin_time                       !begin of CPU time
    real(REALK) end_time                         !end of CPU time
    logical different                            !if result from \fn(carmom_hrr_ket) is different from reference
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on bra center:",            &
                            log_int=orders_hgto_bra, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(log_text="Orders of HGTOs on ket center:",            &
                            log_int=orders_hgto_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    ! gets the begin time
    call xtimer_set(begin_time)
    ! sets the orders of HGTOs on bra center for \fn(carmom_hbra)
    orders_hbra = orders_hgto_bra+orders_hgto_ket
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
    allocate(hket_pints(dim_hgto_bra,dim_hgto_ket), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hket_pints:",        &
                               log_int=dim_hgto_bra*dim_hgto_ket, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      deallocate(hbra_pints)
      test_failed = .true.
      return
    end if
    ! calculates the integrals using \fn(carmom_hrr_ket)
    call carmom_hrr_ket(orders_hgto_bra, EXPONENT_BRA, &
                        orders_hgto_ket, EXPONENT_KET, &
                        dim_hbra, hbra_pints,          &
                        dim_hgto_bra, dim_hgto_ket, hket_pints)
    deallocate(hbra_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(                                             &
            log_text="Time (s) used for carmom_hbra and carmom_hrr_ket:",  &
            log_real=end_time-begin_time, fmt_real="F10.4", io_log=io_log, &
            font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the data used for the module \fn(recur_carmom)
    call recur_carmom_create(COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET, &
                             DIPORG, SCAL_CONST, order_elec)
    ! initializes the address of integrals from \fn(carmom_hrr_ket)
    addr_hket = 0
    ! loops over the orders of HGTOs on ket center
    do order_hket = orders_hgto_ket(1), orders_hgto_ket(2)
      do z_hket = 0, order_hket
        do y_hket = 0, order_hket-z_hket
          x_hket = order_hket-(z_hket+y_hket)
          addr_hket = addr_hket+1
          ! initializes the address of integrals from \fn(carmom_hrr_ket)
          addr_hbra = 0
          ! loops over the orders of HGTOs on bra center
          do order_hbra = orders_hgto_bra(1), orders_hgto_bra(2)
            do z_hbra = 0, order_hbra
              do y_hbra = 0, order_hbra-z_hbra
                x_hbra = order_hbra-(z_hbra+y_hbra)
                ! integral from \fn(recur_carmom_hrr_ket)
                recur_pint = recur_carmom_hrr_ket((/x_hbra,y_hbra,z_hbra/), &
                                                  (/x_hket,y_hket,z_hket/))
                addr_hbra = addr_hbra+1
                call check_difference(recur_pint, &
                                      hket_pints(addr_hbra,addr_hket), different)
                if (different) then
                  call html_log_int_array(log_text="Orders of HGTO on bra:", &
                                          log_int=(/x_hbra,y_hbra,z_hbra/),  &
                                          fmt_int="I3", io_log=io_log)
                  call html_log_int_array(log_text="Orders of HGTO on ket:", &
                                          log_int=(/x_hket,y_hket,z_hket/),  &
                                          fmt_int="I3", io_log=io_log)
                  call html_log_real_number(                              &
                         log_text="Reference from recur_carmom_hrr_ket:", &
                         log_real=recur_pint, fmt_real="Es20.12",         &
                         io_log=io_log, font_color="blue")
                  call html_log_real_number(                       &
                         log_text="Result from carmom_hrr_ket:",   &
                         log_real=hket_pints(addr_hbra,addr_hket), &
                         fmt_real="Es20.12", io_log=io_log, font_color="red")
                  test_failed = .true.
                end if
              end do
            end do
          end do
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_carmom_hrr_ket:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",     &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(hket_pints)
    return
  end subroutine sub_test_carmom_hrr_ket
