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
!!  Tests subroutine carmom_hbra.
!!
!!  2012-02-12, Bin Gao
!!  * first version

  !> \brief tests subroutine carmom_hbra
  !> \author Bin Gao
  !> \date 2012-02-12
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_carmom_hbra(io_log, test_failed)
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
    ! the following tests are for individual cases in \fn(carmom_hbra)
    orders_hgto_bra = (/0,0/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,3/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,6/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,3/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,6/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,3/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/3,7/)
    call sub_test_carmom_hbra(orders_hgto_bra, ORDER_ELEC, io_log, test_failed)
    return
  end subroutine test_carmom_hbra

  !> \brief tests subroutine carmom_hbra with specific orders of HGTOs
  !> \author Bin Gao
  !> \date 2012-02-12
  !> \param orders_hgto_bra contains the orders of HGTOs on bra center
  !> \param order_elec is the order of electronic derivatives
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_carmom_hbra(orders_hgto_bra, order_elec, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions of Cartesian multipole moment integrals
    use recur_carmom
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer dim_hgto_bra                       !dimension of integrals from \fn(carmom_hbra)
    real(REALK), allocatable :: hbra_pints(:)  !integrals from \fn(carmom_hbra)
    integer order_hbra                         !order of HGTO on bra center
    integer x_hbra, y_hbra, z_hbra             !orders of xyz components of HGTO on bra center
    real(REALK) recur_pint                     !results from recursive function \fn(recur_carmom_hbra)
    integer addr_hbra                          !address of \var(hbra_pints)
    integer ierr                               !error information
    real(REALK) begin_time                     !begin of CPU time
    real(REALK) end_time                       !end of CPU time
    logical different                          !if result from \fn(carmom_hbra) is different from reference
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on bra center:",            &
                            log_int=orders_hgto_bra, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    ! gets the begin time
    call xtimer_set(begin_time)
    ! allocates memory for the integrals from \fn(carmom_hbra)
    if (orders_hgto_bra(1)>0) then
      dim_hgto_bra = ((orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                      *(orders_hgto_bra(2)+3)                       &
                   -  orders_hgto_bra(1)*(orders_hgto_bra(1)+1)     &
                      *(orders_hgto_bra(1)+2))/6
    else
      dim_hgto_bra = (orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                   * (orders_hgto_bra(2)+3)/6
    end if
    allocate(hbra_pints(dim_hgto_bra), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hbra_pints:", &
                               log_int=dim_hgto_bra, fmt_int="I12",       &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! calculates the integrals using \fn(carmom_hbra)
    call carmom_hbra(orders_hgto_bra, COORD_BRA, EXPONENT_BRA, &
                     COORD_KET, EXPONENT_KET, order_elec,      &
                     SCAL_CONST, dim_hgto_bra, hbra_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for carmom_hbra:",      &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the data used for the module \fn(recur_carmom)
    call recur_carmom_create(COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET, &
                             DIPORG, SCAL_CONST, order_elec)
    ! initializes the address of integrals from \fn(carmom_hbra)
    addr_hbra = 0
    ! loops over the orders of HGTOs on bra center
    do order_hbra = orders_hgto_bra(1), orders_hgto_bra(2)
      do z_hbra = 0, order_hbra
        do y_hbra = 0, order_hbra-z_hbra
          x_hbra = order_hbra-(z_hbra+y_hbra)
          ! integral from \fn(recur_carmom_hbra)
          recur_pint = recur_carmom_hbra((/x_hbra,y_hbra,z_hbra/))
          addr_hbra = addr_hbra+1
          call check_difference(recur_pint, hbra_pints(addr_hbra), different)
          if (different) then
            call html_log_int_array(log_text="Orders of HGTO on bra:", &
                                    log_int=(/x_hbra,y_hbra,z_hbra/),  &
                                    fmt_int="I3", io_log=io_log)
            call html_log_real_number(log_text="Reference from recur_carmom_hbra:", &
                                      log_real=recur_pint, fmt_real="Es20.12",      &
                                      io_log=io_log, font_color="blue")
            call html_log_real_number(log_text="Result from carmom_hbra:", &
                                      log_real=hbra_pints(addr_hbra),      &
                                      fmt_real="Es20.12", io_log=io_log,   &
                                      font_color="red")
            test_failed = .true.
          end if
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_carmom_hbra:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",  &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(hbra_pints)
    return
  end subroutine sub_test_carmom_hbra
