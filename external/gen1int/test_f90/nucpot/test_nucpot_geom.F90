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
!!  Tests subroutines in nucpot_geom.F90.
!!
!!  2011-09-17, Bin Gao
!!  * first version

  !> \brief tests subroutines in nucpot_geom.F90
  !> \author Bin Gao
  !> \date 2011-09-17
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_nucpot_geom(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer, parameter :: ORDER_ELEC = 2  !order of electronic derivatives
    integer orders_geo_pot(2)             !orders of geometric derivatives on nuclear attraction potential origin
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
    call html_log_real_array(log_text="Coordinates of nuclear potential origin:", &
                             log_real=NUCORG, fmt_real="F16.8", io_log=io_log)
    call html_log_real_number("Scaling constant on operator:", &
                              SCAL_CONST, "F16.8", io_log)
    ! the following tests are for individual cases in \fn(nucpot_geom)
    orders_geo_pot = (/0,0/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/0,1/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/1,1/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/1,7/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/1,8/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/0,2/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/1,2/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/2,2/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/2,7/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/2,8/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/3,7/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/3,8/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/6,15/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    orders_geo_pot = (/6,16/)
    call sub_test_nucpot_geom(orders_geo_pot, ORDER_ELEC, io_log, test_failed)
    return
  end subroutine test_nucpot_geom

  !> \brief tests subroutines in nucpot_geom.F90 with specific orders of geometric derivatives
  !>        on nuclear attraction potential origin
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param orders_geo_pot contains the orders of geometric derivatives
  !> \param order_elec is the order of electronic derivatives
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_nucpot_geom(orders_geo_pot, order_elec, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions of nuclear attraction potential integrals
    use recur_nucpot
    implicit none
    integer, intent(in) :: orders_geo_pot(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer dim_geo_pot          !dimension of primitive integrals with zeroth order HGTOs
    real(REALK), allocatable :: geo_pot_pints(:)
                                 !primitive integrals with zeroth order HGTOs
    integer order_geo            !order of geometric derivatives
    integer x_geo, y_geo, z_geo  !orders of xyz components of geometric derivatives
    real(REALK) recur_pint       !results from recursive function \fn(recur_nucpot_geom)
    integer addr_geo             !address of \var(geo_pot_pints)
    integer ierr                 !error information
    real(REALK) begin_time       !begin of CPU time
    real(REALK) end_time         !end of CPU time
    logical different            !if result from \fn(nucpot_geom) is different from reference
    ! dumps information of tests
    call html_log_int_array(                                                       &
           log_text="Orders of geometric derivative on nuclear potential origin:", &
           log_int=orders_geo_pot, fmt_int="I3", io_log=io_log, separation="&rarr;")
    ! gets the begin time
    call xtimer_set(begin_time)
    ! allocates memory for the integrals from \fn(nucpot_geom)
    if (orders_geo_pot(1)>0) then
      dim_geo_pot = ((orders_geo_pot(2)+1)*(orders_geo_pot(2)+2)*(orders_geo_pot(2)+3) &
                  -  orders_geo_pot(1)*(orders_geo_pot(1)+1)*(orders_geo_pot(1)+2))/6
    else
      dim_geo_pot = (orders_geo_pot(2)+1)*(orders_geo_pot(2)+2)*(orders_geo_pot(2)+3)/6
    end if
    allocate(geo_pot_pints(dim_geo_pot), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate geo_pot_pints:", &
                               log_int=dim_geo_pot, fmt_int="I12",           &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! calculates the integrals using \fn(nucpot_geom)
    call nucpot_geom(COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET, &
                      NUCORG, SCAL_CONST, orders_geo_pot, order_elec,  &
                      dim_geo_pot, geo_pot_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for nucpot_geom:",      &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the data used for the module \fn(recur_nucpot)
    call recur_nucpot_create(coord_bra=COORD_BRA, exponent_bra=EXPONENT_BRA, &
                             coord_ket=COORD_KET, exponent_ket=EXPONENT_KET, &
                             nucpot_origin=NUCORG, dipole_origin=DIPORG,     &
                             scal_const=SCAL_CONST, order_elec=order_elec,   &
                             order_aux=orders_geo_pot(2))
    ! initializes the address of integrals from \fn(nucpot_geom)
    addr_geo = 0
    ! loops over the orders of geometric derivatives
    do order_geo = orders_geo_pot(1), orders_geo_pot(2)
      ! loops over xyz components
      do z_geo = 0, order_geo
        do y_geo = 0, order_geo-z_geo
          x_geo = order_geo-(z_geo+y_geo)
          ! integral from \fn(recur_nucpot_geom)
          recur_pint = recur_nucpot_geom((/x_geo,y_geo,z_geo/), 0)
          addr_geo = addr_geo+1
          call check_difference(recur_pint, geo_pot_pints(addr_geo), &
                                different)
          if (different) then
            call html_log_int_array(                               &
                   log_text="Orders of geometric derivative on "// &
                            "nuclear potential origin:",           &
                   log_int=(/x_geo,y_geo,z_geo/), fmt_int="I3",    &
                   io_log=io_log)
            call html_log_real_number(                           &
                   log_text="Reference from recur_nucpot_geom:", &
                   log_real=recur_pint, fmt_real="Es20.12",      &
                   io_log=io_log, font_color="blue")
            call html_log_real_number(log_text="Result from nucpot_geom:", &
                                      log_real=geo_pot_pints(addr_geo),    &
                                      fmt_real="Es20.12", io_log=io_log,   &
                                      font_color="red")
            test_failed = .true.
          end if
        end do
      end do
    end do
    call recur_nucpot_destroy
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_nucpot_geom:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",  &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(geo_pot_pints)
    return
  end subroutine sub_test_nucpot_geom
