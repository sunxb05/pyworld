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
!!  Tests subroutine prim_hgto_nucpot.
!!
!!  2011-07-27, Bin Gao
!!  * first version

  !> \brief tests subroutine prim_hgto_nucpot
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_prim_hgto_nucpot(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer orders_hgto_bra(2)  !orders of HGTOs on bra center
    integer orders_hgto_ket(2)  !orders of HGTOs on ket center
    integer order_mom           !order of Cartesian multipole moments
    integer order_geo_pot       !order of geometric derivatives on nuclear potential origin
    integer order_elec          !order of electronic derivatives
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
    call html_log_real_array(log_text="Coordinates of nuclear potential origin:", &
                             log_real=NUCORG, fmt_real="F16.8", io_log=io_log)
    call html_log_real_array(log_text="Coordinates of dipole origin:", &
                             log_real=DIPORG, fmt_real="F16.8", io_log=io_log)
    call html_log_real_number("Scaling constant on operator:", &
                              SCAL_CONST, "F16.8", io_log)
    ! the following tests are for individual cases in \fn(prim_hgto_nucpot)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,0/)
    order_mom = 3
    order_geo_pot = 3
    order_elec = 0
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,0/)
    order_mom = 3
    order_geo_pot = 3
    order_elec = 1
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,1/)
    order_mom = 3
    order_geo_pot = 3
    order_elec = 2
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,1/)
    order_mom = 3
    order_geo_pot = 2
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/1,1/)
    order_mom = 2
    order_geo_pot = 3
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/1,1/)
    order_mom = 2
    order_geo_pot = 2
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/1,1/)
    order_mom = 2
    order_geo_pot = 2
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,2/)
    order_mom = 3
    order_geo_pot = 2
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/0,2/)
    order_mom = 1
    order_geo_pot = 2
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/1,2/)
    order_mom = 2
    order_geo_pot = 2
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/1,2/)
    order_mom = 2
    order_geo_pot = 1
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/1,2/)
    order_mom = 1
    order_geo_pot = 1
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/2,2/)
    order_mom = 2
    order_geo_pot = 1
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/2,2/)
    order_mom = 1
    order_geo_pot = 1
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/2,2/)
    order_mom = 1
    order_geo_pot = 1
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    orders_hgto_ket = (/2,2/)
    order_mom = 1
    order_geo_pot = 1
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/3,4/)
    order_mom = 1
    order_geo_pot = 1
    order_elec = 3
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    orders_hgto_bra = (/3,4/)
    orders_hgto_ket = (/3,3/)
    order_mom = 0
    order_geo_pot = 0
    order_elec = 1
    call sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                   order_elec, order_geo_pot, io_log, test_failed)
    return
  end subroutine test_prim_hgto_nucpot

  !> \brief tests subroutine prim_hgto_nucpot with specific orders of
  !>        Hermite Gaussians and Cartesian multipole moments
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param orders_hgto_bra contains the orders of HGTOs on bra center
  !> \param orders_hgto_ket contains the orders of HGTOs on ket center
  !> \param order_mom is the order of Cartesian multipole moments
  !> \param order_geo_pot is the order of geometric derivatives on nuclear potential origin
  !> \param order_elec is the order of electronic derivatives
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_prim_hgto_nucpot(orders_hgto_bra, orders_hgto_ket, order_mom, &
                                       order_elec, order_geo_pot, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions of nuclear attraction potential integrals
    use recur_nucpot
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: order_mom
    integer, intent(in) :: order_geo_pot
    integer, intent(in) :: order_elec
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer max_order_geo           !maximum order of geometric derivatives on nuclear potential origin
    integer dim_hgto_bra            !dimension of HGTOs on bra center
    integer dim_hgto_ket            !dimension of HGTOs on ket center
    integer num_elec                !number of electronic derivatives
    integer num_mom                 !number of Cartesian multipole moments
    integer num_geo_pot             !number of geometric derivatives on nuclear potential origin
    real(REALK), allocatable :: hgto_pints(:,:,:,:,:)
                                    !integrals from \fn(prim_hgto_nucpot)
    integer order_hbra              !order of HGTO on bra center
    integer x_hbra, y_hbra, z_hbra  !orders of xyz components of HGTO on bra center
    integer order_hket              !order of HGTO on ket center
    integer x_hket, y_hket, z_hket  !orders of xyz components of HGTO on ket center
    integer x_elec, y_elec, z_elec  !orders of xyz components of electronic derivative
    integer x_mom, y_mom, z_mom     !orders of xyz components of Cartesian multipole moment
    integer x_geo, y_geo, z_geo     !orders of xyz components of geometric derivatives
    real(REALK) recur_pint          !results from recursive functions
    integer addr_hbra, addr_hket, addr_elec, addr_mom, addr_geo
                                    !addresses of \var(hgto_pints)
    integer ierr                    !error information
    real(REALK) begin_time          !begin of CPU time
    real(REALK) end_time            !end of CPU time
    logical different               !if result from \fn(prim_hgto_nucpot) is different from reference
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on bra center:",            &
                            log_int=orders_hgto_bra, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(log_text="Orders of HGTOs on ket center:",            &
                            log_int=orders_hgto_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_number("Order of electronic derivatives:", &
                             order_elec, "I3", io_log)
    call html_log_int_number("Order of Cartesian multipole moment:", order_mom, &
                             "I3", io_log)
    call html_log_int_number(                                            &
           "Order of geometric derivative on nuclear potential origin:", &
           order_geo_pot, "I3", io_log)
    ! gets the begin time
    call xtimer_set(begin_time)
    ! allocates memory for the integrals from \fn(prim_hgto_nucpot)
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
    num_elec = (order_elec+1)*(order_elec+2)/2
    num_mom = (order_mom+1)*(order_mom+2)/2
    num_geo_pot = (order_geo_pot+1)*(order_geo_pot+2)/2
    allocate(hgto_pints(dim_hgto_bra,dim_hgto_ket,num_elec, &
                        num_mom,num_geo_pot), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hgto_pints:", &
                               log_int=dim_hgto_bra*dim_hgto_ket*num_elec &
                                       *num_mom*num_geo_pot,              &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    call prim_hgto_nucpot(orders_hgto_bra, COORD_BRA, EXPONENT_BRA, &
                          orders_hgto_ket, COORD_KET, EXPONENT_KET, &
                          order_elec, NUCORG, DIPORG, SCAL_CONST,   &
                          order_mom, order_geo_pot, dim_hgto_bra,   &
                          dim_hgto_ket, num_elec, num_mom,          &
                          num_geo_pot, hgto_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for prim_hgto_nucpot:", &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! maximum order of geometric derivatives on nuclear potential origin
    max_order_geo = orders_hgto_bra(2)+orders_hgto_ket(2)+order_elec &
                  + order_mom+order_geo_pot
    ! sets the data used for the module \fn(recur_nucpot)
    call recur_nucpot_create(coord_bra=COORD_BRA, exponent_bra=EXPONENT_BRA, &
                             coord_ket=COORD_KET, exponent_ket=EXPONENT_KET, &
                             nucpot_origin=NUCORG, dipole_origin=DIPORG,     &
                             scal_const=SCAL_CONST, order_elec=order_elec,   &
                             order_aux=max_order_geo)
    ! initializes the address of integrals from \fn(prim_hgto_nucpot)
    addr_geo = 0
    ! loops over the xyz components of geometric derivatives
    do z_geo = 0, order_geo_pot
      do y_geo = 0, order_geo_pot-z_geo
        x_geo = order_geo_pot-(z_geo+y_geo)
        addr_geo = addr_geo+1
        ! initializes the address of integrals from \fn(prim_hgto_nucpot)
        addr_mom = 0
        ! loops over the xyz components of Cartesian multipole moments
        do z_mom = 0, order_mom
          do y_mom = 0, order_mom-z_mom
            x_mom = order_mom-(z_mom+y_mom)
            addr_mom = addr_mom+1
            ! initializes the address of integrals from \fn(prim_hgto_nucpot)
            addr_elec = 0
            ! loops over xyz components of electronic derivative
            do z_elec = 0, order_elec
              do y_elec = 0, order_elec-z_elec
                x_elec = order_elec-(z_elec+y_elec)
                addr_elec = addr_elec+1
                ! initializes the address of integrals from \fn(prim_hgto_nucpot)
                addr_hket = 0
                ! loops over the orders of HGTOs on ket center
                do order_hket = orders_hgto_ket(1), orders_hgto_ket(2)
                  do z_hket = 0, order_hket
                    do y_hket = 0, order_hket-z_hket
                      x_hket = order_hket-(z_hket+y_hket)
                      addr_hket = addr_hket+1
                      ! initializes the address of integrals from \fn(prim_hgto_nucpot)
                      addr_hbra = 0
                      ! loops over the orders of HGTOs on bra center
                      do order_hbra = orders_hgto_bra(1), orders_hgto_bra(2)
                        do z_hbra = 0, order_hbra
                          do y_hbra = 0, order_hbra-z_hbra
                            x_hbra = order_hbra-(z_hbra+y_hbra)
                            ! integral from recursive functions
                            recur_pint = recur_nucpot_moment((/x_hbra,y_hbra,z_hbra/), &
                                                             (/x_hket+x_elec,          &
                                                               y_hket+y_elec,          &
                                                               z_hket+z_elec/),        &
                                                             (/x_mom,y_mom,z_mom/),    &
                                                             (/x_geo,y_geo,z_geo/))
                            addr_hbra = addr_hbra+1
                            call check_difference(recur_pint,                &
                                   hgto_pints(addr_hbra,addr_hket,addr_elec, &
                                              addr_mom,addr_geo),            &
                                   different)
                            if (different) then
                              call html_log_int_array(log_text="Orders of HGTO on bra:", &
                                                      log_int=(/x_hbra,y_hbra,z_hbra/),  &
                                                      fmt_int="I3", io_log=io_log)
                              call html_log_int_array(log_text="Orders of HGTO on ket:", &
                                                      log_int=(/x_hket,y_hket,z_hket/),  &
                                                      fmt_int="I3", io_log=io_log)
                              call html_log_int_array(                               &
                                     log_text="Orders of electronic derivatives:",   &
                                     log_int=(/x_elec,y_elec,z_elec/), fmt_int="I3", &
                                     io_log=io_log)
                              call html_log_int_array(                                 &
                                     log_text="Orders of Cartesian multipole moment ", &
                                     log_int=(/x_mom,y_mom,z_mom/), fmt_int="I3",      &
                                     io_log=io_log)
                              call html_log_int_array(                               &
                                     log_text="Orders of geometric derivative on "// &
                                              "nuclear potential origin:",           &
                                     log_int=(/x_geo,y_geo,z_geo/), fmt_int="I3",    &
                                     io_log=io_log)
                              call html_log_real_number(                             &
                                     log_text="Reference from recur_nucpot_moment:", &
                                     log_real=recur_pint, fmt_real="Es20.12",        &
                                     io_log=io_log, font_color="blue")
                              call html_log_real_number(                                &
                                     log_text="Result from prim_hgto_nucpot:",          &
                                     log_real=hgto_pints(addr_hbra,addr_hket,addr_elec, &
                                                         addr_mom,addr_geo),            &
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
          end do
        end do
      end do
    end do
    call recur_nucpot_destroy
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recursive functions:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",    &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(hgto_pints)
    return
  end subroutine sub_test_prim_hgto_nucpot
