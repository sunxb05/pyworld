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
!!  Tests subroutines in delta_hket.F90.
!!
!!  2012-03-16, Bin Gao
!!  * first version

  !> \brief tests subroutines in delta_hket.F90
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_delta_hket(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer, parameter :: ORDER_ELEC = 2  !order of electronic derivatives
    integer orders_hgto_ket(2)            !orders of HGTOs on ket center
    integer orders_geo_pot(2)             !orders of geometric derivatives
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
    call html_log_real_number("Scaling constant on Dirac delta function:", &
                              SCAL_CONST, "F16.8", io_log)
    ! the following tests are for individual cases in \fn(delta_hket)
    orders_hgto_ket = (/0,0/)
    orders_geo_pot = (/0,3/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,1/)
    orders_geo_pot = (/0,0/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,2/)
    orders_geo_pot = (/0,0/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,3/)
    orders_geo_pot = (/0,0/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,5/)
    orders_geo_pot = (/0,0/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,1/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/1,1/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,2/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,3/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,4/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/4,4/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,5/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/4,5/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/5,5/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/5,6/)
    orders_geo_pot = (/4,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,1/)
    orders_geo_pot = (/0,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/1,1/)
    orders_geo_pot = (/1,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,2/)
    orders_geo_pot = (/0,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,3/)
    orders_geo_pot = (/1,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,4/)
    orders_geo_pot = (/0,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/4,4/)
    orders_geo_pot = (/1,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/0,5/)
    orders_geo_pot = (/0,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/4,5/)
    orders_geo_pot = (/1,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/5,5/)
    orders_geo_pot = (/0,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/5,6/)
    orders_geo_pot = (/1,4/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    orders_hgto_ket = (/5,6/)
    orders_geo_pot = (/1,2/)
    call sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, ORDER_ELEC, &
                             io_log, test_failed)
    return
  end subroutine test_delta_hket

  !> \brief tests subroutines in delta_hket.F90 with specific orders of HGTOs on
  !>        ket center and geometric derivatives on Dirac delta function
  !> \author Bin Gao
  !> \date 2012-03-16
  !> \param orders_hgto_ket contains the orders of HGTOs on ket center
  !> \param orders_geo_pot contains the orders of geometric derivatives on Dirac delta function
  !> \param order_elec is the order of electronic derivatives
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_delta_hket(orders_hgto_ket, orders_geo_pot, order_elec, &
                                 io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions of Dirac delta function integrals
    use recur_delta
    implicit none
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: orders_geo_pot(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer orders_hgto_bra(2)                     !orders of HGTOs on bra center
    integer order_hbra(3)                          !orders of xyz components of HGTOs on bra center
    integer max_dim_geo                            !maximum dimension of geometric derivatives on ...
                                                   !Dirac delta function
    integer min_delta_geom                         !minimum order of geometric derivatives on Dirac ...
                                                   !delta function from \fn(delta_geom)
    integer dim_geo_pot                            !dimension of integrals from \fn(delta_geom)
    real(REALK), allocatable :: geo_pot_pints(:)   !integrals from \fn(delta_geom)
    integer min_geo_hbra                           !minimum order of geometric derivatives on Dirac ...
                                                   !delta function after recovering HGTOs on bra center
    integer dim_hgto_bra                           !dimension of HGTOs on bra center
    integer dim_geo_hbra                           !dimension of geometric derivatives on Dirac delta ...
                                                   !function after recovering HGTOs on bra center
    real(REALK), allocatable :: hbra_pints(:,:)    !integrals after recovering HGTOs on bra center
    integer dim_hgto_ket                           !dimension of HGTOs on ket center
    integer dim_geo_hket                           !dimension of geometric derivatives on Dirac delta ...
                                                   !function after recovering HGTOs on ket center
    real(REALK), allocatable :: hket_pints(:,:,:)  !integrals after recovering HGTOs on ket center
    integer order_geo                              !order of geometric derivatives on Dirac delta function
    integer x_geo, y_geo, z_geo                    !orders of xyz components of geometric derivatives
    integer order_hket                             !order of HGTOs on ket center
    integer x_hket, y_hket, z_hket                 !orders of xyz components of HGTOs on ket center
    real(REALK) recur_pint                         !results from recursive function \fn(rec_delta_hket)
    integer addr_hket, addr_geo                    !addresses of \var(hket_pints)
    integer ierr                                   !error information
    real(REALK) begin_time                         !begin of CPU time
    real(REALK) end_time                           !end of CPU time
    logical different                              !if result from \fn(delta_hket) is different from reference
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on ket center:",            &
                            log_int=orders_hgto_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(log_text="Orders of geometric derivatives:",         &
                            log_int=orders_geo_pot, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    ! sets the orders of HGTOs on bra center
    orders_hgto_bra = (/0,0/)
    order_hbra = (/0,0,0/)
    ! gets the begin time
    call xtimer_set(begin_time)
    ! sets the minimum order of geometric derivatives on Dirac delta function
    ! after recovering HGTOs on bra center
    min_geo_hbra = max(0,orders_geo_pot(1)-orders_hgto_ket(2))
    ! sets the minimum order of geometric derivatives on Dirac delta function from \fn(delta_geom)
    min_delta_geom = max(0,min_geo_hbra-orders_hgto_bra(2))
    ! sets the maximum dimension of geometric derivatives on Dirac delta function
    max_dim_geo = (orders_geo_pot(2)+1)*(orders_geo_pot(2)+2)*(orders_geo_pot(2)+3)/6
    ! allocates memory for integrals from \fn(delta_geom)
    if (min_delta_geom==0) then
      dim_geo_pot = max_dim_geo
    else
      dim_geo_pot = max_dim_geo &
                  - min_delta_geom*(min_delta_geom+1)*(min_delta_geom+2)/6
    end if
    allocate(geo_pot_pints(dim_geo_pot), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate geo_pot_pints:", &
                               log_int=dim_geo_pot, fmt_int="I12",           &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! recovers the geometric derivatives on Dirac delta function using \fn(delta_geom)
    call delta_geom(COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET,         &
                    DELORG, SCAL_CONST, (/min_delta_geom,orders_geo_pot(2)/), &
                    order_elec, dim_geo_pot, geo_pot_pints)
    ! allocates memory for integrals from \fn(delta_hket) on bra center
    if (orders_hgto_bra(1)==0) then
      dim_hgto_bra = (orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                   * (orders_hgto_bra(2)+3)/6
    else
      dim_hgto_bra = ((orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                   *  (orders_hgto_bra(2)+3)                        &
                   -  orders_hgto_bra(1)*(orders_hgto_bra(1)+1)     &
                   *  (orders_hgto_bra(1)+2))/6
    end if
    if (min_geo_hbra==0) then
      dim_geo_hbra = max_dim_geo
    else
      dim_geo_hbra = max_dim_geo-min_geo_hbra*(min_geo_hbra+1)*(min_geo_hbra+2)/6
    end if
    allocate(hbra_pints(dim_hgto_bra,dim_geo_hbra), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hbra_pints:", &
                               log_int=dim_hgto_bra*dim_geo_hbra,         &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(geo_pot_pints)
      test_failed = .true.
      return
    end if
    ! recovers the HGTOs on bra center using \fn(delta_hket)
    call delta_hket(orders_hgto_bra, (/min_geo_hbra,orders_geo_pot(2)/), &
                    COORD_BRA, EXPONENT_BRA, DELORG, 1, dim_geo_pot,     &
                    geo_pot_pints, dim_hgto_bra, dim_geo_hbra, hbra_pints)
    deallocate(geo_pot_pints)
    ! allocates memory for integrals from \fn(delta_hket) on ket center
    if (orders_hgto_ket(1)==0) then
      dim_hgto_ket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                   * (orders_hgto_ket(2)+3)/6
    else
      dim_hgto_ket = ((orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) &
                   *  (orders_hgto_ket(2)+3)                        &
                   -  orders_hgto_ket(1)*(orders_hgto_ket(1)+1)     &
                   *  (orders_hgto_ket(1)+2))/6
    end if
    if (orders_geo_pot(1)==0) then
      dim_geo_hket = max_dim_geo
    else
      dim_geo_hket = max_dim_geo &
                   - orders_geo_pot(1)*(orders_geo_pot(1)+1)*(orders_geo_pot(1)+2)/6
    end if
    allocate(hket_pints(dim_hgto_bra,dim_hgto_ket,dim_geo_hket), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hket_pints:",      &
                               log_int=dim_hgto_bra*dim_hgto_ket*dim_geo_hket, &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(hbra_pints)
      test_failed = .true.
      return
    end if
    ! recovers the HGTOs on ket center using \fn(delta_hket)
    call delta_hket(orders_hgto_ket, orders_geo_pot, COORD_KET, EXPONENT_KET, &
                    DELORG, dim_hgto_bra, dim_geo_hbra, hbra_pints,           &
                    dim_hgto_ket, dim_geo_hket, hket_pints)
    deallocate(hbra_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for delta_geom and delta_hket:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",          &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the data used for the module \fn(recur_delta)
    call recur_delta_create(COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET, &
                            DELORG, DIPORG, SCAL_CONST, order_elec)
    ! initializes the address of integrals from \fn(delta_hket)
    addr_geo = 0
    ! loops over the orders of geometric derivatives
    do order_geo = orders_geo_pot(1), orders_geo_pot(2)
      do z_geo = 0, order_geo
        do y_geo = 0, order_geo-z_geo
          x_geo = order_geo-(z_geo+y_geo)
          addr_geo = addr_geo+1
          ! initializes the address of integrals from \fn(delta_hket)
          addr_hket = 0
          do order_hket = orders_hgto_ket(1), orders_hgto_ket(2)
            do z_hket = 0, order_hket
              do y_hket = 0, order_hket-z_hket
                x_hket = order_hket-(z_hket+y_hket)
                ! integral from \fn(rec_delta_hket)
                recur_pint = recur_delta_hket(order_hbra, (/x_hket,y_hket,z_hket/), &
                                              (/x_geo,y_geo,z_geo/))
                addr_hket = addr_hket+1
                call check_difference(recur_pint, hket_pints(1,addr_hket,addr_geo), &
                                      different)
                if (different) then
                  call html_log_int_array(log_text="Orders of HGTOs on ket center:", &
                                          log_int=(/x_hket,y_hket,z_hket/),          &
                                          fmt_int="I3", io_log=io_log)
                  call html_log_int_array(log_text="Orders of geometric derivatives:", &
                                          log_int=(/x_geo,y_geo,z_geo/), fmt_int="I3", &
                                          io_log=io_log)
                  call html_log_real_number(                          &
                         log_text="Reference from recur_delta_hket:", &
                         log_real=recur_pint, fmt_real="Es20.12",     &
                         io_log=io_log, font_color="blue")
                  call html_log_real_number(log_text="Result from delta_hket:",        &
                                            log_real=hket_pints(1,addr_hket,addr_geo), &
                                            fmt_real="Es20.12", io_log=io_log,         &
                                            font_color="red")
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
    call html_log_real_number(log_text="Time (s) used for recur_delta_hket:", &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(hket_pints)
    return
  end subroutine sub_test_delta_hket
