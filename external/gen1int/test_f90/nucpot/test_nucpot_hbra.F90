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
!!  Tests subroutines in nucpot_hbra.F90.
!!
!!  2011-10-18, Bin Gao
!!  * first version

  !> \brief tests subroutines in nucpot_hbra.F90
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_nucpot_hbra(io_log, test_failed)
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
    integer orders_geo_pot(2)             !orders of geometric derivatives on nuclear potential origin
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
    ! the following tests are for individual cases in \fn(nucpot_hbra)
    orders_hgto_bra = (/0,0/)
    orders_hgto_ket = (/0,0/)
    orders_geo_pot = (/0,6/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    orders_hgto_ket = (/0,1/)
    orders_geo_pot = (/0,6/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    orders_hgto_ket = (/1,1/)
    orders_geo_pot = (/0,6/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    orders_hgto_ket = (/0,2/)
    orders_geo_pot = (/0,6/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    orders_hgto_ket = (/1,2/)
    orders_geo_pot = (/0,0/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    orders_hgto_ket = (/2,2/)
    orders_geo_pot = (/0,6/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/0,3/)
    orders_hgto_ket = (/0,3/)
    orders_geo_pot = (/0,6/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/1,4/)
    orders_hgto_ket = (/1,4/)
    orders_geo_pot = (/0,6/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    orders_hgto_bra = (/3,5/)
    orders_hgto_ket = (/3,5/)
    orders_geo_pot = (/0,3/)
    call sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                              ORDER_ELEC, io_log, test_failed)
    return
  end subroutine test_nucpot_hbra

  !> \brief tests subroutines in nucpot_hbra.F90 with specific orders of
  !>        Hermite Gaussians and geometric derivatives
  !> \author Bin Gao
  !> \date 2011-10-18
  !> \param orders_hgto_bra contains the orders of HGTOs on bra center
  !> \param orders_hgto_ket contains the orders of HGTOs on ket center
  !> \param orders_geo_pot is the order of geometric derivatives on nuclear potential origin
  !> \param order_elec is the order of electronic derivatives
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                                  order_elec, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions of nuclear attraction potential integrals
    use recur_nucpot
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: orders_geo_pot(2)
    integer, intent(in) :: order_elec
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer max_order_nucpot        !maximum order of geometric derivatives for \fn(nucpot_geom)
    integer dim_nucpot              !dimension of geometric derivatives of nuclear potential
                                    !origin with zeroth order HGTOs
    real(REALK), allocatable :: nucpot_pints(:)
                                    !geometric derivatives of nuclear potential origin with zeroth order HGTOs
    integer min_hgto_ket            !minimum order of HGTOs on ket center after recovering HGTOs on ket center
    integer max_geo_hket            !maximum order of geometric derivatives after recovering HGTOs on ket center
    integer dim_hket                !dimension of HGTOs on ket center after recovering HGTOs on ket center
    integer dim_geo_hket            !dimension of geometric derivatives after recovering HGTOs on ket center
    real(REALK), allocatable :: hket_pints(:,:)
                                    !integrals after recovering HGTOs on ket center
    integer dim_hgto_bra            !dimension of HGTOs on bra center after recovering HGTOs on bra center
    integer dim_hgto_ket            !dimension of HGTOs on ket center after recovering HGTOs on bra center
    integer dim_geo_hbra            !dimension of  geometric derivatives after recovering HGTOs on bra center
    real(REALK), allocatable :: hbra_pints(:,:,:)
                                    !integrals after recovering HGTOs on bra center
    integer order_hbra              !order of HGTOs on bra center
    integer x_hbra, y_hbra, z_hbra  !orders of xyz components of HGTOs on bra center
    integer order_hket              !order of HGTOs on ket center
    integer x_hket, y_hket, z_hket  !orders of xyz components of HGTOs on ket center
    integer order_geo               !order of geometric derivatives
    integer x_geo, y_geo, z_geo     !orders of xyz components of geometric derivatives
    real(REALK) recur_pint          !results from recursive function \fn(recur_nucpot_hbra)
    integer addr_hbra, addr_hket, addr_geo
                                    !addresses of \var(hket_pints)
    integer ierr                    !error information
    real(REALK) begin_time          !begin of CPU time
    real(REALK) end_time            !end of CPU time
    logical different               !if result from \fn(nucpot_hbra) is different from reference
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on bra center:",            &
                            log_int=orders_hgto_bra, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(log_text="Orders of HGTOs on ket center:",            &
                            log_int=orders_hgto_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(                                                       &
           log_text="Orders of geometric derivative on nuclear potential origin:", &
           log_int=orders_geo_pot, fmt_int="I3", io_log=io_log, separation="&rarr;")
    ! gets the begin time
    call xtimer_set(begin_time)
    ! sets the maximum order and dimension of geometric derivatives for \fn(nucpot_geom)
    max_order_nucpot = orders_hgto_bra(2)+orders_hgto_ket(2)+orders_geo_pot(2)
    if (orders_geo_pot(1)>0) then
      dim_nucpot = ((max_order_nucpot+1)*(max_order_nucpot+2)*(max_order_nucpot+3) &
                 -  orders_geo_pot(1)*(orders_geo_pot(1)+1)*(orders_geo_pot(1)+2))/6
    else
      dim_nucpot = (max_order_nucpot+1)*(max_order_nucpot+2)*(max_order_nucpot+3)/6
    end if
    ! allocates memory for the integrals from \fn(nucpot_geom)
    allocate(nucpot_pints(dim_nucpot), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate nucpot_pints:", &
                               log_int=dim_nucpot, fmt_int="I12",           &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    call nucpot_geom(COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET,           &
                     NUCORG, SCAL_CONST, (/orders_geo_pot(1),max_order_nucpot/), &
                     order_elec, dim_nucpot, nucpot_pints)
    ! allocates memory for the integrals from \fn(nucpot_hket)
    min_hgto_ket = max(orders_hgto_ket(1)-orders_hgto_bra(2),0)
    if (min_hgto_ket>0) then
      dim_hket = ((orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2)*(orders_hgto_ket(2)+3) &
               -  min_hgto_ket*(min_hgto_ket+1)*(min_hgto_ket+2))/6
    else
      dim_hket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2)*(orders_hgto_ket(2)+3)/6
    end if
    max_geo_hket = orders_hgto_bra(2)+orders_geo_pot(2)
    if (orders_geo_pot(1)>0) then
      dim_geo_hket = ((max_geo_hket+1)*(max_geo_hket+2)*(max_geo_hket+3) &
                   -  orders_geo_pot(1)*(orders_geo_pot(1)+1)*(orders_geo_pot(1)+2))/6
    else
      dim_geo_hket = (max_geo_hket+1)*(max_geo_hket+2)*(max_geo_hket+3)/6
    end if
    allocate(hket_pints(dim_hket,dim_geo_hket), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hket_pints:",    &
                               log_int=dim_hket*dim_geo_hket, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      deallocate(nucpot_pints)
      test_failed = .true.
      return
    end if
    call nucpot_hket((/min_hgto_ket,orders_hgto_ket(2)/),              &
                     (/orders_geo_pot(1),max_geo_hket/),               &
                     COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET, &
                     dim_nucpot, nucpot_pints, dim_hket, dim_geo_hket, &
                     hket_pints)
    deallocate(nucpot_pints)
    ! allocates memory for the integrals from \fn(nucpot_hbra)
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
    if (orders_geo_pot(1)>0) then
      dim_geo_hbra = ((orders_geo_pot(2)+1)*(orders_geo_pot(2)+2) &
                      *(orders_geo_pot(2)+3)                      &
                    - orders_geo_pot(1)*(orders_geo_pot(1)+1)     &
                      *(orders_geo_pot(1)+2))/6
    else
      dim_geo_hbra = (orders_geo_pot(2)+1)*(orders_geo_pot(2)+2) &
                   * (orders_geo_pot(2)+3)/6
    end if
    allocate(hbra_pints(dim_hgto_bra,dim_hgto_ket,dim_geo_hbra), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hbra_pints:",      &
                               log_int=dim_hgto_bra*dim_hgto_ket*dim_geo_hbra, &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(hket_pints)
      test_failed = .true.
      return
    end if
    call nucpot_hbra(orders_hgto_bra, orders_hgto_ket, orders_geo_pot, &
                     COORD_BRA, EXPONENT_BRA, COORD_KET, EXPONENT_KET, &
                     dim_hket, dim_geo_hket, hket_pints, dim_hgto_bra, &
                     dim_hgto_ket, dim_geo_hbra, hbra_pints)
    deallocate(hket_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(                                                     &
           log_text="Time (s) used for nucpot_geom, nucpot_hket and nucpot_hbra:", &
           log_real=end_time-begin_time, fmt_real="F10.4", io_log=io_log,          &
           font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the data used for the module \fn(recur_nucpot)
    call recur_nucpot_create(coord_bra=COORD_BRA, exponent_bra=EXPONENT_BRA, &
                             coord_ket=COORD_KET, exponent_ket=EXPONENT_KET, &
                             nucpot_origin=NUCORG, dipole_origin=DIPORG,     &
                             scal_const=SCAL_CONST, order_elec=order_elec,   &
                             order_aux=max_order_nucpot)
    ! initializes the address of integrals from \fn(nucpot_hbra)
    addr_geo = 0
    ! loops over the orders of geometric derivatives
    do order_geo = orders_geo_pot(1), orders_geo_pot(2)
      ! loops over the xyz components of geometric derivatives
      do z_geo = 0, order_geo
        do y_geo = 0, order_geo-z_geo
          x_geo = order_geo-(z_geo+y_geo)
          addr_geo = addr_geo+1
          ! initializes the address of integrals from \fn(nucpot_hbra)
          addr_hket = 0
          ! loops over the orders of HGTOs on ket center
          do order_hket = orders_hgto_ket(1), orders_hgto_ket(2)
            ! loops over the xyz components of HGTOs on ket center
            do z_hket = 0, order_hket
              do y_hket = 0, order_hket-z_hket
                x_hket = order_hket-(z_hket+y_hket)
                addr_hket = addr_hket+1
                ! initializes the address of integrals from \fn(nucpot_hbra)
                addr_hbra = 0
                ! loops over the orders of HGTOs on bra center
                do order_hbra = orders_hgto_bra(1), orders_hgto_bra(2)
                  ! loops over the xyz components of HGTOs on bra center
                  do z_hbra = 0, order_hbra
                    do y_hbra = 0, order_hbra-z_hbra
                      x_hbra = order_hbra-(z_hbra+y_hbra)
                      ! integral from \fn(recur_nucpot_hbra)
                      recur_pint = recur_nucpot_hbra((/x_hbra,y_hbra,z_hbra/), &
                                                     (/x_hket,y_hket,z_hket/), &
                                                     (/x_geo,y_geo,z_geo/))
                      addr_hbra = addr_hbra+1
                      call check_difference(recur_pint,                &
                             hbra_pints(addr_hbra,addr_hket,addr_geo), &
                             different)
                      if (different) then
                        call html_log_int_array(log_text="Orders of HGTO on bra:", &
                                                log_int=(/x_hbra,y_hbra,z_hbra/),  &
                                                fmt_int="I3", io_log=io_log)
                        call html_log_int_array(log_text="Orders of HGTO on ket:", &
                                                log_int=(/x_hket,y_hket,z_hket/),  &
                                                fmt_int="I3", io_log=io_log)
                        call html_log_int_array(                               &
                               log_text="Orders of geometric derivative on "// &
                                        "nuclear potential origin:",           &
                               log_int=(/x_geo,y_geo,z_geo/), fmt_int="I3",    &
                               io_log=io_log)
                        call html_log_real_number(                           &
                               log_text="Reference from recur_nucpot_hbra:", &
                               log_real=recur_pint, fmt_real="Es20.12",      &
                               io_log=io_log, font_color="blue")
                        call html_log_real_number(log_text="Result from nucpot_hbra:", &
                               log_real=hbra_pints(addr_hbra,addr_hket,addr_geo),      &
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
    call recur_nucpot_destroy
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_nucpot_hbra:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",  &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(hbra_pints)
    return
  end subroutine sub_test_nucpot_hbra
