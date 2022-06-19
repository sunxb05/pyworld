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
!!  This file calculates the partial geometric derivatives using finite difference.
!!
!!  2012-04-15, Bin Gao
!!  * first version

#include "xkind.h"

    ! step length for finite difference. Power of 2 to ensure
    ! addition and subsequent subtraction yields original
    real(8), parameter   :: fdistep = 2d0**(-25)

  !> \brief tests subroutine prim_hgto_nucpot using finite difference with specific centers
  !> \author Bin Gao
  !> \date 2012-04-15
  !> \param coord_bra contains the coordinates of bra center
  !> \param exponent_bra is the exponent of primitive Gaussian of bra center
  !> \param coord_ket contains the coordinates of ket center
  !> \param exponent_ket is the exponent of primitive Gaussian of ket center
  !> \param nucpot_origin contains the coordinates of nuclear potential origin
  !> \param dipole_origin contains the coordinates of dipole origin
  !> \param scal_const is the scale constant for nuclear attraction potential
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine geom_part_fdm(coord_bra, exponent_bra,      &
                                     coord_ket, exponent_ket,      &
                                     nucpot_origin, dipole_origin, &
                                     scal_const, io_log, test_failed)
    ! module of HTML test log routines
    use html_log
    implicit none
    real(REALK), intent(in) :: coord_bra(3)
    real(REALK), intent(in) :: exponent_bra
    real(REALK), intent(in) :: coord_ket(3)
    real(REALK), intent(in) :: exponent_ket
    real(REALK), intent(in) :: nucpot_origin(3)
    real(REALK), intent(in) :: dipole_origin(3)
    real(REALK), intent(in) :: scal_const
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of finite difference
#include "finite_difference.h"
    integer orders_hgto_bra(2)      !orders of HGTOs on bra center
    integer orders_hgto_ket(2)      !orders of HGTOs on ket center
    integer order_geo_pot           !order of geometric derivatives on nuclear potential origin
    integer dim_hgto_bra            !dimension of HGTOs on bra center
    integer dim_hgto_ket            !dimension of HGTOs on ket center
    integer num_elec                !number of electronic derivatives
    integer num_mom                 !number of Cartesian multipole moments
    integer num_geo_pot             !number of geometric derivatives on nuclear potential origin
    real(REALK), allocatable :: low_pints(:,:,:,:,:,:)
    real(REALK), allocatable :: up_pints(:,:,:,:,:)
    real(REALK), allocatable :: fd_pints(:,:,:,:,:)
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
                                    !addresses of \var(low_pints)
    integer ierr                    !error information
    real(REALK) begin_time          !begin of CPU time
    real(REALK) end_time            !end of CPU time
    logical different               !if result from \fn(prim_hgto_nucpot) is different from reference

    ! dumps the orders of electronic derivatives and Cartesian multipole moment
    num_elec = (FD_ORDER_ELEC+1)*(FD_ORDER_ELEC+2)/2
    call html_log_int_number("Order of electronic derivatives:", &
                             FD_ORDER_ELEC, "I3", io_log)
    num_mom = (FD_ORDER_MOM+1)*(FD_ORDER_MOM+2)/2
    call html_log_int_number("Order of Cartesian multipole moment:", FD_ORDER_MOM, &
                             "I3", io_log)
    ! tests HGTOs on bra center using finite difference
    call html_log_text("Tests HGTOs on bra center ...", io_log)
    ! sets the orders of HGTOs on ket center and geometric derivative on nuclear potential origin
    orders_hgto_ket = 0
    dim_hgto_ket = 1
    call html_log_int_array(log_text="Orders of HGTOs on ket center:",            &
                            log_int=orders_hgto_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    order_geo_pot = 0
    num_geo_pot = 1
    call html_log_int_number(                                            &
           "Order of geometric derivative on nuclear potential origin:", &
           order_geo_pot, "I3", io_log)
    ! gets the begin time
    call xtimer_set(begin_time)
    ! loops till the maximum order to test
    do iorder = 0, FD_MAX_ORDER-1
      ! sets the order of HGTOs on bra center
      orders_hgto_bra = iorder
      call html_log_int_array(log_text="Orders of HGTOs on bra center:",            &
                              log_int=orders_hgto_bra, fmt_int="I3", io_log=io_log, &
                              separation="&rarr;")
      ! allocates memory for the integrals from \fn(prim_hgto_nucpot)
      dim_hgto_bra = (iorder+1)*(iorder+2)/2
      allocate(low_pints(dim_hgto_bra,dim_hgto_ket,num_elec, &
                          num_mom,num_geo_pot), stat=ierr)
      if (ierr/=0) then
        call html_log_int_number(log_text="Failed to allocate low_pints:", &
                                 log_int=dim_hgto_bra*dim_hgto_ket*num_elec &
                                         *num_mom*num_geo_pot,              &
                                 fmt_int="I12", io_log=io_log, font_color="red")
        test_failed = .true.
        return
      end if
      call prim_hgto_nucpot(orders_hgto_bra, coord_bra, exponent_bra,    &
                            orders_hgto_ket, coord_ket, exponent_ket,    &
                            FD_ORDER_ELEC, nucpot_origin, dipole_origin, &
                            scal_const, FD_ORDER_MOM, order_geo_pot,     &
                            dim_hgto_bra, dim_hgto_ket, num_elec,        &
                            num_mom, num_geo_pot, low_pints)
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for prim_hgto_nucpot:", &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
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
        do z_mom = 0, FD_ORDER_MOM
          do y_mom = 0, FD_ORDER_MOM-z_mom
            x_mom = FD_ORDER_MOM-(z_mom+y_mom)
            addr_mom = addr_mom+1
            ! initializes the address of integrals from \fn(prim_hgto_nucpot)
            addr_elec = 0
            ! loops over xyz components of electronic derivative
            do z_elec = 0, FD_ORDER_ELEC
              do y_elec = 0, FD_ORDER_ELEC-z_elec
                x_elec = FD_ORDER_ELEC-(z_elec+y_elec)
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
                                   low_pints(addr_hbra,addr_hket,addr_elec, &
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
                                     log_real=low_pints(addr_hbra,addr_hket,addr_elec, &
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
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recursive functions:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",    &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(low_pints)

    ! tests HGTOs on ket center using finite difference
    call html_log_text("Tests HGTOs on ket center ...", io_log)

    ! tests geometric derivatives on nuclear potential origin using finite difference
    call html_log_text("Tests geometric derivatives on nuclear potential origin ...", &
                       io_log)

    return
  end subroutine geom_part_fdm
