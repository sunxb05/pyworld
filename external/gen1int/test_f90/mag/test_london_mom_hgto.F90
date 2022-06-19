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
!!  Tests subroutines in london_mom_hgto.F90.
!!
!!  2011-10-26, Bin Gao
!!  * first version

  !> \brief tests subroutines in london_mom_hgto.F90
  !> \author Bin Gao
  !> \date 2011-10-26
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_london_mom_hgto(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer orders_hgto_ket(2)  !orders of HGTOs on ket center
    integer orders_mom(2)       !orders of Cartesian multipole moments
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_real_array(log_text="Coordinates of ket center:", &
                             log_real=COORD_KET, fmt_real="F16.8", io_log=io_log)
    call html_log_real_number("Orbital exponent of ket center:", &
                              EXPONENT_KET, "F16.8", io_log)
    call html_log_real_array(log_text="Coordinates of dipole origin:", &
                             log_real=DIPORG, fmt_real="F16.8", io_log=io_log)
    call html_log_real_number("Scaling constant on operator:", &
                              SCAL_CONST, "F16.8", io_log)
    ! the following tests are for individual cases in \fn(london_mom_hgto)
    orders_hgto_ket = (/0,6/)
    orders_mom = (/0,0/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/0,6/)
    orders_mom = (/0,1/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/1,6/)
    orders_mom = (/1,1/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/0,5/)
    orders_mom = (/0,2/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/0,5/)
    orders_mom = (/1,2/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/1,5/)
    orders_mom = (/2,2/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/0,4/)
    orders_mom = (/0,3/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/1,4/)
    orders_mom = (/1,3/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/1,3/)
    orders_mom = (/0,6/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    orders_hgto_ket = (/1,3/)
    orders_mom = (/2,6/)
    call sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    return
  end subroutine test_london_mom_hgto

  !> \brief tests subroutines in london_mom_hgto.F90 with specific orders of
  !>        HGTOs on ket center and Cartesian multipole moments
  !> \author Bin Gao
  !> \date 2011-10-26
  !> \param orders_hgto_ket contains the orders of HGTOs on ket center
  !> \param orders_mom contains the orders of Cartesian multipole moments
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_london_mom_hgto(orders_hgto_ket, orders_mom, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions related to magnetic derivatives
    use recur_mag
    implicit none
    integer, intent(in) :: orders_hgto_ket(2)
    integer, intent(in) :: orders_mom(2)
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    ! see src/london_mom_hgto.F90 for the meaning of the following variables
    integer, parameter :: dim_hgto_bra = 3
    integer dim_hket_zero
    integer, parameter :: num_opt = 6
    integer, parameter :: num_geo_derv = 3
    real(REALK), allocatable :: hgto_pints(:,:,:,:)
    integer dim_hgto_ket
    integer dim_cart_mom
    real(REALK), allocatable :: lmom_pints(:,:,:,:,:)
    integer orders_hket_zero(2)     !orders of HGTOs on ket center with zeroth order Cartesian multipole moment
    integer order_hket              !incremental recorder over orders of HGTOs on ket center
    integer x_hket, y_hket, z_hket  !orders of xyz components of HGTOs on ket center
    integer order_mom               !incremental recorder over orders of Cartesian multipole moment
    integer x_mom, y_mom, z_mom     !orders of xyz components of Cartesian multipole moment
    real(REALK) recur_pint          !results from recursive function \fn(recur_london_mom_hgto)
    integer addr_hket, addr_mom     !addresses of \var(lmom_pints)
    integer ibra, iopt, igeo        !incremental recorders
    integer ierr                    !error information
    real(REALK) begin_time          !begin of CPU time
    real(REALK) end_time            !end of CPU time
    logical different               !if result from \fn(london_mom_hgto) is different from reference
    real(REALK) xyz_function        !function returning artificial value for given xyz orders
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on ket center:",            &
                            log_int=orders_hgto_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(log_text="Orders of Cartesian multipole moments:", &
                            log_int=orders_mom, fmt_int="I3", io_log=io_log,   &
                            separation="&rarr;")
    ! gets the begin time
    call xtimer_set(begin_time)
    ! sets the orders and dimension of integrals with zeroth order Cartesian multipole moment
    orders_hket_zero(1) = max(0,orders_hgto_ket(1)-orders_mom(2))
    orders_hket_zero(2) = orders_hgto_ket(2)+orders_mom(2)
    if (orders_hket_zero(1)==0) then
      dim_hket_zero = (orders_hket_zero(2)+1)*(orders_hket_zero(2)+2) &
                    * (orders_hket_zero(2)+3)/6
    else
      dim_hket_zero = ((orders_hket_zero(2)+1)*(orders_hket_zero(2)+2) &
                    *  (orders_hket_zero(2)+3)                         &
                    -  orders_hket_zero(1)*(orders_hket_zero(1)+1)     &
                    *  (orders_hket_zero(1)+2))/6
    end if
    ! allocates the memory for integrals with zeroth order Cartesian multipole moment
    allocate(hgto_pints(dim_hgto_bra,dim_hket_zero,num_opt,num_geo_derv), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hgto_pints:", &
                               log_int=dim_hgto_bra*dim_hket_zero         &
                                       *num_opt*num_geo_derv,             &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! sets the dimensions of integrals with required orders of Cartesian multipole moments
    if (orders_hgto_ket(1)==0) then
      dim_hgto_ket = (orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) & 
                   * (orders_hgto_ket(2)+3)/6
    else
      dim_hgto_ket = ((orders_hgto_ket(2)+1)*(orders_hgto_ket(2)+2) & 
                   *  (orders_hgto_ket(2)+3)                        &
                   -  orders_hgto_ket(1)*(orders_hgto_ket(1)+1)     &
                   *  (orders_hgto_ket(1)+2))/6
    end if
    if (orders_mom(1)==0) then 
      dim_cart_mom = (orders_mom(2)+1)*(orders_mom(2)+2)*(orders_mom(2)+3)/6
    else
      dim_cart_mom = ((orders_mom(2)+1)*(orders_mom(2)+2)*(orders_mom(2)+3) &
                   -  orders_mom(1)*(orders_mom(1)+1)*(orders_mom(1)+2))/6
    end if
    ! allocates memory for the integrals with required orders of Cartesian multipole moments
    allocate(lmom_pints(dim_hgto_bra,dim_hgto_ket,num_opt,dim_cart_mom, &
                        num_geo_derv), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate lmom_pints:", &
                               log_int=dim_hgto_bra*dim_hgto_ket*num_opt  &
                                       *dim_cart_mom*num_geo_derv,        &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(hgto_pints)
      test_failed = .true.
      return
    end if
    ! sets the faked HGTO integrals with zeroth order Cartesian multipole moment
    do igeo = 1, num_geo_derv
      do iopt = 1, num_opt
        addr_hket = 0
        do order_hket = orders_hket_zero(1), orders_hket_zero(2)
          do z_hket = 0, order_hket
            do y_hket = 0, order_hket-z_hket
              x_hket = order_hket-(z_hket+y_hket)
              recur_pint = xyz_function((/x_hket,y_hket,z_hket/))
              addr_hket = addr_hket+1
              do ibra = 1, dim_hgto_bra
                hgto_pints(ibra,addr_hket,iopt,igeo) = recur_pint
              end do
            end do
          end do
        end do
      end do
    end do
    ! calculates the integrals using \fn(london_mom_hgto)
    call london_mom_hgto(orders_hgto_ket, orders_mom, COORD_KET, EXPONENT_KET, &
                         DIPORG, dim_hgto_bra, dim_hket_zero, num_opt,         &
                         num_geo_derv, hgto_pints, dim_hgto_ket, dim_cart_mom, &
                         lmom_pints)
    deallocate(hgto_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for london_mom_hgto:",  &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the data used for the module \fn(recur_mag)
    call recur_mag_create(COORD_KET, EXPONENT_KET, DIPORG)
    do igeo = 1, num_geo_derv
      ! initializes the address of integrals from \fn(london_mom_hgto)
      addr_mom = 0
      ! loops over the orders of Cartesian multipole moments
      do order_mom = orders_mom(1), orders_mom(2)
        ! loops over the xyz components of Cartesian multipole moments
        do z_mom = 0, order_mom
          do y_mom = 0, order_mom-z_mom
            x_mom = order_mom-(z_mom+y_mom)
            addr_mom = addr_mom+1
            do iopt = 1, num_opt
              ! initializes the address of integrals from \fn(london_mom_hgto)
              addr_hket = 0
              ! loops over the orders of HGTOs on ket center
              do order_hket = orders_hgto_ket(1), orders_hgto_ket(2)
                ! loops over the xyz components of HGTOs on ket center
                do z_hket = 0, order_hket
                  do y_hket = 0, order_hket-z_hket
                    x_hket = order_hket-(z_hket+y_hket)
                    addr_hket = addr_hket+1
                    do ibra = 1, dim_hgto_bra
                      ! integral from \fn(recur_london_mom_hgto)
                      recur_pint = recur_london_mom_hgto((/x_mom,y_mom,z_mom/), &
                                                         (/x_hket,y_hket,z_hket/))
                      call check_difference(recur_pint,                     &
                             lmom_pints(ibra,addr_hket,iopt,addr_mom,igeo), &
                             different)
                      if (different) then
                        call html_log_int_array(log_text="Orders of HGTO:",       &
                                                log_int=(/x_hket,y_hket,z_hket/), &
                                                fmt_int="I3", io_log=io_log)
                        call html_log_int_array(                                 &
                               log_text="Orders of Cartesian multipole moment:", &
                               log_int=(/x_mom,y_mom,z_mom/), fmt_int="I3",      &
                               io_log=io_log)
                        call html_log_real_number(                               &
                               log_text="Reference from recur_london_mom_hgto:", &
                               log_real=recur_pint, fmt_real="Es20.12",          &
                               io_log=io_log, font_color="blue")
                        call html_log_real_number(                                     &
                               log_text="Result from london_mom_hgto:",                &
                               log_real=lmom_pints(ibra,addr_hket,iopt,addr_mom,igeo), &
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
    deallocate(lmom_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_london_mom_hgto:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",      &
                              io_log=io_log, font_color="blue")
    return
  end subroutine sub_test_london_mom_hgto
