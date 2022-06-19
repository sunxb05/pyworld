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
!!  Tests subroutine prim_hgto_value.
!!
!!  2012-03-10, Bin Gao
!!  * first version

  !> \brief tests subroutine prim_hgto_value
  !> \author Bin Gao
  !> \date 2012-03-10
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_prim_hgto_value(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer orders_hgto_bra(2)  !orders of HGTOs on bra center
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_real_array(log_text="Coordinates of bra center:", &
                             log_real=COORD_BRA, fmt_real="F16.8", io_log=io_log)
    call html_log_real_number("Orbital exponent of bra center:", &
                              EXPONENT_BRA, "F16.8", io_log)
    ! the following tests are for individual cases in \fn(prim_hgto_value)
    orders_hgto_bra = (/0,0/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/0,1/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/1,1/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/0,2/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/1,2/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/2,2/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/0,4/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/3,4/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/3,5/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/0,6/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/1,6/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    orders_hgto_bra = (/2,6/)
    call sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    return
  end subroutine test_prim_hgto_value

  !> \brief tests subroutine prim_hgto_value with specific orders of HGTOs
  !> \author Bin Gao
  !> \date 2012-03-10
  !> \param orders_hgto_bra contains the orders of HGTOs on bra center
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_prim_hgto_value(orders_hgto_bra, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions related to HGTOs
    use recur_hgto
    implicit none
    integer, intent(in) :: orders_hgto_bra(2)
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer num_points                            !number of grid points
    real(REALK), allocatable :: grid_points(:,:)  !coordinates of grid points
    integer dim_hgto_bra                          !number of Hermite Gaussians on bra center
    real(REALK), allocatable :: hgto_value(:,:)   !overlap distribution using primitive Hermite Gaussians
    integer order_hbra              !order of HGTO on bra center
    integer x_hbra, y_hbra, z_hbra  !orders of xyz components of HGTO on bra center
    real(REALK) elec_coord_bra(3)   !electron coordinates on bra center
    real(REALK) recur_value         !value of HGTO on bra center from \fn(recur_hgto_value)
    integer ipoint                  !incremental recorder over grid points
    integer addr_hbra               !address of \var(hgto_value)
    integer ierr                    !error information
    real(REALK) begin_time          !begin of CPU time
    real(REALK) end_time            !end of CPU time
    logical different               !if result from \fn(prim_hgto_value) is different from reference
    ! dumps information of tests
    call html_log_int_array(log_text="Orders of HGTOs on bra center:",            &
                            log_int=orders_hgto_bra, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    ! gets the begin time
    call xtimer_set(begin_time)
    ! sets grid points
    num_points = 100
    allocate(grid_points(3,num_points), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate grid_points:", &
                               log_int=3*num_points, fmt_int="I12",        &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    call random_number(grid_points)
    ! allocates memory for \fn(prim_hgto_value)
    if (orders_hgto_bra(1)>0) then
      dim_hgto_bra = ((orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                      *(orders_hgto_bra(2)+3)                       &
                   -  orders_hgto_bra(1)*(orders_hgto_bra(1)+1)     &
                      *(orders_hgto_bra(1)+2))/6
    else
      dim_hgto_bra = (orders_hgto_bra(2)+1)*(orders_hgto_bra(2)+2) &
                   * (orders_hgto_bra(2)+3)/6
    end if
    allocate(hgto_value(dim_hgto_bra,num_points), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hgto_value:",      &
                               log_int=dim_hgto_bra*num_points, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      deallocate(grid_points)
      return
    end if
    call prim_hgto_value(orders_hgto_bra, COORD_BRA, EXPONENT_BRA, &
                         num_points, grid_points, dim_hgto_bra,    &
                         hgto_value)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for prim_hgto_value:",  &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! loops over grid points
    do ipoint = 1, num_points
      ! sets the electron coordinates
      elec_coord_bra = grid_points(:,ipoint)-COORD_BRA
      ! initializes the address of overlap distribution from \fn(prim_hgto_value)
      addr_hbra = 0
      ! loops over the orders of HGTOs on bra center
      do order_hbra = orders_hgto_bra(1), orders_hgto_bra(2)
        do z_hbra = 0, order_hbra
          do y_hbra = 0, order_hbra-z_hbra
            x_hbra = order_hbra-(z_hbra+y_hbra)
            addr_hbra = addr_hbra+1
            ! value of HGTO on bra center from \fn(recur_hgto_value)
            recur_value = recur_hgto_value((/x_hbra,y_hbra,z_hbra/), &
                                           -0.5_REALK/EXPONENT_BRA,  &
                                           elec_coord_bra)
            call check_difference(recur_value, hgto_value(addr_hbra,ipoint), &
                                  different)
            if (different) then
              call html_log_int_array(log_text="Orders of HGTO on bra:", &
                                      log_int=(/x_hbra,y_hbra,z_hbra/),  &
                                      fmt_int="I3", io_log=io_log)
              call html_log_real_number(                          &
                     log_text="Reference from recur_hgto_value:", &
                     log_real=recur_value, fmt_real="Es20.12",    &
                     io_log=io_log, font_color="blue")
              call html_log_real_number(                      &
                     log_text="Result from prim_hgto_value:", &
                     log_real=hgto_value(addr_hbra,ipoint),   &
                     fmt_real="Es20.12", io_log=io_log, font_color="red")
              test_failed = .true.
            end if
          end do
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_hgto_value:", &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(grid_points)
    deallocate(hgto_value)
    return
  end subroutine sub_test_prim_hgto_value
