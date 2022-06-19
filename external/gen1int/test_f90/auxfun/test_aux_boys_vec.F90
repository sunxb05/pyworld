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
!!  Tests subroutines in aux_boys_vec.F90.
!!
!!  2010-10-08, Bin Gao
!!  * first version

  !> \brief tests subroutines in aux_boys_vec.F90
  !> \author Bin Gao
  !> \date 2010-10-08
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_aux_boys_vec(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters and referenced results of Boys functions
#include "ref_boys.h"
    real(REALK) arg_boys                    !argument of Boys functions
    real(REALK) val_boys(0:MAX_ORDER_BOYS)  !Boys functions
    integer istep                           !incremental recorder of steps
    integer iorder                          !incremental recorder of orders
    integer iref                            !incremental recorder of referenced Boys functions
    logical different                       !if result from \fn(aux_boys_vec) is different from reference
    real(REALK) begin_time                  !begin of CPU time
    real(REALK) end_time                    !end of CPU time
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_int_number("Maximum order of Boys functions:", &
                             MAX_ORDER_BOYS, "I4", io_log)
    call html_log_real_number("Minimum argument of Boys functions:", &
                              MIN_ARG_BOYS, "F10.6", io_log)
    call html_log_real_number("Step of argument of Boys functions:", &
                              STEP_ARG_BOYS, "F10.6", io_log)
    call html_log_int_number("Number of steps:", NUM_STEPS, "I4", io_log)
    iref = 0
    ! gets the begin time
    call xtimer_set(begin_time)
    ! loops over sampling points
    do istep = 0, NUM_STEPS
      arg_boys = MIN_ARG_BOYS+STEP_ARG_BOYS*real(istep,REALK)
      call aux_boys_vec(0, MAX_ORDER_BOYS, arg_boys, val_boys)
      ! checks the results
      do iorder = 0, MAX_ORDER_BOYS
        iref = iref+1
        call check_difference(REF_BOYS(iref), val_boys(iorder), different)
        if (different) then
          call html_log_int_number("Order of Boys function:", iorder, "I4", io_log)
          call html_log_real_number("Argument of Boys function:", arg_boys, &
                                    "F10.6", io_log)
          call html_log_real_number(log_text="Referenced Boys function:",       &
                                    log_real=REF_BOYS(iref), fmt_real="F20.16", &
                                    io_log=io_log, font_color="red")
          call html_log_real_number(log_text="Result from aux_boys_vec:",         &
                                    log_real=val_boys(iorder), fmt_real="F20.16", &
                                    io_log=io_log, font_color="red")
          test_failed = .true.
        end if
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for aux_boys_vec:",     &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    return
  end subroutine test_aux_boys_vec
