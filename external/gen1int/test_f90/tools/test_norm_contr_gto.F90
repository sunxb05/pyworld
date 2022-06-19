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
!!  Tests normalization of contracted GTOs.
!!
!!  2011-04-04, Bin Gao
!!  * first version

  !> \brief tests normalization of contracted GTOs
  !> \author Bin Gao
  !> \date 2011-04-04
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_norm_contr_gto(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! water molecule using cc-pV6Z basis sets
#include "test_f90/water_ccpv6z.h"
    real(REALK) nrm_contr_coef(TOTAL_NCONTR)  !normalized contracted coefficients of each block
    integer iblock                            !incremental recorder over blocks
    integer icoef                             !incremental recorder over coefficients
    logical different                         !if result from \fn(norm_contr_cgto) or
                                              !\fn(norm_contr_sgto) is different from reference
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_text("Tests water molecule using cc-pV6Z basis sets ...", io_log)
    call html_log_int_number("Number of blocks of sub-shells:", &
                             NUM_SUB_BLOCK, "I6", io_log)
    call html_log_heading(heading="Tests normalization of contracted CGTOs", &
                          io_log=io_log, level=4)
    ! gets the original contraction coefficients of each block
    nrm_contr_coef = CONTR_COEF_BLOCK
    ! loops over blocks of sub-shells
    do iblock = 1, NUM_SUB_BLOCK
      call norm_contr_cgto(ANGULAR_BLOCK(iblock), NUM_PRIM_BLOCK(iblock), &
                           EXPONENT_BLOCK(START_PRIM_BLOCK(iblock)+1      &
                                          :START_PRIM_BLOCK(iblock+1)),   &
                           NUM_CONTR_BLOCK(iblock),                       &
                           nrm_contr_coef(START_CONTR_BLOCK(iblock)+1     &
                                          :START_CONTR_BLOCK(iblock+1)))
      ! checks
      do icoef = START_CONTR_BLOCK(iblock)+1, START_CONTR_BLOCK(iblock+1)
        call check_difference(REF_NRM_CGTO(icoef), nrm_contr_coef(icoef), different)
        if (different) then
          call html_log_int_array(log_text="ID of block and coefficient:", &
                                  log_int=(/iblock,icoef/), fmt_int="I6",  &
                                  io_log=io_log)
          call html_log_real_number(                                   &
                 log_text="Referenced coefficient in water_ccpv6z.h:", &
                 log_real=REF_NRM_CGTO(icoef), fmt_real="Es16.8",      &
                 io_log=io_log, font_color="blue")
          call html_log_real_number(log_text="Result from norm_contr_cgto:", &
                                    log_real=nrm_contr_coef(icoef),          &
                                    fmt_real="Es16.8", io_log=io_log, font_color="red")
          test_failed = .true.
        end if
      end do
    end do
    call html_log_heading(heading="Tests normalization of contracted SGTOs", &
                          io_log=io_log, level=4)
    ! gets the original contraction coefficients of each block
    nrm_contr_coef = CONTR_COEF_BLOCK
    ! loops over blocks of sub-shells
    do iblock = 1, NUM_SUB_BLOCK
      call norm_contr_sgto(ANGULAR_BLOCK(iblock), NUM_PRIM_BLOCK(iblock), &
                           EXPONENT_BLOCK(START_PRIM_BLOCK(iblock)+1      &
                                          :START_PRIM_BLOCK(iblock+1)),   &
                           NUM_CONTR_BLOCK(iblock),                       &
                           nrm_contr_coef(START_CONTR_BLOCK(iblock)+1     &
                                          :START_CONTR_BLOCK(iblock+1)))
      ! checks
      do icoef = START_CONTR_BLOCK(iblock)+1, START_CONTR_BLOCK(iblock+1)
        call check_difference(REF_NRM_SGTO(icoef), nrm_contr_coef(icoef), different)
        if (different) then
          call html_log_int_array(log_text="ID of block and coefficient:", &
                                  log_int=(/iblock,icoef/), fmt_int="I6",  &
                                  io_log=io_log)
          call html_log_real_number(                                   &
                 log_text="Referenced coefficient in water_ccpv6z.h:", &
                 log_real=REF_NRM_SGTO(icoef), fmt_real="Es16.8",      &
                 io_log=io_log, font_color="blue")
          call html_log_real_number(log_text="Result from norm_contr_sgto:", &
                                    log_real=nrm_contr_coef(icoef),          &
                                    fmt_real="Es16.8", io_log=io_log, font_color="red")
          test_failed = .true.
        end if
      end do
    end do
    return
  end subroutine test_norm_contr_gto
