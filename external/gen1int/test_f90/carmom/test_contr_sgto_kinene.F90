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
!!  Tests calculating kinetic energy integrals with SGTOs.
!!
!!  2011-04-03, Bin Gao
!!  * first version

  !> \brief tests calculating kinetic energy integrals with SGTOs
  !> \author Bin Gao
  !> \date 2011-04-03
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_contr_sgto_kinene(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! water molecule using cc-pV6Z basis sets
#include "test_f90/water_ccpv6z.h"
! test parameters and referenced results
#include "test_f90/water_ccpv6z_sgto_kinene.h"
    integer itst                          !incremental recorder over test cases
    integer iblock, jblock                !incremental recorders over blocks of sub-shells
    integer, parameter :: NUM_OPT = 6     !number of operators including derivatives
    real(REALK), allocatable :: contr_ints(:,:,:,:,:)
                                          !contracted integrals
    integer idx_bra, idx_ket              !indices of bra and ket centers
    integer num_sgto_bra, num_sgto_ket    !number of spherical Gaussians of bra and ket centers
    integer ierr                          !error information
    integer addr_ref                      !address of referenced results
    real(REALK) begin_time                !begin of CPU time
    real(REALK) end_time                  !end of CPU time
    logical different                     !if result from \fn(contr_sgto_carmom) is different from reference
    integer icontr, jcontr, isgto, jsgto  !incremental recorders
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_text("Tests water molecule using cc-pV6Z spherical GTOs ...", io_log)
    ! loops over test cases
    addr_ref = 0
    do itst = 1, NUM_TESTS
      ! gets the blocks of sub-shells
      iblock = BRA_BLOCK(itst)
      jblock = KET_BLOCK(itst)
      call html_log_int_number("ID of block of sub-shells on bra center:", &
                               iblock, "I4", io_log)
      call html_log_int_number("ID of block of sub-shells on ket center:", &
                               jblock, "I4", io_log)
      ! gets the begin time
      call xtimer_set(begin_time)
      ! gets the indices of bra and ket centers
      idx_bra = IDX_CENT_BLOCK(iblock)
      idx_ket = IDX_CENT_BLOCK(jblock)
      ! sets the number of spherical Gaussians of bra and ket centers
      num_sgto_bra = 2*ANGULAR_BLOCK(iblock)+1
      num_sgto_ket = 2*ANGULAR_BLOCK(jblock)+1
      ! allocates the contracted integrals
      allocate(contr_ints(num_sgto_bra,NUM_CONTR_BLOCK(iblock), &
                          num_sgto_ket,NUM_CONTR_BLOCK(jblock), &
                          NUM_OPT), stat=ierr)
      if (ierr/=0) then
        call html_log_int_number(log_text="Failed to allocate contr_ints:",    &
                                 log_int=num_sgto_bra*NUM_CONTR_BLOCK(iblock)  &
                                         *num_sgto_ket*NUM_CONTR_BLOCK(jblock) &
                                         *NUM_OPT,                             &
                                 fmt_int="I12", io_log=io_log, font_color="red")
        test_failed = .true.
        cycle
      end if
      ! computes the contracted integrals
      call contr_sgto_carmom(idx_bra, COORDS(3*idx_bra-2:3*idx_bra),        &
                             ANGULAR_BLOCK(iblock), NUM_PRIM_BLOCK(iblock), &
                             EXPONENT_BLOCK(START_PRIM_BLOCK(iblock)+1      &
                                            :START_PRIM_BLOCK(iblock+1)),   &
                             NUM_CONTR_BLOCK(iblock),                       &
                             REF_NRM_SGTO(START_CONTR_BLOCK(iblock)+1       &
                                          :START_CONTR_BLOCK(iblock+1)),    &
                             idx_ket, COORDS(3*idx_ket-2:3*idx_ket),        &
                             ANGULAR_BLOCK(jblock), NUM_PRIM_BLOCK(jblock), &
                             EXPONENT_BLOCK(START_PRIM_BLOCK(jblock)+1      &
                                            :START_PRIM_BLOCK(jblock+1)),   &
                             NUM_CONTR_BLOCK(jblock),                       &
                             REF_NRM_SGTO(START_CONTR_BLOCK(jblock)+1       &
                                          :START_CONTR_BLOCK(jblock+1)),    &
                             2, IDX_DIPORG, DIP_ORIGIN, SCAL_CONST, 0,      &
                             0, 0, 0, 0, IDX_CENT(1:0), ORDER_CENT(1:0),    &
                             num_sgto_bra, num_sgto_ket, NUM_OPT, contr_ints)
      ! sums xx, yy and zz components
      contr_ints(:,:,:,:,1) = contr_ints(:,:,:,:,1) &
                            + contr_ints(:,:,:,:,3)+contr_ints(:,:,:,:,6)
      ! gets the end time
      call xtimer_set(end_time)
      ! prints the CPU elapsed time
      call html_log_real_number(log_text="Time (s) used for contr_sgto_carmom:", &
                                log_real=end_time-begin_time, fmt_real="F10.4",  &
                                io_log=io_log, font_color="blue")
      ! checks the results
      do jcontr = 1, NUM_CONTR_BLOCK(jblock)
        do jsgto = 1, num_sgto_ket
          do icontr = 1, NUM_CONTR_BLOCK(iblock)
            do isgto = 1, num_sgto_bra
              addr_ref = addr_ref+1
              call check_difference(REF_CONTR_INTS(addr_ref),                &
                                    contr_ints(isgto,icontr,jsgto,jcontr,1), &
                                    different)
              if (different) then
                call html_log_int_number("ID of test:", itst, "I3", io_log)
                call html_log_int_array(log_text="ID of contraction:",         &
                                        log_int=(/isgto,icontr,jsgto,jcontr/), &
                                        fmt_int="I3", io_log=io_log)
                call html_log_real_number(                             &
                       log_text="Referenced kinetic energy integral:", &
                       log_real=REF_CONTR_INTS(addr_ref),              &
                       fmt_real="Es16.8", io_log=io_log, font_color="blue")
                call html_log_real_number(                               &
                       log_text="Result from contr_sgto_carmom:",        &
                       log_real=contr_ints(isgto,icontr,jsgto,jcontr,1), &
                       fmt_real="Es16.8", io_log=io_log, font_color="red")
                test_failed = .true.
              end if
            end do
          end do
        end do
      end do
      ! cleans
      deallocate(contr_ints)
    end do
    return
  end subroutine test_contr_sgto_kinene
