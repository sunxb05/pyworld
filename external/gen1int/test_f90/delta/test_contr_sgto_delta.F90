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
!!  Tests subroutines in contr_sgto_delta.F90
!!
!!  2012-03-19, Bin Gao
!!  * first version

  !> \brief tests subroutines in contr_sgto_delta.F90
  !> \author Bin Gao
  !> \date 2012-03-19
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_contr_sgto_delta(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! water molecule using cc-pV6Z basis sets
#include "test_f90/water_ccpv6z.h"
! test parameters and referenced results
#include "water_ccpv6z_sgto_delta.h"
    integer itst, jtst, ktst, ltst      !incremental recorders over test cases
    integer iblock, jblock              !incremental recorders over blocks of sub-shells
    integer num_opt                     !number of operators including derivatives
    real(REALK), allocatable :: contr_ints(:,:,:,:,:)
                                        !contracted integrals
    integer idx_bra, idx_ket            !indices of bra and ket centers
    integer num_sgto_bra, num_sgto_ket  !number of spherical Gaussians of bra and ket centers
    integer ierr                        !error information
    integer addr_ref                    !address of referenced results
    real(REALK) begin_time              !begin of CPU time
    real(REALK) end_time                !end of CPU time
    logical different                   !if result from \fn(contr_sgto_delta) is different from reference
    integer icontr, jcontr, isgto, jsgto, iopt
                                        !incremental recorders
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_text("Tests water molecule using cc-pV6Z spherical GTOs ...", io_log)
    ! loops over test cases
    jtst = 1
    addr_ref = 0
    do itst = 1, NUM_TESTS
      ktst = jtst+2
      call html_log_int_number("Order of electronic derivatives:", &
                               ORDER_ELEC(itst), "I4", io_log)
      if (IDX_DELORG(itst)>0) &
        call html_log_int_number("Index of Dirac delta function origin:", &
                                 IDX_DELORG(itst), "I4", io_log)
      call html_log_real_array("Coordinates of Dirac delta function origin:", &
                               DEL_ORIGIN(jtst:ktst), "F16.8", io_log)
      if (IDX_DIPORG(itst)>0)                               &
        call html_log_int_number("Index of dipole origin:", &
                                 IDX_DIPORG(itst), "I4", io_log)
      call html_log_real_array("Coordinates of dipole origin:", DIP_ORIGIN(jtst:ktst), &
                               "F16.8", io_log)
      call html_log_real_number("Scaling constant:", SCAL_CONST(itst), "F16.8", io_log)
      call html_log_int_number("Order of Cartesian multipole moments:", &
                               ORDER_MOM(itst), "I4", io_log)
      call html_log_int_number("Order of partial geometric derivatives on bra center:", &
                               ORDER_GEO_BRA(itst), "I4", io_log)
      call html_log_int_number("Order of partial geometric derivatives on ket center:", &
                               ORDER_GEO_KET(itst), "I4", io_log)
      call html_log_int_number(                                                &
             "Order of geometric derivatives on Dirac delta function origin:", &
             ORDER_GEO_POT(itst), "I4", io_log)
      call html_log_int_number("Order of geometric derivatives on dipole origin:", &
                               ORDER_GEO_MOM(itst), "I4", io_log)
      ! computes the number of total geometric derivatives
      num_opt = 1
      do ltst = jtst, jtst+NUM_CENTS(itst)-1
        num_opt = num_opt*(ORDER_CENT(ltst)+1)*(ORDER_CENT(ltst)+2)/2
      end do
      ltst = jtst+NUM_CENTS(itst)-1
      if (NUM_CENTS(itst)>0) then
        call html_log_int_number("Numbers of differentiated centers of total "// &
                                 "geometric derivatives:", NUM_CENTS(itst), "I4", io_log)
        call html_log_int_array("Indices of differentiated centers:", &
                                IDX_CENT(jtst:ltst), "I4", io_log)
        call html_log_int_array("Orders of differentiated centers:", &
                                ORDER_CENT(jtst:ltst), "I4", io_log)
        if (NUM_CENTS(itst)>3) then
          call html_log_int_number(log_text="Invalid number of differentiated "//      &
                                            "centers of total geometric derivatives:", &
                                   log_int=NUM_CENTS(itst), fmt_int="I4",              &
                                   io_log=io_log, font_color="red")
          test_failed = .true.
          cycle
        end if
      end if
      ! gets the begin time
      call xtimer_set(begin_time)
      ! computs the number of operators including the derivatives
      num_opt = num_opt                                         &
              * (ORDER_ELEC(itst)+1)*(ORDER_ELEC(itst)+2)       &
              * (ORDER_MOM(itst)+1)*(ORDER_MOM(itst)+2)         &
              * (ORDER_GEO_BRA(itst)+1)*(ORDER_GEO_BRA(itst)+2) &
              * (ORDER_GEO_KET(itst)+1)*(ORDER_GEO_KET(itst)+2) &
              * (ORDER_GEO_POT(itst)+1)*(ORDER_GEO_POT(itst)+2) &
              * (ORDER_GEO_MOM(itst)+1)*(ORDER_GEO_MOM(itst)+2)/64
      ! gets the blocks of sub-shells
      iblock = BRA_BLOCK(itst)
      jblock = KET_BLOCK(itst)
      call html_log_int_number("ID of block of sub-shells on bra center:", &
                               iblock, "I4", io_log)
      call html_log_int_number("ID of block of sub-shells on ket center:", &
                               jblock, "I4", io_log)
      ! gets the indices of bra and ket centers
      idx_bra = IDX_CENT_BLOCK(iblock)
      idx_ket = IDX_CENT_BLOCK(jblock)
      ! sets the number of spherical Gaussians of bra and ket centers
      num_sgto_bra = 2*ANGULAR_BLOCK(iblock)+1
      num_sgto_ket = 2*ANGULAR_BLOCK(jblock)+1
      ! allocates the contracted integrals
      allocate(contr_ints(num_sgto_bra,NUM_CONTR_BLOCK(iblock), &
                          num_sgto_ket,NUM_CONTR_BLOCK(jblock), &
                          num_opt), stat=ierr)
      if (ierr/=0) then
        call html_log_int_number(log_text="Failed to allocate contr_ints:",    &
                                 log_int=num_sgto_bra*NUM_CONTR_BLOCK(iblock)  &
                                         *num_sgto_ket*NUM_CONTR_BLOCK(jblock) &
                                         *num_opt,                             &
                                 fmt_int="I12", io_log=io_log, font_color="red")
        test_failed = .true.
        cycle
      end if
      ! computes the contracted integrals
      call contr_sgto_delta(idx_bra, COORDS(3*idx_bra-2:3*idx_bra),        &
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
                            ORDER_ELEC(itst), IDX_DELORG(itst),            &
                            DEL_ORIGIN(jtst:ktst), IDX_DIPORG(itst),       &
                            DIP_ORIGIN(jtst:ktst), SCAL_CONST(itst),       &
                            ORDER_MOM(itst), ORDER_GEO_BRA(itst),          &
                            ORDER_GEO_KET(itst), ORDER_GEO_POT(itst),      &
                            ORDER_GEO_MOM(itst), NUM_CENTS(itst),          &
                            IDX_CENT(jtst:ltst), ORDER_CENT(jtst:ltst),    &
                            num_sgto_bra, num_sgto_ket, num_opt, contr_ints)
      ! gets the end time
      call xtimer_set(end_time)
      ! prints the CPU elapsed time
      call html_log_real_number(log_text="Time (s) used for contr_sgto_delta:", &
                                log_real=end_time-begin_time, fmt_real="F10.4", &
                                io_log=io_log, font_color="blue")
      ! checks the results
      do iopt = 1, num_opt
        do jcontr = 1, NUM_CONTR_BLOCK(jblock)
          do jsgto = 1, num_sgto_ket
            do icontr = 1, NUM_CONTR_BLOCK(iblock)
              do isgto = 1, num_sgto_bra
                addr_ref = addr_ref+1
                call check_difference(REF_CONTR_INTS(addr_ref),                   &
                                      contr_ints(isgto,icontr,jsgto,jcontr,iopt), &
                                      different)
                if (different) then
                  call html_log_int_number("ID of test:", itst, "I3", io_log)
                  call html_log_int_array(log_text="ID of contraction and operator:", &
                                          log_int=(/isgto,icontr,jsgto,jcontr,iopt/), &
                                          fmt_int="I3", io_log=io_log)
                  call html_log_real_number(                         &
                         log_text="Referenced contracted integral:", &
                         log_real=REF_CONTR_INTS(addr_ref),          &
                         fmt_real="Es16.8", io_log=io_log, font_color="blue")
                  call html_log_real_number(                                  &
                         log_text="Result from contr_sgto_delta:",            &
                         log_real=contr_ints(isgto,icontr,jsgto,jcontr,iopt), &
                         fmt_real="Es16.8", io_log=io_log, font_color="red")
                  test_failed = .true.
                end if
              end do
            end do
          end do
        end do
      end do
      ! cleans
      deallocate(contr_ints)
      ! updates the incremental recorder
      jtst = jtst+3
    end do
    return
  end subroutine test_contr_sgto_delta
