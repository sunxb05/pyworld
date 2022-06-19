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
!!  Tests tests spherical multipole integrals from Fortran 90 modules using SGTOs.
!!
!!  2012-03-20, Bin Gao
!!  * first version

  !> \brief tests spherical multipole integrals from Fortran 90 modules using SGTOs
  !> \author Bin Gao
  !> \date 2012-03-20
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_f90mod_sgto_sphmom(io_log, test_failed)
    use xkind
    ! Fortran 90 module of Gen1Int library
    use gen1int
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! water molecule using cc-pV6Z basis sets
#include "test_f90/water_ccpv6z.h"
! test parameters and referenced results
#include "water_f90mod_sgto_sphmom.h"
    integer itst                        !incremental recorder over test cases
    integer iblock, jblock              !incremental recorders over blocks of sub-shells
    integer idx_bra, idx_ket            !indices of bra and ket centers
    integer num_sgto_bra, num_sgto_ket  !number of spherical Gaussians of bra and ket centers
    type(one_prop_t) one_prop           !one-electron property integral operator
    integer num_prop                    !number of property integral matrices
    real(REALK), allocatable :: contr_ints(:,:,:,:,:)
                                        !contracted integrals
    integer ierr                        !error information
    integer addr_ref                    !address of referenced results
    real(REALK) begin_time              !begin of CPU time
    real(REALK) end_time                !end of CPU time
    logical different                   !if result from \fn(contr_sgto_carmom) is different from reference
    integer iopt                        !incremental recorder over operators
    integer icontr, jcontr              !incremental recorders over contractions
    integer isgto, jsgto                !incremental recorders over GTOs
    ! assumes the test will pass
    test_failed = .false.
    call html_log_text("Tests spherical multipole integrals of water molecule "// &
                       "using cc-pV6Z spherical GTOs ...", io_log)
    ! initializes the information of spherical multipole integrals
    call OnePropCreate(prop_name=INT_SPHER_MULTIPOLE, &
                       one_prop=one_prop,             &
                       info_prop=ierr,                &
                       dipole_origin=DIP_ORIGIN,      &
                       order_mom=ORDER_MOM)
    ! gets the number of property integral matrices
    call OnePropGetNumProp(one_prop=one_prop, num_prop=num_prop)
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
                          num_prop), stat=ierr)
      if (ierr/=0) then
        call html_log_int_number(log_text="Failed to allocate contr_ints:",    &
                                 log_int=num_sgto_bra*NUM_CONTR_BLOCK(iblock)  &
                                         *num_sgto_ket*NUM_CONTR_BLOCK(jblock) &
                                         *num_prop,                            &
                                 fmt_int="I12", io_log=io_log, font_color="red")
        test_failed = .true.
        cycle
      end if
      ! computes the contracted integrals
      call OnePropGetIntegral(idx_bra=idx_bra, coord_bra=COORDS(3*idx_bra-2:3*idx_bra),  &
                              angular_bra=ANGULAR_BLOCK(iblock),                         &
                              num_prim_bra=NUM_PRIM_BLOCK(iblock),                       &
                              exponent_bra=EXPONENT_BLOCK(START_PRIM_BLOCK(iblock)+1     &
                                                          :START_PRIM_BLOCK(iblock+1)),  &
                              num_contr_bra=NUM_CONTR_BLOCK(iblock),                     &
                              contr_coef_bra=REF_NRM_SGTO(START_CONTR_BLOCK(iblock)+1    &
                                                          :START_CONTR_BLOCK(iblock+1)), &
                              idx_ket=idx_ket, coord_ket=COORDS(3*idx_ket-2:3*idx_ket),  &
                              angular_ket=ANGULAR_BLOCK(jblock),                         &
                              num_prim_ket=NUM_PRIM_BLOCK(jblock),                       &
                              exponent_ket=EXPONENT_BLOCK(START_PRIM_BLOCK(jblock)+1     &
                                                          :START_PRIM_BLOCK(jblock+1)),  &
                              num_contr_ket=NUM_CONTR_BLOCK(jblock),                     &
                              contr_coef_ket=REF_NRM_SGTO(START_CONTR_BLOCK(jblock)+1    &
                                                          :START_CONTR_BLOCK(jblock+1)), &
                              spher_gto=.true., one_prop=one_prop,                       &
                              num_gto_bra=num_sgto_bra, num_gto_ket=num_sgto_ket,        &
                              num_opt=num_prop, contr_ints=contr_ints)
      ! gets the end time
      call xtimer_set(end_time)
      ! prints the CPU elapsed time
      call html_log_real_number(log_text="Time (s) used for OnePropGetIntegral:", &
                                log_real=end_time-begin_time, fmt_real="F10.4",   &
                                io_log=io_log, font_color="blue")
      ! checks the results
      do iopt = 1, num_prop
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
                  call html_log_real_number(                                  &
                         log_text="Referenced spherical multipole integral:", &
                         log_real=REF_CONTR_INTS(addr_ref),                   &
                         fmt_real="Es16.8", io_log=io_log, font_color="blue")
                  call html_log_real_number(                                  &
                         log_text="Result from OnePropGetIntegral:",          &
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
    end do
    return
  end subroutine test_f90mod_sgto_sphmom
