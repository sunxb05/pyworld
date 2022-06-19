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
!!  Tests subroutines in reorder_ints.F90.
!!
!!  2011-08-04, Bin Gao
!!  * first version

  !> \brief tests subroutines in reorder_ints.F90
  !> \author Bin Gao
  !> \date 2011-08-04
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_reorder_ints(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
    integer, parameter :: ANG_BRA = 10                   !angular number on bra center
    integer, parameter :: ANG_KET = 12                   !angular number on ket center
    integer, parameter :: NUM_SGTO_BRA = 2*ANG_BRA+1     !number of SGTOs on bra center
    integer, parameter :: NUM_SGTO_KET = 2*ANG_KET+1     !number of SGTOs on ket center
    integer mag_bra(NUM_SGTO_BRA)                        !magnetic numbers on bra center
    integer mag_ket(NUM_SGTO_KET)                        !magnetic numbers on ket center
    integer, parameter :: NUM_CONTR_BRA = 10             !number of contractions on bra center
    integer, parameter :: NUM_CONTR_KET = 12             !number of contractions on ket center
    integer, parameter :: NUM_OPT = 6                    !number of operators
    real(REALK) sgto_ints(NUM_SGTO_BRA,NUM_CONTR_BRA, &  !contracted SGTO integrals from Gen1Int
                          NUM_SGTO_KET,NUM_CONTR_KET,NUM_OPT)
    real(REALK) ro_sints(NUM_SGTO_BRA,NUM_CONTR_BRA, &   !reordered SGTO integrals
                         NUM_SGTO_KET,NUM_CONTR_KET,NUM_OPT)
    integer, parameter :: NUM_CGTO_BRA = (ANG_BRA+1)*(ANG_BRA+2)/2  !number of CGTOs on bra center
    integer, parameter :: NUM_CGTO_KET = (ANG_KET+1)*(ANG_KET+2)/2  !number of CGTOs on ket center
    integer power_bra(3,NUM_CGTO_BRA)                               !Cartesian powers on bra center
    integer power_ket(3,NUM_CGTO_KET)                               !Cartesian powers on ket center
    real(REALK) cgto_ints(NUM_CGTO_BRA,NUM_CONTR_BRA, &  !contracted CGTO integrals from Gen1Int
                          NUM_CGTO_KET,NUM_CONTR_KET,NUM_OPT)
    real(REALK) ro_cints(NUM_CGTO_BRA,NUM_CONTR_BRA, &   !reordered CGTO integrals
                         NUM_CGTO_KET,NUM_CONTR_KET,NUM_OPT)
    integer ibra, iket, icontr, jcontr, iopt             !incremental recorders
    integer ypow, zpow                                   !incremental recorders over yz powers
    integer addr_bra, addr_ket                           !addresses of GTOs on bra and ket centers
    real(REALK) begin_time                               !begin of CPU time
    real(REALK) end_time                                 !end of CPU time
    logical different                                    !if result is different from reference
    ! assumes the test will pass
    test_failed = .false.
    ! sets SGTOs on bra and centers, which are in the reverse order of Gen1Int
    do ibra = 1, NUM_SGTO_BRA
      mag_bra(ibra) = ANG_BRA-ibra+1
    end do
    do iket = 1, NUM_SGTO_KET
      mag_ket(iket) = ANG_KET-iket+1
    end do
    ! sets the integrals
    do iopt = 1, NUM_OPT
      do jcontr = 1, NUM_CONTR_KET
        do iket = 1, NUM_SGTO_KET
          do icontr = 1, NUM_CONTR_BRA
            do ibra = 1, NUM_SGTO_BRA
              sgto_ints(ibra,icontr,iket,jcontr,iopt)                          &
                = real(ibra,REALK)*0.1_REALK+real(icontr,REALK)*0.01_REALK     &
                + real(iket,REALK)*0.001_REALK+real(jcontr,REALK)*0.0001_REALK &
                + real(iopt,REALK)*0.00001_REALK
            end do
          end do
        end do
      end do
    end do
    ! dumps information of tests
    call html_log_heading(heading="Tests subroutine reorder_sgtos", &
                          io_log=io_log, level=4)
    ! gets the begin time
    call xtimer_set(begin_time)
    ! reorders the SGTOs on bra center using \fn(reorder_sgtos)
    call reorder_sgtos(ANG_BRA, NUM_SGTO_BRA, mag_bra, 1, NUM_CONTR_BRA, &
                       NUM_SGTO_KET*NUM_CONTR_KET*NUM_OPT, sgto_ints, ro_sints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for reorder_sgtos:",    &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! checks the results
    do iopt = 1, NUM_OPT
      do jcontr = 1, NUM_CONTR_KET
        do iket = 1, NUM_SGTO_KET
          do icontr = 1, NUM_CONTR_BRA
            do ibra = 1, NUM_SGTO_BRA
              addr_bra = NUM_SGTO_BRA-ibra+1
              call check_difference(sgto_ints(ibra,icontr,iket,jcontr,iopt),    &
                                    ro_sints(addr_bra,icontr,iket,jcontr,iopt), &
                                    different)
              if (different) then
                call html_log_int_array(log_text="ID of contraction and operator:",   &
                                        log_int=(/addr_bra,icontr,iket,jcontr,iopt/), &
                                        fmt_int="I3", io_log=io_log)
                call html_log_real_number(log_text="Referenced result:", &
                       log_real=sgto_ints(ibra,icontr,iket,jcontr,iopt), &
                       fmt_real="Es20.12", io_log=io_log, font_color="blue")
                call html_log_real_number(log_text="Result from reorder_sgtos:", &
                       log_real=ro_sints(addr_bra,icontr,iket,jcontr,iopt),      &
                       fmt_real="Es20.12", io_log=io_log, font_color="red")
                test_failed = .true.
              end if
            end do
          end do
        end do
      end do
    end do
    ! dumps information of tests
    call html_log_heading(heading="Tests subroutine reorder_sgto_ints", &
                          io_log=io_log, level=4)
    ! resets the begin time
    call xtimer_set(begin_time)
    ! reorders the SGTOs on bra and ket centers using \fn(reorder_sgto_ints)
    call reorder_sgto_ints(ANG_BRA, NUM_SGTO_BRA, mag_bra, &
                           ANG_KET, NUM_SGTO_KET, mag_ket, &
                           NUM_CONTR_BRA, NUM_CONTR_KET,   &
                           NUM_OPT, sgto_ints, ro_sints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for reorder_sgto_ints:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",  &
                              io_log=io_log, font_color="blue")
    ! checks the results
    do iopt = 1, NUM_OPT
      do jcontr = 1, NUM_CONTR_KET
        do iket = 1, NUM_SGTO_KET
          addr_ket = NUM_SGTO_KET-iket+1
          do icontr = 1, NUM_CONTR_BRA
            do ibra = 1, NUM_SGTO_BRA
              addr_bra = NUM_SGTO_BRA-ibra+1
              call check_difference(sgto_ints(ibra,icontr,iket,jcontr,iopt),        &
                                    ro_sints(addr_bra,icontr,addr_ket,jcontr,iopt), &
                                    different)
              if (different) then
                call html_log_int_array(log_text="ID of contraction and operator:",    &
                       log_int=(/addr_bra,icontr,addr_ket,jcontr,iopt/), fmt_int="I3", &
                       io_log=io_log)
                call html_log_real_number(log_text="Referenced result:", &
                       log_real=sgto_ints(ibra,icontr,iket,jcontr,iopt), &
                       fmt_real="Es20.12", io_log=io_log, font_color="blue")
                call html_log_real_number(log_text="Result from reorder_sgto_ints:", &
                       log_real=ro_sints(addr_bra,icontr,addr_ket,jcontr,iopt),      &
                       fmt_real="Es20.12", io_log=io_log, font_color="red")
                test_failed = .true.
              end if
            end do
          end do
        end do
      end do
    end do
    ! sets CGTOs on bra and ket centers, which are in the reverse order of Gen1Int
    addr_bra = 0
    do zpow = ANG_BRA, 0, -1
      do ypow = ANG_BRA-zpow, 0, -1
        addr_bra = addr_bra+1
        power_bra(1,addr_bra) = ANG_BRA-(zpow+ypow)
        power_bra(2,addr_bra) = ypow
        power_bra(3,addr_bra) = zpow
      end do
    end do
    addr_ket = 0
    do zpow = ANG_KET, 0, -1
      do ypow = ANG_KET-zpow, 0, -1
        addr_ket = addr_ket+1
        power_ket(1,addr_ket) = ANG_KET-(zpow+ypow)
        power_ket(2,addr_ket) = ypow
        power_ket(3,addr_ket) = zpow
      end do
    end do
    ! sets the integrals
    do iopt = 1, NUM_OPT
      do jcontr = 1, NUM_CONTR_KET
        do iket = 1, NUM_CGTO_KET
          do icontr = 1, NUM_CONTR_BRA
            do ibra = 1, NUM_CGTO_BRA
              cgto_ints(ibra,icontr,iket,jcontr,iopt)                          &
                = real(ibra,REALK)*0.1_REALK+real(icontr,REALK)*0.01_REALK     &
                + real(iket,REALK)*0.001_REALK+real(jcontr,REALK)*0.0001_REALK &
                + real(iopt,REALK)*0.00001_REALK
            end do
          end do
        end do
      end do
    end do
    ! dumps information of tests
    call html_log_heading(heading="Tests subroutine reorder_cgtos", &
                          io_log=io_log, level=4)
    ! resets the begin time
    call xtimer_set(begin_time)
    ! reorders the CGTOs on bra center using \fn(reorder_cgtos)
    call reorder_cgtos(ANG_BRA, NUM_CGTO_BRA, power_bra, 1, NUM_CONTR_BRA, &
                       NUM_CGTO_KET*NUM_CONTR_KET*NUM_OPT, cgto_ints, ro_cints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for reorder_cgtos:",    &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! checks the results
    do iopt = 1, NUM_OPT
      do jcontr = 1, NUM_CONTR_KET
        do iket = 1, NUM_CGTO_KET
          do icontr = 1, NUM_CONTR_BRA
            do ibra = 1, NUM_CGTO_BRA
              addr_bra = NUM_CGTO_BRA+1-ibra
              call check_difference(cgto_ints(ibra,icontr,iket,jcontr,iopt),    &
                                    ro_cints(addr_bra,icontr,iket,jcontr,iopt), &
                                    different)
              if (different) then
                call html_log_int_array(log_text="ID of contraction and operator:",   &
                                        log_int=(/addr_bra,icontr,iket,jcontr,iopt/), &
                                        fmt_int="I3", io_log=io_log)
                call html_log_real_number(log_text="Referenced result:", &
                       log_real=cgto_ints(ibra,icontr,iket,jcontr,iopt), &
                       fmt_real="Es20.12", io_log=io_log, font_color="blue")
                call html_log_real_number(log_text="Result from reorder_cgtos:", &
                       log_real=ro_cints(addr_bra,icontr,iket,jcontr,iopt),      &
                       fmt_real="Es20.12", io_log=io_log, font_color="red")
                test_failed = .true.
              end if
            end do
          end do
        end do
      end do
    end do
    ! dumps information of tests
    call html_log_heading(heading="Tests subroutine reorder_cgto_ints", &
                          io_log=io_log, level=4)
    ! resets the begin time
    call xtimer_set(begin_time)
    ! reorders the CGTOs on bra and ket centers using \fn(reorder_cgto_ints)
    call reorder_cgto_ints(ANG_BRA, NUM_CGTO_BRA, power_bra, &
                           ANG_KET, NUM_CGTO_KET, power_ket, &
                           NUM_CONTR_BRA, NUM_CONTR_KET,     &
                           NUM_OPT, cgto_ints, ro_cints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for reorder_cgto_ints:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",  &
                              io_log=io_log, font_color="blue")
    ! checks the results
    do iopt = 1, NUM_OPT
      do jcontr = 1, NUM_CONTR_KET
        do iket = 1, NUM_CGTO_KET
          addr_ket = NUM_CGTO_KET+1-iket
          do icontr = 1, NUM_CONTR_BRA
            do ibra = 1, NUM_CGTO_BRA
              addr_bra = NUM_CGTO_BRA+1-ibra
              call check_difference(cgto_ints(ibra,icontr,iket,jcontr,iopt),        &
                                    ro_cints(addr_bra,icontr,addr_ket,jcontr,iopt), &
                                    different)
              if (different) then
                call html_log_int_array(log_text="ID of contraction and operator:",    &
                       log_int=(/addr_bra,icontr,addr_ket,jcontr,iopt/), fmt_int="I3", &
                       io_log=io_log)
                call html_log_real_number(log_text="Referenced result:", &
                       log_real=cgto_ints(ibra,icontr,iket,jcontr,iopt), &
                       fmt_real="Es20.12", io_log=io_log, font_color="blue")
                call html_log_real_number(log_text="Result from reorder_cgto_ints:", &
                       log_real=ro_cints(addr_bra,icontr,addr_ket,jcontr,iopt),      &
                       fmt_real="Es20.12", io_log=io_log, font_color="red")
                test_failed = .true.
              end if
            end do
          end do
        end do
      end do
    end do
    return
  end subroutine test_reorder_ints
