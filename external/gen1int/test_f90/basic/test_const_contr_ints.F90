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
!!  Tests subroutine const_contr_ints.
!!
!!  2011-07-27, Bin Gao
!!  * first version

  !> \brief tests subroutine const_contr_ints
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_const_contr_ints(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
    integer, parameter :: NUM_CONTR_BRA = 20                    !number of contractions of bra center
    integer, parameter :: NUM_PRIM_BRA = 16                     !number of primitive Gaussians of bra center
    real(REALK) contr_coef_bra(NUM_CONTR_BRA,NUM_PRIM_BRA)      !contraction coefficients of bra center
    integer, parameter :: NUM_CONTR_KET = 20                    !number of contractions of ket center
    integer, parameter :: NUM_PRIM_KET = 16                     !number of primitive Gaussians of ket center
    real(REALK) contr_coef_ket(NUM_CONTR_KET,NUM_PRIM_KET)      !contraction coefficients of ket center
    integer, parameter :: NUM_CART_BRA = 15                     !number of Cartesian/Hermite GTOs on bra
    integer, parameter :: NUM_CART_KET = 15                     !number of Cartesian/Hermite GTOs on ket
    integer, parameter :: NUM_OPT = 6                           !number of operators including derivatives
    real(REALK) prim_ints(NUM_CART_BRA,NUM_CART_KET,NUM_OPT, &  !primitive Cartesian/Hermite integrals
                          NUM_PRIM_BRA,NUM_PRIM_KET)
    real(REALK) contr_ints(NUM_CART_BRA,NUM_CONTR_BRA, &        !contracted Cartesian/Hermite integrals
                           NUM_CART_KET,NUM_CONTR_KET,NUM_OPT)
    real(REALK) ref_cints(NUM_CART_BRA,NUM_CONTR_BRA, &         !referenced contracted integrals
                          NUM_CART_KET,NUM_CONTR_KET,NUM_OPT)
    integer iprim, jprim    !incremental recorders over primitives
    integer icontr, jcontr  !incremental recorders over contractions
    integer icart, jcart    !incremental recorders over Cartesian/Hermite Gaussians
    integer iopt            !incremental recorder over operators
    logical different       !if result from \fn(const_contr_ints) is different from reference
    real(REALK) begin_time  !begin of CPU time
    real(REALK) end_time    !end of CPU time
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_int_number("Number of contractions of bra center:", &
                             NUM_CONTR_BRA, "I4", io_log)
    call html_log_int_number("Number of primitive Gaussians of bra center:", &
                             NUM_PRIM_BRA, "I4", io_log)
    call html_log_int_number("Number of contractions of ket center:", &
                             NUM_CONTR_KET, "I4", io_log)
    call html_log_int_number("Number of primitive Gaussians of ket center:", &
                             NUM_PRIM_KET, "I4", io_log)
    call html_log_int_number("Number of Cartesian/Hermite GTOs on bra:", &
                             NUM_CART_BRA, "I4", io_log)
    call html_log_int_number("Number of Cartesian/Hermite GTOs on ket:", &
                             NUM_CART_KET, "I4", io_log)
    call html_log_int_number("Number of operators including derivatives:", &
                             NUM_OPT, "I4", io_log)
    ! assigns the contraction coefficients
    do iprim = 1, NUM_PRIM_BRA
      do icontr = 1, NUM_CONTR_BRA
        contr_coef_bra(icontr,iprim) = real(icontr,REALK)*0.2_REALK+real(iprim,REALK)
      end do
    end do
    do iprim = 1, NUM_PRIM_KET
      do icontr = 1, NUM_CONTR_KET
        contr_coef_ket(icontr,iprim) = real(icontr,REALK)*0.1_REALK+real(iprim,REALK)
      end do
    end do
    ! assigns the primitive Cartesian/Hermite integrals
    do jprim = 1, NUM_PRIM_KET
      do iprim = 1, NUM_PRIM_BRA
        do iopt = 1, NUM_OPT
          do jcart = 1, NUM_CART_KET
            do icart = 1, NUM_CART_BRA
              prim_ints(icart,jcart,iopt,iprim,jprim) &
                = real(icart,REALK)*0.0001_REALK      &
                + real(jcart,REALK)*0.001_REALK       &
                + real(iopt,REALK)*0.01_REALK         &
                + real(iprim,REALK)*0.1_REALK         &
                + real(jprim,REALK)
            end do
          end do
        end do
      end do
    end do
    ! gets the begin time
    call xtimer_set(begin_time)
    ! calls \fn(const_contr_ints) to get the const_contr_ints
    call const_contr_ints(NUM_CONTR_BRA, NUM_PRIM_BRA, contr_coef_bra, &
                          NUM_CONTR_KET, NUM_PRIM_KET, contr_coef_ket, &
                          NUM_CART_BRA, NUM_CART_KET, NUM_OPT,         &
                          prim_ints, contr_ints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for const_contr_ints:", &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! assigns the referenced contracted integrals
    ref_cints = 0.0_REALK
    do jprim = 1, NUM_PRIM_KET
      do iprim = 1, NUM_PRIM_BRA
        do iopt = 1, NUM_OPT
          do jcontr = 1, NUM_CONTR_KET
            do jcart = 1, NUM_CART_KET
              do icontr = 1, NUM_CONTR_BRA
                do icart = 1, NUM_CART_BRA
                  ref_cints(icart,icontr,jcart,jcontr,iopt)     &
                    = ref_cints(icart,icontr,jcart,jcontr,iopt) &
                    + contr_coef_bra(icontr,iprim)              &
                    * contr_coef_ket(jcontr,jprim)              &
                    * prim_ints(icart,jcart,iopt,iprim,jprim)
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
    call html_log_real_number(log_text="Time (s) used for explicit loops:",   &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! checks the results using five explicit do loops
    do iopt = 1, NUM_OPT
      do jcontr = 1, NUM_CONTR_KET
        do jcart = 1, NUM_CART_KET
          do icontr = 1, NUM_CONTR_BRA
            do icart = 1, NUM_CART_BRA
              call check_difference(ref_cints(icart,icontr,jcart,jcontr,iopt),  &
                                    contr_ints(icart,icontr,jcart,jcontr,iopt), &
                                    different)
              if (different) then
                call html_log_int_array(log_text="ID of contraction and operator:", &
                                        log_int=(/icart,icontr,jcart,jcontr,iopt/), &
                                        fmt_int="I3", io_log=io_log)
                call html_log_real_number(log_text="Referenced contracted integral:", &
                       log_real=ref_cints(icart,icontr,jcart,jcontr,iopt),            &
                       fmt_real="Es16.8", io_log=io_log, font_color="blue")
                call html_log_real_number(                                  &
                       log_text="Result from const_contr_ints:",            &
                       log_real=contr_ints(icart,icontr,jcart,jcontr,iopt), &
                       fmt_real="Es16.8", io_log=io_log, font_color="red")
                test_failed = .true.
              end if
            end do
          end do
        end do
      end do
    end do
    return
  end subroutine test_const_contr_ints
