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
!!  Tests subroutines in binom_coeff.F90.
!!
!!  2011-07-27, Bin Gao
!!  * first version

  !> \brief tests subroutines in binom_coeff.F90
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_binom_coeff(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
    integer, parameter :: MAX_ORDER_BINOM = 8                      !maximum order of the binomial coefficients
    integer ref_pascal(0:MAX_ORDER_BINOM,0:MAX_ORDER_BINOM)        !referenced Pascal's triangle
    real(REALK) binom_pascal(0:MAX_ORDER_BINOM,0:MAX_ORDER_BINOM)  !Pascal's triangle from \fn(pascal_triangle)
    integer binom_coef                                             !binomial coefficient from \fn(dbinom_coeff)
    integer num_n, num_k                                           !incremental recorders
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_int_number("Maximum order of the binomial coefficients:", &
                             MAX_ORDER_BINOM, "I4", io_log)
    ! assigns the referenced Pascal's triangle
    ref_pascal(:,0) = (/1,  0,  0,  0,  0,  0,  0,  0,  0/)
    ref_pascal(:,1) = (/1,  1,  0,  0,  0,  0,  0,  0,  0/)
    ref_pascal(:,2) = (/1,  2,  1,  0,  0,  0,  0,  0,  0/)
    ref_pascal(:,3) = (/1,  3,  3,  1,  0,  0,  0,  0,  0/)
    ref_pascal(:,4) = (/1,  4,  6,  4,  1,  0,  0,  0,  0/)
    ref_pascal(:,5) = (/1,  5, 10, 10,  5,  1,  0,  0,  0/)
    ref_pascal(:,6) = (/1,  6, 15, 20, 15,  6,  1,  0,  0/)
    ref_pascal(:,7) = (/1,  7, 21, 35, 35, 21,  7,  1,  0/)
    ref_pascal(:,8) = (/1,  8, 28, 56, 70, 56, 28,  8,  1/)
    ! calculates the Pascal's triangle from \fn(pascal_triangle)
    call pascal_triangle(MAX_ORDER_BINOM, binom_pascal)
    do num_n = 0, MAX_ORDER_BINOM
      do num_k = 0, num_n
        ! tests \fn(dbinom_coeff)
        call dbinom_coeff(num_n, num_k, binom_coef)
        if (binom_coef/=ref_pascal(num_k,num_n)) then
          call html_log_int_array(log_text="Numbers n and k:", &
                                  log_int=(/num_n,num_k/),     &
                                  fmt_int="I3", io_log=io_log)
          call html_log_int_number(log_text="Referenced result:",   &
                                   log_int=ref_pascal(num_k,num_n), &
                                   fmt_int="I3", io_log=io_log,     &
                                   font_color="blue")
          call html_log_int_number(log_text="Result from dbinom_coeff:", &
                                   log_int=binom_coef, fmt_int="I3",     &
                                   io_log=io_log, font_color="red")
          test_failed = .true.
        end if
        ! tests \fn(pascal_triangle)
        if (int(binom_pascal(num_k,num_n))/=ref_pascal(num_k,num_n)) then
          call html_log_int_array(log_text="Numbers n and k:", &
                                  log_int=(/num_n,num_k/),     &
                                  fmt_int="I3", io_log=io_log)
          call html_log_int_number(log_text="Referenced result:",   &
                                   log_int=ref_pascal(num_k,num_n), &
                                   fmt_int="I3", io_log=io_log,     &
                                   font_color="blue")
          call html_log_int_number(log_text="Result from pascal_triangle:", &
                                   log_int=int(binom_pascal(num_k,num_n)),  &
                                   fmt_int="I3", io_log=io_log, font_color="red")
          test_failed = .true.
        end if
      end do
    end do
    return
  end subroutine test_binom_coeff
