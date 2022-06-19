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
!!  This subroutine converts an integer to string.
!!
!!  2009-10-25, Bin Gao:
!!  * first version

  !> \brief converts an integer to string, with prefix zeros if required
  !> \author Bin Gao
  !> \date 2009-10-25
  !> \param the_int is the integer
  !> \param len_str is the length of string
  !> \param left_align indicates if the string is left align
  !> \param zero_prefix indicats if the string is with prefix zeros when right align
  !> \return the_str is the converted string, with prefix zeros if required
  subroutine int_to_str(the_int, len_str, left_align, zero_prefix, the_str)
    use xkind
    implicit none
    integer, intent(in) :: the_int
    integer, intent(in) :: len_str
    logical, intent(in) :: left_align
    logical, intent(in) :: zero_prefix
    character(len_str), intent(out) :: the_str
! digit from 0 to 9
#include "digit.h"
    integer abs_int    !absolute value of the integer
    integer remainder  !remainder
    integer ic         !incremental recorder over characters
    ! initializes the string
    the_str = ""
    ! assigns the absolute value of the integer at first
    abs_int = abs(the_int)
    ! returns string with prefix zeros
    if (zero_prefix) then
      ! adds "-" for negative integer
      if (the_int<0) the_str(1:1) = "-"
      ! loops over the characters in the string
      do ic = len_str, 2, -1
        ! remainder
        remainder = mod(abs_int,10)
        ! converts the remainder to a character
        the_str(ic:ic) = DIGIT(remainder)
        ! updates the absolute value of the integer
        abs_int = (abs_int-remainder)/10
      end do
      ! checks if the length of string is enough
      if (abs_int>0)                  &
        call error_stop("int_to_str", &
                        "failed to convert the integer ("//the_str//")", the_int)
    ! returns string without prefix zeros
    else
      ! loops over the characters in the string
      do ic = len_str, 2, -1
        ! remainder
        remainder = mod(abs_int,10)
        ! converts the remainder to a character
        the_str(ic:ic) = DIGIT(remainder)
        ! updates the absolute value of the integer
        abs_int = (abs_int-remainder)/10
        if (abs_int<=0) then
          ! adds "-" for negative integer
          if (the_int<0) the_str(ic-1:ic-1) = "-"
          exit
        end if
      end do
      ! checks if the length of string is enough
      if (abs_int>0)                  &
        call error_stop("int_to_str", &
                        "failed to convert the integer ("//the_str//")", the_int)
      if (left_align) the_str = adjustl(the_str)
    end if
    return
  end subroutine int_to_str
