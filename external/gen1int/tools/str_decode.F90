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
!!  This file converts the input string into integers.
!!
!!  2012-03-09, Bin Gao:
!!  * first version

#include "stdout.h"

!> \brief module of strind decoding
!> \author Bin Gao
!> \date 2012-03-09
module str_decode

  use xkind
  implicit none

  ! information of decoding string
  type, public :: decode_str_t
    private
    ! length of the string +1
    integer :: len_str = 0
    ! decoding characters
    character, allocatable :: decode_char(:)
    ! connection character
    character :: conn_char = "-"
    ! separation character
    character :: sep_char = ","
  end type decode_str_t

  public :: StrDecodeCreate
  public :: StrDecodeGetInts
  public :: StrDecodeDestroy

  contains

  !> \brief finds the number of integers from a string and intializes the information
  !>        of the decoding string
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param the_str is the string to decode
  !> \param conn_char is the connection character
  !> \param sep_char is the separation character
  !> \return num_ints is the number of found integers
  !> \return decode_str contains the information of string to decode
  subroutine StrDecodeCreate(the_str, conn_char, sep_char, num_ints, decode_str)
    character*(*), intent(in) :: the_str
    character, intent(in) :: conn_char
    character, intent(in) :: sep_char
    integer, intent(out) :: num_ints
    type(decode_str_t), intent(inout) :: decode_str
! digit from 0 to 9
#include "digit.h"
    integer sep_range(2)   !minimum and maximum integers between separations
    integer id_update      !pointer to updating \var(sep_range)
    integer num_sep_ints   !number of converted integers in between separations
    logical illegal_char   !if an illegal character
    integer istr           !incremental recorder over the string
    integer idx            !incremental recorder over digits
    integer ierr           !error information
#ifdef XTIME
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! gets the length of string
    decode_str%len_str = len_trim(the_str)+1
    ! sets the decoding characters
    allocate(decode_str%decode_char(decode_str%len_str), stat=ierr)
    if (ierr/=0)                                        &
      call error_stop("StrDecodeCreate",                &
                      "failed to allocate decode_char", &
                      decode_str%len_str)
    do istr = 1, decode_str%len_str-1
      decode_str%decode_char(istr) = the_str(istr:istr)
    end do
    decode_str%decode_char(decode_str%len_str) = sep_char
    ! initializes the number of converted integers
    num_ints = 0
    ! initializes the minimum and maximum integers between separations, and the pointer
    sep_range = 0
    id_update = 1
    ! initializes the recorder looping over the string
    istr = 1
    ! loops over the string
    do while (istr<=decode_str%len_str)
      ! connection character, we get the first integer, we then
      ! (1) change to sep_range(2) and,
      ! (2) increase the recorder
      if (decode_str%decode_char(istr)==conn_char) then
        id_update = 2
        istr = istr+1
      ! separation character or we meet the end of the string, we get the second figure, and
      ! (1) update the number of converted integers
      ! (2) change to sep_range(1)
      ! (3) clear the sep_range
      ! (4) increase the recorder
      else if (decode_str%decode_char(istr)==sep_char .or. &
               istr==decode_str%len_str) then
        ! number of converted integers before current separation
        if (sep_range(2)==0) sep_range(2)=sep_range(1)  !for the case of ...,5,...
        num_sep_ints = sep_range(2)-sep_range(1)+1
        num_ints = num_ints+num_sep_ints
        id_update = 1
        sep_range = 0
        istr = istr+1
      ! digit or illegal character, we then
      ! (1) update the sep_range
      ! (2) increase the recorder
      else
        illegal_char = .true.
        ! checks if it is a digit
        do idx = 0, 9
          if (DIGIT(idx)==decode_str%decode_char(istr)) then
            sep_range(id_update) = 10*sep_range(id_update)+idx
            illegal_char = .false.
            exit
          end if
        end do
        ! illegal character
        if (illegal_char) then
          call error_stop("StrDecodeCreate",                                  &
                          "illegal character "//decode_str%decode_char(istr), &
                          istr)
          ! or we just skip the illegal character
          !-istr = istr+1
          !-cycle
        end if
        istr = istr+1
      end if
    end do
    ! stores the connection and separation characters
    decode_str%conn_char = conn_char
    decode_str%sep_char = sep_char
#ifdef XTIME
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "StrDecodeCreate", STDOUT)
#endif
  end subroutine StrDecodeCreate

  !> \brief converts a string into integers
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param decode_str contains the information of the string to convert
  !> \param num_ints is the number of converted integers
  !> \return convert_ints contains the converted integers
  subroutine StrDecodeGetInts(decode_str, num_ints, convert_ints)
    type(decode_str_t), intent(in) :: decode_str
    integer, intent(in) :: num_ints
    integer, intent(out) :: convert_ints(num_ints)
! digit from 0 to 9
#include "digit.h"
    integer num_cnvrt_ints  !current number of converted integers
    integer sep_range(2)    !minimum and maximum integers between separations
    integer id_update       !pointer to updating \var(sep_range)
    integer num_sep_ints    !number of converted integers in between separations
    logical illegal_char    !if an illegal character
    integer istr            !incremental recorder over the string
    integer idx             !incremental recorder over digits
#ifdef XTIME
    real(REALK) curr_time   !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! initializes current number of converted integers
    num_cnvrt_ints = 0
    ! initializes the minimum and maximum integers between separations, and the pointer
    sep_range = 0
    id_update = 1
    ! initializes the recorder looping over the string
    istr = 1
    ! loops over the string
    do while (istr<=decode_str%len_str)
      ! connection character, we get the first integer, we then
      ! (1) change to sep_range(2) and,
      ! (2) increase the recorder
      if (decode_str%decode_char(istr)==decode_str%conn_char) then
        id_update = 2
        istr = istr+1
      ! separation character or we meet the end of the string, we get the second figure, and
      ! (1) gest the converted integers
      ! (2) change to sep_range(1)
      ! (3) clear the sep_range
      ! (4) increase the recorder
      else if (decode_str%decode_char(istr)==decode_str%sep_char .or. &
               istr==decode_str%len_str) then
        ! number of converted integers before current separation
        if (sep_range(2)==0) sep_range(2)=sep_range(1)  !for the case of ...,5,...
        num_sep_ints = sep_range(2)-sep_range(1)+1
        ! checks the size of converted integers
        if (num_cnvrt_ints+num_sep_ints>num_ints)                     &
          call error_stop("StrDecodeGetInts",                         &
                          "size of converted integers is not enough", &
                          num_cnvrt_ints+num_sep_ints)
        ! gets the converted integers before current separation
        do idx = 0, num_sep_ints-1
          num_cnvrt_ints = num_cnvrt_ints + 1
          convert_ints(num_cnvrt_ints) = sep_range(1)+idx
        end do
        id_update = 1
        sep_range = 0
        istr = istr+1
      ! digit or illegal character, we then
      ! (1) update the sep_range
      ! (2) increase the recorder
      else
        illegal_char = .true.
        ! checks if it is a digit
        do idx = 0, 9
          if (DIGIT(idx)==decode_str%decode_char(istr)) then
            sep_range(id_update) = 10*sep_range(id_update)+idx
            illegal_char = .false.
            exit
          end if
        end do
        ! illegal character
        if (illegal_char) then
          call error_stop("StrDecodeGetInts",                                 &
                          "illegal character "//decode_str%decode_char(istr), &
                          istr)
          ! or we just skip the illegal character
          !-istr = istr+1
          !-cycle
        end if
        istr = istr+1
      end if
    end do
#ifdef XTIME
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "StrDecodeGetInts", STDOUT)
#endif
  end subroutine StrDecodeGetInts

  !> \brief frees the space by decoding string
  !> \author Bin Gao
  !> \date 2012-03-09
  !> \param decode_str contains the information of decoding string
  subroutine StrDecodeDestroy(decode_str)
    type(decode_str_t), intent(inout) :: decode_str
#ifdef XTIME
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    decode_str%len_str = 0
    if (allocated(decode_str%decode_char)) deallocate(decode_str%decode_char)
    decode_str%conn_char = "-"
    decode_str%sep_char = ","
#ifdef XTIME
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "StrDecodeDestroy", STDOUT)
#endif
  end subroutine StrDecodeDestroy

end module str_decode
