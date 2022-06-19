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
!!  This subroutine checks the difference of two real numbers.
!!
!!  2012-02-07, Bin Gao
!!  * first version

  !> \brief checks the difference of two real numbers
  !> \author Bin Gao
  !> \date 2012-02-07
  !> \param ref_value is the referenced number
  !> \param chk_value is the number to check
  !> \return different indicates if these two numbers are different
  subroutine check_difference(ref_value, chk_value, different)
    use xkind
    implicit none
    real(REALK), intent(in) :: ref_value
    real(REALK), intent(in) :: chk_value
    logical, intent(out) :: different
!nearly negligible number relative to 1
#include "zero_thrsh.h"
! threshold of error
#include "err_thrsh.h"
    real(REALK) ratio_value  !ratio between referenced and checked numbers
    ! referenced number ~ 0.0
    if (abs(ref_value)<ZERO_THRSH) then
      ! checked number ~ 0.0
      if (abs(chk_value)<ZERO_THRSH) then
        different = .false.
      else
        ratio_value = ref_value/chk_value
        different = ratio_value<RATIO_THRSH(1) .or. ratio_value>RATIO_THRSH(2)
      end if
    else
      ratio_value = chk_value/ref_value
      different = ratio_value<RATIO_THRSH(1) .or. ratio_value>RATIO_THRSH(2)
    end if
    return
  end subroutine check_difference
