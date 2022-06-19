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
!!  This file computes the scaled modified spherical Bessel function of the first kind.
!!
!!  2012-03-20, Bin Gao:
!!  * first version

#include "xkind.h"
#include "stdout.h"

         (c) dauxfun_msphi: scaled modified spherical Bessel function of the first kind,
                            used for effective core potential

  !> \brief computes the scaled modified spherical Bessel function of the first kind
  !> \author Bin Gao
  !> \date 2012-03-20
  !> \param 
  !> \param 
  !> \param 
  !> \return 
  subroutine aux_msphi_vec()
    implicit none

#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif

#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "aux_msphi_vec", STDOUT)
#endif
  end subroutine aux_msphi_vec
