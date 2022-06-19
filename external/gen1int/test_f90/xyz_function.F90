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
!!  Returns artificial value for given orders of xyz components.
!!
!!  2011-07-27, Bin Gao
!!  * first version

  !> \brief returns artificial value for given orders of xyz components
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param order_xyz contains the orders of xyz components
  !> \return art_value is the artificial value
  function xyz_function(order_xyz) result(art_value)
    use xkind
    implicit none
    real(REALK) art_value
    integer, intent(in) :: order_xyz(3)
    art_value = real(order_xyz(1),REALK)*0.0006_REALK  &
              + real(order_xyz(2),REALK)*0.00002_REALK &
              + real(order_xyz(3),REALK)*0.00004_REALK
  end function xyz_function
