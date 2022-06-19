!!  gen1int: compute derivatives of one-electron integrals using Hermite Gaussians
!!  Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen
!!
!!  This file is part of gen1int.
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
!!  This file contains subroutines related to partial geometric derivatives
!!  of integrals of two-center operators.
!!
!!  2011-07-02, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief sets the orders of partial geometric derivatives and the parameters
  !>        for operators with two centers
  !> \author Bin Gao
  !> \date 2011-07-02
  !> \param num_cents is the number of differentiated centers, could be 1, 2, 3 or 4
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the order of derivatives of the differentiated centers
  !> \param sum_max_order is the sum of two maximum orders of derivatives
  !> \param idx_part_cent contains the indices of operator, bra and ket centers; it will
  !>        be sorted in ascending order on exit
  !> \return order_part_cent contains the orders of partial geometric derivatives of
  !>         operator, bra and ket centers as the same sequence on entry
  !> \return zero_ints indicates if the integrals are zero
  !> \return neg_one indicates if the integrals will be multiplied by -1
  !> \return binom_exp indicates if binomial expansion of total geometric derivatives is performed
  !> \return num_binom is the number of terms in binomial expansion except for
  !>         the first and last terms
  !> \return binom_cent is the operator, bra and ket centers involved in binomial expansion
  !> \return binom_coef contains the binomial coefficients up to the order of binomial expansion,
  !>         will be used later for generating the total geometric derivatives from partial ones
  !> \return num_part_geo is the number of resulted partial geometric derivatives
  !> \return seq_part_geo contains the sequence of operator, bra and ket centers
  !>         for partial geometric derivatives
  subroutine geom_part_two_param(num_cents, idx_cent, order_cent, sum_max_order,     &
                                 idx_part_cent, order_part_cent, zero_ints, neg_one, &
                                 binom_exp, num_binom, binom_cent, binom_coef,       &
                                 num_part_geo, seq_part_geo)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: sum_max_order
    integer, intent(inout) :: idx_part_cent(4)
    integer, intent(inout) :: order_part_cent(4)
    logical, intent(out) :: zero_ints
    logical, intent(out) :: neg_one
    logical, intent(out) :: binom_exp
    integer, intent(out) :: num_binom
    integer, intent(out) :: binom_cent(2)
    real(REALK), intent(out) :: binom_coef(0:sum_max_order,0:sum_max_order)
    integer, intent(out) :: num_part_geo(3)
    integer, intent(out) :: seq_part_geo(4,3)
!f2py intent(hide) :: num_cents
!f2py intent(in) :: idx_cent
!f2py intent(in) :: order_cent
!f2py depend(num_cents) :: order_cent
!f2py intent(in) :: sum_max_order
!f2py intent(inout) :: idx_part_cent
!f2py intent(inout) :: order_part_cent
!f2py intent(out) :: zero_ints
!f2py intent(out) :: neg_one
!f2py intent(out) :: binom_exp
!f2py intent(out) :: num_binom
!f2py intent(out) :: binom_cent
!f2py intent(out) :: binom_coef
!f2py depend(sum_max_order) :: binom_coef
!f2py intent(out) :: num_part_geo
!f2py intent(out) :: seq_part_geo
    integer, parameter :: NUM_PART_CENTS = 4  !number of operator, bra and ket centers
    integer tag_part_cent(NUM_PART_CENTS)     !marking the indices before sorting
    integer num_non_cents                     !number of unique non-atomic centers
    integer id_non_cent(NUM_PART_CENTS)       !position of unique non-atomic centers
    integer num_non_iden(NUM_PART_CENTS)      !number of identical non-atomic centers
    integer num_atom_cents                    !number of unique atomic centers
    integer id_atom_cent(NUM_PART_CENTS)      !position of unique atomic centers
    integer num_atom_iden(NUM_PART_CENTS)     !number of identical atomic centers
#if defined(XTIME)
    real(REALK) curr_time                     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! sorts the centers, and finds the unique ones and corresponding number of identical ones
    tag_part_cent = (/1,2,3,4/)
    call sort_gen_cents(NUM_PART_CENTS, idx_part_cent, tag_part_cent, &
                        num_non_cents, id_non_cent, num_non_iden,     &
                        num_atom_cents, id_atom_cent, num_atom_iden)
    ! too many differentiated centers
    if (num_cents>num_atom_cents) then
      zero_ints = .true.
    else
      ! initializes
      order_part_cent = 0
      select case(num_non_cents)
      ! there does not exist any non-atomic center
      case(0)
        select case(num_atom_cents)
        ! zeros due to translational invariance, no need to check the differentiated centers
        case(1)
          zero_ints = .true.
        ! there are two unique atomic centers
        case(2)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              select case(num_atom_iden(1))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! binomial expansion
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              ! translation invariance by performing partial geometric derivatives
              ! on another unique atomic center
              case(3)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = mod(order_cent(1),2)==1
                binom_exp = .false.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(1) (0/2/1)!"
              end select
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              select case(num_atom_iden(2))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! binomial expansion
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(2)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              ! translation invariance by performing partial geometric derivatives
              ! on another unique atomic center
              case(3)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = mod(order_cent(1),2)==1
                binom_exp = .false.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(2) (0/2/1)!"
              end select
            else
              zero_ints = .true.
            end if
          ! two-center total geometric derivatives
          case(2)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
              select case(num_atom_iden(1))
              ! translational invairance by performing partial geometric derivatives
              ! on the first unique atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)+order_cent(2)
                zero_ints = .false.
                neg_one = mod(order_cent(2),2)==1
                binom_exp = .false.
              ! translational invairance by performing partial geometric derivatives
              ! on the first unique atomic center, but with binomial expansion
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)+order_cent(2)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_part_cent(seq_part_geo(1,2))-1  !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)                 !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_part_cent(seq_part_geo(1,2)),              &
                                     binom_coef(0:order_part_cent(seq_part_geo(1,2)), &
                                                0:order_part_cent(seq_part_geo(1,2))))
                zero_ints = .false.
                neg_one = mod(order_cent(2),2)==1
                binom_exp = .true.
              ! translational invairance by performing partial geometric derivatives
              ! on the second unique atomic center
              case(3)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)+order_cent(2)
                zero_ints = .false.
                neg_one = mod(order_cent(1),2)==1
                binom_exp = .false.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(1) (0/2/2)!"
              end select
            else
              zero_ints = .true.
            end if
          case default
            stop "geom_part_two_param>> invalid num_cents (0/2/X)!"
          end select
        ! there are three unique atomic centers, the only possible choice is there are
        ! two single unique atomic centers, and one double unique atomic center
        case(3)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              select case(num_atom_iden(1))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! binomial expansion
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(1) (0/3/1)!"
              end select
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              select case(num_atom_iden(2))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! binomial expansion
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(2)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(2) (0/3/1)!"
              end select
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(3))) then
              select case(num_atom_iden(3))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(3))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! binomial expansion
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(3)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(3)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(3) (0/3/1)!"
              end select
            else
              zero_ints = .true.
            end if
          ! two-center total geometric derivatives
          case(2)
            ! notice that both centers are sorted in ascending order, so that
            ! we divide them into two parts according to the third atomic center
            if (idx_cent(2)<idx_part_cent(id_atom_cent(3))) then
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                  idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
                ! binomial expansion on the first atomic center
                if (num_atom_iden(1)==2) then
                  num_part_geo(1) = 2  !the first and last monomial terms have only one
                  num_part_geo(2) = 3  !differentiated center, while others have two,
                  num_part_geo(3) = 2  !we then add another atomic center
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                  seq_part_geo(3,2) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                  seq_part_geo(2,3) = tag_part_cent(id_atom_cent(2))
                  ! this gives the first monomial term, we loops over others outside
                  order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,2)) = 0
                  order_part_cent(seq_part_geo(3,2)) = order_cent(2)
                  num_binom = order_cent(1)-1        !number of binomial terms
                  binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                  binom_cent(2) = seq_part_geo(2,2)
                  ! prepares binomial coefficients for later use
                  call pascal_triangle(order_cent(1), &
                                       binom_coef(0:order_cent(1),0:order_cent(1)))
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .true.
                ! binomial expansion on the second atomic center
                else if (num_atom_iden(2)==2) then
                  num_part_geo(1) = 2  !the first and last monomial terms have only one
                  num_part_geo(2) = 3  !differentiated center, while others have two,
                  num_part_geo(3) = 2  !we then add another atomic center
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(3,2) = tag_part_cent(id_atom_cent(2)+1)
                  seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,3) = tag_part_cent(id_atom_cent(2)+1)
                  order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                  ! this gives the first monomial term, we loops over others outside
                  order_part_cent(seq_part_geo(2,2)) = order_cent(2)
                  order_part_cent(seq_part_geo(3,2)) = 0
                  num_binom = order_cent(2)-1        !number of binomial terms
                  binom_cent(1) = seq_part_geo(2,2)  !centers of binomial terms
                  binom_cent(2) = seq_part_geo(3,2)
                  ! prepares binomial coefficients for later use
                  call pascal_triangle(order_cent(2), &
                                       binom_coef(0:order_cent(2),0:order_cent(2)))
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .true.
                ! only one partial geometric derivative term on the first and second atomic centers
                else
                  num_part_geo(1) = 2
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                  order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .false.
                end if
              else
                zero_ints = .true.
              end if
            else if (idx_cent(2)==idx_part_cent(id_atom_cent(3))) then
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
                ! binomial expansion on the first atomic center
                if (num_atom_iden(1)==2) then
                  num_part_geo(1) = 2  !the first and last monomial terms have only one
                  num_part_geo(2) = 3  !differentiated center, while others have two,
                  num_part_geo(3) = 2  !we then add another atomic center
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                  seq_part_geo(3,2) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                  seq_part_geo(2,3) = tag_part_cent(id_atom_cent(3))
                  ! this gives the first monomial term, we loops over others outside
                  order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,2)) = 0
                  order_part_cent(seq_part_geo(3,2)) = order_cent(2)
                  num_binom = order_cent(1)-1        !number of binomial terms
                  binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                  binom_cent(2) = seq_part_geo(2,2)
                  ! prepares binomial coefficients for later use
                  call pascal_triangle(order_cent(1), &
                                       binom_coef(0:order_cent(1),0:order_cent(1)))
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .true.
                ! binomial expansion on the third atomic center
                else if (num_atom_iden(3)==2) then
                  num_part_geo(1) = 2  !the first and last monomial terms have only one
                  num_part_geo(2) = 3  !differentiated center, while others have two,
                  num_part_geo(3) = 2  !we then add another atomic center
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,2) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(3,2) = tag_part_cent(id_atom_cent(3)+1)
                  seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,3) = tag_part_cent(id_atom_cent(3)+1)
                  order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                  ! this gives the first monomial term, we loops over others outside
                  order_part_cent(seq_part_geo(2,2)) = order_cent(2)
                  order_part_cent(seq_part_geo(3,2)) = 0
                  num_binom = order_cent(2)-1        !number of binomial terms
                  binom_cent(1) = seq_part_geo(2,2)  !centers of binomial terms
                  binom_cent(2) = seq_part_geo(3,2)
                  ! prepares binomial coefficients for later use
                  call pascal_triangle(order_cent(2), &
                                       binom_coef(0:order_cent(2),0:order_cent(2)))
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .true.
                ! only one partial geometric derivative term on the first and third atomic centers
                else
                  num_part_geo(1) = 2
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                  order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .false.
                end if
              else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
                ! binomial expansion on the second atomic center
                if (num_atom_iden(2)==2) then
                  num_part_geo(1) = 2  !the first and last monomial terms have only one
                  num_part_geo(2) = 3  !differentiated center, while others have two,
                  num_part_geo(3) = 2  !we then add another atomic center
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(1,2) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2)+1)
                  seq_part_geo(3,2) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(1,3) = tag_part_cent(id_atom_cent(2)+1)
                  seq_part_geo(2,3) = tag_part_cent(id_atom_cent(3))
                  ! this gives the first monomial term, we loops over others outside
                  order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,2)) = 0
                  order_part_cent(seq_part_geo(3,2)) = order_cent(2)
                  num_binom = order_cent(1)-1        !number of binomial terms
                  binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                  binom_cent(2) = seq_part_geo(2,2)
                  ! prepares binomial coefficients for later use
                  call pascal_triangle(order_cent(1), &
                                       binom_coef(0:order_cent(1),0:order_cent(1)))
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .true.
                ! binomial expansion on the third atomic center
                else if (num_atom_iden(3)==2) then
                  num_part_geo(1) = 2  !the first and last monomial terms have only one
                  num_part_geo(2) = 3  !differentiated center, while others have two,
                  num_part_geo(3) = 2  !we then add another atomic center
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(1,2) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(2,2) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(3,2) = tag_part_cent(id_atom_cent(3)+1)
                  seq_part_geo(1,3) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(2,3) = tag_part_cent(id_atom_cent(3)+1)
                  order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                  ! this gives the first monomial term, we loops over others outside
                  order_part_cent(seq_part_geo(2,2)) = order_cent(2)
                  order_part_cent(seq_part_geo(3,2)) = 0
                  num_binom = order_cent(2)-1        !number of binomial terms
                  binom_cent(1) = seq_part_geo(2,2)  !centers of binomial terms
                  binom_cent(2) = seq_part_geo(3,2)
                  ! prepares binomial coefficients for later use
                  call pascal_triangle(order_cent(2), &
                                       binom_coef(0:order_cent(2),0:order_cent(2)))
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .true.
                ! only one partial geometric derivative term on the second and third atomic centers
                else
                  num_part_geo(1) = 2
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                  order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .false.
                end if
              else
                zero_ints = .true.
              end if
            else
              zero_ints = .true.
            end if
          ! three-center total geometric derivatives
          case(3)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2)) .and. &
                idx_cent(3)==idx_part_cent(id_atom_cent(3))) then
              ! binomial expansion on the first atomic center
              if (num_atom_iden(1)==2) then
                num_part_geo(1) = 3  !the first and last monomial terms have only one
                num_part_geo(2) = 4  !differentiated center, while others have two,
                num_part_geo(3) = 3  !we then add another two atomic centers
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,1) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(3,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(4,2) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(2,3) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,3) = tag_part_cent(id_atom_cent(3))
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                order_part_cent(seq_part_geo(3,2)) = order_cent(2)
                order_part_cent(seq_part_geo(4,2)) = order_cent(3)
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              ! binomial expansion on the second atomic center
              else if (num_atom_iden(2)==2) then
                num_part_geo(1) = 3  !the first and last monomial terms have only one
                num_part_geo(2) = 4  !differentiated center, while others have two,
                num_part_geo(3) = 3  !we then add another two atomic centers
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,1) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,2) = tag_part_cent(id_atom_cent(2)+1)
                seq_part_geo(4,2) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,3) = tag_part_cent(id_atom_cent(2)+1)
                seq_part_geo(3,3) = tag_part_cent(id_atom_cent(3))
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(2,2)) = order_cent(2)
                order_part_cent(seq_part_geo(3,2)) = 0
                order_part_cent(seq_part_geo(4,2)) = order_cent(3)
                num_binom = order_cent(2)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(2,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(3,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(2), &
                                     binom_coef(0:order_cent(2),0:order_cent(2)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              ! binomial expansion on the third atomic center
              else if (num_atom_iden(3)==2) then
                num_part_geo(1) = 3  !the first and last monomial terms have only one
                num_part_geo(2) = 4  !differentiated center, while others have two,
                num_part_geo(3) = 3  !we then add another two atomic centers
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,1) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,2) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(4,2) = tag_part_cent(id_atom_cent(3)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,3) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,3) = tag_part_cent(id_atom_cent(3)+1)
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = order_cent(2)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(3,2)) = order_cent(3)
                order_part_cent(seq_part_geo(4,2)) = 0
                num_binom = order_cent(3)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(3,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(4,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(3), &
                                     binom_coef(0:order_cent(3),0:order_cent(3)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              else
                stop "geom_part_two_param>> invalid num_atom_iden (0/3/3)!"
              end if
            else
              zero_ints = .true.
            end if
          case default
            stop "geom_part_two_param>> invalid num_cents (0/3/X)!"
          end select
        ! there are four unique atomic centers
        case(4)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            ! notice that atomic centers are sorted in ascending order, so that
            ! we divide them into two parts according to the third atomic center
            if (idx_cent(1)<idx_part_cent(id_atom_cent(3))) then
              ! partial geometric derivatives on the first atomic center
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! partial geometric derivatives on the second atomic center
              else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              else
                zero_ints = .true.
              end if
            else
              ! partial geometric derivatives on the third atomic center
              if (idx_cent(1)==idx_part_cent(id_atom_cent(3))) then
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(3))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! partial geometric derivatives on the fourth atomic center
              else if (idx_cent(1)==idx_part_cent(id_atom_cent(4))) then
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(4))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              else
                zero_ints = .true.
              end if
            end if
          ! two-center total geometric derivatives
          case(2)
            ! notice that both centers are sorted in ascending order, so that
            ! we divide them into two parts according to the third atomic center
            if (idx_cent(2)<idx_part_cent(id_atom_cent(3))) then
              ! partial geometric derivatives on the first and second atomic centers
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                  idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
                num_part_geo(1) = 2
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              else
                zero_ints = .true.
              end if
            else
              if (idx_cent(1)>idx_part_cent(id_atom_cent(2))) then
                ! partial geometric derivatives on the third and fourth atomic centers
                if (idx_cent(1)==idx_part_cent(id_atom_cent(3)) .and. &
                    idx_cent(2)==idx_part_cent(id_atom_cent(4))) then
                  num_part_geo(1) = 2
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(3))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(4))
                  order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .false.
                else
                  zero_ints = .true.
                end if
              else
                if (idx_cent(2)==idx_part_cent(id_atom_cent(3))) then
                  ! partial geometric derivatives on the first and third atomic centers
                  if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
                    num_part_geo(1) = 2
                    seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                    seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                    order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                    order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                    zero_ints = .false.
                    neg_one = .false.
                    binom_exp = .false.
                  ! partial geometric derivatives on the second and thrid atomic centers
                  else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
                    num_part_geo(1) = 2
                    seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                    seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                    order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                    order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                    zero_ints = .false.
                    neg_one = .false.
                    binom_exp = .false.
                  else
                    zero_ints = .true.
                  end if
                else if (idx_cent(2)==idx_part_cent(id_atom_cent(4))) then
                  ! partial geometric derivatives on the first and fourth atomic centers
                  if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
                    num_part_geo(1) = 2
                    seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                    seq_part_geo(2,1) = tag_part_cent(id_atom_cent(4))
                    order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                    order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                    zero_ints = .false.
                    neg_one = .false.
                    binom_exp = .false.
                  ! partial geometric derivatives on the second and fourth atomic centers
                  else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
                    num_part_geo(1) = 2
                    seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                    seq_part_geo(2,1) = tag_part_cent(id_atom_cent(4))
                    order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                    order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                    zero_ints = .false.
                    neg_one = .false.
                    binom_exp = .false.
                  else
                    zero_ints = .true.
                  end if
                else
                  zero_ints = .true.
                end if
              end if
            end if
          ! three-center total geometric derivatives
          case(3)
            ! notice that both centers are sorted in ascending order, so that
            ! we divide them into two parts according to the third atomic center
            if (idx_cent(2)<idx_part_cent(id_atom_cent(3))) then
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                  idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
                ! partial geometric derivatives on the first, second and third atomic centers
                if (idx_cent(3)==idx_part_cent(id_atom_cent(3))) then
                  num_part_geo(1) = 3
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(3,1) = tag_part_cent(id_atom_cent(3))
                  order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                  order_part_cent(seq_part_geo(3,1)) = order_cent(3)
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .false.
                ! partial geometric derivatives on the first, second and fourth atomic centers
                else if (idx_cent(3)==idx_part_cent(id_atom_cent(4))) then
                  num_part_geo(1) = 3
                  seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                  seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                  seq_part_geo(3,1) = tag_part_cent(id_atom_cent(4))
                  order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                  order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                  order_part_cent(seq_part_geo(3,1)) = order_cent(3)
                  zero_ints = .false.
                  neg_one = .false.
                  binom_exp = .false.
                else
                  zero_ints = .true.
                end if
              else
                zero_ints = .true.
              end if
            else if (idx_cent(2)==idx_part_cent(id_atom_cent(3)) .and. &
                     idx_cent(3)==idx_part_cent(id_atom_cent(4))) then
              ! partial geometric derivatives on the first, third and fourth atomic centers
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
                num_part_geo(1) = 3
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(3,1) = tag_part_cent(id_atom_cent(4))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                order_part_cent(seq_part_geo(3,1)) = order_cent(3)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! partial geometric derivatives on the second, third and fourth atomic centers
              else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
                num_part_geo(1) = 3 
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(3,1) = tag_part_cent(id_atom_cent(4))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                order_part_cent(seq_part_geo(3,1)) = order_cent(3)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              else
                zero_ints = .true.
              end if
            else
              zero_ints = .true.
            end if
          ! four-center total geometric derivatives
          case(4)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2)) .and. &
                idx_cent(3)==idx_part_cent(id_atom_cent(3)) .and. &
                idx_cent(4)==idx_part_cent(id_atom_cent(4))) then
                num_part_geo(1) = 4
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,1) = tag_part_cent(id_atom_cent(3))
                seq_part_geo(4,1) = tag_part_cent(id_atom_cent(4))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                order_part_cent(seq_part_geo(3,1)) = order_cent(3)
                order_part_cent(seq_part_geo(4,1)) = order_cent(4)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
            else
              zero_ints = .true.
            end if
          case default
            stop "geom_part_two_param>> invalid num_cents (0/3/X)!"
          end select
        case default
          stop "geom_part_two_param>> invalid num_atom_cents (0/X)!"
        end select
      ! one unique non-atomic center
      case(1)
        select case(num_atom_cents)
        ! there is one unique atomic center
        case(1)
          ! the only non-zero one-center total geometric derivatives
          if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
            select case(num_atom_iden(1))
            ! partial geometric derivative is performed only on this atomic center
            case(1)
              num_part_geo(1) = 1
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .false.
            ! binomial expansion, since there are also two identical non-atomic centers
            case(2)
              num_part_geo(1) = 1  !the first and last monomial terms have only one
              num_part_geo(2) = 2  !differentiated center, while others have two
              num_part_geo(3) = 1
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
              seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
              seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
              seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
              ! this gives the first monomial term, we loops over others outside
              order_part_cent(seq_part_geo(1,2)) = order_cent(1)
              order_part_cent(seq_part_geo(2,2)) = 0
              num_binom = order_cent(1)-1        !number of binomial terms
              binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
              binom_cent(2) = seq_part_geo(2,2)
              ! prepares binomial coefficients for later use
              call pascal_triangle(order_cent(1), &
                                   binom_coef(0:order_cent(1),0:order_cent(1)))
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .true.
            ! translation invariance by performing partial geometric derivatives
            ! on the unique non-atomic center
            case(3)
              num_part_geo(1) = 1
              seq_part_geo(1,1) = tag_part_cent(id_non_cent(1))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              zero_ints = .false.
              neg_one = mod(order_cent(1),2)==1
              binom_exp = .false.
            case default
              stop "geom_part_two_param>> invalid num_atom_iden(1) (1/1/1)!"
            end select
          else
            zero_ints = .true.
          end if
        ! there are two unique atomic centers
        case(2)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              select case(num_atom_iden(1))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! binomial expansion, since there are another unique atomic center and non-atomic center
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(1) (1/2/1)!"
              end select
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              select case(num_atom_iden(2))
              ! partial geometric derivative is performed only on this atomic center
              case(1)
                num_part_geo(1) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! binomial expansion, since there are another unique atomic center and non-atomic center
              case(2)
                num_part_geo(1) = 1  !the first and last monomial terms have only one
                num_part_geo(2) = 2  !differentiated center, while others have two
                num_part_geo(3) = 1
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(2)+1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              case default
                stop "geom_part_two_param>> invalid num_atom_iden(2) (1/2/1)!"
              end select
            else
              zero_ints = .true.
            end if
          ! two-center total geometric derivatives
          case(2)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
              ! binomial expansion on the first atomic center
              if (num_atom_iden(1)==2) then
                num_part_geo(1) = 2  !the first and last monomial terms have only one
                num_part_geo(2) = 3  !differentiated center, while others have two,
                num_part_geo(3) = 2  !we then add another atomic center
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(3,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
                seq_part_geo(2,3) = tag_part_cent(id_atom_cent(2))
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                order_part_cent(seq_part_geo(2,2)) = 0
                order_part_cent(seq_part_geo(3,2)) = order_cent(2)
                num_binom = order_cent(1)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(2,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(1), &
                                     binom_coef(0:order_cent(1),0:order_cent(1)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              ! binomial expansion on the second atomic center
              else if (num_atom_iden(2)==2) then
                num_part_geo(1) = 2  !the first and last monomial terms have only one
                num_part_geo(2) = 3  !differentiated center, while others have two,
                num_part_geo(3) = 2  !we then add another atomic center
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,2) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,2) = tag_part_cent(id_atom_cent(2)+1)
                seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,3) = tag_part_cent(id_atom_cent(2)+1)
                order_part_cent(seq_part_geo(1,2)) = order_cent(1)
                ! this gives the first monomial term, we loops over others outside
                order_part_cent(seq_part_geo(2,2)) = order_cent(2)
                order_part_cent(seq_part_geo(3,2)) = 0
                num_binom = order_cent(2)-1        !number of binomial terms
                binom_cent(1) = seq_part_geo(2,2)  !centers of binomial terms
                binom_cent(2) = seq_part_geo(3,2)
                ! prepares binomial coefficients for later use
                call pascal_triangle(order_cent(2), &
                                     binom_coef(0:order_cent(2),0:order_cent(2)))
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .true.
              ! only one partial geometric derivative term on the first and second atomic centers
              else
                num_part_geo(1) = 2
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              end if
            else
              zero_ints = .true.
            end if
          case default
            stop "geom_part_two_param>> invalid num_cents (1/2/X)!"
          end select
        ! there are three unique atomic centers, so only one non-atomic center
        case(3)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            ! partial geometric derivatives on the first atomic center
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              num_part_geo(1) = 1
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .false.
            ! partial geometric derivatives on the second atomic center
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              num_part_geo(1) = 1 
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .false.
            ! partial geometric derivatives on the third atomic center
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(3))) then
              num_part_geo(1) = 1 
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(3))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .false.
            else
              zero_ints = .true.
            end if
          ! two-center total geometric derivatives
          case(2)
            ! notice that both centers are sorted in ascending order, so that
            ! we divide them into two parts according to the third atomic center
            if (idx_cent(2)<idx_part_cent(id_atom_cent(3))) then
              ! partial geometric derivatives on the first and second atomic centers
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                  idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
                num_part_geo(1) = 2
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              else
                zero_ints = .true.
              end if
            else if (idx_cent(2)==idx_part_cent(id_atom_cent(3))) then
              ! partial geometric derivatives on the first and third atomic centers
              if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
                num_part_geo(1) = 2
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              ! partial geometric derivatives on the second and thrid atomic centers
              else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
                num_part_geo(1) = 2
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(3))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
              else
                zero_ints = .true.
              end if
            else
              zero_ints = .true.
            end if
          ! three-center total geometric derivatives
          case(3)
            ! notice that both center sequences are sorted in ascending order
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2)) .and. &
                idx_cent(3)==idx_part_cent(id_atom_cent(3))) then
                num_part_geo(1) = 3
                seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
                seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
                seq_part_geo(3,1) = tag_part_cent(id_atom_cent(3))
                order_part_cent(seq_part_geo(1,1)) = order_cent(1)
                order_part_cent(seq_part_geo(2,1)) = order_cent(2)
                order_part_cent(seq_part_geo(3,1)) = order_cent(3)
                zero_ints = .false.
                neg_one = .false.
                binom_exp = .false.
            else
              zero_ints = .true.
            end if
          case default
            stop "geom_part_two_param>> invalid num_cents (1/3/X)!"
          end select
        case default
          stop "geom_part_two_param>> invalid num_atom_cents (1/X)!"
        end select
      ! two unique non-atomic centers, so there are two atomic centers (might be identical)
      case(2)
        select case(num_atom_cents)
        ! there is one unique atomic center (two identical atomic centers)
        case(1)
          ! the only non-zero one-center total geometric derivatives, by binomial expansion
          if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
            num_part_geo(1) = 1  !the first and last monomial terms have only one
            num_part_geo(2) = 2  !differentiated center, while others have two
            num_part_geo(3) = 1
            seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
            seq_part_geo(1,2) = tag_part_cent(id_atom_cent(1))
            seq_part_geo(2,2) = tag_part_cent(id_atom_cent(1)+1)
            seq_part_geo(1,3) = tag_part_cent(id_atom_cent(1)+1)
            ! this gives the first monomial term, we loops over others outside
            order_part_cent(seq_part_geo(1,2)) = order_cent(1)
            order_part_cent(seq_part_geo(2,2)) = 0
            num_binom = order_cent(1)-1        !number of binomial terms
            binom_cent(1) = seq_part_geo(1,2)  !centers of binomial terms
            binom_cent(2) = seq_part_geo(2,2)
            ! prepares binomial coefficients for later use
            call pascal_triangle(order_cent(1), &
                                 binom_coef(0:order_cent(1),0:order_cent(1)))
            zero_ints = .false.
            neg_one = .false.
            binom_exp = .true.
          else
            zero_ints = .true.
          end if
        ! there are two unique atomic centers
        case(2)
          select case(num_cents)
          ! one-center total geometric derivatives
          case(1)
            ! partial geometric derivatives on the first atomic center
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1))) then
              num_part_geo(1) = 1
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .false.
            ! partial geometric derivatives on the second atomic center
            else if (idx_cent(1)==idx_part_cent(id_atom_cent(2))) then
              num_part_geo(1) = 1
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(2))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .false.
            else
              zero_ints = .true.
            end if
          ! two-center total geometric derivatives
          case(2)
            ! partial geometric derivatives on these two atomic centers
            if (idx_cent(1)==idx_part_cent(id_atom_cent(1)) .and. &
                idx_cent(2)==idx_part_cent(id_atom_cent(2))) then
              num_part_geo(1) = 2
              seq_part_geo(1,1) = tag_part_cent(id_atom_cent(1))
              seq_part_geo(2,1) = tag_part_cent(id_atom_cent(2))
              order_part_cent(seq_part_geo(1,1)) = order_cent(1)
              order_part_cent(seq_part_geo(2,1)) = order_cent(2)
              zero_ints = .false.
              neg_one = .false.
              binom_exp = .false.
            else
              zero_ints = .true.
            end if
          case default
            stop "geom_part_two_param>> invalid num_cents (2/2/X)!"
          end select
        case default
          stop "geom_part_two_param>> invalid num_atom_cents (2/X)!"
        end select
      case default
        stop "geom_part_two_param>> invalid num_non_cents!"
      end select
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "geom_part_two_param", STDOUT)
#endif
    return
  end subroutine geom_part_two_param

!  !> \brief 
!  !> \author Bin Gao
!  !> \date 2011-07-02
!  !> \param 
!  !> \param 
!  !> \param 
!  !> \return 
!  subroutine geom_part_two_assign()
!    use xkind
    implicit none
!    integer, intent(in) :: dim_cints
!    integer, intent(in) :: dim_geo_op1
!    integer, intent(in) :: dim_geo_op2
!    integer, intent(in) :: dim_mag
!    integer, intent(in) :: dim_geo_bra
!    integer, intent(in) :: dim_geo_ket
!    real(REALK), intent(in) :: part_geo(dim_cints,dim_geo_op1,dim_geo_op2,dim_mag, &
!                                    dim_geo_bra,dim_geo_ket)
!    integer, intent(in) :: num_geo_op1
!    integer, intent(in) :: num_geo_op2
!    integer, intent(in) :: num_geo_bra
!    integer, intent(in) :: num_geo_ket
!    integer, intent(in) :: num_geo_total
!    real(REALK), intent(inout) :: total_geo(dim_cints,num_geo_op1,num_geo_op2,dim_mag, &
!                                        num_geo_bra,num_geo_ket,num_geo_total)
!    return
!  end subroutine geom_part_two_assign
!
!  subroutine geom_part_two_binom
!    use xkind
    implicit none
!    logical, intent(in) :: neg_one
!    integer, intent(in) :: binom_order
!    integer, intent(in) :: monom_order
!    real(REALK), intent(in) :: binom_coef(0:binom_order,0:binom_order)
!    return
!  end subroutine geom_part_two_binom
!
!  subroutine geom_part_two_add
!    use xkind
    implicit none
!    return
!  end subroutine geom_part_two_add
