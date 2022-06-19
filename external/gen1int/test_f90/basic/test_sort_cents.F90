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
!!  Tests subroutines in sort_cents.F90.
!!
!!  2011-07-27, Bin Gao
!!  * first version

#include "max_idx_non.h"

  !> \brief tests subroutines in sort_cents.F90
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_sort_cents(io_log, test_failed)
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
    ! see \fn(sort_gen_cents) and \fn(sort_atom_cents) for the details of the following variables
    integer, parameter :: NUM_CENTS = 10
    integer idx_sort_cent(NUM_CENTS)
    integer tag_cent(NUM_CENTS)
    integer num_non_cents
    integer id_non_cent(NUM_CENTS)
    integer num_non_iden(NUM_CENTS)
    integer num_atom_cents
    integer id_atom_cent(NUM_CENTS)
    integer num_atom_iden(NUM_CENTS)
    integer idx_ref_cent(NUM_CENTS)  !indices of centers before sorting
    integer num_ref_cents            !referenced number of unique centers
    integer id_ref_cent(NUM_CENTS)   !referenced places of unique centers
    integer num_ref_iden(NUM_CENTS)  !referenced number of identical centers
    integer icent, jcent, kcent      !incremental recorder
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_heading(heading="Tests subroutine sort_gen_cents", &
                          io_log=io_log, level=4)
    call html_log_int_number("Number of centers:", NUM_CENTS, "I4", io_log)
    ! assigns the indices of the centers
    do icent = 1, NUM_CENTS
      tag_cent(icent) = icent
      idx_sort_cent(icent) = min(icent,NUM_CENTS-icent)
    end do
    idx_sort_cent(1) = MAX_IDX_NON
    idx_sort_cent(NUM_CENTS) = MAX_IDX_NON
    idx_ref_cent = idx_sort_cent
    call html_log_int_array(log_text="Indices of the centers before sorting:", &
                            log_int=idx_sort_cent, fmt_int="I3", io_log=io_log)
    ! tests \fn(sort_gen_cents)
    call sort_gen_cents(NUM_CENTS, idx_sort_cent, tag_cent,       &
                        num_non_cents, id_non_cent, num_non_iden, &
                        num_atom_cents, id_atom_cent, num_atom_iden)
    call html_log_int_array(log_text="Indices of the centers after sorting:", &
                            log_int=idx_sort_cent, fmt_int="I3", io_log=io_log)
    call html_log_int_number("Number of unique non-atomic centers:", &
                             num_non_cents, "I4", io_log)
    call html_log_int_array(log_text="Places of unique non-atomic centers:",    &
                            log_int=id_non_cent(1:num_non_cents), fmt_int="I3", &
                            io_log=io_log)
    call html_log_int_array(log_text="Number of identical non-atomic centers:",  &
                            log_int=num_non_iden(1:num_non_cents), fmt_int="I3", &
                            io_log=io_log)
    call html_log_int_number("Number of unique atomic centers:", &
                             num_atom_cents, "I4", io_log)
    call html_log_int_array(log_text="Places of unique atomic centers:", &
                            log_int=id_atom_cent(1:num_atom_cents),      &
                            fmt_int="I3", io_log=io_log)
    call html_log_int_array(log_text="Number of identical atomic centers:", &
                            log_int=num_atom_iden(1:num_atom_cents),        &
                            fmt_int="I3", io_log=io_log)
    ! checks the results
    do icent = 1, NUM_CENTS
      ! we should be able to find the corresponding center before sorting
      if (idx_sort_cent(icent)/=idx_ref_cent(tag_cent(icent))) then
        call html_log_int_number("Center:", icent, "I4", io_log)
        call html_log_int_number(log_text="Referenced index:",          &
                                 log_int=idx_ref_cent(tag_cent(icent)), &
                                 fmt_int="I3", io_log=io_log, font_color="blue")
        call html_log_int_number(log_text="Result from sort_gen_cents:",     &
                                 log_int=idx_sort_cent(icent), fmt_int="I3", &
                                 io_log=io_log, font_color="red")
        test_failed = .true.
      end if
      ! sorting centers should be in ascending order
      if (icent<NUM_CENTS) then
        if (idx_sort_cent(icent)>idx_sort_cent(icent+1)) then
          call html_log_int_number("Center:", icent, "I4", io_log)
          call html_log_int_number(log_text="Result from sort_gen_cents:",     &
                                   log_int=idx_sort_cent(icent), fmt_int="I3", &
                                   io_log=io_log, font_color="red")
          call html_log_int_number("Center:", icent+1, "I4", io_log)
          call html_log_int_number(log_text="Result from sort_gen_cents:",       &
                                   log_int=idx_sort_cent(icent+1), fmt_int="I3", &
                                   io_log=io_log, font_color="red")
          call html_log_text(log_text="Index of center from sort_gen_cents "// &
                                      "is less than its following center.",    &
                             io_log=io_log, font_color="red")
          test_failed = .true.
        end if
      end if
    end do
    if (.not.test_failed) then
      ! gets the referenced unique non-atomic centers
      if (idx_sort_cent(1)<=MAX_IDX_NON) then
        jcent = 1
        num_ref_cents = 1
        id_ref_cent(jcent) = 1
        num_ref_iden(jcent) = 1
        do icent = 2, NUM_CENTS
          ! exists for atomic center
          if (idx_sort_cent(icent)>MAX_IDX_NON) then
            exit
          else
            ! found another unique atomic center
            if (idx_sort_cent(icent)/=idx_sort_cent(id_ref_cent(jcent))) then
              jcent = jcent+1
              id_ref_cent(jcent) = icent
              num_ref_iden(jcent) = 1
              num_ref_cents = num_ref_cents+1
            else
              num_ref_iden(jcent) = num_ref_iden(jcent)+1
            end if
          end if
        end do
        ! checks the unique non-atomic centers
        if (num_non_cents/=num_ref_cents) then
          call html_log_int_number(                                          &
                 log_text="Referenced number of unique non-atomic centers:", &
                 log_int=num_ref_cents, fmt_int="I3", io_log=io_log, font_color="blue")
          call html_log_int_number(log_text="Result from sort_gen_cents:", &
                                   log_int=num_non_cents, fmt_int="I3",    &
                                   io_log=io_log, font_color="red")
          test_failed = .true.
        else
          do kcent = 1, num_non_cents
            if (id_non_cent(kcent)/=id_ref_cent(kcent)) then
              call html_log_int_number("Unique non-atomic center:", kcent, &
                                       "I4", io_log)
              call html_log_int_number(                                        &
                     log_text="Referenced place of unique non-atomic center:", &
                     log_int=id_ref_cent(kcent), fmt_int="I3", io_log=io_log,  &
                     font_color="blue")
              call html_log_int_number(log_text="Result from sort_gen_cents:",   &
                                       log_int=id_non_cent(kcent), fmt_int="I3", &
                                       io_log=io_log, font_color="red")
              test_failed = .true.
            end if
            if (num_non_iden(kcent)/=num_ref_iden(kcent)) then
              call html_log_int_number("Identical non-atomic centers:", kcent, &
                                       "I4", io_log)
              call html_log_int_number(                                             &
                     log_text="Referenced number of identical non-atomic centers:", &
                     log_int=num_ref_iden(kcent), fmt_int="I3", io_log=io_log,      &
                     font_color="blue")
              call html_log_int_number(log_text="Result from sort_gen_cents:",    &
                                       log_int=num_non_iden(kcent), fmt_int="I3", &
                                       io_log=io_log, font_color="red")
              test_failed = .true.
            end if
          end do
        end if
      else
        icent = 1
      end if
      ! gets the referenced unique atomic centers
      jcent = 1
      num_ref_cents = 1
      id_ref_cent(jcent) = icent
      num_ref_iden(jcent) = 1
      do kcent = icent+1, NUM_CENTS
        ! found another unique atomic center
        if (idx_sort_cent(kcent)/=idx_sort_cent(id_ref_cent(jcent))) then
          jcent = jcent+1
          id_ref_cent(jcent) = kcent
          num_ref_iden(jcent) = 1
          num_ref_cents = num_ref_cents+1
        else
          num_ref_iden(jcent) = num_ref_iden(jcent)+1
        end if
      end do
      ! checks the unique atomic centers
      if (num_atom_cents/=num_ref_cents) then
        call html_log_int_number(                                      &
               log_text="Referenced number of unique atomic centers:", &
               log_int=num_ref_cents, fmt_int="I3", io_log=io_log, font_color="blue")
        call html_log_int_number(log_text="Result from sort_gen_cents:", &
                                 log_int=num_atom_cents, fmt_int="I3",   &
                                 io_log=io_log, font_color="red")
        test_failed = .true.
      else
        do kcent = 1, num_atom_cents
          if (id_atom_cent(kcent)/=id_ref_cent(kcent)) then
            call html_log_int_number("Unique atomic center:", kcent, &
                                     "I4", io_log)
            call html_log_int_number(                                       &
                   log_text="Referenced place of unique atomic center:",    &
                   log_int=id_ref_cent(kcent), fmt_int="I3", io_log=io_log, &
                   font_color="blue")
            call html_log_int_number(log_text="Result from sort_gen_cents:",    &
                                     log_int=id_atom_cent(kcent), fmt_int="I3", &
                                     io_log=io_log, font_color="red")
            test_failed = .true.
          end if
          if (num_atom_iden(kcent)/=num_ref_iden(kcent)) then
            call html_log_int_number("Identical atomic centers:", kcent, &
                                     "I4", io_log)
            call html_log_int_number(                                         &
                   log_text="Referenced number of identical atomic centers:", &
                   log_int=num_ref_iden(kcent), fmt_int="I3", io_log=io_log,  &
                   font_color="blue")
            call html_log_int_number(log_text="Result from sort_gen_cents:",     &
                                     log_int=num_atom_iden(kcent), fmt_int="I3", &
                                     io_log=io_log, font_color="red")
            test_failed = .true.
          end if
        end do
      end if
    end if
    ! dumps information of tests
    call html_log_heading(heading="Tests subroutine sort_atom_cents", &
                          io_log=io_log, level=4)
    call html_log_int_number("Number of centers:", NUM_CENTS, "I4", io_log)
    ! assigns the indices of the atomic centers
    do icent = 1, NUM_CENTS
      tag_cent(icent) = icent
      idx_sort_cent(icent) = min(icent,NUM_CENTS-icent+1)
    end do
    idx_ref_cent = idx_sort_cent
    call html_log_int_array(log_text="Indices of the centers before sorting:", &
                            log_int=idx_sort_cent, fmt_int="I3", io_log=io_log)
    ! test \fn(sort_atom_cents)
    call sort_atom_cents(NUM_CENTS, idx_sort_cent, tag_cent, &
                         num_atom_cents, id_atom_cent, num_atom_iden)
    call html_log_int_array(log_text="Indices of the centers after sorting:", &
                            log_int=idx_sort_cent, fmt_int="I3", io_log=io_log)
    call html_log_int_number("Number of unique atomic centers:", &
                             num_atom_cents, "I4", io_log)
    call html_log_int_array(log_text="Places of unique atomic centers:", &
                            log_int=id_atom_cent(1:num_atom_cents),      &
                            fmt_int="I3", io_log=io_log)
    call html_log_int_array(log_text="Number of identical atomic centers:", &
                            log_int=num_atom_iden(1:num_atom_cents),        &
                            fmt_int="I3", io_log=io_log)
    ! checks the results
    do icent = 1, NUM_CENTS
      ! we should be able to find the corresponding center before sorting
      if (idx_sort_cent(icent)/=idx_ref_cent(tag_cent(icent))) then
        call html_log_int_number("Center:", icent, "I4", io_log)
        call html_log_int_number(log_text="Referenced index:",          &
                                 log_int=idx_ref_cent(tag_cent(icent)), &
                                 fmt_int="I3", io_log=io_log, font_color="blue")
        call html_log_int_number(log_text="Result from sort_atom_cents:",    &
                                 log_int=idx_sort_cent(icent), fmt_int="I3", &
                                 io_log=io_log, font_color="red")
        test_failed = .true.
      end if
      ! sorting centers should be in ascending order
      if (icent<NUM_CENTS) then
        if (idx_sort_cent(icent)>idx_sort_cent(icent+1)) then
          call html_log_int_number("Center:", icent, "I4", io_log)
          call html_log_int_number(log_text="Result from sort_atom_cents:",    &
                                   log_int=idx_sort_cent(icent), fmt_int="I3", &
                                   io_log=io_log, font_color="red")
          call html_log_int_number("Center:", icent+1, "I4", io_log)
          call html_log_int_number(log_text="Result from sort_atom_cents:",      &
                                   log_int=idx_sort_cent(icent+1), fmt_int="I3", &
                                   io_log=io_log, font_color="red")
          call html_log_text(log_text="Index of center from sort_atom_cents "// &
                                      "is less than its following center.",     &
                             io_log=io_log, font_color="red")
          test_failed = .true.
        end if
      end if
    end do
    if (.not.test_failed) then
      ! gets the referenced unique atomic centers
      jcent = 1
      num_ref_cents = 1
      id_ref_cent(jcent) = 1
      num_ref_iden(jcent) = 1
      do icent = 2, NUM_CENTS
        ! found another unique atomic center
        if (idx_sort_cent(icent)/=idx_sort_cent(id_ref_cent(jcent))) then
          jcent = jcent+1
          id_ref_cent(jcent) = icent
          num_ref_iden(jcent) = 1
          num_ref_cents = num_ref_cents+1
        else
          num_ref_iden(jcent) = num_ref_iden(jcent)+1
        end if
      end do
      ! checks the unique atomic centers
      if (num_atom_cents/=num_ref_cents) then
        call html_log_int_number(                                      &
               log_text="Referenced number of unique atomic centers:", &
               log_int=num_ref_cents, fmt_int="I3", io_log=io_log, font_color="blue")
        call html_log_int_number(log_text="Result from sort_atom_cents:", &
                                 log_int=num_atom_cents, fmt_int="I3",    &
                                 io_log=io_log, font_color="red")
        test_failed = .true.
      else
        do kcent = 1, num_atom_cents
          if (id_atom_cent(kcent)/=id_ref_cent(kcent)) then
            call html_log_int_number("Unique atomic center:", kcent, &
                                     "I4", io_log)
            call html_log_int_number(                                       &
                   log_text="Referenced place of unique atomic center:",    &
                   log_int=id_ref_cent(kcent), fmt_int="I3", io_log=io_log, &
                   font_color="blue")
            call html_log_int_number(log_text="Result from sort_atom_cents:",   &
                                     log_int=id_atom_cent(kcent), fmt_int="I3", &
                                     io_log=io_log, font_color="red")
            test_failed = .true.
          end if
          if (num_atom_iden(kcent)/=num_ref_iden(kcent)) then
            call html_log_int_number("Identical atomic centers:", kcent, &
                                     "I4", io_log)
            call html_log_int_number(                                         &
                   log_text="Referenced number of identical atomic centers:", &
                   log_int=num_ref_iden(kcent), fmt_int="I3", io_log=io_log,  &
                   font_color="blue")
            call html_log_int_number(log_text="Result from sort_atom_cents:",    &
                                     log_int=num_atom_iden(kcent), fmt_int="I3", &
                                     io_log=io_log, font_color="red")
            test_failed = .true.
          end if
        end do
      end if
    end if
    return
  end subroutine test_sort_cents
