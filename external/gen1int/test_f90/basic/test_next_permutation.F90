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
!!  Tests subroutine \fn(next_permutation).
!!
!!  2011-11-07, Bin Gao
!!  * first version

  !> \brief tests subroutine \fn(next_permutation)
  !> \author Bin Gao
  !> \date 2011-11-07
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_next_permutation(io_log, test_failed)
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
    integer, parameter :: NUM_IDX = 10     !number of indices
    integer, parameter :: FIRST = 3        !index of the first number in permutation
    integer, parameter :: LAST = 8         !index of the last number in permutation
    integer num_perm                       !number of permutations
    integer, allocatable :: ref_perm(:,:)  !referenced permutations
    integer perm_num(NUM_IDX)              !permutation from \fn(next_permutation)
    logical do_perm                        !indicates if there exits next permutation
    integer iperm, idx                     !incremental recorders
    integer ierr                           !error information
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_int_number("Number of indices:", NUM_IDX, "I3", io_log)
    ! assigns the first permutation and referenced permutations
    perm_num = (/1,2,3,3,3,6,7,7,9,10/)
    num_perm = 6*5*4*3*2/(3*2)/2
    allocate(ref_perm(NUM_IDX,num_perm), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate ref_perm:", &
                               log_int=NUM_IDX*num_perm, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ref_perm(:, 1) = (/1, 2, 3, 3, 3, 7, 6, 7, 9, 10/)
    ref_perm(:, 2) = (/1, 2, 3, 3, 3, 7, 7, 6, 9, 10/)
    ref_perm(:, 3) = (/1, 2, 3, 3, 6, 3, 7, 7, 9, 10/)
    ref_perm(:, 4) = (/1, 2, 3, 3, 6, 7, 3, 7, 9, 10/)
    ref_perm(:, 5) = (/1, 2, 3, 3, 6, 7, 7, 3, 9, 10/)
    ref_perm(:, 6) = (/1, 2, 3, 3, 7, 3, 6, 7, 9, 10/)
    ref_perm(:, 7) = (/1, 2, 3, 3, 7, 3, 7, 6, 9, 10/)
    ref_perm(:, 8) = (/1, 2, 3, 3, 7, 6, 3, 7, 9, 10/)
    ref_perm(:, 9) = (/1, 2, 3, 3, 7, 6, 7, 3, 9, 10/)
    ref_perm(:,10) = (/1, 2, 3, 3, 7, 7, 3, 6, 9, 10/)
    ref_perm(:,11) = (/1, 2, 3, 3, 7, 7, 6, 3, 9, 10/)
    ref_perm(:,12) = (/1, 2, 3, 6, 3, 3, 7, 7, 9, 10/)
    ref_perm(:,13) = (/1, 2, 3, 6, 3, 7, 3, 7, 9, 10/)
    ref_perm(:,14) = (/1, 2, 3, 6, 3, 7, 7, 3, 9, 10/)
    ref_perm(:,15) = (/1, 2, 3, 6, 7, 3, 3, 7, 9, 10/)
    ref_perm(:,16) = (/1, 2, 3, 6, 7, 3, 7, 3, 9, 10/)
    ref_perm(:,17) = (/1, 2, 3, 6, 7, 7, 3, 3, 9, 10/)
    ref_perm(:,18) = (/1, 2, 3, 7, 3, 3, 6, 7, 9, 10/)
    ref_perm(:,19) = (/1, 2, 3, 7, 3, 3, 7, 6, 9, 10/)
    ref_perm(:,20) = (/1, 2, 3, 7, 3, 6, 3, 7, 9, 10/)
    ref_perm(:,21) = (/1, 2, 3, 7, 3, 6, 7, 3, 9, 10/)
    ref_perm(:,22) = (/1, 2, 3, 7, 3, 7, 3, 6, 9, 10/)
    ref_perm(:,23) = (/1, 2, 3, 7, 3, 7, 6, 3, 9, 10/)
    ref_perm(:,24) = (/1, 2, 3, 7, 6, 3, 3, 7, 9, 10/)
    ref_perm(:,25) = (/1, 2, 3, 7, 6, 3, 7, 3, 9, 10/)
    ref_perm(:,26) = (/1, 2, 3, 7, 6, 7, 3, 3, 9, 10/)
    ref_perm(:,27) = (/1, 2, 3, 7, 7, 3, 3, 6, 9, 10/)
    ref_perm(:,28) = (/1, 2, 3, 7, 7, 3, 6, 3, 9, 10/)
    ref_perm(:,29) = (/1, 2, 3, 7, 7, 6, 3, 3, 9, 10/)
    ref_perm(:,30) = (/1, 2, 6, 3, 3, 3, 7, 7, 9, 10/)
    ref_perm(:,31) = (/1, 2, 6, 3, 3, 7, 3, 7, 9, 10/)
    ref_perm(:,32) = (/1, 2, 6, 3, 3, 7, 7, 3, 9, 10/)
    ref_perm(:,33) = (/1, 2, 6, 3, 7, 3, 3, 7, 9, 10/)
    ref_perm(:,34) = (/1, 2, 6, 3, 7, 3, 7, 3, 9, 10/)
    ref_perm(:,35) = (/1, 2, 6, 3, 7, 7, 3, 3, 9, 10/)
    ref_perm(:,36) = (/1, 2, 6, 7, 3, 3, 3, 7, 9, 10/)
    ref_perm(:,37) = (/1, 2, 6, 7, 3, 3, 7, 3, 9, 10/)
    ref_perm(:,38) = (/1, 2, 6, 7, 3, 7, 3, 3, 9, 10/)
    ref_perm(:,39) = (/1, 2, 6, 7, 7, 3, 3, 3, 9, 10/)
    ref_perm(:,40) = (/1, 2, 7, 3, 3, 3, 6, 7, 9, 10/)
    ref_perm(:,41) = (/1, 2, 7, 3, 3, 3, 7, 6, 9, 10/)
    ref_perm(:,42) = (/1, 2, 7, 3, 3, 6, 3, 7, 9, 10/)
    ref_perm(:,43) = (/1, 2, 7, 3, 3, 6, 7, 3, 9, 10/)
    ref_perm(:,44) = (/1, 2, 7, 3, 3, 7, 3, 6, 9, 10/)
    ref_perm(:,45) = (/1, 2, 7, 3, 3, 7, 6, 3, 9, 10/)
    ref_perm(:,46) = (/1, 2, 7, 3, 6, 3, 3, 7, 9, 10/)
    ref_perm(:,47) = (/1, 2, 7, 3, 6, 3, 7, 3, 9, 10/)
    ref_perm(:,48) = (/1, 2, 7, 3, 6, 7, 3, 3, 9, 10/)
    ref_perm(:,49) = (/1, 2, 7, 3, 7, 3, 3, 6, 9, 10/)
    ref_perm(:,50) = (/1, 2, 7, 3, 7, 3, 6, 3, 9, 10/)
    ref_perm(:,51) = (/1, 2, 7, 3, 7, 6, 3, 3, 9, 10/)
    ref_perm(:,52) = (/1, 2, 7, 6, 3, 3, 3, 7, 9, 10/)
    ref_perm(:,53) = (/1, 2, 7, 6, 3, 3, 7, 3, 9, 10/)
    ref_perm(:,54) = (/1, 2, 7, 6, 3, 7, 3, 3, 9, 10/)
    ref_perm(:,55) = (/1, 2, 7, 6, 7, 3, 3, 3, 9, 10/)
    ref_perm(:,56) = (/1, 2, 7, 7, 3, 3, 3, 6, 9, 10/)
    ref_perm(:,57) = (/1, 2, 7, 7, 3, 3, 6, 3, 9, 10/)
    ref_perm(:,58) = (/1, 2, 7, 7, 3, 6, 3, 3, 9, 10/)
    ref_perm(:,59) = (/1, 2, 7, 7, 6, 3, 3, 3, 9, 10/)
    ref_perm(:,60) = (/1, 2, 3, 3, 3, 6, 7, 7, 9, 10/)
    ! loops till we get the last permutation
    iperm = 0
    do_perm = .true.
    do while (do_perm .and. iperm<num_perm)
      call next_permutation(FIRST, LAST, perm_num, do_perm)
      iperm = iperm+1
      ! checks the results
      do idx = 1, NUM_IDX
        if (perm_num(idx)/=ref_perm(idx,iperm)) then
          call html_log_int_array(log_text="Index and permutation:", &
                                  log_int=(/idx, iperm/), fmt_int="I3", io_log=io_log)
          call html_log_int_number(log_text="Referenced permutation:",        &
                                   log_int=ref_perm(idx,iperm), fmt_int="I3", &
                                   io_log=io_log, font_color="blue")
          call html_log_int_number(log_text="Permutation from next_permutation:",      &
                                   log_int=perm_num(idx), fmt_int="I3", io_log=io_log, &
                                   font_color="red")
          test_failed = .true.
        end if
      end do
    end do
    if (iperm/=num_perm) then
      call html_log_int_number(log_text="Referenced number of permutations:", &
                               log_int=num_perm, fmt_int="I4", io_log=io_log, &
                               font_color="blue")
      call html_log_int_number(                                        &
             log_text="Number of permutations from next_permutation:", &
             log_int=iperm, fmt_int="I4", io_log=io_log, font_color="red")

      test_failed = .true.
    end if
    deallocate(ref_perm)
    return
  end subroutine test_next_permutation
