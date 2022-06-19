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
!!  Tests subroutines in geom_part_two.F90.
!!
!!  2011-08-02, Bin Gao
!!  * first version

#include "private/tag_cent.h"

  !> \brief tests subroutines in geom_part_two.F90
  !> \author Bin Gao
  !> \date 2011-08-02
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_geom_part_two(io_log, test_failed)
    use xkind
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters and referenced results
#include "geom_part_two.h"
    integer itst, jtst, ktst  !incremental recorders over test cases
    ! we first test shell_scatter.F90 which will be used for testing \fn(geom_part_two_scatter)
    write(io_log,100) "we first test subroutine shell_scatter ..."
    call test_shell_scatter(io_log, test_failed)
    if (test_failed) then
      write(io_log,999) "tests subroutine shell_scatter failed! no further test will be performed!"
      return
    else
      write(io_log,100) "tests subroutine shell_scatter passed ..."
    end if
    ! loops over test cases
    jtst = 1
    do itst = 1, NUM_TOT_GEO
      ! dumps information of tests
      write(io_log,110) "tests subroutine geom_part_two_param with number of centers", &
                        NUM_CENTS(itst)
      write(io_log,120) "indices of bra, ket and operator centers", IDX_BK_OPT(jtst:jtst+2)
      ! tests \fn(geom_part_two_param)
      ktst = jtst+NUM_CENTS(itst)-1
      call test_geom_part_two_param(NUM_CENTS(itst), IDX_GEO_CENT(jtst:ktst),           &
                                    ORDER_GEO_CENT(jtst:ktst), IDX_BK_OPT(jtst:jtst+2), &
                                    ORDER_GEO_BRA(itst), ORDER_GEO_KET(itst),           &
                                    ORDER_GEO_OPT(itst), REF_ORDER_PART(jtst:jtst+2),   &
                                    REF_ZERO_INTS(itst), REF_NEG_ONE(itst),             &
                                    REF_SCATTER_DER(itst), REF_SEQ_PART(jtst:ktst),     &
                                    io_log, test_failed)
      jtst = jtst+3
    end do
    return
100 format("<em>",A,"</em><br>")
110 format("<em><b>",A,"</b></em>: ",I3,"<br>")
120 format("<em>",A,"</em>: ",2(I3,","),I3,"<br>")
999 format("<font color=""red"">",A,"!</font><br>")
  end subroutine test_geom_part_two

  !> \brief tests subroutine geom_part_two_param
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param num_cents is the number of differentiated centers, could be 1, 2 or 3
  !> \param idx_cent contains the indices of different differentiated centers
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param idx_part_cent contains the indices of bra, ket and operator centers
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param order_geo_opt is the order of partial geometric derivatives on operator center
  !> \param ref_order_part contains referenced final orders of partial geometric derivatives
  !>        by adding those from total geometric derivatives
  !> \param ref_zero_ints indicates if the total geometric derivatives are zero
  !> \param ref_neg_one indicates if the integrals will be multiplied by -1
  !> \param ref_scatter_der indicates if scattering the geometric derivatives later on
  !> \param ref_seq_part contains the sequence of bra, ket and operator centers for partial derivative terms
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_geom_part_two_param(num_cents, idx_cent, order_cent, idx_part_cent, &
                                      order_geo_bra, order_geo_ket, order_geo_opt,    &
                                      ref_order_part,  ref_zero_ints, ref_neg_one,    &
                                      ref_scatter_deriv, ref_seq_part, io_log, test_failed)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: idx_part_cent(3)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_opt
    integer, intent(in) :: ref_order_part(3)
    logical, intent(in) :: ref_zero_ints
    logical, intent(in) :: ref_neg_one
    logical, intent(in) :: ref_scatter_deriv
    integer, intent(in) :: ref_seq_part(num_cents)
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
    ! see \fn(geom_part_two_param) for the details of the following variables
    integer order_part_cent(3)
    logical zero_ints
    logical neg_one
    logical scatter_deriv
    integer seq_part_geo(3)
    ! dumps information of tests
    select case(num_cents)
    case(1)
      write(io_log,100) "indices of differentiated centers", idx_cent
      write(io_log,100) "orders of differentiated centers", order_cent
    case(2)
      write(io_log,110) "indices of differentiated centers", idx_cent
      write(io_log,110) "orders of differentiated centers", order_cent
    case(3)
      write(io_log,120) "indices of differentiated centers", idx_cent
      write(io_log,120) "orders of differentiated centers", order_cent
    case default
      write(io_log,100) "zero derivatives due to too many differentiated centers", num_cents
      return
    end select
    ! sets the orders of partial geometric derivatives
    order_part_cent = (/order_geo_bra,order_geo_ket,order_geo_opt/)
    write(io_log,120) "initial orders of partial geometric derivatives", order_part_cent
    ! tests \fn(geom_part_two_param)
    call geom_part_two_param(num_cents, idx_cent, order_cent, idx_part_cent,     &
                             order_part_cent, zero_ints, neg_one, scatter_deriv, &
                             seq_part_geo)
    ! zero integrals
    if (ref_zero_ints) then
      if (zero_ints) then
        write(io_log,130) "total geometric derivatives are zero"
      else
        write(io_log,130) "total geometric derivatives should be zero"
        write(io_log,999) "zero_ints from geom_part_two_param is wrong"
        test_failed = .true.
      end if
    ! non-zero integrals
    else
      if (zero_ints) then
        write(io_log,130) "total geometric derivatives should not be zero"
        write(io_log,999) "zero_ints from geom_part_two_param is wrong"
        test_failed = .true.
      else
        ! checks the resulted orders of partial geometric derivatives
        if (any(order_part_cent/=ref_order_part)) then
          write(io_log,120) "correct final orders of partial geometric derivatives", &
                            ref_order_part
          write(io_log,998) "results from geom_part_two_param are wrong", order_part_cent
          test_failed = .true.
        else
          write(io_log,120) "final orders of partial geometric derivatives", order_part_cent
        end if
        ! integrals will be multiplied by -1
        if (ref_neg_one) then
          write(io_log,130) "integrals will be multiplied by -1"
          if (.not.neg_one) then
            write(io_log,999) "neg_one from geom_part_two_param is wrong"
            test_failed = .true.
          end if
        else
          write(io_log,130) "integrals will not be multiplied by -1"
          if (neg_one) then
            write(io_log,999) "neg_one from geom_part_two_param is wrong"
            test_failed = .true.
          end if
        end if
        ! needs to scatter integrals
        if (ref_scatter_deriv) then
          if (scatter_deriv) then
            if (any(seq_part_geo(1:num_cents)/=ref_seq_part)) then
              select case(num_cents)
              case(1)
                write(io_log,100) "correct place of final differentiated center", &
                                  ref_seq_part
                write(io_log,996) "place of final differentiated center is wrong", &
                                  seq_part_geo(1:num_cents)
              case(2)
                write(io_log,110) "correct places of final differentiated centers", &
                                  ref_seq_part
                write(io_log,997) "places of final differentiated centers are wrong", &
                                  seq_part_geo(1:num_cents)
              case(3)
                write(io_log,120) "correct places of final differentiated centers", &
                                  ref_seq_part
                write(io_log,998) "places of final differentiated centers are wrong", &
                                  seq_part_geo(1:num_cents)
              case default
                write(io_log,996) "invalid number of total geometric derivatives", num_cents
              end select
              test_failed = .true.
            else
              select case(num_cents)
              case(1)
                write(io_log,100) "place of final differentiated center", &
                                  seq_part_geo(1:num_cents)
              case(2)
                write(io_log,110) "places of final differentiated centers", &
                                  seq_part_geo(1:num_cents)
              case(3)
                write(io_log,120) "places of final differentiated centers", &
                                  seq_part_geo(1:num_cents)
              case default
                write(io_log,996) "invalid number of total geometric derivatives", num_cents
                test_failed = .true.
              end select
            end if
            ! tests \fn(geom_part_two_scatter)
            if (.not.test_failed)                                             &
              call test_geom_part_two_scatter(num_cents, order_cent,          &
                                              seq_part_geo(1:num_cents),      &
                                              order_part_cent, order_geo_bra, &
                                              order_geo_ket, order_geo_opt,   &
                                              io_log, test_failed)
          else
            write(io_log,999) "scatter_deriv from geom_part_two_param is wrong"
            test_failed = .true.
          end if
        else
          write(io_log,130) "scattering integrals is not needed"
          if (scatter_deriv) then
            write(io_log,999) "scatter_deriv from geom_part_two_param is wrong"
            test_failed = .true.
          end if
        end if
      end if
    end if
    return
100 format("<em>",A,"</em>: ",I3,"<br>")
110 format("<em>",A,"</em>: ",I3,",",I3,"<br>")
120 format("<em>",A,"</em>: ",2(I3,","),I3,"<br>")
130 format("<em>",A,"</em><br>")
996 format("<font color=""red"">",A," (",I3,")!</font><br>")
997 format("<font color=""red"">",A," (",I3,","I3,")!</font><br>")
998 format("<font color=""red"">",A," (",2(I3,","),I3,")!</font><br>")
999 format("<font color=""red"">",A,"!</font><br>")
  end subroutine test_geom_part_two_param

  !> \brief tests subroutine geom_part_two_scatter
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param num_cents is the number of differentiated centers, could be 1, 2 or 3
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param seq_part_geo contains the sequence of bra, ket and operator centers for
  !>        partial geometric derivatives, from \fn(geom_part_two_param)
  !> \param order_part_cent contains the final orders of partial geometric derivatives by adding
  !>        those from total geometric derivatives
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param order_geo_opt is the order of partial geometric derivatives on operator center
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_geom_part_two_scatter(num_cents, order_cent, seq_part_geo,  &
                                        order_part_cent, order_geo_bra,       &
                                        order_geo_ket, order_geo_opt, io_log, &
                                        test_failed)
    use xkind
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: seq_part_geo(num_cents)
    integer, intent(in) :: order_part_cent(3)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: order_geo_opt
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! threshold of error
#include "err_thrsh.h"
    ! see \fn(geom_part_two_scatter) for the details of the following variables
    integer, parameter :: DIM_CINTS = 21
    integer dim_geo_bra
    integer dim_geo_ket
    integer dim_geo_opt
    integer, parameter :: NUM_OPT = 21
    real(REALK), allocatable :: part_cints(:,:,:,:,:)
    integer num_geo_bra
    integer num_geo_ket
    integer num_geo_opt
    integer num_tot_geo
    real(REALK), allocatable :: total_cints(:,:,:,:,:,:)
    real(REALK), allocatable :: ref_cints(:,:,:,:,:,:)
                                                 !referenced total geometric derivatives
    integer num_tot_bra                          !number of total geometric derivatives on bra center
    integer num_tot_ket                          !number of total geometric derivatives on ket center
    integer num_tot_opt                          !number of total geometric derivatives on operator center
    integer num_tot_cent2                        !number of total geometric derivatives on two centers
    real(REALK), allocatable :: scatt1_cints(:,:,:,:,:,:)
    real(REALK), allocatable :: scatt2_cints(:,:,:,:,:,:)
                                                 !integrals after scattering derivatives
    integer icint, ibra, iket, iopt, itot, nopt  !incremental recorders
    integer ierr                                 !error information
    real(REALK) curr_time                        !current CPU time
    ! dumps information of tests
    select case(num_cents)
    ! one-center total geometric derivatives
    case(1)
      write(io_log,100) "tests subroutine geom_part_two_scatter "// &
                        "with one-center total geometric derivatives"
      write(io_log,110) "order of one-center total geometric derivatives", order_cent
      write(io_log,110) "place of one-center total geometric derivative center", seq_part_geo
    ! two-center total geometric derivatives
    case(2)
      write(io_log,100) "tests subroutine geom_part_two_scatter "// &
                        "with two-center total geometric derivatives"
      write(io_log,120) "orders of two-center total geometric derivatives", order_cent
      write(io_log,120) "places of two-center total geometric derivative centers", seq_part_geo
    ! three-center total geometric derivatives
    case(3)
      write(io_log,100) "tests subroutine geom_part_two_scatter "// &
                        "with two-center total geometric derivatives"
      write(io_log,130) "orders of two-center total geometric derivatives", order_cent
      write(io_log,130) "places of two-center total geometric derivative centers", seq_part_geo
    ! error
    case default
      write(io_log,999) "invalid number of total geometric derivatives", num_cents
      test_failed = .true.
      return
    end select
    write(io_log,130) "final orders of partial geometric derivatives", order_part_cent
    write(io_log,110) "order of partial geometric derivatives on bra center", order_geo_bra
    write(io_log,110) "order of partial geometric derivatives on ket center", order_geo_ket
    write(io_log,110) "order of partial geometric derivatives on operator center", order_geo_opt
    ! assigns the faked partial geometric derivatives
    dim_geo_bra = (order_part_cent(1)+1)*(order_part_cent(1)+2)/2
    dim_geo_ket = (order_part_cent(2)+1)*(order_part_cent(2)+2)/2
    dim_geo_opt = (order_part_cent(3)+1)*(order_part_cent(3)+2)/2
    allocate(part_cints(DIM_CINTS,dim_geo_bra,dim_geo_ket,dim_geo_opt,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      write(io_log,999) "failed to allocate part_cints", ierr
      test_failed = .true.
      return
    end if
    do nopt = 1, NUM_OPT
      do iopt = 1, dim_geo_opt
        do iket = 1, dim_geo_ket
          do ibra = 1, dim_geo_bra
            do icint = 1, DIM_CINTS
              part_cints(icint,ibra,iket,iopt,nopt)                        &
                = real(icint,REALK)+real(ibra,REALK)*0.1_REALK             &
                + real(iket,REALK)*0.01_REALK+real(iopt,REALK)*0.001_REALK &
                + real(nopt,REALK)*0.0001_REALK
            end do
          end do
        end do
      end do
    end do
    ! allocates memory for total geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
    num_geo_opt = (order_geo_opt+1)*(order_geo_opt+2)/2
    num_tot_geo = 1
    do itot = 1, num_cents
      num_tot_geo = num_tot_geo*(order_cent(itot)+1)*(order_cent(itot)+2)/2
    end do
    allocate(total_cints(DIM_CINTS,num_geo_bra,num_geo_ket,num_geo_opt, &
                         num_tot_geo,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      write(io_log,999) "failed to allocate total_cints", ierr
      deallocate(part_cints)
      test_failed = .true.
      return
    end if
    ! allocates memory for referenced total geometric derivatives
    allocate(ref_cints(DIM_CINTS,num_geo_bra,num_geo_ket,num_geo_opt, &
                       num_tot_geo,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      write(io_log,999) "failed to allocate ref_cints", ierr
      deallocate(part_cints)
      deallocate(total_cints)
      test_failed = .true.
      return
    end if
    ! resets the time
    call xtimer_set(curr_time)
    ! gets the total and partial geometric derivatives using \fn(geom_part_two_scatter)
    call geom_part_two_scatter(num_cents, order_cent, seq_part_geo,         &
                               order_geo_bra, order_geo_ket, order_geo_opt, &
                               DIM_CINTS, dim_geo_bra, dim_geo_ket,         &
                               dim_geo_opt, NUM_OPT, part_cints,            &
                               num_geo_bra, num_geo_ket, num_geo_opt,       &
                               num_tot_geo, total_cints)
    ! prints the CPU elapsed time
    write(io_log,"(A)") "<font color=""blue""><em>time used for geom_part_two_scatter</em></font><br>"
    call xtimer_print(curr_time, "geom_part_two_scatter", io_log)
    write(io_log,"(A)") "<br>"
    ! resets the time
    call xtimer_set(curr_time)
    ! gets the total and partial geometric derivatives by calling \fn(shell_scatter)
    select case(num_cents)
    ! one-center total geometric derivatives
    case(1)
      ! scatters on bra center
      if (seq_part_geo(1)==TAG_BRA) then
        call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*dim_geo_opt,      &
                           NUM_OPT, part_cints, (/order_geo_bra,order_geo_bra/), &
                           order_cent(1), num_geo_bra, num_tot_geo, ref_cints)
      ! scatters on ket center
      else if (seq_part_geo(1)==TAG_KET) then
        call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, dim_geo_opt,      &
                           NUM_OPT, part_cints, (/order_geo_ket,order_geo_ket/), &
                           order_cent(1), num_geo_ket, num_tot_geo, ref_cints)
      ! scatters on operator center
      else
        call shell_scatter(DIM_CINTS*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1,    &
                           NUM_OPT, part_cints, (/order_geo_opt,order_geo_opt/), &
                           order_cent(1), num_geo_opt, num_tot_geo, ref_cints)
      end if
    ! two-center total geometric derivatives
    case(2)
      ! first differentiated center is bra center
      if (seq_part_geo(1)==TAG_BRA) then
        ! scatters on bra centers
        if (seq_part_geo(2)==TAG_BRA) then
          itot = order_geo_bra+order_cent(1)
          num_tot_ket = (itot+1)*(itot+2)/2
          num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on bra center
          allocate(scatt1_cints(DIM_CINTS,num_tot_ket,dim_geo_ket,dim_geo_opt, &
                                num_tot_opt,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on bra", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on bra center
          call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*dim_geo_opt,   &
                             NUM_OPT, part_cints, (/itot,itot/), order_cent(2), &
                             num_tot_ket, num_tot_opt, scatt1_cints)
          ! scatters on bra center
          num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS, num_tot_ket, dim_geo_ket*dim_geo_opt, &
                             num_tot_opt*NUM_OPT, scatt1_cints,               &
                             (/order_geo_bra,order_geo_bra/), order_cent(1),  &
                             num_geo_bra, num_tot_bra, ref_cints)
        ! scatters on ket center first
        else if (seq_part_geo(2)==TAG_KET) then
          num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on ket center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,num_geo_ket,dim_geo_opt, &
                                num_tot_ket,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on ket", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on ket center
          call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, dim_geo_opt,      &
                             NUM_OPT, part_cints, (/order_geo_ket,order_geo_ket/), &
                             order_cent(2), num_geo_ket, num_tot_ket, scatt1_cints)
          ! scatters on bra center
          num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS, dim_geo_bra, num_geo_ket*dim_geo_opt, &
                             num_tot_ket*NUM_OPT, scatt1_cints,               &
                             (/order_geo_bra,order_geo_bra/), order_cent(1),  &
                             num_geo_bra, num_tot_bra, ref_cints)
        ! scatters on operator center first
        else
          num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on operator center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,dim_geo_ket,num_geo_opt, &
                                num_tot_opt,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on operator", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on operator center
          call shell_scatter(DIM_CINTS*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1,    &
                             NUM_OPT, part_cints, (/order_geo_opt,order_geo_opt/), &
                             order_cent(2), num_geo_opt, num_tot_opt, scatt1_cints)
          ! scatters on bra center
          num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*num_geo_opt, &
                             num_tot_opt*NUM_OPT, scatt1_cints,               &
                             (/order_geo_bra,order_geo_bra/), order_cent(1),  &
                             num_geo_bra, num_tot_bra, ref_cints)
        end if
      ! first differentiated center is ket center
      else if (seq_part_geo(1)==TAG_KET) then
        ! scatters on bra center first
        if (seq_part_geo(2)==TAG_BRA) then
          num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on bra center
          allocate(scatt1_cints(DIM_CINTS,num_geo_bra,dim_geo_ket,dim_geo_opt, &
                                num_tot_bra,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on bra", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on bra center
          call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*dim_geo_opt,       &
                             NUM_OPT, part_cints, (/order_geo_bra,order_geo_bra/),  &
                             order_cent(2), num_geo_bra, num_tot_bra, scatt1_cints)
          ! scatters on ket center
          num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS*num_geo_bra, dim_geo_ket, dim_geo_opt, &
                             num_tot_bra*NUM_OPT, scatt1_cints,               &
                             (/order_geo_ket,order_geo_ket/), order_cent(1),  &
                             num_geo_ket, num_tot_ket, ref_cints)
        ! scatters on ket center
        else if (seq_part_geo(2)==TAG_KET) then
          itot = order_geo_ket+order_cent(1)
          num_tot_ket = (itot+1)*(itot+2)/2
          num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on ket center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,num_tot_ket,dim_geo_opt, &
                                num_tot_opt,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on ket", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on ket center
          call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, dim_geo_opt,   &
                             NUM_OPT, part_cints, (/itot,itot/), order_cent(2), &
                             num_tot_ket, num_tot_opt, scatt1_cints)
          ! scatters on ket center
          num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS*dim_geo_bra, num_tot_ket, dim_geo_opt, &
                             num_tot_opt*NUM_OPT, scatt1_cints,               &
                             (/order_geo_ket,order_geo_ket/), order_cent(1),  &
                             num_geo_ket, num_tot_bra, ref_cints)
        ! scatters on operator center first
        else
          num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on operator center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,dim_geo_ket,num_geo_opt, &
                                num_tot_opt,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on operator", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on operator center
          call shell_scatter(DIM_CINTS*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1,     &
                             NUM_OPT, part_cints, (/order_geo_opt,order_geo_opt/),  &
                             order_cent(2), num_geo_opt, num_tot_opt, scatt1_cints)
          ! scatters on ket center
          num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, num_geo_opt, &
                             num_tot_opt*NUM_OPT, scatt1_cints,               &
                             (/order_geo_ket,order_geo_ket/), order_cent(1),  &
                             num_geo_ket, num_tot_ket, ref_cints)
        end if
      ! first differentiated center is operator center
      else
        ! scatters on bra center first
        if (seq_part_geo(2)==TAG_BRA) then
          num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on bra center
          allocate(scatt1_cints(DIM_CINTS,num_geo_bra,dim_geo_ket,dim_geo_opt, &
                                num_tot_bra,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on bra", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on bra center
          call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*dim_geo_opt,       &
                             NUM_OPT, part_cints, (/order_geo_bra,order_geo_bra/),  &
                             order_cent(2), num_geo_bra, num_tot_bra, scatt1_cints)
          ! scatters on operator center
          num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS*num_geo_bra*dim_geo_ket, dim_geo_opt,  &
                             1, num_tot_bra*NUM_OPT, scatt1_cints,            &
                             (/order_geo_opt,order_geo_opt/), order_cent(1),  &
                             num_geo_opt, num_tot_opt, ref_cints)
        ! scatters on ket center first
        else if (seq_part_geo(2)==TAG_KET) then
          num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on ket center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,num_geo_ket,dim_geo_opt, &
                                num_tot_ket,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on ket", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on ket center
          call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, dim_geo_opt,      &
                             NUM_OPT, part_cints, (/order_geo_ket,order_geo_ket/), &
                             order_cent(2), num_geo_ket, num_tot_ket, scatt1_cints)
          ! scatters on operator center
          num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS*dim_geo_bra*num_geo_ket, dim_geo_opt,  &
                             1, num_tot_ket*NUM_OPT, scatt1_cints,            &
                             (/order_geo_opt,order_geo_opt/), order_cent(1),  &
                             num_geo_opt, num_tot_opt, ref_cints)
        ! scatters on operator centers
        else
          itot = order_geo_opt+order_cent(1)
          num_tot_ket = (itot+1)*(itot+2)/2
          num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
          ! allocates memory for integrals after scattering derivatives on operator center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,dim_geo_ket,num_tot_ket, &
                                num_tot_opt,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on operator", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on operator center
          call shell_scatter(DIM_CINTS*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1, &
                             NUM_OPT, part_cints, (/itot,itot/), order_cent(2), &
                             num_tot_ket, num_tot_opt, scatt1_cints)
          ! scatters on operator center
          num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
          call shell_scatter(DIM_CINTS*dim_geo_bra*dim_geo_ket, num_tot_ket, &
                             1, num_tot_opt*NUM_OPT, scatt1_cints,           &
                             (/order_geo_opt,order_geo_opt/), order_cent(1), &
                             num_geo_opt, num_tot_bra, ref_cints)
        end if
      end if
      ! cleans
      deallocate(scatt1_cints)
    ! three-center total geometric derivatives
    case(3)
      ! first differentiated center is bra center
      if (seq_part_geo(1)==TAG_BRA) then
        ! scatters on operator center first
        if (seq_part_geo(2)==TAG_KET) then
          num_tot_opt = (order_cent(3)+1)*(order_cent(3)+2)/2
          ! allocates memory for integrals after scattering derivatives on operator center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,dim_geo_ket,num_geo_opt, &
                                num_tot_opt,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on operator", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on operator center
          call shell_scatter(DIM_CINTS*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1,    &
                             NUM_OPT, part_cints, (/order_geo_opt,order_geo_opt/), &
                             order_cent(3), num_geo_opt, num_tot_opt, scatt1_cints)
          num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
          num_tot_cent2 = num_tot_ket*num_tot_opt
          ! allocates memory for integrals after scattering derivatives on ket center
          allocate(scatt2_cints(DIM_CINTS,dim_geo_bra,num_geo_ket,num_geo_opt, &
                                num_tot_cent2,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt2_cints on ket", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            deallocate(scatt1_cints)
            test_failed = .true.
            return
          end if
          ! scatters on ket center
          call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, num_geo_opt, &
                             num_tot_opt*NUM_OPT, scatt1_cints,               &
                             (/order_geo_ket,order_geo_ket/), order_cent(2),  &
                             num_geo_ket, num_tot_ket, scatt2_cints)
        ! scatters on ket center first
        else
          num_tot_ket = (order_cent(3)+1)*(order_cent(3)+2)/2
          ! allocates memory for integrals after scattering derivatives on ket center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,num_geo_ket,dim_geo_opt, &
                                num_tot_ket,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on ket", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on ket center
          call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, dim_geo_opt,      &
                             NUM_OPT, part_cints, (/order_geo_ket,order_geo_ket/), &
                             order_cent(3), num_geo_ket, num_tot_ket, scatt1_cints)
          num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
          num_tot_cent2 = num_tot_opt*num_tot_ket
          ! allocates memory for integrals after scattering derivatives on operator center
          allocate(scatt2_cints(DIM_CINTS,dim_geo_bra,num_geo_ket,num_geo_opt, &
                                num_tot_cent2,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt2_cints on operator", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            deallocate(scatt1_cints)
            test_failed = .true.
            return
          end if
          ! scatters on operator center
          call shell_scatter(DIM_CINTS*dim_geo_bra*num_geo_ket, dim_geo_opt, &
                             1, num_tot_ket*NUM_OPT, scatt1_cints,           &
                             (/order_geo_opt,order_geo_opt/), order_cent(2), &
                             num_geo_opt, num_tot_opt, scatt2_cints)
        end if
        ! cleans
        deallocate(scatt1_cints)
        ! scatters on bra center
        num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
        call shell_scatter(DIM_CINTS, dim_geo_bra, num_geo_ket*num_geo_opt, &
                           num_tot_cent2*NUM_OPT, scatt2_cints,             &
                           (/order_geo_bra,order_geo_bra/), order_cent(1),  &
                           num_geo_bra, num_tot_bra, ref_cints)
      ! first differentiated center is ket center
      else if (seq_part_geo(1)==TAG_KET) then
        ! scatters on operator center first
        if (seq_part_geo(2)==TAG_BRA) then
          num_tot_opt = (order_cent(3)+1)*(order_cent(3)+2)/2
          ! allocates memory for integrals after scattering derivatives on operator center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,dim_geo_ket,num_geo_opt, &
                                num_tot_opt,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on opreator", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on operator center
          call shell_scatter(DIM_CINTS*dim_geo_bra*dim_geo_ket, dim_geo_opt, 1,    &
                             NUM_OPT, part_cints, (/order_geo_opt,order_geo_opt/), &
                             order_cent(3), num_geo_opt, num_tot_opt, scatt1_cints)
          num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
          num_tot_cent2 = num_tot_bra*num_tot_opt
          ! allocates memory for integrals after scattering derivatives on bra center
          allocate(scatt2_cints(DIM_CINTS,num_geo_bra,dim_geo_ket,num_geo_opt, &
                                num_tot_cent2,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt2_cints on bra", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            deallocate(scatt1_cints)
            test_failed = .true.
            return
          end if
          ! scatters on bra center
          call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*num_geo_opt, &
                             num_tot_opt*NUM_OPT, scatt1_cints,               &
                             (/order_geo_bra,order_geo_bra/), order_cent(2),  &
                             num_geo_bra, num_tot_bra, scatt2_cints)
        ! scatters on bra center first
        else
          num_tot_bra = (order_cent(3)+1)*(order_cent(3)+2)/2
          ! allocates memory for integrals after scattering derivatives on bra center
          allocate(scatt1_cints(DIM_CINTS,num_geo_bra,dim_geo_ket,dim_geo_opt, &
                                num_tot_bra,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on bra", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on bra center
          call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*dim_geo_opt,      &
                             NUM_OPT, part_cints, (/order_geo_bra,order_geo_bra/), &
                             order_cent(3), num_geo_bra, num_tot_bra, scatt1_cints)
          num_tot_opt = (order_cent(2)+1)*(order_cent(2)+2)/2
          num_tot_cent2 = num_tot_opt*num_tot_bra
          ! allocates memory for integrals after scattering derivatives on operator center
          allocate(scatt2_cints(DIM_CINTS,num_geo_bra,dim_geo_ket,num_geo_opt, &
                                num_tot_cent2,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt2_cints on operator", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            deallocate(scatt1_cints)
            test_failed = .true.
            return
          end if
          ! scatters on operator center
          call shell_scatter(DIM_CINTS*num_geo_bra*dim_geo_ket, dim_geo_opt,  &
                             1, num_tot_bra*NUM_OPT, scatt1_cints,            &
                             (/order_geo_opt,order_geo_opt/), order_cent(2),  &
                             num_geo_opt, num_tot_opt, scatt2_cints)
        end if
        ! cleans
        deallocate(scatt1_cints)
        ! scatters on ket center
        num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
        call shell_scatter(DIM_CINTS*num_geo_bra, dim_geo_ket, num_geo_opt, &
                           num_tot_cent2*NUM_OPT, scatt2_cints,             &
                           (/order_geo_ket,order_geo_ket/), order_cent(1),  &
                           num_geo_ket, num_tot_ket, ref_cints)
      ! first differentiated center is operator center
      else
        ! scatters on ket center first
        if (seq_part_geo(2)==TAG_BRA) then
          num_tot_ket = (order_cent(3)+1)*(order_cent(3)+2)/2
          ! allocates memory for integrals after scattering derivatives on ket center
          allocate(scatt1_cints(DIM_CINTS,dim_geo_bra,num_geo_ket,dim_geo_opt, &
                                num_tot_ket,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on ket", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on ket center
          call shell_scatter(DIM_CINTS*dim_geo_bra, dim_geo_ket, dim_geo_opt,      &
                             NUM_OPT, part_cints, (/order_geo_ket,order_geo_ket/), &
                             order_cent(3), num_geo_ket, num_tot_ket, scatt1_cints)
          num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
          num_tot_cent2 = num_tot_bra*num_tot_ket
          ! allocates memory for integrals after scattering derivatives on bra center
          allocate(scatt2_cints(DIM_CINTS,num_geo_bra,num_geo_ket,dim_geo_opt, &
                                num_tot_cent2,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt2_cints on bra", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            deallocate(scatt1_cints)
            test_failed = .true.
            return
          end if
          ! scatters on bra center
          call shell_scatter(DIM_CINTS, dim_geo_bra, num_geo_ket*dim_geo_opt, &
                             num_tot_ket*NUM_OPT, scatt1_cints,               &
                             (/order_geo_bra,order_geo_bra/), order_cent(2),  &
                             num_geo_bra, num_tot_bra, scatt2_cints)
        ! scatters on bra center first
        else
          num_tot_bra = (order_cent(3)+1)*(order_cent(3)+2)/2
          ! allocates memory for integrals after scattering derivatives on bra center
          allocate(scatt1_cints(DIM_CINTS,num_geo_bra,dim_geo_ket,dim_geo_opt, &
                                num_tot_bra,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt1_cints on bra", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            test_failed = .true.
            return
          end if
          ! scatters on bra center
          call shell_scatter(DIM_CINTS, dim_geo_bra, dim_geo_ket*dim_geo_opt,      &
                             NUM_OPT, part_cints, (/order_geo_bra,order_geo_bra/), &
                             order_cent(3), num_geo_bra, num_tot_bra, scatt1_cints)
          num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
          num_tot_cent2 = num_tot_ket*num_tot_bra
          ! allocates memory for integrals after scattering derivatives on ket center
          allocate(scatt2_cints(DIM_CINTS,num_geo_bra,num_geo_ket,dim_geo_opt, &
                                num_tot_cent2,NUM_OPT), stat=ierr)
          if (ierr/=0) then
            write(io_log,999) "failed to allocate scatt2_cints on ket", ierr
            deallocate(part_cints)
            deallocate(total_cints)
            deallocate(ref_cints)
            deallocate(scatt1_cints)
            test_failed = .true.
            return
          end if
          ! scatters on ket center
          call shell_scatter(DIM_CINTS*num_geo_bra, dim_geo_ket, dim_geo_opt, &
                             num_tot_bra*NUM_OPT, scatt1_cints,               &
                             (/order_geo_ket,order_geo_ket/), order_cent(2),  &
                             num_geo_ket, num_tot_ket, scatt2_cints)
        end if
        ! cleans
        deallocate(scatt1_cints)
        ! scatters on operator center
        num_tot_opt = (order_cent(1)+1)*(order_cent(1)+2)/2
        call shell_scatter(DIM_CINTS*num_geo_bra*num_geo_ket, dim_geo_opt, &
                           1, num_tot_cent2*NUM_OPT, scatt2_cints,         &
                           (/order_geo_opt,order_geo_opt/), order_cent(1), &
                           num_geo_opt, num_tot_opt, ref_cints)
      end if
      ! cleans
      deallocate(scatt2_cints)
    end select
    ! prints the CPU elapsed time
    write(io_log,"(A)") "<font color=""blue""><em>time used for calling shell_scatter</em></font><br>"
    call xtimer_print(curr_time, "calling shell_scatter", io_log)
    write(io_log,"(A)") "<br>"
    ! cleans
    deallocate(part_cints)
    ! checks the results
    do nopt = 1, NUM_OPT
      do itot = 1, num_tot_geo
        do iopt = 1, num_geo_opt
          do iket = 1, num_geo_ket
            do ibra = 1, num_geo_bra
              do icint = 1, DIM_CINTS
                if (abs(total_cints(icint,ibra,iket,iopt,itot,nopt) &
                        -ref_cints(icint,ibra,iket,iopt,itot,nopt))>ERR_THRSH) then
                  write(io_log,140) "correct total geometric derivatives", &
                                    icint, ibra, iket, iopt, itot, nopt,   &
                                    ref_cints(icint,ibra,iket,iopt,itot,nopt)
                  write(io_log,998) "result from geom_part_two_scatter is wrong", &
                                    total_cints(icint,ibra,iket,iopt,itot,nopt)
                  test_failed = .true.
                end if
              end do
            end do
          end do
        end do
      end do
    end do
    ! cleans
    deallocate(total_cints)
    deallocate(ref_cints)
    return
100 format("<em>",A,"</em><br>")
110 format("<em>",A,"</em>: ",I3,"<br>")
120 format("<em>",A,"</em>: ",I3,",",I3,"<br>")
130 format("<em>",A,"</em>: ",2(I3,","),I3,"<br>")
140 format("<em>",A," (",5(I3,","),I3,"): ",F16.6,"<br>")
998 format("<font color=""red"">",A," (",F16.6,")!</font><br>")
999 format("<font color=""red"">",A," (",I8,")!</font><br>")
  end subroutine test_geom_part_two_scatter
