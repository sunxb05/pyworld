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
!!  Tests subroutines in geom_part_zero.F90.
!!
!!  2011-08-02, Bin Gao
!!  * first version

#include "private/tag_cent.h"

  !> \brief tests subroutines in geom_part_zero.F90
  !> \author Bin Gao
  !> \date 2011-08-02
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_geom_part_zero(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters and referenced results
#include "geom_part_zero.h"
    integer itst, jtst, ktst  !incremental recorders over test cases
    ! we first test shell_scatter.F90 which will be used for testing \fn(geom_part_zero_scatter)
    call html_log_text("We first test subroutines in shell_scatter.F90 ...", io_log)
    call test_shell_scatter(io_log, test_failed)
    if (test_failed) then
      call html_log_text(log_text="Tests subroutines in shell_scatter.F90 failed! "// &
                                  "No further test will be performed!",               &
                         io_log=io_log, font_color="red")
      return
    else
      call html_log_text("Tests subroutines in shell_scatter.F90 passed!", io_log)
    end if
    ! loops over test cases
    jtst = 1
    do itst = 1, NUM_TOT_GEO
      ! dumps information of tests
      call html_log_heading(heading="Tests subroutine geom_part_zero_param", &
                            io_log=io_log, level=4)
      call html_log_int_number("Number of differentiated centers:", &
                               NUM_CENTS(itst), "I4", io_log)
      call html_log_int_array("Indices of bra and ket centers:", &
                               IDX_BRA_KET(jtst:jtst+1), "I4", io_log)
      ! tests \fn(geom_part_zero_param)
      ktst = jtst+NUM_CENTS(itst)-1
      call test_geom_part_zero_param(NUM_CENTS(itst), IDX_GEO_CENT(jtst:ktst),    &
                                     ORDER_GEO_CENT(jtst:ktst),                   &
                                     IDX_BRA_KET(jtst:jtst+1),                    &
                                     ORDER_GEO_BRA(itst), ORDER_GEO_KET(itst),    &
                                     REF_ORDER_PART(jtst:jtst+1),                 &
                                     REF_ZERO_INTS(itst), REF_SCATTER_DER(itst),  &
                                     REF_SEQ_PART(jtst:ktst), io_log, test_failed)
      jtst = jtst+2
    end do
    return
  end subroutine test_geom_part_zero

  !> \brief tests subroutine geom_part_zero_param
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param num_cents is the number of differentiated centers, could be 1 or 2
  !> \param idx_cent contains the indices of different differentiated centers
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param idx_part_cent contains the indices of bra and ket centers
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param ref_order_part contains referenced final orders of partial geometric derivatives
  !>        by adding those from total geometric derivatives
  !> \param ref_zero_ints indicates if the total geometric derivatives are zero
  !> \param ref_scatter_der indicates if scattering the geometric derivatives later on
  !> \param ref_seq_part contains the sequence of bra and ket centers for partial derivative terms
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_geom_part_zero_param(num_cents, idx_cent, order_cent, idx_part_cent, &
                                       order_geo_bra, order_geo_ket, ref_order_part,   &
                                       ref_zero_ints, ref_scatter_der, ref_seq_part,   &
                                       io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: idx_cent(num_cents)
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: idx_part_cent(2)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: ref_order_part(2)
    logical, intent(in) :: ref_zero_ints
    logical, intent(in) :: ref_scatter_der
    integer, intent(in) :: ref_seq_part(num_cents)
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
    ! see \fn(geom_part_zero_param) for the details of the following variables
    integer order_part_cent(2)
    logical zero_ints
    logical scatter_deriv
    integer seq_part_geo(2)
    ! dumps information of tests
    select case(num_cents)
    case(1,2)
      call html_log_int_array("Indices of differentiated centers:", &
                              idx_cent, "I4", io_log)
      call html_log_int_array("Orders of differentiated centers:", &
                              order_cent, "I4", io_log)
    case default
      call html_log_int_number("Zero derivatives due to too many differentiated centers:", &
                               num_cents, "I4", io_log)
      return
    end select
    ! sets the orders of partial geometric derivatives
    order_part_cent = (/order_geo_bra,order_geo_ket/)
    call html_log_int_array("Initial orders of partial geometric derivatives:", &
                            order_part_cent, "I4", io_log)
    ! tests \fn(geom_part_zero_param)
    call geom_part_zero_param(num_cents, idx_cent, order_cent, idx_part_cent, &
                              order_part_cent, zero_ints, scatter_deriv, seq_part_geo)
    ! zero integrals
    if (ref_zero_ints) then
      if (zero_ints) then
        call html_log_text("Total geometric derivatives are zero.", io_log)
      else
        call html_log_text("Total geometric derivatives shoud be zero.", io_log)
        call html_log_text(log_text="zero_ints from geom_part_zero_param is wrong!", &
                           io_log=io_log, font_color="red")
        test_failed = .true.
      end if
    ! non-zero integrals
    else
      if (zero_ints) then
        call html_log_text("Total geometric derivatives shoud not be zero.", io_log)
        call html_log_text(log_text="zero_ints from geom_part_zero_param is wrong!", &
                           io_log=io_log, font_color="red")
        test_failed = .true.
      else
        ! checks the resulted orders of partial geometric derivatives
        if (any(order_part_cent/=ref_order_part)) then
          call html_log_int_array(                                                     &
                 log_text="Referenced final orders of partial geometric derivatives:", &
                 log_int=ref_order_part, fmt_int="I3", io_log=io_log, font_color="blue")
          call html_log_int_array(                                        &
                 log_text="Results from geom_part_zero_param are wrong:", &
                 log_int=order_part_cent, fmt_int="I3", io_log=io_log, font_color="red")
          test_failed = .true.
        else
          call html_log_int_array(                                          &
                 log_text="Final orders of partial geometric derivatives:", &
                 log_int=order_part_cent, fmt_int="I3", io_log=io_log)
        end if
        ! needs to scatter integrals
        if (ref_scatter_der) then
          if (scatter_deriv) then
            if (any(seq_part_geo(1:num_cents)/=ref_seq_part)) then
              select case(num_cents)
              case(1,2)
                call html_log_int_array(                                             &
                       log_text="Referenced places of final differentiated center:", &
                       log_int=ref_seq_part, fmt_int="I4", io_log=io_log,            &
                       font_color="blue")
                call html_log_int_array(                                             &
                       log_text="Places of final differentiated centers are wrong:", &
                       log_int=seq_part_geo(1:num_cents), fmt_int="I4",              &
                       io_log=io_log, font_color="red")
              case default
                call html_log_int_number(                                         &
                       log_text="Invalid number of total geometric derivatives:", &
                       log_int=num_cents, fmt_int="I4", io_log=io_log, font_color="red")
              end select
              test_failed = .true.
            else
              select case(num_cents)
              case(1,2)
                call html_log_int_array(                                   &
                       log_text="Places of final differentiated centers:", &
                       log_int=seq_part_geo(1:num_cents), fmt_int="I4",    &
                       io_log=io_log)
              case default
                call html_log_int_number(                                         &
                       log_text="Invalid number of total geometric derivatives:", &
                       log_int=num_cents, fmt_int="I4", io_log=io_log, font_color="red")
                test_failed = .true.
              end select
            end if
            ! tests \fn(geom_part_zero_scatter)
            if (.not.test_failed) then
              call html_log_heading(heading="Tests subroutine geom_part_zero_scatter", &
                                    io_log=io_log, level=4)
              call test_geom_part_zero_scatter(num_cents, order_cent,          &
                                               seq_part_geo(1:num_cents),      &
                                               order_part_cent, order_geo_bra, &
                                               order_geo_ket, io_log, test_failed)
            end if
          else
            call html_log_text(                                                  &
                   log_text="scatter_deriv from geom_part_zero_param is wrong!", &
                   io_log=io_log, font_color="red")
            test_failed = .true.
          end if
        else
          call html_log_text("Scattering integrals is not needed.", io_log)
          if (scatter_deriv) then
            call html_log_text(                                                  &
                   log_text="scatter_deriv from geom_part_zero_param is wrong!", &
                   io_log=io_log, font_color="red")
            test_failed = .true.
          end if
        end if
      end if
    end if
    return
  end subroutine test_geom_part_zero_param

  !> \brief tests subroutine geom_part_zero_scatter
  !> \author Bin Gao
  !> \date 2011-08-03
  !> \param num_cents is the number of differentiated centers, could be 1 or 2
  !> \param order_cent contains the order of total geometric derivatives of differentiated centers
  !> \param seq_part_geo contains the sequence of bra and ket centers for
  !>        partial geometric derivatives, from \fn(geom_part_zero_param)
  !> \param order_part_cent contains the final orders of partial geometric derivatives by adding
  !>        those from total geometric derivatives
  !> \param order_geo_bra is the order of partial geometric derivatives on bra center
  !> \param order_geo_ket is the order of partial geometric derivatives on ket center
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_geom_part_zero_scatter(num_cents, order_cent, seq_part_geo, &
                                         order_part_cent, order_geo_bra,      &
                                         order_geo_ket, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: num_cents
    integer, intent(in) :: order_cent(num_cents)
    integer, intent(in) :: seq_part_geo(num_cents)
    integer, intent(in) :: order_part_cent(2)
    integer, intent(in) :: order_geo_bra
    integer, intent(in) :: order_geo_ket
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
    ! see \fn(geom_part_zero_scatter) for the details of the following variables
    integer, parameter :: DIM_CINTS = 21
    integer dim_geo_bra
    integer dim_geo_ket
    integer, parameter :: NUM_OPT = 21
    real(REALK), allocatable :: part_cints(:,:,:,:)
    integer num_geo_bra
    integer num_geo_ket
    integer num_tot_geo
    real(REALK), allocatable :: total_cints(:,:,:,:,:)
    real(REALK), allocatable :: ref_cints(:,:,:,:,:)
                            !referenced total geometric derivatives
    integer num_tot_bra     !number of total geometric derivatives on bra center
    integer num_tot_ket     !number of total geometric derivatives on ket center
    real(REALK), allocatable :: scatter_cints(:,:,:,:,:)
                            !integrals after scattering derivatives on bra or ket
    integer icint, ibra, iket, itot, nopt
                            !incremental recorders
    integer ierr            !error information
    real(REALK) begin_time  !begin of CPU time
    real(REALK) end_time    !end of CPU time
    logical different       !if result from \fn(geom_part_zero_scatter) is different from reference
    ! dumps information of tests
    select case(num_cents)
    case(1,2)
      call html_log_int_array("Orders of total geometric derivatives:", &
                              order_cent, "I4", io_log)
      call html_log_int_array("Places of total geometric derivative centers:", &
                              seq_part_geo, "I4", io_log)
    case default
      call html_log_int_number(                                         &
             log_text="Invalid number of total geometric derivatives:", &
             log_int=num_cents, fmt_int="I4", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end select
    call html_log_int_array("Final orders of partial geometric derivatives:", &
                            order_part_cent, "I4", io_log)
    call html_log_int_number("Order of partial geometric derivatives on bra center:", &
                             order_geo_bra, "I4", io_log)
    call html_log_int_number("Order of partial geometric derivatives on ket center:", &
                             order_geo_ket, "I4", io_log)
    ! assigns the faked partial geometric derivatives
    dim_geo_bra = (order_part_cent(1)+1)*(order_part_cent(1)+2)/2
    dim_geo_ket = (order_part_cent(2)+1)*(order_part_cent(2)+2)/2
    allocate(part_cints(DIM_CINTS,dim_geo_bra,dim_geo_ket,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate part_ints:",          &
                               log_int=DIM_CINTS*dim_geo_bra*dim_geo_ket*NUM_OPT, &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    do nopt = 1, NUM_OPT
      do iket = 1, dim_geo_ket
        do ibra = 1, dim_geo_bra
          do icint = 1, DIM_CINTS
            part_cints(icint,ibra,iket,nopt)                            &
              = real(icint,REALK)*0.1_REALK+real(ibra,REALK)*0.01_REALK &
              + real(iket,REALK)*0.001_REALK+real(nopt,REALK)*0.0001_REALK
          end do
        end do
      end do
    end do
    ! allocates memory for total geometric derivatives
    num_geo_bra = (order_geo_bra+1)*(order_geo_bra+2)/2
    num_geo_ket = (order_geo_ket+1)*(order_geo_ket+2)/2
    num_tot_geo = 1
    do itot = 1, num_cents
      num_tot_geo = num_tot_geo*(order_cent(itot)+1)*(order_cent(itot)+2)/2
    end do
    allocate(total_cints(DIM_CINTS,num_geo_bra,num_geo_ket,num_tot_geo,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate total_ints:", &
                               log_int=DIM_CINTS*num_geo_bra*num_geo_ket  &
                                       *num_tot_geo*NUM_OPT,              &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(part_cints)
      test_failed = .true.
      return
    end if
    ! allocates memory for referenced total geometric derivatives
    allocate(ref_cints(DIM_CINTS,num_geo_bra,num_geo_ket,num_tot_geo,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate ref_ints:",   &
                               log_int=DIM_CINTS*num_geo_bra*num_geo_ket  &
                                       *num_tot_geo*NUM_OPT,              &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(part_cints)
      deallocate(total_cints)
      test_failed = .true.
      return
    end if
    ! gets the begin time
    call xtimer_set(begin_time)
    ! gets the total and partial geometric derivatives using \fn(geom_part_zero_scatter)
    call geom_part_zero_scatter(num_cents, order_cent, seq_part_geo,     &
                                order_geo_bra, order_geo_ket, DIM_CINTS, &
                                dim_geo_bra, dim_geo_ket, NUM_OPT,       &
                                part_cints, num_geo_bra, num_geo_ket,    &
                                num_tot_geo, total_cints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for geom_part_zero_scatter:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",       &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! gets the total and partial geometric derivatives by calling \fn(scatter_single)
    select case(num_cents)
    ! one-center total geometric derivatives
    case(1)
      ! scatters on bra center
      if (seq_part_geo(1)==TAG_BRA) then
        call scatter_single(DIM_CINTS, dim_geo_bra, dim_geo_ket, NUM_OPT, &
                            part_cints, order_geo_bra, order_cent(1),     &
                            num_geo_bra, num_tot_geo, ref_cints)
      ! scatters on ket center
      else
        call scatter_single(DIM_CINTS*dim_geo_bra, dim_geo_ket, 1, NUM_OPT, &
                            part_cints, order_geo_ket, order_cent(1),       &
                            num_geo_ket, num_tot_geo, ref_cints)
      end if
    ! two-center total geometric derivatives
    case(2)
      ! scatters on ket center first
      if (seq_part_geo(1)==TAG_BRA) then
        num_tot_ket = (order_cent(2)+1)*(order_cent(2)+2)/2
        ! allocates memory for integrals after scattering derivatives on ket center
        allocate(scatter_cints(DIM_CINTS,dim_geo_bra,num_geo_ket,num_tot_ket,NUM_OPT), stat=ierr)
        if (ierr/=0) then
          call html_log_int_number(                                             &
                 log_text="Failed to allocate scatter_cints on ket:",           &
                 log_int=DIM_CINTS*dim_geo_bra*num_geo_ket*num_tot_ket*NUM_OPT, &
                 fmt_int="I12", io_log=io_log, font_color="red")
          deallocate(part_cints)
          deallocate(total_cints)
          deallocate(ref_cints)
          test_failed = .true.
          return
        end if
        ! scatters on ket center
        call scatter_single(DIM_CINTS*dim_geo_bra, dim_geo_ket, 1, NUM_OPT, &
                            part_cints, order_geo_ket, order_cent(2),       &
                            num_geo_ket, num_tot_ket, scatter_cints)
        ! scatters on bra center
        num_tot_bra = (order_cent(1)+1)*(order_cent(1)+2)/2
        call scatter_single(DIM_CINTS, dim_geo_bra, num_geo_ket, num_tot_ket*NUM_OPT, &
                            scatter_cints, order_geo_bra, order_cent(1), num_geo_bra, &
                            num_tot_bra, ref_cints)
        ! cleans
        deallocate(scatter_cints)
      ! scatters on bra center first
      else
        num_tot_bra = (order_cent(2)+1)*(order_cent(2)+2)/2
        ! allocates memory for integrals after scattering derivatives on bra center
        allocate(scatter_cints(DIM_CINTS,num_geo_bra,dim_geo_ket,num_tot_bra,NUM_OPT), stat=ierr)
        if (ierr/=0) then
          call html_log_int_number(                                             &
                 log_text="Failed to allocate scatter_cints on bra:",           &
                 log_int=DIM_CINTS*num_geo_bra*dim_geo_ket*num_tot_bra*NUM_OPT, &
                 fmt_int="I12", io_log=io_log, font_color="red") 
          deallocate(part_cints)
          deallocate(total_cints)
          deallocate(ref_cints)
          test_failed = .true.
          return
        end if
        ! scatters on bra center
        call scatter_single(DIM_CINTS, dim_geo_bra, dim_geo_ket, NUM_OPT, &
                            part_cints, order_geo_bra, order_cent(2),     &
                            num_geo_bra, num_tot_bra, scatter_cints)
        ! scatters on ket center
        num_tot_ket = (order_cent(1)+1)*(order_cent(1)+2)/2
        call scatter_single(DIM_CINTS*num_geo_bra, dim_geo_ket, 1, &
                            num_tot_bra*NUM_OPT, scatter_cints,    &
                            order_geo_ket, order_cent(1),          &
                            num_geo_ket, num_tot_ket, ref_cints)
        ! cleans
        deallocate(scatter_cints)
      end if
    end select
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for calling scatter_single:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",       &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(part_cints)
    ! checks the results
    do nopt = 1, NUM_OPT
      do itot = 1, num_tot_geo
        do iket = 1, num_geo_ket
          do ibra = 1, num_geo_bra
            do icint = 1, DIM_CINTS
              call check_difference(ref_cints(icint,ibra,iket,itot,nopt),   &
                                    total_cints(icint,ibra,iket,itot,nopt), &
                                    different)
              if (different) then
                call html_log_int_array(log_text="ID of total geometric derivatives:", &
                                        log_int=(/icint,ibra,iket,itot,nopt/),         &
                                        fmt_int="I3", io_log=io_log)
                call html_log_real_number(                                 &
                       log_text="Referenced total geometric derivatives:", &
                       log_real=ref_cints(icint,ibra,iket,itot,nopt),      &
                       fmt_real="Es16.8", io_log=io_log, font_color="blue")
                call html_log_real_number(                              &
                       log_text="Result from geom_part_zero_scatter:",  &
                       log_real=total_cints(icint,ibra,iket,itot,nopt), &
                       fmt_real="Es16.8", io_log=io_log, font_color="red")
                test_failed = .true.
              end if
            end do
          end do
        end do
      end do
    end do
    ! cleans
    deallocate(total_cints)
    deallocate(ref_cints)
    return
  end subroutine test_geom_part_zero_scatter
