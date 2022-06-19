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
!!  Tests subroutines in shell_scatter.F90.
!!
!!  2011-07-27, Bin Gao
!!  * first version

  !> \brief tests subroutines in shell_scatter.F90
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_shell_scatter(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
    ! see shell_scater.F90 for the details of the following variables
    integer, parameter :: DIM_INTS = 6
    integer dim_gather
    integer, parameter :: DIM_BETWEEN = 3
    integer, parameter :: DIM_OUTMOST = 3
    real(REALK), allocatable :: gather_ints(:,:,:,:)
    integer :: orders_inner(2) = (/0,4/)
    integer :: orders_outer(2) = (/0,4/)
    integer dim_inner
    integer dim_outer
    real(REALK), allocatable :: scatter_ints(:,:,:,:,:)
    integer orders_gather(2)                 !orders of gathered integrals
    integer iint, igat, ibet, iout           !incremental recorders
    integer xout, yout, zout, xin, yin, zin  !incremental recorders over xyz components
    integer start_gat                        !start address of gathered integrals
    integer addr_gat                         !address of gathered integrals
    integer addr_out, addr_in                !addresses of outer and inner integrals
    integer ierr                             !error information
    real(REALK) begin_time                   !begin of CPU time
    real(REALK) end_time                     !end of CPU time
    logical different                        !if result is different from reference
    ! assumes the test will pass
    test_failed = .false.
    ! (1) dumps information of tests
    call html_log_heading(heading="Tests subroutine scatter_single", &
                          io_log=io_log, level=4)
    call html_log_int_number("Order of inner shell:", orders_inner(2), "I3", io_log)
    call html_log_int_number("Order of outer shell:", orders_outer(2), "I3", io_log)
    ! sets the order of gathered integrals
    orders_gather(2) = orders_inner(2)+orders_outer(2)
    ! allocates memory for the gathered integrals
    dim_gather = (orders_gather(2)+1)*(orders_gather(2)+2)/2
    allocate(gather_ints(DIM_INTS,dim_gather,DIM_BETWEEN,DIM_OUTMOST), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate gather_ints:",          &
                               log_int=DIM_INTS*dim_gather*DIM_BETWEEN*DIM_OUTMOST, &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! allocates memory for the scattered integrals
    dim_inner = (orders_inner(2)+1)*(orders_inner(2)+2)/2
    dim_outer = (orders_outer(2)+1)*(orders_outer(2)+2)/2
    allocate(scatter_ints(DIM_INTS,dim_inner,DIM_BETWEEN,dim_outer,DIM_OUTMOST), stat=ierr)
    if (ierr/=0) then
      deallocate(gather_ints)
      call html_log_int_number(log_text="Failed to allocate scatter_ints:", &
                               log_int=DIM_INTS*dim_inner*DIM_BETWEEN       &
                                       *dim_outer*DIM_OUTMOST,              &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! assigns the gathered integrals
    do iout = 1, DIM_OUTMOST
      do ibet = 1, DIM_BETWEEN
        do igat = 1, dim_gather
          do iint = 1, DIM_INTS
            gather_ints(iint,igat,ibet,iout)                           &
              = real(iint,REALK)*0.1_REALK+real(igat,REALK)*0.01_REALK &
              + real(ibet,REALK)*0.001_REALK+real(iout,REALK)*0.0001_REALK
          end do
        end do
      end do
    end do
    ! gets the begin time
    call xtimer_set(begin_time)
    ! tests \fn(scatter_multi_inner)
    call scatter_single(DIM_INTS, dim_gather, DIM_BETWEEN, DIM_OUTMOST, &
                        gather_ints, orders_inner(2), orders_outer(2),  &
                        dim_inner, dim_outer, scatter_ints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for scatter_single:",   &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! checks the results
    do iout = 1, DIM_OUTMOST
      addr_out = 0
      do zout = 0, orders_outer(2)
        do yout = 0, orders_outer(2)-zout
          xout = orders_outer(2)-(yout+zout)
          addr_out = addr_out+1
          addr_in = 0
          ! initializes the start address of value in the gathered integrals
          start_gat = 0
          do zin = 0, orders_inner(2)
            do yin = 0, orders_inner(2)-zin
              xin = orders_inner(2)-(yin+zin)
              addr_in = addr_in+1
              ! index of x^{l}y^{m}z^{n} with l+m+n=angm is 1+m+(2*angm+3-n)*n/2
              addr_gat = start_gat+1+(yin+yout) &
                       + (2*(orders_outer(2)+orders_inner(2))+3-(zin+zout))*(zin+zout)/2
              do ibet = 1, DIM_BETWEEN
                do iint = 1, DIM_INTS
                  call check_difference(gather_ints(iint,addr_gat,ibet,iout),          &
                                        scatter_ints(iint,addr_in,ibet,addr_out,iout), &
                                        different)
                  if (different) then
                    call html_log_int_array(log_text="Orders of inner shell::",    &
                                            log_int=(/xin,yin,zin/), fmt_int="I3", &
                                            io_log=io_log)
                    call html_log_int_array(log_text="Orders of outer shell::",       &
                                            log_int=(/xout,yout,zout/), fmt_int="I3", &
                                            io_log=io_log)
                    call html_log_real_number(log_text="Referenced result:", &
                           log_real=gather_ints(iint,addr_gat,ibet,iout),    &
                           fmt_real="F12.8", io_log=io_log, font_color="blue")
                    call html_log_real_number(log_text="Result from scatter_single:", &
                           log_real=scatter_ints(iint,addr_in,ibet,addr_out,iout),    &
                           fmt_real="F12.8", io_log=io_log, font_color="red")
                    test_failed = .true.
                  end if
                end do
              end do
            end do
          end do
          ! updates the start address of value in the gathered integrals
          start_gat = start_gat                           &
                    + (orders_outer(2)+orders_inner(2)+1) &
                    * (orders_outer(2)+orders_inner(2)+2)/2
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for explicit address:", &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(gather_ints)
    deallocate(scatter_ints)
    ! (2) dumps information of tests
    call html_log_heading(heading="Tests subroutine scatter_multi_inner", &
                          io_log=io_log, level=4)
    call html_log_int_array(log_text="Orders of inner shells:", &
                            log_int=orders_inner, fmt_int="I3", &
                            io_log=io_log, separation="&rarr;")
    call html_log_int_number("Order of outer shell:", orders_outer(2), "I3", io_log)
    ! sets the orders of gathered integrals
    orders_gather(1) = orders_inner(1)+orders_outer(2)
    orders_gather(2) = orders_inner(2)+orders_outer(2)
    ! allocates memory for the gathered integrals
    if (orders_gather(1)>0) then
      dim_gather = ((orders_gather(2)+1)*(orders_gather(2)+2)*(orders_gather(2)+3) &
                 -  orders_gather(1)*(orders_gather(1)+1)*(orders_gather(1)+2))/6
    else
      dim_gather = (orders_gather(2)+1)*(orders_gather(2)+2)*(orders_gather(2)+3)/6
    end if
    allocate(gather_ints(DIM_INTS,dim_gather,DIM_BETWEEN,DIM_OUTMOST), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate gather_ints:",          &
                               log_int=DIM_INTS*dim_gather*DIM_BETWEEN*DIM_OUTMOST, &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! allocates memory for the scattered integrals
    if (orders_inner(1)>0) then
      dim_inner = ((orders_inner(2)+1)*(orders_inner(2)+2)*(orders_inner(2)+3) &
                -  orders_inner(1)*(orders_inner(1)+1)*(orders_inner(1)+2))/6
    else
      dim_inner = (orders_inner(2)+1)*(orders_inner(2)+2)*(orders_inner(2)+3)/6
    end if
    dim_outer = (orders_outer(2)+1)*(orders_outer(2)+2)/2
    allocate(scatter_ints(DIM_INTS,dim_inner,DIM_BETWEEN,dim_outer,DIM_OUTMOST), stat=ierr)
    if (ierr/=0) then
      deallocate(gather_ints)
      call html_log_int_number(log_text="Failed to allocate scatter_ints:", &
                               log_int=DIM_INTS*dim_inner*DIM_BETWEEN       &
                                       *dim_outer*DIM_OUTMOST,              &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! assigns the gathered integrals
    do iout = 1, DIM_OUTMOST
      do ibet = 1, DIM_BETWEEN
        do igat = 1, dim_gather
          do iint = 1, DIM_INTS
            gather_ints(iint,igat,ibet,iout)                           &
              = real(iint,REALK)*0.1_REALK+real(igat,REALK)*0.01_REALK &
              + real(ibet,REALK)*0.001_REALK+real(iout,REALK)*0.0001_REALK
          end do
        end do
      end do
    end do
    ! gets the begin time
    call xtimer_set(begin_time)
    ! tests \fn(scatter_multi_inner)
    call scatter_multi_inner(DIM_INTS, dim_gather, DIM_BETWEEN, DIM_OUTMOST, &
                             gather_ints, orders_inner, orders_outer(2),     &
                             dim_inner, dim_outer, scatter_ints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for scatter_multi_inner:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",    &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! checks the results
    do iout = 1, DIM_OUTMOST
      addr_out = 0
      do zout = 0, orders_outer(2)
        do yout = 0, orders_outer(2)-zout
          xout = orders_outer(2)-(yout+zout)
          addr_out = addr_out+1
          addr_in = 0
          ! initializes the start address of value in the gathered integrals
          start_gat = 0
          do igat = orders_inner(1), orders_inner(2)
            do zin = 0, igat
              do yin = 0, igat-zin
                xin = igat-(yin+zin)
                addr_in = addr_in+1
                ! index of x^{l}y^{m}z^{n} with l+m+n=angm is 1+m+(2*angm+3-n)*n/2
                addr_gat = start_gat+1+(yin+yout) &
                         + (2*(orders_outer(2)+igat)+3-(zin+zout))*(zin+zout)/2
                do ibet = 1, DIM_BETWEEN
                  do iint = 1, DIM_INTS
                    call check_difference(gather_ints(iint,addr_gat,ibet,iout),          &
                                          scatter_ints(iint,addr_in,ibet,addr_out,iout), &
                                          different)
                    if (different) then
                      call html_log_int_array(log_text="Orders of inner shell::",    &
                                              log_int=(/xin,yin,zin/), fmt_int="I3", &
                                              io_log=io_log)
                      call html_log_int_array(log_text="Orders of outer shell::",       &
                                              log_int=(/xout,yout,zout/), fmt_int="I3", &
                                              io_log=io_log)
                      call html_log_real_number(log_text="Referenced result:", &
                             log_real=gather_ints(iint,addr_gat,ibet,iout),    &
                             fmt_real="F12.8", io_log=io_log, font_color="blue")
                      call html_log_real_number(                                     &
                             log_text="Result from scatter_multi_inner:",            &
                             log_real=scatter_ints(iint,addr_in,ibet,addr_out,iout), &
                             fmt_real="F12.8", io_log=io_log, font_color="red")
                      test_failed = .true.
                    end if
                  end do
                end do
              end do
            end do
            ! updates the start address of value in the gathered integrals
            start_gat = start_gat+(orders_outer(2)+igat+1)*(orders_outer(2)+igat+2)/2
          end do
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for explicit address:", &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(gather_ints)
    deallocate(scatter_ints)
    return
  end subroutine test_shell_scatter
