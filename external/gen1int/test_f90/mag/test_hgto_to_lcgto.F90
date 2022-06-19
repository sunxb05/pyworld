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
!!  Tests subroutine hgto_to_lcgto.
!!
!!  2011-07-27, Bin Gao
!!  * first version

  !> \brief tests subroutine hgto_to_lcgto
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_hgto_to_lcgto(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer angulars_ket(2)    !angular numbers of ket center
    integer orders_geo_ket(2)  !orders of geometric derivatives of ket center
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_real_number("Orbital exponent of ket center:", &
                              EXPONENT_KET, "F16.8", io_log)
    ! the following tests are for individual cases in \fn(hgto_to_lcgto)
    angulars_ket = (/0,0/)
    orders_geo_ket = (/6,8/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/0,1/)
    orders_geo_ket = (/6,8/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/1,1/)
    orders_geo_ket = (/6,8/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/0,2/)
    orders_geo_ket = (/4,7/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/1,2/)
    orders_geo_ket = (/4,7/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/2,2/)
    orders_geo_ket = (/4,7/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/3,6/)
    orders_geo_ket = (/4,7/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/5,5/)
    orders_geo_ket = (/4,7/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/4,6/)
    orders_geo_ket = (/4,7/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    angulars_ket = (/0,5/)
    orders_geo_ket = (/4,10/)
    call sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    return
  end subroutine test_hgto_to_lcgto

  !> \brief tests subroutine hgto_to_lcgto with specific orders of
  !>        CGTO and geometric derivative
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param angulars_ket is the angular number
  !> \param orders_geo_ket is the order of geometric derivative
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine sub_test_hgto_to_lcgto(angulars_ket, orders_geo_ket, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    ! recursive functions related to HGTOs
    use recur_hgto
    implicit none
    integer, intent(in) :: angulars_ket(2)
    integer, intent(in) :: orders_geo_ket(2)
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
! parameters of testing primitive integrals
#include "test_f90/prim_gto.h"
    integer, parameter :: DIM_CGTO_BRA = 3  !dimensions of CGTOs on bra center
    integer dim_hgto_ket                    !dimension of HGTOs on ket center
    integer, parameter :: DIM_GEO_BRA = 3   !dimensions of geometric derivative on bra center
    integer, parameter :: NUM_OPT = 1       !number of operators
    real(REALK), allocatable :: hgto_pints(:,:,:,:)
                                            !primitive HGTO integrals
    integer dim_cgto_ket                    !dimension of CGTOs on ket center
    integer dim_geo_ket                     !dimension of geometric derivative on ket center
    real(REALK), allocatable :: cgto_pints(:,:,:,:,:)
                                            !primitive CGTO integrals
    integer order_hgto_ket(2)               !orders of HGTOs on ket center
    integer iorder                          !incremental recorder over orders of HGTO on ket center
    integer x_hgto, y_hgto, z_hgto          !orders of xyz components of HGTO on ket center
    integer x_cgto, y_cgto, z_cgto          !orders of xyz components of CGTO on ket center
    integer addr_hgto, addr_cgto            !addresses of HGTO and CGTO
    integer iopt, igeo, ibra                !incremental recorders
    integer angular_ket                     !angular number of ket center
    integer order_geo_ket                   !order of geometric derivatives of ket center
    real(REALK) hrp_ket_expnt               !\f$1/2b\f$, where \f$b\f$ is the exponent on ket center
    real(REALK) recur_pint                  !results from recursive functions
    integer ierr                            !error information
    real(REALK) begin_time                  !begin of CPU time
    real(REALK) end_time                    !end of CPU time
    logical different                       !if result from \fn(hgto_to_lcgto) is different from reference
    real(REALK) xyz_function                !function returning artificial value for given xyz orders
    ! dumps information of tests
    call html_log_int_array(log_text="Angular momentums on ket center:",       &
                            log_int=angulars_ket, fmt_int="I3", io_log=io_log, &
                            separation="&rarr;")
    call html_log_int_array(log_text="Orders of geometric derivative on ket center:", &
                            log_int=orders_geo_ket, fmt_int="I3", io_log=io_log,      &
                            separation="&rarr;")
    ! gets the begin time
    call xtimer_set(begin_time)
    ! computes the minimum and maximum orders of HGTOs
    order_hgto_ket(1) = orders_geo_ket(1)
    order_hgto_ket(2) = orders_geo_ket(2)+angulars_ket(2)
    ! allocates memory for Hermite integrals
    if (order_hgto_ket(1)>0) then
      dim_hgto_ket = ((order_hgto_ket(2)+1)*(order_hgto_ket(2)+2)*(order_hgto_ket(2)+3) &
                   -  order_hgto_ket(1)*(order_hgto_ket(1)+1)*(order_hgto_ket(1)+2))/6
    else
      dim_hgto_ket = (order_hgto_ket(2)+1)*(order_hgto_ket(2)+2)*(order_hgto_ket(2)+3)/6
    end if
    allocate(hgto_pints(DIM_CGTO_BRA,dim_hgto_ket,DIM_GEO_BRA,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate hgto_pints:",             &
                               log_int=DIM_CGTO_BRA*dim_hgto_ket*DIM_GEO_BRA*NUM_OPT, &
                               fmt_int="I12", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! assigns the faked Hermite integrals
    do iopt = 1, NUM_OPT
      do igeo = 1, DIM_GEO_BRA
        addr_hgto = 0
        do iorder = order_hgto_ket(1), order_hgto_ket(2)
          do z_hgto = 0, iorder
            do y_hgto = 0, iorder-z_hgto
              x_hgto = iorder-(z_hgto+y_hgto)
              addr_hgto = addr_hgto+1
              hgto_pints(:,addr_hgto,igeo,iopt) = xyz_function((/x_hgto,y_hgto,z_hgto/))
            end do
          end do
        end do
      end do
    end do
    ! allocate memory for Cartesian integrals
    if (angulars_ket(1)>0) then
      dim_cgto_ket = ((angulars_ket(2)+1)*(angulars_ket(2)+2)*(angulars_ket(2)+3) &
                   -  angulars_ket(1)*(angulars_ket(1)+1)*(angulars_ket(1)+2))/6
    else
      dim_cgto_ket = (angulars_ket(2)+1)*(angulars_ket(2)+2)*(angulars_ket(2)+3)/6
    end if
    if (orders_geo_ket(1)>0) then
      dim_geo_ket = ((orders_geo_ket(2)+1)*(orders_geo_ket(2)+2)*(orders_geo_ket(2)+3) &
                  -  orders_geo_ket(1)*(orders_geo_ket(1)+1)*(orders_geo_ket(1)+2))/6
    else
      dim_geo_ket = (orders_geo_ket(2)+1)*(orders_geo_ket(2)+2)*(orders_geo_ket(2)+3)/6
    end if
    allocate(cgto_pints(DIM_CGTO_BRA,dim_cgto_ket,DIM_GEO_BRA,dim_geo_ket,NUM_OPT), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate cgto_pints:",    &
                               log_int=DIM_CGTO_BRA*dim_cgto_ket*DIM_GEO_BRA &
                                       *dim_geo_ket*NUM_OPT,                 &
                               fmt_int="I12", io_log=io_log, font_color="red")
      deallocate(hgto_pints)
      test_failed = .true.
      return
    end if
    ! calls \fn(hgto_to_lcgto) to get the integrals using CGTOs
    call hgto_to_lcgto(angulars_ket, orders_geo_ket, EXPONENT_KET, &
                       DIM_CGTO_BRA, dim_hgto_ket, DIM_GEO_BRA,    &
                       NUM_OPT, hgto_pints, dim_cgto_ket, dim_geo_ket, cgto_pints)
    ! cleans
    deallocate(hgto_pints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for hgto_to_lcgto:",    &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! \f$1/(2b)\f$
    hrp_ket_expnt = 0.5_REALK/EXPONENT_KET
    ! checks the results
    do iopt = 1, NUM_OPT
      do igeo = 1, DIM_GEO_BRA
        ! initializes the address of integrals from \fn(hgto_to_lcgto)
        addr_hgto = 0
        ! loops over different orders of geometric derivatives
        do order_geo_ket = orders_geo_ket(1), orders_geo_ket(2)
          ! loops over the xyz components of geometric derivative on ket center
          do z_hgto = 0, order_geo_ket
            do y_hgto = 0, order_geo_ket-z_hgto
              x_hgto = order_geo_ket-(z_hgto+y_hgto)
              addr_hgto = addr_hgto+1
              ! initializes the address of integrals from \fn(hgto_to_lcgto)
              addr_cgto = 0
              ! loops over angular numbers
              do angular_ket = angulars_ket(1), angulars_ket(2)
                ! loops over the xyz components of CGTO on ket center
                do z_cgto = 0, angular_ket
                  do y_cgto = 0, angular_ket-z_cgto
                    x_cgto = angular_ket-(z_cgto+y_cgto)
                    addr_cgto = addr_cgto+1
                    do ibra = 1, DIM_CGTO_BRA
                      ! integral from \fn(recur_hgto_to_cgto)
                      recur_pint = recur_hgto_to_cgto(hrp_ket_expnt,            &
                                                      (/x_cgto,y_cgto,z_cgto/), &
                                                      (/x_hgto,y_hgto,z_hgto/))
                      call check_difference(recur_pint,                      &
                             cgto_pints(ibra,addr_cgto,igeo,addr_hgto,iopt), &
                             different)
                      if (different) then
                        call html_log_int_array(log_text="Orders of CGTO:",       &
                                                log_int=(/x_cgto,y_cgto,z_cgto/), &
                                                fmt_int="I3", io_log=io_log)
                        call html_log_int_array(log_text="Orders of HGTO:",       &
                                                log_int=(/x_hgto,y_hgto,z_hgto/), &
                                                fmt_int="I3", io_log=io_log)
                        call html_log_real_number(                             &
                               log_text="Reference from recur_hgto_to_lcgto:", &
                               log_real=recur_pint, fmt_real="Es20.12",        &
                               io_log=io_log, font_color="blue")
                        call html_log_real_number(                                      &
                               log_text="Result from hgto_to_lcgto:",                   &
                               log_real=cgto_pints(ibra,addr_cgto,igeo,addr_hgto,iopt), &
                               fmt_real="Es20.12", io_log=io_log, font_color="red")
                        test_failed = .true.
                      end if
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for recur_hgto_to_lcgto:", &
                              log_real=end_time-begin_time, fmt_real="F10.4",    &
                              io_log=io_log, font_color="blue")
    ! cleans
    deallocate(cgto_pints)
    return
  end subroutine sub_test_hgto_to_lcgto
