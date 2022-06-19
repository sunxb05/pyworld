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
!!  Tests subroutine carmom_deriv.
!!
!!  2011-07-27, Bin Gao
!!  * first version

  !> \brief tests subroutine carmom_deriv
  !> \author Bin Gao
  !> \date 2011-07-27
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_carmom_deriv(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
    ! see \fn(carmom_deriv) for the details of the following variables
    integer, parameter :: ORDER_LOW_MOM = 12
    integer, parameter :: DIM_CINTS = 21
    integer, parameter :: NUM_LOW_MOM = (ORDER_LOW_MOM+1)*(ORDER_LOW_MOM+2)/2
    integer, parameter :: NUM_OPT = 21
    real(REALK) low_mom_cints(DIM_CINTS,NUM_LOW_MOM,NUM_OPT)
    integer, parameter :: ORDER_MOM = 24
    integer, parameter :: ORDER_GEO_MOM = ORDER_MOM-ORDER_LOW_MOM
    integer, parameter :: NUM_MOM = (ORDER_MOM+1)*(ORDER_MOM+2)/2
    integer, parameter :: NUM_GEO_MOM = (ORDER_GEO_MOM+1)*(ORDER_GEO_MOM+2)/2
    real(REALK) geo_mom_cints(DIM_CINTS,NUM_MOM,NUM_OPT,NUM_GEO_MOM)
    integer x_geo, y_geo, z_geo  !order of geometric derivatives
    integer x_mom, y_mom, z_mom  !order of Cartesian multipole moments
    integer x_low, y_low, z_low  !order of lower order Cartesian multipole moments
    integer addr_geo             !address of geometric derivatives
    integer addr_mom             !address of Cartesian multipole moment
    integer addr_low             !address of lower order Cartesian multipole moment
    real(REALK) prod_coeff       !prodcut coefficient from xyz components
    real(REALK) neg_one          !prefactor of -1 or +1 of the prodcut coefficient
    integer iopt, icint, ilow    !incremental recorders
    logical different            !if result from \fn(carmom_deriv) is different from reference
    real(REALK) begin_time       !begin of CPU time
    real(REALK) end_time         !end of CPU time
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_int_number("Order of lower order Cartesian multipole moment:", &
                             ORDER_LOW_MOM, "I4", io_log)
    call html_log_int_number("Order of Cartesian multipole moment:", &
                             ORDER_MOM, "I4", io_log)
    call html_log_int_number("Order of geometric derivatives:", &
                             ORDER_GEO_MOM, "I4", io_log)
    ! assigns the lower order Cartesian multipole moment integrals
    do iopt = 1, NUM_OPT
      do ilow = 1, NUM_LOW_MOM
        do icint = 1, DIM_CINTS
          low_mom_cints(icint,ilow,iopt) = real(icint,REALK)*0.01_REALK &
                                         + real(ilow,REALK)*0.001_REALK &
                                         + real(iopt,REALK)*0.0001_REALK
        end do
      end do
    end do
    ! gets the begin time
    call xtimer_set(begin_time)
    ! gets the geometric derivatives of higher order Cartesian multipole moment integrals
    ! from \fn(carmom_deriv)
    call carmom_deriv(ORDER_LOW_MOM, DIM_CINTS, NUM_LOW_MOM, NUM_OPT, low_mom_cints, &
                      ORDER_MOM, ORDER_GEO_MOM, NUM_MOM, NUM_GEO_MOM, geo_mom_cints)
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for carmom_deriv:",     &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    ! resets the begin time
    call xtimer_set(begin_time)
    ! sets the prefactor of -1 or +1 of the prodcut coefficient
    if (mod(ORDER_GEO_MOM,2)==0) then
      neg_one = 1.0_REALK
    else
      neg_one = -1.0_REALK
    end if
    ! checks the results, the following do loops are the same as \fn(carmom_deriv)
    !
    ! loops over the xyz componets of geometric derivative
    addr_geo = 0
    do z_geo = 0, ORDER_GEO_MOM
      do y_geo = 0, ORDER_GEO_MOM-z_geo
        x_geo = ORDER_GEO_MOM-(z_geo+y_geo)
        addr_geo = addr_geo+1
        ! loops over different operators
        do iopt = 1, NUM_OPT
          ! loops over the xyz components of Cartesian multipole moment
          addr_mom = 0
          do z_mom = 0, ORDER_MOM
            do y_mom = 0, ORDER_MOM-z_mom
              x_mom = ORDER_MOM-(z_mom+y_mom)
              addr_mom = addr_mom+1
              ! the corresponding orders of lower order Cartesian multipole moment
              x_low = x_mom-x_geo
              y_low = y_mom-y_geo
              z_low = z_mom-z_geo
              if (x_low<0 .or. y_low<0 .or. z_low<0) then
                do icint = 1, DIM_CINTS
                  if (geo_mom_cints(icint,addr_mom,iopt,addr_geo)/=0.0_REALK) then
                    call html_log_int_array(                                &
                           log_text="Order of Cartesian multipole moment:", &
                           log_int=(/x_mom,y_mom,z_mom/), fmt_int="I3", io_log=io_log)
                    call html_log_int_array(                           &
                           log_text="Order of geometric derivatives:", &
                           log_int=(/x_geo,y_geo,z_geo/), fmt_int="I3", io_log=io_log)
                    call html_log_real_number(log_text="Referenced result:",         &
                                              log_real=0.0_REALK, fmt_real="Es16.8", &
                                              io_log=io_log, font_color="blue")
                    call html_log_real_number(                                   &
                           log_text="Result from carmom_deriv:",                 &
                           log_real=geo_mom_cints(icint,addr_mom,iopt,addr_geo), &
                           fmt_real="Es16.8", io_log=io_log, font_color="red")
                    test_failed = .true.
                  end if
                end do
              else
                ! sets the product coefficient
                prod_coeff = neg_one
                do ilow = x_low+1, x_mom
                  prod_coeff = prod_coeff*real(ilow,REALK)
                end do
                do ilow = y_low+1, y_mom
                  prod_coeff = prod_coeff*real(ilow,REALK)
                end do
                do ilow = z_low+1, z_mom
                  prod_coeff = prod_coeff*real(ilow,REALK)
                end do
                ! index of x^{l}y^{m}z^{n} with l+m+n=angm is 1+m+(2*angm+3-n)*n/2
                addr_low = 1+y_low+(2*ORDER_LOW_MOM+3-z_low)*z_low/2
                do icint = 1, DIM_CINTS
                  call check_difference(geo_mom_cints(icint,addr_mom,iopt,addr_geo),   &
                                        prod_coeff*low_mom_cints(icint,addr_low,iopt), &
                                        different)
                  if (different) then
                    call html_log_int_array(                                &
                           log_text="Order of Cartesian multipole moment:", &
                           log_int=(/x_mom,y_mom,z_mom/), fmt_int="I3", io_log=io_log)
                    call html_log_int_array(                           &
                           log_text="Order of geometric derivatives:", &
                           log_int=(/x_geo,y_geo,z_geo/), fmt_int="I3", io_log=io_log)
                    call html_log_real_number(log_text="Referenced result:",       &
                           log_real=prod_coeff*low_mom_cints(icint,addr_low,iopt), &
                           fmt_real="Es16.8", io_log=io_log, font_color="blue")
                    call html_log_real_number(                                   &
                           log_text="Result from carmom_deriv:",                 &
                           log_real=geo_mom_cints(icint,addr_mom,iopt,addr_geo), &
                           fmt_real="Es16.8", io_log=io_log, font_color="red")
                    test_failed = .true.
                  end if
                end do
              end if
            end do
          end do
        end do
      end do
    end do
    ! gets the end time
    call xtimer_set(end_time)
    ! prints the CPU elapsed time
    call html_log_real_number(log_text="Time (s) used for local loops:",      &
                              log_real=end_time-begin_time, fmt_real="F10.4", &
                              io_log=io_log, font_color="blue")
    return
  end subroutine test_carmom_deriv
