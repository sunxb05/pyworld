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
!!  Tests subroutines in geom_total.F90.
!!
!!  2011-03-04, Bin Gao
!!  * first version

  !> \brief tests subroutines in geom_total.F90
  !> \author Bin Gao
  !> \date 2011-03-04
  !> \param io_log is the IO unit of log file
  !> \return test_failed indicates if the test is failed
  subroutine test_geom_total(io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: io_log
    logical, intent(out) :: test_failed
! parameters and referenced results of total geometric derivatives
#include "geom_total.h"
    integer num_paths                    !number of different paths from \fn{geom_total_tree_search}
    integer visit_height                 !height of atom node to visit
    integer, allocatable :: idx_node(:)  !selected atom nodes
    integer, allocatable :: wt_node(:)   !weights of the selected atom nodes
    integer idx_cent(MAX_NCENT(1))       !indices of differentiated centers
    integer order_cent(MAX_NCENT(1))     !orders of geometric derivatives of differentiated centers
    integer curr_num_geo                 !number of geometric derivatives of differentiated centers
    integer tot_num_geo                  !total number of geometric derivatives
    integer ipath                        !incremental recorder over different paths
    integer inode                        !incremental recorder over nodes
    integer base_ref                     !base address of referenced results
    integer itst                         !incremental recorder over tests
    integer, parameter :: NUM_COORD = 3*NUM_ATOMS(1)
                                         !number of coordinates
    integer power_nc(ORDER_GEO(1))       !powers of \var(NUM_COORD)
    integer base_nc(ORDER_GEO(1))        !base-\var(NUM_COORD)
    integer, parameter :: NUM_DENS = 3   !number of AO density matrices
    integer, parameter :: DIM_REDUNT = NUM_COORD**ORDER_GEO(1)
                                         !dimensions of redundant total geometric derivatives
    real(REALK), allocatable :: redunt_expt(:,:)
                                         !expectation values of redundant total geometric derivatives
    real(REALK), allocatable :: ref_rexpt(:,:)
                                         !referenced expectation values of redundant total geometric derivatives
    integer iderv, jderv                 !incremental recorder over derivatives
    integer iorder                       !incremental recorder over orders of derivatives
    integer idens                        !incremental recorder over AO density matrices
    integer ierr                         !error information
    logical different                    !if result from \fn(geom_total_redunt_expectation)
                                         !is different from reference
    ! assumes the test will pass
    test_failed = .false.
    ! dumps information of tests
    call html_log_heading(heading="Tests subroutine geom_total_num_derv", &
                          io_log=io_log, level=4)
    do itst = 1, NUM_TEST
      call html_log_int_number("Number of atoms:", NUM_ATOMS(itst), "I4", io_log)
      if (ORDER_GEO(1)<=0) then
        call html_log_int_number(                                           &
               log_text="Order of geometric derivatives is less than 1, "// &
                        "skips current test:",                              &
               log_int=ORDER_GEO(itst), fmt_int="I4", io_log=io_log, font_color="red")
        test_failed = .true.
        cycle
      else
        call html_log_int_number("Order of geometric derivatives:", &
                                 ORDER_GEO(itst), "I4", io_log)
      end if
      call html_log_int_number("Maximum number of differentiated centers:", &
                               MAX_NCENT(itst), "I4", io_log)
      call geom_total_num_derv(ORDER_GEO(itst), MAX_NCENT(itst), &
                               NUM_ATOMS(itst), tot_num_geo)
      ! checks the result
      if (tot_num_geo/=REF_TOT_NGEO(itst)) then
        call html_log_int_number(                                            &
               log_text="Referenced number of total geometric derivatives:", &
               log_int=REF_TOT_NGEO(itst), fmt_int="I8", io_log=io_log,      &
               font_color="blue")
        call html_log_int_number(log_text="Number of total geometric derivatives "// &
                                          "from geom_total_num_derv:",               &
                                 log_int=tot_num_geo, fmt_int="I8", io_log=io_log,   &
                                 font_color="red")
        test_failed = .true.
        cycle
      else
        call html_log_int_number("Number of total geometric derivatives "// &
                                 "from geom_total_num_derv:", tot_num_geo,  &
                                 "I8", io_log)
      end if
    end do
    ! dumps information of tests
    call html_log_heading(heading="Tests other subroutines in geom_total.F90", &
                          io_log=io_log, level=4)
    call html_log_int_number("Number of atoms:", NUM_ATOMS(1), "I4", io_log)
    if (ORDER_GEO(1)<=0) then

      call html_log_int_number(                                           &
             log_text="Order of geometric derivatives is less than 1, "// &
                      "no test will be performed:",                       &
             log_int=ORDER_GEO(1), fmt_int="I4", io_log=io_log, font_color="red")
      test_failed = .true.
      return
    else
      call html_log_int_number("Order of geometric derivatives:", &
                               ORDER_GEO(1), "I4", io_log)
    end if
    call html_log_int_number("Maximum number of differentiated centers:", &
                             MAX_NCENT(1), "I4", io_log)
    ! allocates memory for the selected atom nodes and their weights
    allocate(idx_node(ORDER_GEO(1)), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate idx_node:", &
                               log_int=ORDER_GEO(1), fmt_int="I12",     &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    allocate(wt_node(ORDER_GEO(1)), stat=ierr)
    if (ierr/=0) then
      deallocate(idx_node)
      call html_log_int_number(log_text="Failed to allocate wt_node:", & 
                               log_int=ORDER_GEO(1), fmt_int="I12",    & 
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! allocates memory for the expectation values of redundant total geometric derivatives
    allocate(redunt_expt(DIM_REDUNT,NUM_DENS), stat=ierr)
    if (ierr/=0) then
      deallocate(idx_node)
      deallocate(wt_node)
      call html_log_int_number(log_text="Failed to allocate redunt_expt:", & 
                               log_int=DIM_REDUNT*NUM_DENS, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    allocate(ref_rexpt(DIM_REDUNT,NUM_DENS), stat=ierr)
    if (ierr/=0) then
      deallocate(idx_node)
      deallocate(wt_node)
      deallocate(redunt_expt)
      call html_log_int_number(log_text="Failed to allocate ref_rexpt:",   &
                               log_int=DIM_REDUNT*NUM_DENS, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      test_failed = .true.
      return
    end if
    ! sets the powers of \var(NUM_COORD)
    power_nc(1) = 1
    do iorder = 2, ORDER_GEO(1)
      power_nc(iorder) = NUM_COORD*power_nc(iorder-1)
    end do
    ! sets the expectation values of redundant total geometric derivatives
    redunt_expt = 0.0_REALK
    do iderv = 1, DIM_REDUNT
      ref_rexpt(iderv,1) = 1.0_REALK
      jderv = iderv-1
      do iorder = ORDER_GEO(1), 2, -1
        if (jderv>=power_nc(iorder)) then
          base_nc(iorder) = (jderv-mod(jderv,power_nc(iorder)))/power_nc(iorder)
          ref_rexpt(iderv,1) = ref_rexpt(iderv,1)*real(base_nc(iorder)+1,REALK)*0.1_REALK
          jderv = jderv-base_nc(iorder)*power_nc(iorder)
        else
          ref_rexpt(iderv,1) = ref_rexpt(iderv,1)*0.1_REALK
        end if
      end do
      ref_rexpt(iderv,1) = ref_rexpt(iderv,1)*real(jderv+1,REALK)*0.1_REALK
    end do
    do idens = 2, NUM_DENS
      do iderv = 1, DIM_REDUNT
        ref_rexpt(iderv,idens) = ref_rexpt(iderv,1)+real(idens-1,REALK)
      end do
    end do
    ! initializes the total number of geometric derivatives
    tot_num_geo = 0
    ! computes the number of different paths and generates the first compositions of centers
    call geom_total_tree_init(NUM_ATOMS(1), ORDER_GEO(1), MAX_NCENT(1),   &
                              num_paths, visit_height, idx_node, wt_node, &
                              idx_cent, order_cent, curr_num_geo)
    ! checks the number of different paths
    if (num_paths/=REF_NUM_PATHS) then
      call html_log_int_number(log_text="Referenced number of different paths:",   &
                               log_int=REF_NUM_PATHS, fmt_int="I4", io_log=io_log, &
                               font_color="blue")
      call html_log_int_number(                                                 &
               log_text="Number of different paths from geom_total_tree_init:", &
               log_int=num_paths, fmt_int="I4", io_log=io_log, font_color="red")
      deallocate(idx_node)
      deallocate(wt_node)
      deallocate(ref_rexpt)
      deallocate(redunt_expt)
      test_failed = .true.
      return
    else
      call html_log_int_number("Number of different paths:", num_paths, "I4", io_log)
    end if
    ! begins the table in log file
    call html_log_table_begin(io_log)
    ! writes the header cells of the table
    call html_log_table_header((/"Path            ", "Start           ",   &
                                 "N<sub>cent</sub>", "N<sub>geo</sub> ",   &
                                 "Centers         ", "Orders          "/), &
                               io_log)
    ipath = 1
    base_ref = 0
    ! checks the first path
    if (wt_node(ORDER_GEO(1))/=REF_NUM_CENTS(ipath)) test_failed = .true.
    if (curr_num_geo/=REF_NGEO_CENT(ipath)) test_failed = .true.
    do inode = 1, wt_node(ORDER_GEO(1))
      if (idx_cent(inode)/=REF_IDX_CENT(base_ref+inode) .or. &
          order_cent(inode)/=REF_ORDER_CENT(base_ref+inode)) then
        test_failed = .true.
        exit
      end if
    end do
    if (test_failed) then
      call html_log_table_row(                                                           &
             cell_int=(/ipath,tot_num_geo+1, REF_NUM_CENTS(ipath), REF_NGEO_CENT(ipath), &
                        REF_IDX_CENT(base_ref+1:base_ref+REF_NUM_CENTS(ipath)),          &
                        REF_ORDER_CENT(base_ref+1:base_ref+REF_NUM_CENTS(ipath))/),      &
             num_int=(/1,1,1,1,REF_NUM_CENTS(ipath),REF_NUM_CENTS(ipath)/),              &
             fmt_int=(/"I4","I6","I3","I4","I3","I2"/), io_log=io_log, font_color="blue")
      call html_log_table_row(cell_int=(/ipath,tot_num_geo+1,                   &
                                         wt_node(ORDER_GEO(1)),                 &
                                         curr_num_geo,                          &
                                         idx_cent(1:wt_node(ORDER_GEO(1))),     &
                                         order_cent(1:wt_node(ORDER_GEO(1)))/), &
                              num_int=(/1,1,1,1,wt_node(ORDER_GEO(1)),          &
                                        wt_node(ORDER_GEO(1))/),                &
                              fmt_int=(/"I4","I6","I3","I4","I3","I2"/),        &
                              io_log=io_log, font_color="red")
      call html_log_table_end(io_log)
      deallocate(idx_node)
      deallocate(wt_node)
      deallocate(ref_rexpt)
      deallocate(redunt_expt)
      ! backs without any further test
      return
    else
      ! dumps the first path
      call html_log_table_row(cell_int=(/ipath,tot_num_geo+1,                   &
                                         wt_node(ORDER_GEO(1)),                 &
                                         curr_num_geo,                          &
                                         idx_cent(1:wt_node(ORDER_GEO(1))),     &
                                         order_cent(1:wt_node(ORDER_GEO(1)))/), &
                              num_int=(/1,1,1,1,wt_node(ORDER_GEO(1)),          &
                                        wt_node(ORDER_GEO(1))/),                &
                              fmt_int=(/"I4","I6","I3","I4","I3","I2"/),        &
                              io_log=io_log)
      ! recovers expectation values of redundant total geometric derivatives
      call test_redundant_expt(NUM_ATOMS(1), wt_node(ORDER_GEO(1)), &
                               idx_cent(1:wt_node(ORDER_GEO(1))),   &
                               order_cent(1:wt_node(ORDER_GEO(1))), &
                               curr_num_geo, DIM_REDUNT, NUM_DENS,  &
                               redunt_expt, io_log, test_failed)
    end if
    ! updates the base address of referenced results
    base_ref = base_ref+MAX_NCENT(1)
    ! updates the total number of geometric derivatives
    tot_num_geo = tot_num_geo+curr_num_geo
    ! the left compositions of centers
    do ipath = 2, num_paths  !or do while (idx_node(1)/=NUM_ATOMS(1))
      call geom_total_tree_search(NUM_ATOMS(1), ORDER_GEO(1), MAX_NCENT(1),  &
                                  visit_height, idx_node, wt_node, idx_cent, &
                                  order_cent, curr_num_geo)
      ! checks current path
      if (wt_node(ORDER_GEO(1))/=REF_NUM_CENTS(ipath)) test_failed = .true.
      if (curr_num_geo/=REF_NGEO_CENT(ipath)) test_failed = .true.
      do inode = 1, wt_node(ORDER_GEO(1))
        if (idx_cent(inode)/=REF_IDX_CENT(base_ref+inode) .or. &
            order_cent(inode)/=REF_ORDER_CENT(base_ref+inode)) then
          test_failed = .true.
          exit
        end if
      end do
      if (test_failed) then
        call html_log_table_row(                                                           &
               cell_int=(/ipath,tot_num_geo+1, REF_NUM_CENTS(ipath), REF_NGEO_CENT(ipath), &
                          REF_IDX_CENT(base_ref+1:base_ref+REF_NUM_CENTS(ipath)),          &
                          REF_ORDER_CENT(base_ref+1:base_ref+REF_NUM_CENTS(ipath))/),      &
               num_int=(/1,1,1,1,REF_NUM_CENTS(ipath),REF_NUM_CENTS(ipath)/),              &
               fmt_int=(/"I4","I6","I3","I4","I3","I2"/), io_log=io_log, font_color="blue")
        call html_log_table_row(cell_int=(/ipath,tot_num_geo+1,                   &
                                           wt_node(ORDER_GEO(1)),                 &
                                           curr_num_geo,                          &
                                           idx_cent(1:wt_node(ORDER_GEO(1))),     &
                                           order_cent(1:wt_node(ORDER_GEO(1)))/), &
                                num_int=(/1,1,1,1,wt_node(ORDER_GEO(1)),          &
                                          wt_node(ORDER_GEO(1))/),                &
                                fmt_int=(/"I4","I6","I3","I4","I3","I2"/),        &
                                io_log=io_log, font_color="red")
        call html_log_table_end(io_log)
        deallocate(idx_node)
        deallocate(wt_node)
        deallocate(ref_rexpt)
        deallocate(redunt_expt)
        ! backs without any further test
        return
      else
        ! dumps current path
        call html_log_table_row(cell_int=(/ipath,tot_num_geo+1,                   &
                                           wt_node(ORDER_GEO(1)),                 &
                                           curr_num_geo,                          &
                                           idx_cent(1:wt_node(ORDER_GEO(1))),     &
                                           order_cent(1:wt_node(ORDER_GEO(1)))/), &
                                num_int=(/1,1,1,1,wt_node(ORDER_GEO(1)),          &
                                          wt_node(ORDER_GEO(1))/),                &
                                fmt_int=(/"I4","I6","I3","I4","I3","I2"/),        &
                                io_log=io_log)
        ! recovers expectation values of redundant total geometric derivatives
        call test_redundant_expt(NUM_ATOMS(1), wt_node(ORDER_GEO(1)), &
                                 idx_cent(1:wt_node(ORDER_GEO(1))),   &
                                 order_cent(1:wt_node(ORDER_GEO(1))), &
                                 curr_num_geo, DIM_REDUNT, NUM_DENS,  &
                                 redunt_expt, io_log, test_failed)
      end if
      ! updates the base address of referenced results
      base_ref = base_ref+MAX_NCENT(1)
      ! updates the total number of geometric derivatives
      tot_num_geo = tot_num_geo+curr_num_geo
    end do
    call html_log_table_end(io_log)
    ! cleans
    deallocate(idx_node)
    deallocate(wt_node)
    ! checks
    if (tot_num_geo/=REF_TOT_NGEO(1)) then
      call html_log_int_number(                                            &
             log_text="Referenced total number of geometric derivatives:", &
             log_int=REF_TOT_NGEO(1), fmt_int="I6", io_log=io_log, font_color="blue")
      call html_log_int_number(log_text="Total number of geometric derivatives:", &
                               log_int=tot_num_geo, fmt_int="I6", io_log=io_log,  &
                               font_color="red")
      test_failed = .true.
      return
    else
      call html_log_int_number("Total number of geometric derivatives:", &
                               tot_num_geo, "I6", io_log)
    end if
    ! checks the expectation values of redundant total geometric derivatives
    do idens = 1, NUM_DENS
      do iderv = 1, DIM_REDUNT
        call check_difference(ref_rexpt(iderv,idens), &
                              redunt_expt(iderv,idens), different)
        if (different) then
          call html_log_int_array(                                          &
                 log_text="ID of geometric derivative and density matrix:", &
                 log_int=(/iderv,idens/), fmt_int="I6", io_log=io_log)
          call html_log_real_number(                                 &
                 log_text="Referenced redunt_expt(iderv,idens):",    &
                 log_real=ref_rexpt(iderv,idens), fmt_real="Es16.8", &
                 io_log=io_log, font_color="blue")
          call html_log_real_number(                                    &
                 log_text="Result from geom_total_redunt_expectation:", &
                 log_real=redunt_expt(iderv,idens), fmt_real="Es16.8",  &
                 io_log=io_log, font_color="red")
          test_failed = .true.
        end if
      end do
    end do
    deallocate(ref_rexpt)
    deallocate(redunt_expt)
    return
  end subroutine test_geom_total

  !> \brief tests subroutine \fn(geom_total_redunt_expectation)
  !> \author Bin Gao
  !> \date 2011-11-07
  !> \param num_atoms is the number of atoms
  !> \param num_cent is the number of differentiated centers
  !> \param idx_cent contains the indices of differentiated centers
  !> \param order_cent contains the orders of differentiated centers
  !> \param num_unique_geo is the number of unique total geometric derivatives
  !> \param dim_redunt is the dimensions of redundant total geometric derivatives,
  !>        equals to \var(3*num_atoms)^sum(\var(order_cent))
  !> \param num_dens is the number of AO density matrices
  !> \param io_log is the IO unit of log file
  !> \return redunt_expt contains the updated expectation values of redundant total
  !>         geometric derivatives on exit
  !> \return test_failed indicates if the test is failed
  subroutine test_redundant_expt(num_atoms, num_cent, idx_cent, order_cent, &
                                 num_unique_geo, dim_redunt, num_dens,      &
                                 redunt_expt, io_log, test_failed)
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, intent(in) :: num_atoms
    integer, intent(in) :: num_cent
    integer, intent(in) :: idx_cent(num_cent)
    integer, intent(in) :: order_cent(num_cent)
    integer, intent(in) :: num_unique_geo
    integer, intent(in) :: dim_redunt
    integer, intent(in) :: num_dens
    real(REALK), intent(inout) :: redunt_expt(dim_redunt,num_dens)
    integer, intent(in) :: io_log
    logical, intent(inout) :: test_failed
    real(REALK), allocatable :: unique_expt(:,:)
                                    !expectation values of unique total geometric derivatives
    integer iderv, jderv            !incremental recorder over derivatives
    integer idens                   !incremental recorder over AO density matrices
    integer last_prod_nderv         !product of number of derivatives up to last differentiated center
    integer curr_prod_nderv         !product of number of derivatives up to current differentiated center
    real(REALK) addr_cent_derv(3)   !addresses of first order derivatives for a given center
    real(REALK) expt_cent           !faked expectation values for a given center
    integer num_cent_derv           !number of derivatives for a given center
    integer x_derv, y_derv, z_derv  !xyz components of derivatives
    integer addr_uniq_geo           !address of unique total geometric derivatives
    integer inode                   !incremental recorder over nodes
    integer ierr                    !error information
    ! sets the expectation values of unique total geometric derivatives
    allocate(unique_expt(num_unique_geo,num_dens), stat=ierr)
    if (ierr/=0) then
      call html_log_int_number(log_text="Failed to allocate unique_expt:",     &
                               log_int=num_unique_geo*num_dens, fmt_int="I12", &
                               io_log=io_log, font_color="red")
      test_failed = .true.
    else
      last_prod_nderv = 1
      unique_expt(:,1) = 1.0_REALK
      do inode = 1, num_cent
        addr_cent_derv(1) = 0.3_REALK*real(idx_cent(inode),REALK)-0.2_REALK  !px
        addr_cent_derv(2) = addr_cent_derv(1)+0.1_REALK                      !py
        addr_cent_derv(3) = addr_cent_derv(2)+0.1_REALK                      !pz
        num_cent_derv = (order_cent(inode)+1)*(order_cent(inode)+2)/2
        curr_prod_nderv = num_cent_derv*last_prod_nderv
        ! generates the addresses of different total geometric derivatives
        addr_uniq_geo = 0
        do iderv = 1, num_unique_geo/curr_prod_nderv
          ! xyz components for the current differentitated center
          do z_derv = 0, order_cent(inode)
            do y_derv = 0, order_cent(inode)-z_derv
              x_derv = order_cent(inode)-(y_derv+z_derv)
              expt_cent = addr_cent_derv(1)**x_derv &
                        * addr_cent_derv(2)**y_derv &
                        * addr_cent_derv(3)**z_derv
              ! geometric derivatives of preceding differentitated centers
              do jderv = 1, last_prod_nderv
                addr_uniq_geo = addr_uniq_geo+1
                unique_expt(addr_uniq_geo,1) = unique_expt(addr_uniq_geo,1)*expt_cent
              end do
            end do
          end do
        end do
        last_prod_nderv = curr_prod_nderv
      end do
      do idens = 2, num_dens
        do iderv = 1, num_unique_geo
          unique_expt(iderv,idens) = unique_expt(iderv,1)+real(idens-1,REALK)
        end do
      end do
      ! calls subroutine \fn(geom_total_redunt_expectation) to recover expectation values
      ! of redundant total geometric derivatives
      call geom_total_redunt_expectation(num_atoms, num_cent, idx_cent, order_cent, &
                                         1, num_unique_geo, num_dens, unique_expt,  &
                                         dim_redunt, redunt_expt)
      deallocate(unique_expt)
    end if
    return
  end subroutine test_redundant_expt
