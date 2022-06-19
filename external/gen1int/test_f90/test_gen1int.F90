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
!!  This file is the FORTRAN test suite of Gen1Int.
!!
!!  2012-02-05, Bin Gao
!!  * uses module of HTML test log routines
!!  * adds option to test different recurrence relations
!!
!!  2011-07-27, Bin Gao
!!  * dumps the results in HTML file
!!
!!  2010-07-28, Bin Gao
!!  * first version

#include "stdout.h"

  !> \brief FORTRAN test suite of Gen1Int
  !> \author Bin Gao
  !> \date 2010-07-28
  program test_gen1int
    use xkind
    ! module of HTML test log routines
    use html_log
    implicit none
    integer, parameter :: IO_LOG = 7                            !IO of test log file
    character(*), parameter :: FNAME_LOG = "test_gen1int.html"  !file name of test log
    integer num_options               !number of user given options
    integer iopt                      !incremental recorder over different options
    character(7) user_option          !user given options
    logical :: run_test = .true.      !if running test suite
    logical :: test_tools = .false.   !if testing tool subroutines
    logical :: test_basic = .false.   !if testing some basic subroutines
    logical :: test_geo = .false.     !if testing subroutines related to geometric derivatives
    logical :: test_mag = .false.     !if testing subroutines related to magnetic derivatives
    logical :: test_carmom = .false.  !if testing Cartesian multipole moment integrals
    logical :: test_delta = .false.   !if testing Dirac delta function integrals
    logical :: test_nucpot = .false.  !if testing nuclear attraction potential integrals
    logical :: test_gaupot = .false.  !if testing Gaussian charge potential integrals
    logical :: test_odist = .false.   !if testing subroutines related to overlap distribution
#if defined(BUILD_F90_MODULE)
    logical :: test_f90mod = .false.  !if testing Fortran 90 modules
#endif
    integer iargc                     !function to get the number of user given options
    integer num_tests                 !number of tests
    logical test_failed               !indicates if the individual test is failed
    integer num_failed                !number of failed tests
    real(REALK) begin_time            !begin of CPU time
    real(REALK) end_time              !end of CPU time
    write(STDOUT,100) "test suite of Gen1Int ;-)"
    ! gets the number of user given options
    num_options = iargc()
    ! no option given, prints help menu
    if (num_options==0) then
      call test_help_menu()
      run_test = .false.
    else
      ! reads different user given options
      do iopt = 1, num_options
        call getarg(iopt, user_option)
        select case(trim(user_option))
        case("-all")
          test_tools = .true.
          test_basic = .true.
          test_geo = .true.
          test_mag = .true.
          test_carmom = .true.
          test_delta = .true.
          test_nucpot = .true.
          test_gaupot = .true.
          test_odist = .true.
#if defined(BUILD_F90_MODULE)
          test_f90mod = .true.
#endif
        case("-tools")
          test_tools = .true.
        case("-basic")
          test_basic = .true.
        case("-geo")
          test_geo = .true.
        case("-mag")
          test_mag = .true.
        case("-carmom")
          test_carmom = .true.
        case("-delta")
          test_delta = .true.
        case("-nucpot")
          test_nucpot = .true.
        case("-gaupot")
          test_gaupot = .true.
        case("-odist")
          test_odist = .true.
#if defined(BUILD_F90_MODULE)
        case("-f90mod")
          test_f90mod = .true.
#endif
        case("-help")
          call test_help_menu()
          run_test = .false.
        case default
          write(STDOUT,100) "Unkown option: "//trim(user_option)
          call test_help_menu()
          run_test = .false.
        end select
      end do
    end if
    ! test suite begins ...
    if (run_test) then
      ! creates an HTML test log file
      write(STDOUT,100) "test log is written into file "//FNAME_LOG
      open(IO_LOG, file=FNAME_LOG, status="unknown")
      call html_log_begin(IO_LOG)
      call html_log_head_begin(IO_LOG)
      call html_log_title("Test log of Gen1Int", IO_LOG)
      call html_log_head_end(IO_LOG)
      call html_log_body_begin(IO_LOG)
      call html_log_heading("Test log of Gen1Int", IO_LOG, 2)
      call html_log_date_and_time(IO_LOG)
      call html_log_text("Note: errors are marked in red.", IO_LOG, "4", "red")
      call html_log_text("Note: CPU usage and reference results are marked in blue.", &
                         IO_LOG, "4", "blue")
      call html_log_horizontal_line(IO_LOG)
      ! gets the begin time
      call xtimer_set(begin_time)
      ! tests different subroutines
      num_tests = 0
      num_failed = 0
      ! (1) tests tool subroutines
      if (test_tools) then
        ! tests subroutines \fn(norm_contr_cgto) and \fn(norm_contr_sgto)
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines norm_contr_cgto and norm_contr_sgto ..."
        call html_log_heading("Tests subroutines norm_contr_cgto and norm_contr_sgto", &
                              IO_LOG, 3)
        call test_norm_contr_gto(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines norm_contr_cgto and norm_contr_sgto failed!"
          num_failed = num_failed+1
        end if
        ! test subroutines in reorder_ints.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in reorder_ints.F90 ..."
        call html_log_heading("Tests subroutines in reorder_ints.F90", IO_LOG, 3)
        call test_reorder_ints(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in reorder_ints.F90 failed!"
          num_failed = num_failed+1
        end if
        ! test subroutines in trace_ints.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in trace_ints.F90 ..."
        call html_log_heading("Tests subroutines in trace_ints.F90", IO_LOG, 3)
        call test_trace_ints(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in trace_ints.F90 failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (2) tests some basic subroutines
      if (test_basic) then
        ! test subroutines in binom_coeff.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in binom_coeff.F90 ..."
        call html_log_heading("Tests subroutines in binom_coeff.F90", IO_LOG, 3)
        call test_binom_coeff(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in binom_coeff.F90 failed!"
          num_failed = num_failed+1
        end if
        ! test subroutine next_permutation
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine next_permutation ..."
        call html_log_heading("Tests subroutine next_permutation", IO_LOG, 3)
        call test_next_permutation(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine next_permutation failed!"
          num_failed = num_failed+1
        end if
        ! tests sort_cents.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in sort_cents.F90 ..."
        call html_log_heading("Tests subroutines in sort_cents.F90", IO_LOG, 3)
        call test_sort_cents(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in sort_cents.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests shell_scatter.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in shell_scatter.F90 ..."
        call html_log_heading("Tests subroutines in shell_scatter.F90", IO_LOG, 3)
        call test_shell_scatter(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in shell_scatter.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in hgto_to_cgto.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in hgto_to_cgto.F90 ..."
        call html_log_heading("Tests subroutines in hgto_to_cgto.F90", IO_LOG, 3)
        call test_hgto_to_cgto(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in hgto_to_cgto.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutine const_contr_ints
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine const_contr_ints ..."
        call html_log_heading("Tests subroutine const_contr_ints", IO_LOG, 3)
        call test_const_contr_ints(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine const_contr_ints failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (3) tests subroutines related to geometric derivatives
      if (test_geo) then
        ! tests geom_total.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in geom_total.F90 ..."
        call html_log_heading("Tests subroutines in geom_total.F90", IO_LOG, 3)
        call test_geom_total(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in geom_total.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in geom_part_zero.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in geom_part_zero.F90 ..."
        call html_log_heading("Tests subroutines in geom_part_zero.F90", IO_LOG, 3)
        call test_geom_part_zero(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in geom_part_zero.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in geom_part_one.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in geom_part_one.F90 ..."
        call html_log_heading("Tests subroutines in geom_part_one.F90", IO_LOG, 3)
        call test_geom_part_one(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in geom_part_one.F90 failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (4) tests subroutines related to magnetic derivatives
      if (test_mag) then  
        ! tests subroutines in london_mom_hgto.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in london_mom_hgto.F90 ..."
        call html_log_heading("Tests subroutines in london_mom_hgto.F90", IO_LOG, 3)
        call test_london_mom_hgto(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in london_mom_hgto.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutine in hgto_to_lcgto
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine hgto_to_lcgto ..."
        call html_log_heading("Tests subroutine hgto_to_lcgto", IO_LOG, 3)
        call test_hgto_to_lcgto(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine hgto_to_lcgto failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (5) tests Cartesian multipole moment integrals
      if (test_carmom) then 
        ! tests subroutine carmom_deriv
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine carmom_deriv ..."
        call html_log_heading("Tests subroutine carmom_deriv", IO_LOG, 3)
        call test_carmom_deriv(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine carmom_deriv failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutine carmom_hbra
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine carmom_hbra ..."
        call html_log_heading("Tests subroutine carmom_hbra", IO_LOG, 3)
        call test_carmom_hbra(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine carmom_hbra failed!"
          num_failed = num_failed+1
        end if
        ! tests carmom_hrr_ket.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in carmom_hrr_ket.F90 ..."
        call html_log_heading("Tests subroutines in carmom_hrr_ket.F90", IO_LOG, 3)
        call test_carmom_hrr_ket(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in carmom_hrr_ket.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests carmom_moment.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in carmom_moment.F90 ..."
        call html_log_heading("Tests subroutines in carmom_moment.F90", IO_LOG, 3)
        call test_carmom_moment(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in carmom_moment.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests prim_hgto_carmom.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine prim_hgto_carmom ..."
        call html_log_heading("Tests subroutine prim_hgto_carmom", IO_LOG, 3)
        call test_prim_hgto_carmom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine prim_hgto_carmom failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in contr_cgto_carmom.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in contr_cgto_carmom.F90 ..."
        call html_log_heading("Tests subroutines in contr_cgto_carmom.F90", IO_LOG, 3)
        call test_contr_cgto_carmom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in contr_cgto_carmom.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in contr_sgto_carmom.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in contr_sgto_carmom.F90 ..."
        call html_log_heading("Tests subroutines in contr_sgto_carmom.F90", IO_LOG, 3)
        call test_contr_sgto_carmom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in contr_sgto_carmom.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in contr_csgto_carmom.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in contr_csgto_carmom.F90 ..."
        call html_log_heading("Tests subroutines in contr_csgto_carmom.F90", IO_LOG, 3)
        call test_contr_csgto_carmom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in contr_csgto_carmom.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests calculating kinetic energy integrals with SGTOs
        num_tests = num_tests+1
        write(STDOUT,100) "tests calculating kinetic energy integrals with SGTOs ..."
        call html_log_heading("Tests calculating kinetic energy integrals with SGTOs", &
                              IO_LOG, 3)
        call test_contr_sgto_kinene(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests calculating kinetic energy integrals with SGTOs failed!"
          num_failed = num_failed+1
        end if
        ! tests calculating kinetic energy integrals with SGTOs
        num_tests = num_tests+1
        write(STDOUT,100) "tests calculating kinetic energy integrals with SGTOs ..."
        call html_log_heading("Tests calculating kinetic energy integrals with SGTOs", &
                              IO_LOG, 3)
        call test_contr_csgto_kinene(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests calculating kinetic energy integrals with SGTOs failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (6) tests Dirac delta function integrals
      if (test_delta) then
        ! tests subroutine delta_geom
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine delta_geom ..."
        call html_log_heading("Tests subroutine delta_geom", IO_LOG, 3)
        call test_delta_geom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine delta_geom failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in delta_hket.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in delta_hket.F90 ..."
        call html_log_heading("Tests subroutines in delta_hket.F90", IO_LOG, 3)
        call test_delta_hket(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in delta_hket.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in delta_moment.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in delta_moment.F90 ..."
        call html_log_heading("Tests subroutines in delta_moment.F90", IO_LOG, 3)
        call test_delta_moment(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in delta_moment.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests prim_hgto_delta.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine prim_hgto_delta ..."
        call html_log_heading("Tests subroutine prim_hgto_delta", IO_LOG, 3)
        call test_prim_hgto_delta(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine prim_hgto_delta failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in contr_cgto_delta.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in contr_cgto_delta.F90 ..."
        call html_log_heading("Tests subroutines in contr_cgto_delta.F90", IO_LOG, 3)
        call test_contr_cgto_delta(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in contr_cgto_delta.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in contr_sgto_delta.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in contr_sgto_delta.F90 ..."
        call html_log_heading("Tests subroutines in contr_sgto_delta.F90", IO_LOG, 3)
        call test_contr_sgto_delta(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in contr_sgto_delta.F90 failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (7) tests nuclear attraction potential integrals
      if (test_nucpot) then
        ! tests subroutines in aux_boys_vec.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in aux_boys_vec.F90 ..."
        call html_log_heading("Tests subroutines in aux_boys_vec.F90", IO_LOG, 3)
        call test_aux_boys_vec(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in aux_boys_vec.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests nucpot_geom.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in nucpot_geom.F90 ..."
        call html_log_heading("Tests subroutines in nucpot_geom.F90", IO_LOG, 3)
        call test_nucpot_geom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in nucpot_geom.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests nucpot_hket.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in nucpot_hket.F90 ..."
        call html_log_heading("Tests subroutines in nucpot_hket.F90", IO_LOG, 3)
        call test_nucpot_hket(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in nucpot_hket.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests nucpot_hbra.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in nucpot_hbra.F90 ..."
        call html_log_heading("Tests subroutines in nucpot_hbra.F90", IO_LOG, 3)
        call test_nucpot_hbra(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in nucpot_hbra.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests prim_hgto_nucpot.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine prim_hgto_nucpot ..."
        call html_log_heading("Tests subroutine prim_hgto_nucpot", IO_LOG, 3)
        call test_prim_hgto_nucpot(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine prim_hgto_nucpot failed!"
          num_failed = num_failed+1
        !! tests prim_hgto_nucpot.F90 using finite difference
        !else
        !  num_tests = num_tests+1
        !  write(STDOUT,100) "tests subroutine prim_hgto_nucpot using "// &
        !                    "finite difference ..."
        !  call html_log_heading("Tests subroutine prim_hgto_nucpot using "// &
        !                        "finite difference", IO_LOG, 3)
        !  call test_fd_prim_hgto_nucpot(IO_LOG, test_failed)
        !  if (test_failed) then
        !    write(STDOUT,100) "tests subroutine prim_hgto_nucpot using "// &
        !                      "finite difference failed!"
        !    num_failed = num_failed+1
        !  end if
        end if
        ! tests subroutines in contr_cgto_nucpot.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in contr_cgto_nucpot.F90 ..."
        call html_log_heading("Tests subroutines in contr_cgto_nucpot.F90", IO_LOG, 3)
        call test_contr_cgto_nucpot(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in contr_cgto_nucpot.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutines in contr_sgto_nucpot.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in contr_sgto_nucpot.F90 ..."
        call html_log_heading("Tests subroutines in contr_sgto_nucpot.F90", IO_LOG, 3)
        call test_contr_sgto_nucpot(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in contr_sgto_nucpot.F90 failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (8) tests Gaussian charge potential integrals
      if (test_gaupot) then
        ! tests subroutines in aux_boys_vec.F90
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutines in aux_boys_vec.F90 ..."
        call html_log_heading("Tests subroutines in aux_boys_vec.F90", IO_LOG, 3)
        call test_aux_boys_vec(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutines in aux_boys_vec.F90 failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutine gaupot_geom
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine gaupot_geom ..."
        call html_log_heading("Tests subroutine gaupot_geom", IO_LOG, 3)
        call test_gaupot_geom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine gaupot_geom failed!"
          num_failed = num_failed+1
        end if
      end if
      ! (9) tests subroutines related to overlap distribution
      if (test_odist) then
        ! tests subroutine prim_hgto_odist
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine prim_hgto_odist ..."
        call html_log_heading("Tests subroutine prim_hgto_odist", IO_LOG, 3)
        call test_prim_hgto_odist(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine prim_hgto_odist failed!"
          num_failed = num_failed+1
        end if
        ! tests subroutine prim_hgto_value
        num_tests = num_tests+1
        write(STDOUT,100) "tests subroutine prim_hgto_value ..."
        call html_log_heading("Tests subroutine prim_hgto_value", IO_LOG, 3)
        call test_prim_hgto_value(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests subroutine prim_hgto_value failed!"
          num_failed = num_failed+1
        end if
      end if
#if defined(BUILD_F90_MODULE)
      ! tests Fortran 90 modules
      if (test_f90mod) then
        write(STDOUT,100) "tests Fortran 90 modules ..."
        call html_log_heading("Tests Fortran 90 modules", IO_LOG, 3)
        ! tests kinetic energy integrals using SGTOs
        num_tests = num_tests+1
        write(STDOUT,100) "tests kinetic energy integrals using SGTOs ..."
        call html_log_heading("Tests kinetic energy integrals using SGTOs", &
                              IO_LOG, 4)
        call test_f90mod_sgto_kinene(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests kinetic energy integrals using SGTOs failed!"
          num_failed = num_failed+1
        end if
        ! tests Cartesian multipole integrals using SGTOs
        num_tests = num_tests+1
        write(STDOUT,100) "tests Cartesian multipole integrals using SGTOs ..."
        call html_log_heading("Tests Cartesian multipole integrals using SGTOs", &
                              IO_LOG, 4)
        call test_f90mod_sgto_carmom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests Cartesian multipole integrals using SGTOs failed!"
          num_failed = num_failed+1
        end if
        ! tests spherical multipole integrals using SGTOs
        num_tests = num_tests+1
        write(STDOUT,100) "tests spherical multipole integrals using SGTOs ..."
        call html_log_heading("Tests spherical multipole integrals using SGTOs", &
                              IO_LOG, 4)
        call test_f90mod_sgto_sphmom(IO_LOG, test_failed)
        if (test_failed) then
          write(STDOUT,100) "tests spherical multipole integrals using SGTOs failed!"
          num_failed = num_failed+1
        end if
      end if
#endif
      ! alll tests done, writes summary
      call html_log_horizontal_line(IO_LOG)
      call html_log_date_and_time(IO_LOG)
      ! gets the end time
      call xtimer_set(end_time)
      ! prints the CPU elapsed time
      call html_log_real_number(log_text="Total time (s) used by test suite :", &
                                log_real=end_time-begin_time, fmt_real="F10.4", &
                                io_log=IO_LOG, font_color="blue")
      write(STDOUT,100) "number of tests:", num_tests
      write(STDOUT,100) "number of failed tests:", num_failed
      call html_log_int_number("Number of tests:", num_tests, "I4", IO_LOG)
      if (num_failed==0) then
        call html_log_text("Congratulations! All tests passed! ;-)", IO_LOG)
        write(STDOUT,100) "all tests passed!"
        write(STDOUT,100) "you may check the log file for further details"
      else
        call html_log_int_number(log_text="Number of failed tests:", &
                                 log_int=num_failed, fmt_int="I4",   &
                                 io_log=IO_LOG, font_color="red")
        call html_log_href("Please report the errors, and  the information of "//  &
                           "your operating system and compilers to the author",    &
                           "mailto:bin.gao@uit.no", "Bin Gao", IO_LOG)
        write(STDOUT,100) "please check the log file for the detailed error message!"
        write(STDOUT,100) "error message is marked in red color!"
      end if
      call html_log_body_end(IO_LOG)
      call html_log_end(IO_LOG)
      close(IO_LOG)
    end if
100 format("test_gen1int>> ",A,I6)
  end program test_gen1int

  !> \brief prints help menu of test suite
  !> \author Bin Gao
  !> \date 2012-02-06
  subroutine test_help_menu()
    use xkind
    implicit none
    write(STDOUT,100) "Usage: ./test_gen1int -option"
    write(STDOUT,100) "Valid options are:"
    write(STDOUT,100) "  -all:    run all tests"
    write(STDOUT,100) "  -tools:  test tool subroutines"
    write(STDOUT,100) "  -basic:  test some basic subroutines"
    write(STDOUT,100) "  -geo:    test subroutines related to geometric derivatives"
    write(STDOUT,100) "  -mag:    test subroutines related to magnetic derivatives"
    write(STDOUT,100) "  -carmom: test Cartesian multipole moment integrals"
    write(STDOUT,100) "  -delta:  test Dirac delta function integrals"
    write(STDOUT,100) "  -nucpot: test nuclear attraction potential integrals"
    write(STDOUT,100) "  -gaupot: test Gaussian charge potential integrals"
    write(STDOUT,100) "  -odist:  test subroutines related to overlap distribution"
#if defined(BUILD_F90_MODULE)
    write(STDOUT,100) "  -f90mod: test Fortran 90 modules"
#endif
    write(STDOUT,100) "  -help:   print help menu of test suite"
    return
100 format(4X,A)
  end subroutine test_help_menu
