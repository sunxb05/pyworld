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
!!  This file contains subroutines for auxiliary functions.
!!
!!  2009-09-22, Bin Gao:
!!  * change subroutine "gen1int_aux_order" due to the modifications of types "g1prop",
!!    "g1int1" and "g1opt"
!!
!!  2009-08-24, Bin Gao:
!!  * finish subroutine of calculations of G function; however, our previous defined
!!    G function ( \int_0^1\exp(-T(1-u^2))(1-u^2)^{n}du ) does not have a stable upward
!!    recursion relation, so we have to give it up, and choose another G function
!!    \int_0^1\exp(-T(1-u^2))u^{2n}du, which is much easier to evaluate compared with
!!    our previous G function; but it requires different recursion relations for
!!    calculating integrals
!!  * finish the tests of the accuray of auxiliary functions, *NOTE* that the referenced
!!    results may not be true sometimes, like for G function, its downward recursion
!!    relation is not stable for large argument, so if the tests involve higher order
!!    G function, the referenced result of 0th order is not always accurate
!!
!!  2009-08-23, Bin Gao:
!!  * change subroutine "gen1int_aux_order" with operators and geometric derivatives
!!    as arguments
!!  * add "max_order_aux" (maximum order of auxiliary functions)
!!  * finish subroutine of calculations of Boys function
!!  * add tests of the accuray of auxiliary functions (not finished)
!!
!!  2009-08-19, Bin Gao:
!!  * add determination of orders of auxiliary functions from the *input* generalized
!!    one-electron operator and/or *calculated* property integrals in subroutine
!!    gen1int_aux_order
!!
!!  2009-08-17, Bin Gao:
!!  * make changes because the modification of *input* generalized one-electron operator
!!
!!  2009-08-10, Bin Gao:
!!  * move the definitions of variables related to auxiliary functions from gen1int_data.F90.
!!
!!  2009-08-03, Bin Gao:
!!  * make changes according to the modifications of operators.
!!
!!  2009-07-23, Bin Gao:
!!  * remove the determination of orders of auxiliary functions from the
!!    property integrals in subroutine gen1int_aux_order, since this
!!    is taken care in gen1int_input_process.
!!
!!  2009-06-27, Bin Gao:
!!  * first version

module gen1int_rfunction
  ! precision
  use xprecision
  ! some constants
  use xconst
  ! some basic variables & subroutines
  use xsystem
  ! basis sets
  use xbasis, only : max_angular_num
  ! some basic variables of GEN1INT
  use gen1int_data
  implicit none
  ! maximum order of auxiliary functions, for safety, and also our code may break down
  ! for large order and argument's auxiliary functions
  integer, parameter :: max_order_aux = 24
  ! number of terms of Taylor expansion for the tabulated auxiliary functions
  integer, parameter :: num_taylor_boys = 6, num_taylor_gfun = 6
  ! NB: Please do not change the following variables, otherwise, you have to
  !     update gen1int_aux_boys & gen1int_aux_gfun !!
  ! number of sampling points for the tabulated auxiliary functions
  integer, parameter :: num_sample_boys = 120, num_sample_gfun = 120
  ! we now use sampling points from 0 to 12, interval is 0.1
  real(xp), parameter :: sample0_boys = zero, interval_boys = tenth, &
                         sample0_gfun = zero, interval_gfun = tenth
  ! the last sampling points
  real(xp), parameter :: sample1_boys = sample0_boys+interval_boys*num_sample_boys, &
                         sample1_gfun = sample0_gfun+interval_gfun*num_sample_gfun
  ! tabulated auxiliary functions
  real(xp), allocatable :: g1_boys_tab(:,:), g1_gfun_tab(:,:)
  ! maximum argument for power series expansions
  !
  ! NOTE!!!!
  ! 1. For Boys function F_n(T), according to our tests, the current asymptotic
  !    series expansions (gen1int_boys_asymp) can achieve a 10^-12 relative
  !    accuracy after T>25 (BUT for quite smaller n if T is not large enough),
  !    we may need to provide a new asymptotic series expansions later ...;
  !    the modified asymptotic series expansions do not have a stable upward
  !    recursion relations when n>=2T (n>=T/2 for 45.0<=T<=100.0, but the
  !    value of Boys function is quite small) according to our tests (negtive
  !    values may occur); but for the time being, it seems to be enough;
  !    last, for the power series expansions, we have used downward recursion
  !    relation, which seems to be very stable according to our tests (0<=T<=400,
  !    n=10, 1000), even the value of largest order is not correct (for T<24.0
  !    it accurate enough for n=1000; and it is accurate enough for T<28.0, n=100;
  !    according to our tests, for larger T, it may become reasonable after n<2T)
  ! 2. For G function G_n(T), since we have used downward recursion relation after
  !    the power series expansions of the largest order N, the downward recursion
  !    relation is not stable for large T; according to our tests, it is stable for
  !    even 1000th order with T<20.0 (indeed, it becomes better for larger N), but
  !    it is unstable when T>=20.0, even the value of the largest order is still
  !    correct (with an appropriate "max_pterms_gfun", it is accurate enough even
  !    for n=1000, T=400.0); the asymptotic series expansions seems to only work
  !    for T>30.0, but it, as well as the rational Chebyshev approximation use
  !    upward recursion relations, which is not stable for large n, according to
  !    our tests, it seems that n must satisfy n<2T+a (a=10 for T<=40.0; a may become
  !    smaller for larger T); for the time being, it seems to be enough, since we use
  !    rational Chebyshev approximation after T>=12, that means we can arrive at least
  !    24th order for T==12; however, if higher orders are required, a non-efficient
  !    way may use power series expansions for each order (oh my god ;-) ), well,
  !    we may need to improve it larter ...
  real(xp), parameter :: max_arg_pow_boys = 30.0D0
  real(xp), parameter :: max_arg_pow_gfun = 30.0D0
  ! maximum terms for the power series expansions
!fixme: how to determine max_pterms, which depends on the order & argument
!fixme: 200 terms are enough for arguments smaller than, like 20
    integer, parameter :: max_pterms_boys = 200
    integer, parameter :: max_pterms_gfun = 200
  ! cut of for auxiliary function
  real(xp), parameter :: g1_cut_aux = 1.0D-15

  ! At present, this module contains the following public subroutines,
  ! 1. gen1int_aux_tabulate: tabulate the auxiliary functions
  ! 2. gen1int_aux_test: test the accuracy of calcualtions of auxiliary functions
  ! 3. gen1int_aux_order: calculate the maximum order of auxiliary functions
  !    from the given property, and/or from *calcualted* pre-defined properties
  ! 4. gen1int_aux_boys: calculate Boys function F_N(T) using (1) Taylor expansions
  !    and downward recursion relation for small T, and (2) modified asymptotic series
  !    and upward recursion relation for large T
  ! 5. gen1int_aux_gfun: calculate G function G_N(T) for operator 1/r^2, using
  !    (1) Taylor expansions and downward recursion relation for small T, and
  !    (2) rational Chebyshev approximation and upward recursion relation for
  !    large T
  !
  ! and the following private subroutines,
  ! 1. gen1int_boys_power: power series expansions of Boys function
  ! 2. gen1int_boys_asymp: asymptotic series expansions of Boys funciton
  ! 3. gen1int_gfun_power: power series expansions of G function
  ! 4. gen1int_gfun_asymp: asymptotic series expansions of G function
  public :: gen1int_aux_tabulate, &
            gen1int_aux_test,     &
            gen1int_aux_order,    &
            gen1int_aux_boys,     &
            gen1int_aux_gfun
  private :: gen1int_boys_power, &
             gen1int_boys_asymp, &
             gen1int_gfun_power, &
             gen1int_gfun_asymp

  contains

  ! tabulate the auxiliary functions
  subroutine gen1int_aux_tabulate
    implicit none
    ! push current subroutine into the stack
    call xsub_enter('gen1int_aux_tabulate')
    ! check if we got max_angular_num
    if ( max_angular_num.lt.0 ) then
      write(x_lupri,1000)'Maximum angular number: ',max_angular_num
      call xstop('Wrong maximum angular number!')
    end if
    ! determine the maximum orders of the auxiliary functions from
    ! *calculated* property integrals, if they have not been done
    if ( (g1_boys_used .and. op1_order_boys.lt.0) .or. &
         (g1_gfun_used .and. op1_order_gfun.lt.0) )    &
      call gen1int_aux_order( check_prop = .true. )
    ! now we get the maximum order of the auxiliary functions
!fixme: maximum orders depend on op1_orders, max_angular_num, and the initial accuracy of F_N(T);
!fixme: it is used to increase the accuracy (should be less than g1_cut_aux) after downward recursions
!fixme: so, if max(2,2*max_angular_num) ok??
    if (g1_boys_used) then
      max_order_boys = op1_order_boys + max(2,2*max_angular_num)
      ! check order
      if ( max_order_boys.gt.max_order_aux ) then
        write(x_lupri,1000)'Maximum order of auxiliary functions: ',max_order_aux
        write(x_lupri,1000)'Order of Boys functions:              ',max_order_boys
        write(x_lupri,1000)'Be careful! The results of high order may be unstable!!'
        call xstop('Increase max_order_aux in gen1int_aux.F90')
      end if
    end if
    if (g1_gfun_used) then
      max_order_gfun = op1_order_gfun + max(2,2*max_angular_num)
      ! check order
      if ( max_order_gfun.gt.max_order_aux ) then
        write(x_lupri,1000)'Maximum order of auxiliary functions: ',max_order_aux
        write(x_lupri,1000)'Order of G functions:                 ',max_order_gfun
        write(x_lupri,1000)'Be careful! The results of high order may be unstable!!'
        call xstop('Increase max_order_aux in gen1int_aux.F90')
      end if
    end if
    ! dump
    if (xlib_lprint.ge.10) then
      if (g1_boys_used) write(x_lupri,1000)'Maximum order of Boys function: ',max_order_boys
      if (g1_gfun_used) write(x_lupri,1000)'Maximum order of G function:    ',max_order_gfun
      write(x_lupri,'()')
    end if
    ! tabulate Boys function
    if (g1_boys_used) then
      allocate( g1_boys_tab(0:max_order_boys+num_taylor_boys, &
                            0:num_sample_boys), stat=x_ierr)
      if ( x_ierr.ne.0 ) call xstop('allocate g1_boys_tab!')
      call gen1int_boys_power( 0, max_order_boys+num_taylor_boys, &
                               sample0_boys, interval_boys, num_sample_boys, g1_boys_tab )
    end if
    ! tabulate G function
    if (g1_gfun_used) then
      allocate( g1_gfun_tab(0:max_order_gfun+num_taylor_gfun, &
                            0:num_sample_gfun), stat=x_ierr)
      if ( x_ierr.ne.0 ) call xstop('allocate g1_gfun_tab!')
      call gen1int_gfun_power( 0, max_order_gfun+num_taylor_gfun, &
                               sample0_gfun, interval_gfun, num_sample_gfun, g1_gfun_tab )
    end if
    ! pop the stack
    call xsub_exit
    return
1000  format('ATAB ',A,I6)
  end subroutine gen1int_aux_tabulate

  ! test the accuracy of calcualtions of auxiliary functions
  subroutine gen1int_aux_test
!fixme: think a much better referenced results, some software could provide infinite accuracy
    implicit none
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! check orders
    integer, parameter :: my_order = 50
    ! number of arguments
    integer, parameter :: my_num_arg = 1000
    ! incremental recorder for orders
    integer iord
    ! incremental recorder for arguments
    integer iarg
    ! arguments
    real(xp) my_arg
    ! results from power series expansions or asymptotic series expansions
    real(xp), allocatable :: my_refs(:,:)
    ! results from subroutines used in the calculations
    real(xp), allocatable :: my_chks(:)
    ! error
    real(xp) my_err
    ! push current subroutine into the stack
    call xsub_enter('gen1int_aux_test')
    write(x_lupri,1010)'Threshold of test: ',g1_cut_aux
    ! results from power series expansions or asymptotic series expansions
    allocate(my_refs(0:my_order,0:my_num_arg), stat=x_ierr)
    if ( x_ierr.ne.0 ) call xstop('allocate my_refs!')
    ! results from subroutines used in the calculations
    allocate(my_chks(0:my_order), stat=x_ierr)
    if ( x_ierr.ne.0 ) call xstop('allocate my_chks!')
    ! tabulate Boys function
    if ( allocated(g1_boys_tab) ) deallocate(g1_boys_tab)
    allocate(g1_boys_tab(0:my_order+num_taylor_boys,0:num_sample_boys), stat=x_ierr)
    if ( x_ierr.ne.0 ) call xstop('allocate g1_boys_tab!')
    call gen1int_boys_power( 0, my_order+num_taylor_boys, &
                             sample0_boys, interval_boys, num_sample_boys, g1_boys_tab )
    ! test Boys functions
    write(x_lupri,1000)
    write(x_lupri,1000)'Test the accuracy of Boys functions ...'
    write(x_lupri,1000)'Order  Argument      Result             Reference      log10(relative error)'
    ! we use power series expansions and asymptotic series expansions as the referenced
    ! results, but we should know that they may be not accurate enough sometimes!!
    call gen1int_boys_power( 0, my_order, zero, one, my_num_arg, my_refs )
    ! initialize argument
    my_arg = zero
    ! loop over arguments
    do iarg = 0, my_num_arg
      ! results from subroutines used in the calculations
      call gen1int_aux_boys( my_order, my_arg, my_chks )
      ! compare them
      do iord = 0, my_order
        my_err = abs( (my_chks(iord)-my_refs(iord,iarg))/my_refs(iord,iarg) )
        ! Boys function is the monotonically decreasing function of T
        ! F_n(0) = 1/(2n+1) <= 1
        if ( my_chks(iord).lt.zero .or. my_refs(iord,iarg).lt.zero .or. &
             my_chks(iord).gt.one  .or. my_refs(iord,iarg).gt.one  .or. &
             my_err.gt.g1_cut_aux ) then
          write(x_lupri,1020) iord,my_arg,my_chks(iord),my_refs(iord,iarg),log10(my_err)
        end if
      end do
      ! increase argument
      my_arg = my_arg + one
    end do
    ! clean tabulated functions
    deallocate(g1_boys_tab)
    ! tabulate G function
    if ( allocated(g1_gfun_tab) ) deallocate(g1_gfun_tab)
    allocate(g1_gfun_tab(0:my_order+num_taylor_gfun,0:num_sample_gfun), stat=x_ierr)
    if ( x_ierr.ne.0 ) call xstop('allocate g1_gfun_tab!')
    call gen1int_gfun_power( 0, my_order+num_taylor_gfun, &
                             sample0_gfun, interval_gfun, num_sample_gfun, g1_gfun_tab )
    ! test G functions
    write(x_lupri,1000)
    write(x_lupri,1000)'Test the accuracy of G functions ...'
    write(x_lupri,1000)'Order  Argument      Result             Reference      log10(relative error)'
    ! we use power series expansions and asymptotic series expansions as the referenced
    ! results, but we should know that they may be not accurate enough sometimes!!
    call gen1int_gfun_power( 0, my_order, zero, one, my_num_arg, my_refs )
    ! initialize argument
    my_arg = zero
    ! loop over arguments
    do iarg = 0, my_num_arg
      ! results from subroutines used in the calculations
      call gen1int_aux_gfun( my_order, my_arg, my_chks )
      ! compare them
      do iord = 0, my_order
        my_err = abs( (my_chks(iord)-my_refs(iord,iarg))/my_refs(iord,iarg) )
        ! G function is the monotonically decreasing function of T
        ! G_n(0) = 1/(2n+1) <= 1
        if ( my_chks(iord).lt.zero .or. my_refs(iord,iarg).lt.zero .or. &
             my_chks(iord).gt.one  .or. my_refs(iord,iarg).gt.one  .or. &
             my_err.gt.g1_cut_aux ) then
          write(x_lupri,1020) iord,my_arg,my_chks(iord),my_refs(iord,iarg),log10(my_err)
        end if
      end do
      ! increase argument
      my_arg = my_arg + one
    end do
    ! clean test results
    deallocate(my_refs)
    deallocate(my_chks)
    ! clean tabulated functions
    deallocate(g1_gfun_tab)
    write(x_lupri,'()')
    ! tabulate auxiliary functions again
    call gen1int_aux_tabulate
    ! pop the stack
    call xsub_exit
    return
1000  format('AXTT ',A)
1010  format('AXTT ',A,E20.12)
1020  format('AXTT ',I4,F10.3,2E20.12,F10.3)
  end subroutine gen1int_aux_test

  ! calculate the maximum order of auxiliary functions from the given property,
  ! and/or from *calcualted* pre-defined properties
  subroutine gen1int_aux_order( aux_prop, check_prop )
    ! property integrals
    use gen1int_property
    ! the maximum order of auxiliary function should be the sum of
    ! (1) powers of multipole moments,
    ! (2) orders of derivatives of 1/r_C^m (m=1,2) function,
    ! (3) orders of derivatives with respect to electrons,
    ! (4) orders of geometric derivatives
    implicit none
    ! property
    type(g1prop), optional, intent(in) :: aux_prop
    ! calcualte the order of auxiliary functions from *calcualted* pre-defined properties
    logical, optional, intent(in) :: check_prop
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! temporary order of auxiliary functions
    integer tmp_order
    ! incremental recorder over integrals
    integer int1
    ! incremental recorder over operators and *calculated* pre-defined properties
    integer iop1
    ! incremental recorder
    integer ior1
    ! temporary stuff for the type of r_C function
    character chk_prop_rfun
    ! push current subroutine into the stack
    call xsub_enter('gen1int_aux_order')
    ! calcualte the order from given operators
    if ( present(aux_prop) ) then
      if ( allocated(aux_prop%range_g1int) ) then
        ! loop over integrals
        do int1 = aux_prop%range_g1int(1), aux_prop%range_g1int(2)
          ! loop over operators
          do iop1 = 1, aux_prop%g1prop_int1(int1)%num_g1opt
!fixme: consider (1/r_K)*(1/r_L)
            select case ( aux_prop%g1prop_int1(int1)%g1int_opt(iop1)%r_func(1) )
            ! 1/r_C
            case ('1')
              g1_boys_used = .true.
              ! orders of geometric derivatives
              if ( allocated(aux_prop%order_gdcent) ) then
                tmp_order = max( sum(aux_prop%order_gdcent), &
                                 aux_prop%order_gderv )
              else
                tmp_order = aux_prop%order_gderv
              end if
              ! orders of derivatives with respect to electrons
              tmp_order = tmp_order &
                        + aux_prop%g1prop_int1(int1)%g1int_opt(iop1)%ederv_order(4)
              do ior1 = 1, 3
                ! powers of multipole moments
                tmp_order = tmp_order &
                          + aux_prop%g1prop_int1(int1)%g1int_opt(iop1)%mult_power(ior1)
                ! orders of derivatives of 1/r_C^m (m=1,2) function
!fixme: consider (1/r_K)*(1/r_L)
                tmp_order = tmp_order &
                          + aux_prop%g1prop_int1(int1)%g1int_opt(iop1)%rderv_order(ior1,1)
              end do
              ! update order of auxiliary funtions
              if ( tmp_order.gt.op1_order_boys ) op1_order_boys = tmp_order
            ! 1/r_C^2
            case ('2')
              g1_gfun_used = .true.
              ! orders of geometric derivatives
              if ( allocated(aux_prop%order_gdcent) ) then
                tmp_order = max( sum(aux_prop%order_gdcent), &
                                 aux_prop%order_gderv )
              else
                tmp_order = aux_prop%order_gderv
              end if
              ! orders of derivatives with respect to electrons
              tmp_order = tmp_order &
                        + aux_prop%g1prop_int1(int1)%g1int_opt(iop1)%ederv_order(4)
              do ior1 = 1, 3
                ! powers of multipole moments
                tmp_order = tmp_order &
                          + aux_prop%g1prop_int1(int1)%g1int_opt(iop1)%mult_power(ior1)
                ! orders of derivatives of 1/r_C^m (m=1,2) function
!fixme: consider (1/r_K)*(1/r_L)
                tmp_order = tmp_order &
                          + aux_prop%g1prop_int1(int1)%g1int_opt(iop1)%rderv_order(ior1,1)
              end do
              ! update order of auxiliary funtions
              if ( tmp_order.gt.op1_order_gfun ) op1_order_gfun = tmp_order
            end select
          end do ! loop over operators
        end do ! loop over integrals
      end if
    end if
    ! calcualte the order from *calcualted* pre-defined properties
    if ( present(check_prop) ) then
      if ( check_prop ) then
        do iop1 = 1, num_calc_prop
          ! get the order from the operators of this *calcualted* pre-defined property
          call gen1int_proper_create( prop_name = def_prop(iop1)%prop_name, &
                                      prop_rfun = chk_prop_rfun,            &
                                      prop_aux_order = tmp_order )
          ! type of r_C function
          select case (chk_prop_rfun)
          case ('1')
            g1_boys_used = .true.
            ! consider the order of geometric derivatives
            if ( allocated(def_prop(iop1)%order_gdcent) ) then
              tmp_order = tmp_order &
                + max( sum(def_prop(iop1)%order_gdcent), def_prop(iop1)%order_gderv )
            else
              tmp_order = tmp_order + def_prop(iop1)%order_gderv
            end if
            ! update order of auxiliary funtions
            if ( tmp_order.gt.op1_order_boys ) op1_order_boys = tmp_order
          case ('2')
            g1_gfun_used = .true.
            ! consider the order of geometric derivatives
            if ( allocated(def_prop(iop1)%order_gdcent) ) then
              tmp_order = tmp_order &
                + max( sum(def_prop(iop1)%order_gdcent), def_prop(iop1)%order_gderv )
            else
              tmp_order = tmp_order + def_prop(iop1)%order_gderv
            end if
            ! update order of auxiliary funtions
            if ( tmp_order.gt.op1_order_gfun ) op1_order_gfun = tmp_order
          end select
        end do
      end if
    end if
    ! pop the stack
    call xsub_exit
    return
1000  format('AXOR ',A)
  end subroutine gen1int_aux_order

  ! calculate Boys function F_N(T) using (1) Taylor expansions and downward
  ! recursion relation for small T, and (2) modified asymptotic series and
  ! upward recursion relation for large T
  subroutine gen1int_aux_boys(order_boys,arg_boys,val_boys)
    implicit none
    ! order of Boys function
    integer, intent(in) :: order_boys
    ! argument of Boys function
    real(xp), intent(in) :: arg_boys
    ! value of Boys function
    real(xp), intent(out) :: val_boys(0:order_boys)
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! recorder for the argument
    integer iarg
    ! recorder for the orders
    integer iord
    ! (T*-T)
    real(xp) my_dif_arg
    ! twice of the argument
    real(xp) my_arg2
    ! F_J(T)
    real(xp) val_boys_max
    ! exp(-T)
    real(xp) my_pre_arg
    ! 1/T
    real(xp) my_iarg
    ! 1/(2T)
    real(xp) my_i2arg
    ! divisor in the downward recursion relations
    real(xp) div_down
    ! four term modified asymptotic series
    real(xp), parameter :: g_asym3(0:3) = (/0.4999489092D0,-0.2473631686D0, &
                                            0.321180909D0,-0.3811559346D0/)
    ! check if we got the tabulated auxiliary functions
    if ( g1_boys_used .and. .not. allocated(g1_boys_tab) ) &
      call xstop('Boys function is not tabulated!')
    ! 0 <= T < 12: using Taylor expansion and downward recursion relation
    if ( arg_boys.lt.sample1_boys ) then
      ! the nearest T*, NB: since sample0_boys == 0!
      iarg = nint(arg_boys/interval_boys)
      ! T*-T
      my_dif_arg = iarg*interval_boys - arg_boys
      ! (T*-T)^k/k!, k==0
      my_pre_arg = one
      ! the 1st Taylor's term
!N      val_boys_max = g1_boys_tab(max_order_boys,iarg)
      val_boys(order_boys) = g1_boys_tab(order_boys,iarg)
      ! the left Taylor's terms
      do iord = 1, num_taylor_boys
        my_pre_arg = my_pre_arg*my_dif_arg/xfloat(iord)
!N        val_boys_max = val_boys_max + my_pre_arg*g1_boys_tab(iord+max_order_boys,iarg)
        val_boys(order_boys) = val_boys(order_boys) + my_pre_arg*g1_boys_tab(iord+order_boys,iarg)
      end do
      ! exp(-T)
      my_pre_arg = exp(-arg_boys)
      ! 2*T
      my_arg2 = arg_boys + arg_boys
!N      ! use downward recursion relations to get F_j(T) from F_J(T) (val_boys_max)
!N      div_down = xfloat(max_order_boys+max_order_boys+1)
      div_down = xfloat(order_boys+order_boys+1)
!N      do iord = max_order_boys-1, order_boys, -1
!N        div_down = div_down - two
!N        val_boys(order_boys) = (my_arg2*val_boys_max+my_pre_arg)/div_down
!N        val_boys_max = val_boys(order_boys)
!N      end do
      ! use downward recursion relations for other F_j(T)
      do iord = order_boys-1, 0, -1
        div_down = div_down - two
        val_boys(iord) = (my_arg2*val_boys(iord+1)+my_pre_arg)/div_down
      end do ! loop over orders
    ! T >= 12: using modified asymptotic series expansions and upward recursion relation
    ! for T <= 2J+36, we use four asymptotic terms, and exact upward recursion relation
!fixme: this modified asymptotic series expansion seems to be bad for larger order, like 50th,
!fixme: with smaller arguments, like 12.0, 13.0, ...
    else if ( arg_boys.le.xfloat(max_order_boys+max_order_boys+36) ) then
      if ( order_boys.ge.floor(arg_boys+arg_boys) ) then
        write(x_lupri,1000)'The maximum order of Boys functions: ',order_boys
        write(x_lupri,1010)'The argument of Boys functions:      ',arg_boys
        write(x_lupri,1000)'The upward recursion relation is unstable for them!'
#ifndef TEST_AUX
        call xstop('Write to author asking them to improve it!')
#endif
      end if
      ! 1/T
      my_iarg = one/arg_boys
      ! exp(-T)
      my_pre_arg = exp(-arg_boys)
      ! g
      val_boys_max = g_asym3(0) + my_iarg*( g_asym3(1) + my_iarg*( g_asym3(2) + my_iarg*g_asym3(3) ) )
      ! F_0(T) = sqrt(pi/T)/2 - exp(-T)*g/T
      val_boys(0) = half*sqrt(pi*my_iarg) - my_pre_arg*val_boys_max*my_iarg
      ! 1/(2T)
      my_i2arg = half*my_iarg
      ! exp(-T)/(2T)
      my_pre_arg = my_i2arg*my_pre_arg
      ! upward recursion relation for other F_j(T)
      do iord = 1, order_boys
        val_boys(iord) = my_i2arg*val_boys(iord-1) - my_pre_arg
        my_i2arg = my_iarg + my_i2arg
      end do
    ! for T > 2J+36, we use F_0(T) = sqrt(pi/T)/2, and approximate upward recursion relation
    ! F_{j+1}(T) = (2T)^{-1}(2j+1)F_j(T)
    else
      if ( order_boys.ge.floor(half*arg_boys) ) then
        write(x_lupri,1000)'The maximum order of Boys functions: ',order_boys
        write(x_lupri,1010)'The argument of Boys functions:      ',arg_boys
        write(x_lupri,1000)'The upward recursion relation is unstable for them!'
#ifndef TEST_AUX
        call xstop('Write to author asking them to improve it!')
#endif
      end if
      ! 1/T
      my_iarg = one/arg_boys
      val_boys(0) = half*sqrt(pi*my_iarg)
      ! 1/(2T)
      my_i2arg = half*my_iarg
      ! upward approximate recursion relation for other F_j(T)
      do iord = 1, order_boys
        val_boys(iord) = my_i2arg*val_boys(iord-1)
        my_i2arg = my_iarg + my_i2arg
      end do
    end if
    return
1000  format('BOYS ',A,I6)
1010  format('BOYS ',A,F14.8)
  end subroutine gen1int_aux_boys

  ! calculate G function G_N(T) for operator 1/r^2, using (1) Taylor expansions
  ! and downward recursion relation for small T, and (2) rational Chebyshev
  ! approximation and upward recursion relation for large T
  subroutine gen1int_aux_gfun( order_gfun, arg_gfun, val_gfun )
    implicit none
    ! order of G function
    integer, intent(in) :: order_gfun
    ! argument of G function
    real(xp), intent(in) :: arg_gfun
    ! values of G function
    real(xp), intent(out) :: val_gfun(0:order_gfun)
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! recorder for the argument
    integer iarg
    ! recorder for the orders
    integer iord
    ! (T*-T)
    real(xp) my_dif_arg
    ! twice of the argument
    real(xp) my_arg2
    ! 1/T
    real(xp) my_iarg
    ! 1/(2T)
    real(xp) my_i2arg
    ! k/T
    real(xp) my_pre_arg
    ! divisor in the downward recursion relations
    real(xp) div_down
    ! coefficients of rational Chebyshev approximation of 0th order G function
    real(xp), parameter :: cheb0p1(0:5) = (/ 0.1828591449264953D0, &
                                            -0.0208111908526560D0, &
                                             0.1333867391625493D0, &
                                            -0.0069563805854603D0, &
                                            -0.0006399406464297D0, &
                                             0.0002262208851365D0/)
    real(xp), parameter :: cheb0q1(1:6) = (/-0.5001558049033900D0, &
                                            -0.0422332795713542D0, &
                                             0.2634740589566963D0, &
                                            -0.0130363366131583D0, &
                                            -0.0015170455558345D0, &
                                             0.0004525493560858D0/)
    real(xp), parameter :: cheb0p2(0:4) = (/-0.0237209201842862D0, & 
                                             0.0369293462259353D0, &
                                            -0.1782265603862903D0, &
                                             0.0305190498240735D0, &
                                            -0.0019473305442354D0/)
    real(xp), parameter :: cheb0q2(1:5) = (/-0.0929707538344341D0, &
                                             0.2362992390268955D0, &
                                            -0.3854264829228933D0, &
                                             0.0629939171730423D0, &
                                            -0.0038947347910983D0/)
    real(xp), parameter :: cheb0p3(0:2) = (/-0.7973613445962394D0, &
                                            -0.1763813954302433D0, &
                                             0.0562314564328704D0/)
    real(xp), parameter :: cheb0q3(1:4) = (/-1.4886748949591277D0, &
                                            -0.4084032309951818D0, &
                                             0.1124504768032380D0, &
                                             0.0000001053836247D0/)
    real(xp), parameter :: cheb0p4(0:2) = (/-1.6901358291748010D0, &
                                             0.9095465039538543D0, &
                                            -0.0821332073172411D0/)
    real(xp), parameter :: cheb0q4(1:3) = (/-4.2082066572381995D0, &
                                             1.9012356348172637D0, &
                                            -0.1642664818611426D0/)
    real(xp), parameter :: cheb0p5(0:2) = (/-1.7168354719492460D0, &
                                             0.9510390010759011D0, &
                                            -0.0888365521252258D0/)
    real(xp), parameter :: cheb0q5(1:3) = (/-4.2963007179996406D0, &
                                             1.9909219459253031D0, &
                                            -0.1776731550443674D0/)
    real(xp), parameter :: cheb0p6(0:2) = (/-1.8429493331531948D0, &
                                             1.1598159456555421D0, &
                                            -0.1246514775701103D0/)
    real(xp), parameter :: cheb0q6(1:3) = (/-4.7212032569199138D0, &
                                             2.4442852443824989D0, &
                                            -0.2493029652188083D0/)
    real(xp), parameter :: cheb0p7(0:2) = (/-1.9744562443420302D0, &
                                             1.3986614326037992D0, &
                                            -0.1691251258024091D0/)
    real(xp), parameter :: cheb0q7(1:3) = (/-5.1784724951222785D0, &
                                             2.9664481966513980D0, &
                                            -0.3382502523273182D0/)
    real(xp), parameter :: cheb0p8(0:2) = (/-2.1099021081579972D0, &
                                             1.6701190095543932D0, &
                                            -0.2240542865988810D0/)
    real(xp), parameter :: cheb0q8(1:3) = (/-5.6658694074754621D0, &
                                             3.5642923069732153D0, &
                                            -0.4481085731988854D0/)
    real(xp), parameter :: cheb0p9(0:1) = (/-1.6616791513685247D0, &
                                             0.6616886439852467D0/)
    real(xp), parameter :: cheb0q9(1:2) = (/-3.9850469515266829D0, &
                                             1.3233772879713688D0/)
    ! check if we got the tabulated auxiliary functions
    if ( g1_gfun_used .and. .not. allocated(g1_gfun_tab) ) &
      call xstop('G function is not tabulated!')
    ! 0 <= T < 12: using Taylor expansions for the largest order and downward
    ! recursion relations for others
    if ( arg_gfun.lt.sample1_gfun ) then
      ! the nearest T*, NB: since sample0_gfun == 0!
      iarg = nint(arg_gfun/interval_gfun)
      ! T*-T
      my_dif_arg = iarg*interval_gfun - arg_gfun
      ! (T*-T)^k/k!, k==0
      my_iarg = one
      ! the 1st Taylor's term
      val_gfun(order_gfun) = g1_gfun_tab(order_gfun,iarg)
      ! the left Taylor's terms
      do iord = 1, num_taylor_gfun
        my_iarg = my_iarg*my_dif_arg/xfloat(iord)
        val_gfun(order_gfun) = val_gfun(order_gfun) + my_iarg*g1_gfun_tab(iord+order_gfun,iarg)
      end do
      ! left orders using downward recursions
      ! 2*T
      my_arg2 = arg_gfun + arg_gfun
      ! use downward recursion relations for other orders
      div_down = xfloat(order_gfun+order_gfun+1)
      do iord = order_gfun-1, 0, -1
        div_down = div_down - two
        val_gfun(iord) = (one-my_arg2*val_gfun(iord+1))/div_down
      end do ! loop over orders
    ! T >= 12: using rational Chebyshev approximation for the 0th order
    ! and upward recursion relations for others
!fixme: the upward recursion relation is unstable for larger order (seems to be >=2T+10)
!fixme: well, it may not be exactly 2T+10, but depends on T (larger T smaller n)
    else
      if ( order_gfun.ge.floor(arg_gfun+arg_gfun)+10 ) then
        write(x_lupri,1000)'The maximum order of G functions: ',order_gfun
        write(x_lupri,1010)'The argument of G functions:      ',arg_gfun
        write(x_lupri,1000)'The upward recursion relation is unstable for them!'
#ifndef TEST_AUX
        call xstop('Write to author asking them to improve it!')
#endif
      end if
!fixme: we may use asymptotic series expansions for very large T, like > 10,000
      !-----------
      ! 0th order
      !-----------
      ! 12 <= T <= 15.8 (fitted from 12 <= T <= 16)
      if ( arg_gfun.le. 15.8D0 ) then
        my_iarg = cheb0p1(0) + arg_gfun*( cheb0p1(1) + arg_gfun*( cheb0p1(2) &
                + arg_gfun*( cheb0p1(3) + arg_gfun*( cheb0p1(4) + arg_gfun*cheb0p1(5) ) ) ) )
        my_i2arg = one + arg_gfun*( cheb0q1(1) + arg_gfun*( cheb0q1(2) &
                 + arg_gfun*( cheb0q1(3) + arg_gfun*( cheb0q1(4) + arg_gfun*( cheb0q1(5) &
                 + arg_gfun*cheb0q1(6) ) ) ) ) )
      ! 15.8 < T <= 18.9 (fitted from 16 <= T <= 19)
      else if ( arg_gfun.le.18.9D0 ) then
        my_iarg = cheb0p2(0) + arg_gfun*( cheb0p2(1) + arg_gfun*( cheb0p2(2) &
                + arg_gfun*( cheb0p2(3) + arg_gfun*cheb0p2(4) ) ) )
        my_i2arg = one + arg_gfun*( cheb0q2(1) + arg_gfun*( cheb0q2(2) &
                 + arg_gfun*( cheb0q2(3) + arg_gfun*( cheb0q2(4) + arg_gfun*cheb0q2(5) ) ) ) )
      ! 18.9 < T <= 21.9 (fitted from 19 <= T <= 22)
      else if ( arg_gfun.le.21.9D0 ) then
        my_iarg = cheb0p3(0) + arg_gfun*( cheb0p3(1) + arg_gfun*cheb0p3(2) )
        my_i2arg = one + arg_gfun*( cheb0q3(1) + arg_gfun*( cheb0q3(2) &
                 + arg_gfun*( cheb0q3(3) + arg_gfun*cheb0q3(4) ) ) )
      ! 21.9 < T <= 25 (fitted from 22 <= T <= 25)
      else if ( arg_gfun.le.25.0D0) then
        my_iarg = cheb0p4(0) + arg_gfun*( cheb0p4(1) + arg_gfun*cheb0p4(2) )
        my_i2arg = one + arg_gfun*( cheb0q4(1) + arg_gfun*( cheb0q4(2) + arg_gfun*cheb0q4(3) ) )
      ! 25 < T <= 30 (fitted from 25 <= T <= 30)
      else if ( arg_gfun.le.30.0D0 ) then
        my_iarg = cheb0p5(0) + arg_gfun*( cheb0p5(1) + arg_gfun*cheb0p5(2) )
        my_i2arg = one + arg_gfun*( cheb0q5(1) + arg_gfun*( cheb0q5(2) + arg_gfun*cheb0q5(3) ) )
      ! 30 < T <= 39 (fitted from 30 <= T <= 39)
      else if ( arg_gfun.le.39.0D0 ) then
        my_iarg = cheb0p6(0) + arg_gfun*( cheb0p6(1) + arg_gfun*cheb0p6(2) )
        my_i2arg = one + arg_gfun*( cheb0q6(1) + arg_gfun*( cheb0q6(2) + arg_gfun*cheb0q6(3) ) )
      ! 39 < T <= 65.3 (fitted from 39 <= T <= 65)
      else if ( arg_gfun.le.65.3D0 ) then
        my_iarg = cheb0p7(0) + arg_gfun*( cheb0p7(1) + arg_gfun*cheb0p7(2) )
        my_i2arg = one + arg_gfun*( cheb0q7(1) + arg_gfun*( cheb0q7(2) + arg_gfun*cheb0q7(3) ) )
      ! 65.3 < T <= 614.9 (fitted from 65 <= T <= 580)
      else if ( arg_gfun.le.614.9D0 ) then
        my_iarg = cheb0p8(0) + arg_gfun*( cheb0p8(1) + arg_gfun*cheb0p8(2) )
        my_i2arg = one + arg_gfun*( cheb0q8(1) + arg_gfun*( cheb0q8(2) + arg_gfun*cheb0q8(3) ) )
      ! T > 614.9 (fitted from 600 <= T <= 10000)
      else
        my_iarg = cheb0p9(0) + arg_gfun*cheb0p9(1)
        my_i2arg = one + arg_gfun*( cheb0q9(1) + arg_gfun*cheb0q9(2) )
      end if
      val_gfun(0) = my_iarg/my_i2arg
      ! 1/T
      my_i2arg = half/arg_gfun
      my_iarg = my_i2arg + my_i2arg
      ! recorder of (2n+1)/2T
      div_down = -my_i2arg
      ! loop over other orders using upward recursion relations
      do iord = 0, order_gfun-1
        ! plus 1/T
        div_down = div_down + my_iarg
        ! perform recursion relation
        val_gfun(iord+1) = my_i2arg - div_down*val_gfun(iord)
      end do ! loop over other orders
    end if
    return
1000  format('GFUN ',A,I6)
1010  format('GFUN ',A,F14.8)
  end subroutine gen1int_aux_gfun

  ! power series expansions of Boys function (for small argument!!)
  subroutine gen1int_boys_power( strt_order, end_order, strt_arg, step_arg, num_arg, val_boys )
    ! the power series expansions can be found, for example, in:
    ! V. R. Saunders. An introduction to molecular integral evaluation.
    ! In G.H.F. Diercksen, B.T. Sutcliffe, and A. Veillard, editors,
    ! Computational Techniques in Quantum Chemistry and Molecular Physics,
    ! Eq. (39), page 347, 1975.
    implicit none
    ! start order of Boys function
    integer, intent(in) :: strt_order
    ! end order of Boys function
    integer, intent(in) :: end_order
    ! start argument of Boys function
    real(xp), intent(in) :: strt_arg
    ! step of the argument
    real(xp), intent(in) :: step_arg
    ! number of arguments
    integer, intent(in) :: num_arg
    ! values of Boys function
    real(xp), intent(out) :: val_boys(strt_order:end_order,0:num_arg)
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! argument of Boys function at each step
    real(xp) my_arg
    ! twice of the argument
    real(xp) my_arg2
    ! incremental recorder for the orders
    integer iord
    ! incremental recorder for the arguments
    integer iarg
    ! value of power series expansions' term
    real(xp) val_power
    ! divisor in the term of power series expansions
    real(xp) div_power
    ! incremental recorder for the summation
    integer ipower
    ! exp(-T)
    real(xp) my_pre_arg
    ! push current subroutine into the stack
    call xsub_enter('gen1int_boys_power')
    ! the minimum argument
    my_arg = min( strt_arg+step_arg*xfloat(num_arg), strt_arg )
    ! the maximum argument
    my_arg2 = max( strt_arg+step_arg*xfloat(num_arg), strt_arg )
    ! check if the value of largest order is reasonable
    if ( my_arg2.gt.twenty3 .and. end_order.ge.nint(my_arg+my_arg) ) then
      write(x_lupri,1000)'Your minimum argument of Boys function: ',my_arg
      write(x_lupri,1010)'Your maximum order of Boys function:    ',end_order
      write(x_lupri,1000)'Warning! Some values of the largest order may not be accurate or correct!!'
    end if
    ! loop over number of arguments
    do iarg = 0, num_arg
      ! for accuracy, this is much better than add "step_arg" to "my_arg" at each time
      my_arg = strt_arg + step_arg*xfloat(iarg)
      ! if the argument exceeds "max_arg_pow_boys", we change to asymptotic series expansions
      if ( my_arg.gt.max_arg_pow_boys ) then
        if (xlib_lprint.ge.200) then
          write(x_lupri,1000)'Your argument of Boys function: ',my_arg
          write(x_lupri,1000)'It is too large, we change to asymptotic series expansions'
        end if
        call gen1int_boys_asymp( strt_order, end_order, my_arg, step_arg, 0, val_boys(:,iarg) )
      else
        ! 2*T
        my_arg2 = my_arg + my_arg
        ! we use power series expansions for the largest order, and use downward recursion
        ! relations for others
        div_power = xfloat(end_order+end_order+1)
        val_power = one/div_power
        ! fisrt term
        val_boys(end_order,iarg) = val_power
        ! loop over power series expansions
        do ipower = 1, max_pterms_boys
          div_power = div_power + two
          val_power = val_power*my_arg2/div_power
          val_boys(end_order,iarg) = val_boys(end_order,iarg) + val_power
          if (val_power.le.g1_cut_aux) exit
        end do
        ! exp(-T)
        my_pre_arg = exp(-my_arg)
        val_boys(end_order,iarg) = my_pre_arg*val_boys(end_order,iarg)
        ! use downward recursion relations for other orders
        div_power = xfloat(end_order+end_order+1)
        do iord = end_order-1, strt_order, -1
          div_power = div_power - two
          val_boys(iord,iarg) = (my_arg2*val_boys(iord+1,iarg)+my_pre_arg)/div_power
        end do ! loop over orders
      end if
    end do ! loop over number of arguments
    ! pop the stack
    call xsub_exit
    return
1000  format('BPOW ',A,F14.8)
1010  format('BPOW ',A,I6)
  end subroutine gen1int_boys_power

  ! asymptotic series expansions of Boys funciton (for large argument!!)
!fixme: we may think about another asymptotic series expansions
  subroutine gen1int_boys_asymp( strt_order, end_order, strt_arg, step_arg, num_arg, val_boys )
    ! we use F_n(T) = (2n-1)!!/2^(n+1)*sqrt(pi/T^(2n+1)), see, for example,
    ! Trygve Helgaker, Poul Jorgensen, Jeppe Olsen, Molecular Electronic Structure Theory,
    ! Eq. (9.8.9), page 365
    implicit none
    ! start order of Boys function
    integer, intent(in) :: strt_order
    ! end order of Boys function
    integer, intent(in) :: end_order
    ! start argument of Boys function
    real(xp), intent(in) :: strt_arg
    ! step of the argument
    real(xp), intent(in) :: step_arg
    ! number of arguments
    integer, intent(in) :: num_arg
    ! values of Boys function
    real(xp), intent(out) :: val_boys(strt_order:end_order,0:num_arg)
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! argument of Boys function at each step
    real(xp) my_arg
    ! 1/T
    real(xp) my_iarg
    ! incremental recorder for the orders
    integer iord
    ! incremental recorder for the arguments
    integer iarg
    ! dividend in the asymptotic term
    real(xp) div_asymp
    ! incremental recorder for the asymptotic series expansions
    integer iasymp
    ! sqrt(pi/T)/2
    real(xp) my_pre_arg
    ! push current subroutine into the stack
    call xsub_enter('gen1int_boys_asymp')
    ! the minimum argument
    my_arg = min( strt_arg+step_arg*xfloat(num_arg), strt_arg )
    ! check if the minimum argument is large enough
    if ( my_arg.le.max_arg_pow_boys ) then
      write(x_lupri,1000)'Your minimum argument of Boys function: ',my_arg
      call xstop('Too small argument for asymptotic series expansions!')
    end if
    ! loop over number of arguments
    do iarg = 0, num_arg
      my_arg = strt_arg + step_arg*xfloat(iarg)
      ! 1/T
      my_iarg = one/my_arg
      ! -1/(2T)
      div_asymp = -half*my_iarg
      ! sqrt(pi/T)/2
      val_boys(strt_order,iarg) = half*sqrt(pi/my_arg)
      ! the smallest order
      do iasymp = 1, strt_order
        ! plus 1/T
        div_asymp = div_asymp + my_iarg
        ! times (2n-1)/(2T)
        val_boys(strt_order,iarg) = val_boys(strt_order,iarg)*div_asymp
      end do
      ! the left orders
      do iord = strt_order, end_order-1
        ! plus 1/T
        div_asymp = div_asymp + my_iarg
        val_boys(iord+1,iarg) = val_boys(iord,iarg)*div_asymp
      end do ! loop over orders
    end do ! loop over number of arguments
    ! pop the stack
    call xsub_exit
    return
1000  format('BASY ',A,F14.8)
  end subroutine gen1int_boys_asymp

  ! power series expansions of G function (for small argument!!)
  ! NB: the downward recursion relation is very bad for large arguments!!
  subroutine gen1int_gfun_power( strt_order, end_order, strt_arg, step_arg, num_arg, val_gfun )
!fixme: give the reference
    implicit none
    ! start order of G function
    integer, intent(in) :: strt_order
    ! end order of G function
    integer, intent(in) :: end_order
    ! start argument of G function
    real(xp), intent(in) :: strt_arg
    ! step of the argument
    real(xp), intent(in) :: step_arg
    ! number of arguments
    integer, intent(in) :: num_arg
    ! values of G function
    real(xp), intent(out) :: val_gfun(strt_order:end_order,0:num_arg)
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! argument of G function at each step
    real(xp) my_arg
    ! twice of the argument
    real(xp) my_arg2
    ! incremental recorder for the orders
    integer iord
    ! incremental recorder for the arguments
    integer iarg
    ! value of power series expansions' term
    real(xp) val_power
    ! divisor in the term of power series expansions
    real(xp) div_power
    ! incremental recorder for the summation
    integer ipower
    ! exp(-T)
    real(xp) my_pre_arg
    ! push current subroutine into the stack
    call xsub_enter('gen1int_gfun_power')
    ! the maximum argument
    my_arg = max( strt_arg+step_arg*xfloat(num_arg), strt_arg )
    ! check if the value of largest order is reasonable
    if ( my_arg.ge.twenty ) then
      write(x_lupri,1000)'Your maximum argument of G function: ',my_arg
      write(x_lupri,1010)'Your maximum order of G function:    ',end_order
      write(x_lupri,1000)'Warning! The downward recursion relation may be unstable!!'
    end if
    ! loop over number of arguments
    do iarg = 0, num_arg
      ! for accuracy, this is much better than add "step_arg" to "my_arg" at each time
      my_arg = strt_arg + step_arg*xfloat(iarg)
      ! if the argument exceeds "max_arg_pow_gfun", we change to asymptotic series expansions
      if ( my_arg.gt.max_arg_pow_gfun ) then
        if (xlib_lprint.ge.200) then
          write(x_lupri,1000)'Your argument of G function: ',my_arg
          write(x_lupri,1000)'It is too large, we change to asymptotic series expansions'
        end if
        call gen1int_gfun_asymp( strt_order, end_order, my_arg, step_arg, 0, val_gfun(:,iarg) )
      else
        ! we use power series expansions for the largest order, and use downward recursion
        ! relations for others
        ! recorder of 2n+2k+1
        div_power = xfloat(end_order+end_order+1)
        ! recorder for T^k/k!
        my_arg2 = one
        ! fisrt term 1/(2n+1)
        val_gfun(end_order,iarg) = my_arg2/div_power
        ! loop over power series expansions
        do ipower = 1, max_pterms_gfun
          div_power = div_power + two
          ! time T/k
          my_arg2 = my_arg2*my_arg/xfloat(ipower)
          ! divided by 2n+2k+1
          val_power = my_arg2/div_power
          val_gfun(end_order,iarg) = val_gfun(end_order,iarg) + val_power
          if (val_power.le.g1_cut_aux) exit
        end do
        ! exp(-T)
        my_pre_arg = exp(-my_arg)
        val_gfun(end_order,iarg) = val_gfun(end_order,iarg)*my_pre_arg
        ! 2*T
        my_arg2 = my_arg + my_arg
        ! use downward recursion relations for other orders
        div_power = xfloat(end_order+end_order+1)
        do iord = end_order-1, strt_order, -1
          div_power = div_power - two
          val_gfun(iord,iarg) = (one-my_arg2*val_gfun(iord+1,iarg))/div_power
        end do ! loop over orders
      end if
    end do ! loop over number of arguments
    ! pop the stack
    call xsub_exit
    return
1000  format('GPOW ',A,F14.8)
1010  format('GPOW ',A,I6)
  end subroutine gen1int_gfun_power

  ! asymptotic series expansions of G function (for large argument!! > 30)
  subroutine gen1int_gfun_asymp( strt_order, end_order, strt_arg, step_arg, num_arg, val_gfun )
    ! we use G_n(T) = 0.5*sum_{j=0}^{J} (0.5-n)_{j} (1/T)^{j+1} = sum_{j=0}^{J} V_{j}
    ! V_{0} = 1/(2T), V_{j} = (0.5-n+j-1)*V_{j-1}/T
!fixme: give the reference
    implicit none
    ! start order of G function
    integer, intent(in) :: strt_order
    ! end order of G function
    integer, intent(in) :: end_order
    ! start argument of G function
    real(xp), intent(in) :: strt_arg
    ! step of the argument
    real(xp), intent(in) :: step_arg
    ! number of arguments
    integer, intent(in) :: num_arg
    ! values of G function
    real(xp), intent(out) :: val_gfun(strt_order:end_order,0:num_arg)
    !////////////////////////////////////////////////////////////////
    !                       local variables
    !////////////////////////////////////////////////////////////////
    ! argument of G function at each step
    real(xp) my_arg
    ! 1/T
    real(xp) my_iarg
    ! 1/(2T)
    real(xp) my_i2arg
    ! dividend in the asymptotic term
    real(xp) div_asymp
    ! incremental recorder for the orders
    integer iord
    ! incremental recorder for the arguments
    integer iarg
    ! value of asymptotic term
    real(xp) val_asymp
    ! maximum terms for asymptotic series expansions, for safety
!fixme: how to determine it, depends on the order and argument
    integer, parameter :: max_aterms_gfun = 100
    ! incremental recorder for the asymptotic series expansions
    integer iasymp
    ! push current subroutine into the stack
    call xsub_enter('gen1int_gfun_asymp')
    ! the minimum argument
    my_arg = min( strt_arg+step_arg*xfloat(num_arg), strt_arg )
    ! check if the minimum argument is large enough
    if ( my_arg.le.max_arg_pow_gfun ) then
      write(x_lupri,1000)'Your minimum argument of G function: ',my_arg
      call xstop('Too small argument for asymptotic series expansions!')
    end if
    ! check if the minimum argument and maximum order are good for a stable upward recursion relation
    if ( end_order.ge.floor(my_arg+my_arg)+10 ) then
      write(x_lupri,1010)'The maximum order of G functions:    ',end_order
      write(x_lupri,1000)'The minimum argument of G functions: ',my_arg
      write(x_lupri,1000)'The upward recursion relation is unstable for them!'
#ifndef TEST_AUX
      call xstop('Write to author asking them to improve it!')
#endif
    end if
    ! loop over number of arguments
    do iarg = 0, num_arg
      my_arg = strt_arg + step_arg*xfloat(iarg)
      ! we use asymptotic series expansions for the smallest order, and upward recursion relations
      ! for others (if there are)
      ! 1/(2T)
      my_i2arg = half/my_arg
      ! 1/T
      my_iarg = my_i2arg + my_i2arg
      ! (0.5-n+j-1)/T == -(2n+1)/(2T), j == 0
      div_asymp = -xfloat(strt_order+strt_order+1)*my_i2arg
      ! first term in asymptotic series expansions: 1/T
      val_asymp = my_i2arg
      val_gfun(strt_order,iarg) = val_asymp
      ! left terms
      do iasymp = 1, max_aterms_gfun
        ! plus 1/T
        div_asymp = div_asymp + my_iarg
        ! times (0.5-n+j-1)/T
        val_asymp = val_asymp*div_asymp
        val_gfun(strt_order,iarg) = val_gfun(strt_order,iarg) + val_asymp
        if ( abs(val_asymp).lt.g1_cut_aux ) exit
      end do
      ! recorder of (2n+1)/2T
      div_asymp = xfloat(strt_order+strt_order-1)*my_i2arg
      ! loop over other orders using upward recursion relations
      do iord = strt_order, end_order-1
        ! plus 1/T
        div_asymp = div_asymp + my_iarg
        ! perform recursion relation
        val_gfun(iord+1,iarg) = my_i2arg - div_asymp*val_gfun(iord,iarg)
      end do ! loop over other orders
    end do ! loop over number of arguments
    ! pop the stack
    call xsub_exit
    return
1000  format('GASY ',A,F14.8)
1010  format('GASY ',A,I6)
  end subroutine gen1int_gfun_asymp

end module gen1int_rfunction

