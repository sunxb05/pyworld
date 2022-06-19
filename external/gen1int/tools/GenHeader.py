#!/usr/bin/env python
#
#  gen1int: compute one-electron integrals using rotational London atomic-orbitals
#  Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen
#
#  gen1int is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  gen1int is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
#
#  This file generates the header files used in Gen1Int.
#
#  2011-03-29, Bin Gao:
#  * first version

# generates the header file containing the number of xyz-components for a given order
def XYZTriangle(max_gen_order, fl_header):
    # leading spaces
    leading_spaces = 4*" "
    # how many terms per line
    num_term = 6
    # left orders in the last line (+1 for considering zeroth order)
    left_order = (max_gen_order+1)%num_term
    if (left_order==0):
        left_order = num_term
    # writes the accumulated number of xyz-components for a given order
    fl_header.write("  ! accumulated number of xyz-components for a given order\n")
    fl_header.write("  integer, parameter :: TRI_ACC_SIZE(0:MAX_GEN_ORDER) = (/ &\n")
    # processes the lines before the last one
    for iorder in xrange(0,max_gen_order+1-left_order,num_term):
        dump_str = leading_spaces
        for jorder in xrange(iorder,iorder+num_term):
            tri_acc_size = (jorder+1)*(jorder+2)*(jorder+3)/6
            dump_str = dump_str+repr(tri_acc_size).rjust(8)+","
        fl_header.write(dump_str+" &\n")
    # left orders in the last line
    iorder = max_gen_order+1-left_order
    dump_str = leading_spaces
    for jorder in xrange(iorder,max_gen_order):
        tri_acc_size = (jorder+1)*(jorder+2)*(jorder+3)/6
        dump_str = dump_str+repr(tri_acc_size).rjust(8)+","
    tri_acc_size = (max_gen_order+1)*(max_gen_order+2)*(max_gen_order+3)/6
    dump_str = dump_str+repr(tri_acc_size).rjust(8)+"/)\n"
    fl_header.write(dump_str)
    return

# generates the header file of pretabulated Boys function
def TabBoys(max_gen_order, min_arg_tab, interv_tab, nsteps_tab, fl_header):
    # leading spaces
    leading_spaces = 4*" "
    fl_header.write(leading_spaces+"! pre-tabulated Boys function\n")
    fl_header.write(leading_spaces+"real(REALK), parameter :: PRE_TAB_BOYS(" \
                    +repr((max_gen_order+1)*(nsteps_tab+1))+") = (/ &\n")
    # computes the tabulated Boys function using incomplete Gamma function
    try:
        val_boys = BoysGamma(0, max_gen_order, min_arg_tab, interv_tab, nsteps_tab)
    # computes the tabulated Boys function using power series expansion and/or asymptotic series expansion
    except:
        val_boys = BoysPower(0, max_gen_order, min_arg_tab, interv_tab, nsteps_tab)
    # writes the Boys function to header file
    #
    # how many terms per line
    num_term = 2
    # left orders in the last line
    left_order = (max_gen_order+1)%num_term
    if (left_order==0):
        left_order = num_term
#FIXME: users may give the number of digits to use?
    # format of dumping
    fmt_per_line = leading_spaces+num_term*"%29.25F_REALK,"+2*" "+"&\n"
    fmt_last_line = leading_spaces+left_order*"%29.25F_REALK," \
                    +36*(num_term-left_order)*" "+2*" "+"&\n"
    for iarg in xrange(nsteps_tab):
        for iord in xrange(0,max_gen_order-left_order+1,num_term):
            fl_header.write(fmt_per_line %tuple(val_boys[iarg][iord:iord+num_term]))
        # left orders
        fl_header.write(fmt_last_line %tuple(val_boys[iarg][max_gen_order-left_order+1:max_gen_order+1]))
    # the last argument
    fmt_last_line = leading_spaces+(left_order-1)*"%29.25F_REALK,"+"%29.25F_REALK/)\n"
    for iord in xrange(0,max_gen_order-left_order+1,num_term):
        fl_header.write(fmt_per_line %tuple(val_boys[nsteps_tab][iord:iord+num_term]))
    # left orders
    fl_header.write(fmt_last_line %tuple(val_boys[nsteps_tab][max_gen_order-left_order+1:max_gen_order+1]))
    return

# computes the Boys function using incomplete Gamma function
def BoysGamma(strt_order, end_order, strt_arg, step_arg, num_arg):
    # Gamma function and incomplete Gamma integral
    from scipy.special import gamma,gammainc
    # power function
    from math import pow
    # initializes the Boys function
    val_boys = [[0 for iord in xrange(end_order-strt_order+1)] for iarg in xrange(num_arg+1)]
    # zero argument
    if (strt_arg==0):
        val_power = BoysPower(strt_order, end_order, strt_arg, step_arg, 0)
        val_boys[0] = val_power[0]
    else:
        for iord in xrange(strt_order,end_order+1):
            order_func = iord+0.5
            val_boys[0][iord-strt_order] = 0.5*gamma(order_func) \
                                         * gammainc(order_func,strt_arg)/pow(strt_arg,order_func)
    # loops over sampling points and orders
    for iarg in xrange(1,num_arg+1):
        arg_boys = strt_arg+step_arg*iarg  #argument of Boys function
        for iord in xrange(strt_order,end_order+1):
            order_func = iord+0.5
            val_boys[iarg][iord-strt_order] = 0.5*gamma(order_func) \
                                            * gammainc(order_func,arg_boys)/pow(arg_boys,order_func)
    return val_boys

## \brief power series expansion of Boys function (for small argument)
## \detail the power series expansion can be found, for example, in:
##         V. R. Saunders. An introduction to molecular integral evaluation.
##         In G.H.F. Diercksen, B.T. Sutcliffe, and A. Veillard, editors,
##         Computational Techniques in Quantum Chemistry and Molecular Physics,
##         Eq. (39), page 347, 1975.
## \author Bin Gao
## \date 2009-06-27
## \param strt_order is the start order of Boys function
## \param end_order is the end order of Boys function
## \param strt_arg is the start argument of Boys function
## \param step_arg is the step of the argument
## \param num_arg is the number of arguments
## \return val_boys contains the values of Boys function
def BoysPower(strt_order, end_order, strt_arg, step_arg, num_arg):
    # power and exponential functions
    from math import pow,exp
    # error stop
    from sys import exit
    # maximum argument for power series expansion
    max_arg_pow = 30.0
    # maximum number of terms for the power series expansion
#FIXME how to determine?? which depends on the order & argument
#FIXME 200 terms are enough for arguments smaller than, like 20
    max_nterm_pow = 200
    # cut off for the terms in power series expansion
    cut_off_pow = pow(10.0,-18)
    # the minimum argument
    curr_arg = min(strt_arg+step_arg*num_arg, strt_arg)
    # the maximum argument
    twice_arg = max(strt_arg+step_arg*num_arg, strt_arg)
    # checks if the value of largest order is reasonable
    if (twice_arg>23 and end_order>=2*curr_arg):
        print "BoysPower>> minimum argument of Boys function:", curr_arg
        print "BoysPower>> maximum order of Boys function:", end_order
        print "BoysPower>> warning! values of the largest order may not be accurate or correct!"
    # initializes the Boys function
    val_boys = [[0 for iord in xrange(end_order-strt_order+1)] for iarg in xrange(num_arg+1)]
    # loops over arguments
    for iarg in xrange(num_arg+1):
        # for accuracy, this is much better than adding \var(step_arg) to \var(curr_arg) at each time
        curr_arg = strt_arg+step_arg*iarg
        # if the argument exceeds \var(max_arg_pow), we change to asymptotic series expansion
        if (curr_arg>max_arg_pow):
            print "BoysPower>> argument of Boys function:", curr_arg
            print "BoysPower>> it is too large, we change to asymptotic series expansion ..."
            val_asymp = BoysAsymp(strt_order, end_order, curr_arg, step_arg, 0)
            val_boys[iarg] = val_asymp[0]
        else:
            twice_arg = curr_arg+curr_arg
            # we use power series expansion for the largest order,
            # and downward recurrence relations for others
            div_power = end_order+end_order+1
            val_power = 1.0/div_power
            # the first term
            val_boys[iarg][end_order-strt_order] = val_power
            # loops over power series expansion
            not_converged = True
            for ipower in xrange(max_nterm_pow):
                div_power = div_power+2.0
                val_power = val_power*twice_arg/div_power
                val_boys[iarg][end_order-strt_order] \
                    = val_boys[iarg][end_order-strt_order]+val_power
                if (val_power<=cut_off_pow):
                    not_converged = False
                    break
            if (not_converged):
                print "BoysPower>> power series expansion is not converged after step:", \
                      max_nterm_pow, "with argument:", curr_arg, "and order:", end_order
                print "BoysPower>> try to increase max_nterm_pow in GenHeader.py!"
                exit(1)
            curr_exp_arg = exp(-curr_arg)
            val_boys[iarg][end_order-strt_order] \
                = curr_exp_arg*val_boys[iarg][end_order-strt_order]
            # using downward recurrence relations for other orders
            div_power = end_order+end_order+1
            for iord in xrange(end_order-1,strt_order-1,-1):
                div_power = div_power-2.0
                val_boys[iarg][iord-strt_order] \
                    = (twice_arg*val_boys[iarg][iord+1-strt_order]+curr_exp_arg)/div_power
    return val_boys

## \brief asymptotic series expansion of Boys funciton (for large argument!!)
## \detail We use \f$F_n(T) = \frac{(2n-1)!!}{2^{n+1}}\sqrt{\frac{\pi}{T^{2n+1}}}\f$,
##         see, for example, Trygve Helgaker, Poul Jorgensen, Jeppe Olsen,
##         Molecular Electronic Structure Theory, Eq. (9.8.9), page 365.
## \author Bin Gao
## \date 2009-06-27
## \param strt_order is the start order of Boys function
## \param end_order is the end order of Boys function
## \param strt_arg is the start argument of Boys function
## \param step_arg is the step of the argument
## \param num_arg is the number of arguments
## \return val_boys contains the values of Boys function
def BoysAsymp(strt_order, end_order, strt_arg, step_arg, num_arg):
    # square root function function and pi constant
    from math import sqrt,pi
    # error stop
    from sys import exit
    # minimum argument for asymptotic series expansion
    min_arg_asym = 30.0
    # the minimum argument
    curr_arg = min(strt_arg+step_arg*num_arg, strt_arg)
    # checks if the minimum argument is large enough
    if (curr_arg<=min_arg_asym):
        print "BoysAsymp>> minimum argument of Boys function:", curr_arg
        print "BoysAsymp>> minimum argument for asymptotic series expansion:", min_arg_asym
        print "BoysAsymp>> too small argument for asymptotic series expansion!"
        exit(1)
    # initializes the Boys function
    val_boys = [[0 for iord in xrange(end_order-strt_order+1)] for iarg in xrange(num_arg+1)]
    # loops over number of arguments
    for iarg in xrange(num_arg+1):
        curr_arg = strt_arg+step_arg*iarg
        curr_iarg = 1.0/curr_arg
        div_asymp = -curr_iarg*0.5
        # \f$\sqrt{\frac{\pi}{T}}/2\f$
        val_boys[iarg][0] = sqrt(pi/curr_arg)*0.5
        # the smallest order
        for iasymp in xrange(strt_order):
            # adds \f$\frac{1}{T}\f$
            div_asymp = div_asymp+curr_iarg
            # multiplied by \f$\frac{2n-1}{2T}\f$
            val_boys[iarg][0] = val_boys[iarg][0]*div_asymp
        # the left orders
        for iord in xrange(strt_order,end_order):
            # adds \f$\frac{1}{T}\f$
            div_asymp = div_asymp+curr_iarg
            val_boys[iarg][iord+1-strt_order] = val_boys[iarg][iord-strt_order]*div_asymp
    return val_boys

if __name__ == '__main__':
    #-# prints creating date
    #-from datetime import date
    # reads the maximum order for Gen1Int
    max_gen_order = input(">> Please input the maximum order for Gen1Int: ")
    # header file of maximum order for Gen1Int
    print ">> Generating the header file of maximum order for Gen1Int ..."
    fl_header = open("max_gen_order.h", "w")
    fl_header.write("#define MAX_GEN_ORDER "+repr(max_gen_order)+"\n")
    fl_header.close()
    #-# header file containing the number of xyz-components for a given order
    #-print ">> Generating the header file of accumulated number of xyz-components ..."
    #-fl_header = open("xyz_triangle.h", "w")
    #-fl_header.write("  ! created by GenHeader.XYZTriangle on "+date.today().isoformat()+"\n")
    #-fl_header.write('#include "max_gen_order.h"\n')
    #-XYZTriangle(max_gen_order, fl_header)
    #-fl_header.close()
    # header file used by power series expansion of Boys function
    print ">> Generating the header file used by power series expansion of Boys function ..."
    fl_header = open("boys_power.h", "w")
    min_arg_tab = 0.0                    #minimum argument for the tabulated Boys function
    max_arg_tab = 12.0                   #maximum argument for the tabulated Boys function
    nsteps_tab = 120                     #number of sampling points for the tabulated Boys function
    interv_tab = max_arg_tab/nsteps_tab  #interval of sampling points for the tabulated Boys function
    fl_header.write("#define MIN_ARG_TAB "+repr(min_arg_tab)+"_REALK\n")
    fl_header.write("#define MAX_ARG_TAB "+repr(max_arg_tab)+"_REALK\n")
    fl_header.write("#define NSTEPS_TAB "+repr(nsteps_tab)+"\n")
    fl_header.write("#define INTERV_TAB %3.1F_REALK\n" %(interv_tab))
    fl_header.close()
    # header file of pretabulated Boys function
    print ">> Generating the header file of pretabulated Boys function ..."
    fl_header = open("tab_boys.h", "w")
    TabBoys(max_gen_order, min_arg_tab, interv_tab, nsteps_tab, fl_header)
    fl_header.close()
    print ">> Well done ;-)"
    print ">> Please update the following header files:"
    print "   (*) src/max_gen_order.h"
    print "   (*) src/boys_power.h"
    #-print "   (*) src/xyz_triangle.h"
    print "   (*) src/tab_boys.h"
