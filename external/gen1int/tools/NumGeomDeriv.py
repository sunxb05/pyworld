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
#  This file computes the number of geometric derivative integrals of given
#  number of (1) atoms, (2) operators, and (3) total order of geometric derivatives.
#
#  2009-09-01, Bin Gao:
#  * use reverse lexical ordering for all the geometric derivatives, i.e., increments
#    first from x, then y, and z, the orders with smaller component in z comes first,
#    and then y, x
#
#  2009-08-27, Bin Gao:
#  * add three- and four-center geometric derivatives
#
#  2009-08-26, Bin Gao:
#  * first version

# read single characters, unbuffered, from standard input
# from Python cookbook, Av Alex Martelli, Anna Martelli Ravenscroft, David Ascher
try:
    from msvcrt import getch
except ImportError:
    """We're not on Windows, so we try the Unix-like approach"""
    def getch():
        import sys, tty, termios
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch

# one-center geometric derivatives
def one_cent_gderiv(tot_order):
    # index
    index_derv = 0
    # derivative along z direction
    for k in xrange(tot_order+1):
        # derivative along y direction
        for j in xrange(tot_order-k+1):
            # index of derivatives
            # index_derv = 1+j+(2*tot_order+3-k)*k/2
            # derivative along x direction
            i = tot_order-(j+k)
            # dump
            index_derv += 1
            print '    ',index_derv,'  (',i,',',j,',',k,')'

# two-center geometric derivatives
def two_cent_gderiv(tot_order):
    # index
    index_derv = 0
    # total order of derivatives acting on the 1st center
    for ijk in xrange(1,tot_order):
        # total order of derivatives acting on the 2nd center
        lmn = tot_order-ijk
        # derivative along z direction (1st center)
        for k in xrange(ijk+1):
            # derivative along y direction (1st center)
            for j in xrange(ijk-k+1):
                # derivative along x direction (1st center)
                i = ijk-(j+k)
                # derivative along z direction (2nd center)
                for n in range(lmn+1):
                    # derivative along y direction (2nd center)
                    for m in range(lmn-n+1):
                        # derivative along x direction (2nd center)
                        l = lmn-(m+n)
                        # dump
                        index_derv += 1
                        print '    ',index_derv,       \
                              '  (',i,',',j,',',k,')', \
                              ' (',l,',',m,',',n,')'

# three-center geometric derivatives
def three_cent_gderiv(tot_order):
    # index
    index_derv = 0
    # total order of derivatives acting on the 1st center
    for ijk in xrange(1,tot_order-1):
        # total order of derivatives acting on the 2nd center
        for lmn in xrange(1,tot_order-ijk):
            # total order of derivatives acting on the 3rd center
            rst = tot_order-(ijk+lmn)
            # derivative along z direction (1st center)
            for k in xrange(ijk+1):
                # derivative along y direction (1st center)
                for j in xrange(ijk-k+1):
                    # derivative along x direction (1st center)
                    i = ijk-(j+k)
                    # derivative along z direction (2nd center)
                    for n in range(lmn+1):
                        # derivative along y direction (2nd center)
                        for m in range(lmn-n+1):
                            # derivative along x direction (2nd center)
                            l = lmn-(m+n)
                            # derivative along z direction (3rd center)
                            for t in xrange(rst+1):
                                # derivative along y direction (3rd center)
                                for s in xrange(rst-t+1):
                                    # derivative along x direction (3rd center)
                                    r = rst-(s+t)
                                    # dump
                                    index_derv += 1
                                    print '    ',index_derv,       \
                                          '  (',i,',',j,',',k,')', \
                                          ' (',l,',',m,',',n,')',  \
                                          ' (',r,',',s,',',t,')'

# four-center geometric derivatives
def four_cent_gderiv(tot_order):
    # index
    index_derv = 0
    # total order of derivatives acting on the 1st center
    for ijk in xrange(1,tot_order-2):
        # total order of derivatives acting on the 2nd center
        for lmn in xrange(1,tot_order-(ijk+1)):
            # total order of derivatives acting on the 3rd center
            for rst in xrange(1,tot_order-(ijk+lmn)):
                # total order of derivatives acting on the 4th center
                uvw = tot_order-(ijk+lmn+rst)
                # derivative along z direction (1st center)
                for k in xrange(ijk+1):
                    # derivative along y direction (1st center)
                    for j in xrange(ijk-k+1):
                        # derivative along x direction (1st center)
                        i = ijk-(j+k)
                        # derivative along z direction (2nd center)
                        for n in range(lmn+1):
                            # derivative along y direction (2nd center)
                            for m in range(lmn-n+1):
                                # derivative along x direction (2nd center)
                                l = lmn-(m+n)
                                # derivative along z direction (3rd center)
                                for t in xrange(rst+1):
                                    # derivative along y direction (3rd center)
                                    for s in xrange(rst-t+1):
                                        # derivative along x direction (3rd center)
                                        r = rst-(s+t)
                                        # derivative along z direction (4th center)
                                        for w in xrange(uvw+1):
                                            # derivative along y direction (4th center)
                                            for v in xrange(uvw-w+1):
                                                # derivative along x direction (4th center)
                                                u = uvw-(v+w)
                                                # dump
                                                index_derv += 1
                                                print '    ',index_derv,       \
                                                      '  (',i,',',j,',',k,')', \
                                                      ' (',l,',',m,',',n,')',  \
                                                      ' (',r,',',s,',',t,')',  \
                                                      ' (',u,',',v,',',w,')'

# binomial function
def binomial(n,k):
    bcoef = 1
    n_inc = n+1
    if k<n-k:
        for i in xrange(n,n_inc-k-1,-1):
            bcoef = bcoef*i/(n_inc-i)
    else:
        for i in xrange(n,k,-1):
            bcoef = bcoef*i/(n_inc-i)
    return bcoef

############################################################
#                  main program
############################################################
#
# reads the number of (1) atoms, (2) operators, and (3) total order of derivatives
num_atom = input(' >> Number of atoms: ')
num_op1  = input(' >> Number of operators: ')
tot_derv = input(' >> Total order of geometric derivatives: ')
print ''

# power function
from math import pow

# initialization
num_gderv1 = 0
num_gderv2 = 0
num_gderv3 = 0
num_gderv4 = 0

# considers the different total order
if tot_derv>=1:
    # number of one-center geometric derivatives
    num_gderv1 = (tot_derv+1)*(tot_derv+2)/2
    print ' >> Number of unique one-center geometric derivatives: ',num_gderv1
    # multiplied by the number of atoms and operators
    num_gderv1 *= num_atom*num_op1
    print ' >> Number of unique one-center geometric derivative integrals: ',num_gderv1
    print ''
    # dumps the geometric derivatives
    print ' >> Do you want to show the unique one-center geometric derivatives? (y/n) '
    if getch() == 'y':
        one_cent_gderiv(tot_derv)
    print ''
    # 2nd order geometric derivatives
    if tot_derv>=2:
        # number of two-center geometric derivatives
        # got from Maxima:
        #   sum((i+1)*(i+2)*(n-i+1)*(n-i+2)/4,i,1,n-1),simpsum;ratsimp(%);factor(%);
        num_gderv2 = (tot_derv-1)*(tot_derv+1)*(tot_derv+2)*((tot_derv+13)*tot_derv+60)/120
        print ' >> Number of unique two-center geometric derivatives: ',num_gderv2
        # multiplied by the number of atoms and operators
        # number of atoms from Maxima: f(n):=sum(n-i,i,1,n);
        num_gderv2 *= (num_atom-1)*num_atom*num_op1/2
        print ' >> Number of unique two-center geometric derivative integrals: ',num_gderv2
        print ''
        # dumps the two-center geometric derivatives
        print ' >> Do you want to show the unique two-center geometric derivatives? (y/n) '
        if getch() == 'y':
            two_cent_gderiv(tot_derv)
        print ''
        # 3rd order geometric derivatives
        if tot_derv>=3:
            # number of three-center geometric derivatives
            # got from Maxima:
            #   sum(sum((i+1)*(i+2)*(j+1)*(j+2)*(n-i-j+1)*(n-i-j+2)/8,j,1,n-i-1),i,1,n-2),simpsum;
            #   ratsimp(%);factor(%);
            num_gderv3 = (tot_derv-2)*(tot_derv-1)*(tot_derv+1)*(tot_derv+2) \
                         *((((tot_derv+36)*tot_derv+551)*tot_derv+3708)*tot_derv+10080)/40320
            print ' >> Number of unique three-center geometric derivatives: ',num_gderv3
            # multiplied by the number of atoms and operators
            # number of atoms from Maxima: f(n):=sum(sum(n-j,j,i+1,n),i,1,n);
            num_gderv3 *= ((num_atom-3)*num_atom+2)*num_atom*num_op1/6
            print ' >> Number of unique three-center geometric derivative integrals: ',num_gderv3
            print ''
            # dumps the three-center geometric derivatives
            print ' >> Do you want to show the unique three-center geometric derivatives? (y/n) '
            if getch() == 'y':
                three_cent_gderiv(tot_derv)
            print ''
            # 4th order or higher
            if tot_derv >= 4:
                # number of four-center geometric derivatives
                # got from Maxima:
                #   sum(sum(sum((i+1)*(i+2)*(j+1)*(j+2)*(k+1)*(k+2)*(n-i-j-k+1)*(n-i-j-k+2)/16,
                #     k,1,n-i-j-1),j,1,n-i-2),i,1,n-3),simpsum;ratsimp(%);factor(%);
                num_gderv4 = (tot_derv-3)*(tot_derv-2)*(tot_derv-1)*(tot_derv+1)*(tot_derv+2) \
                             *((((((tot_derv+69)*tot_derv+2137)*tot_derv+35451)*tot_derv \
                             +330862)*tot_derv+1612920)*tot_derv+3326400)/39916800
                print ' >> Number of unique four-center geometric derivatives: ',num_gderv4
                # multiplied by the number of atoms and operators
                # number of atoms from Maxima: f(n):=sum(sum(sum(k-1,k,1,j-1),j,1,i-1),i,1,n);
                num_gderv4 *= (((num_atom-6)*num_atom+11)*num_atom-6)*num_atom*num_op1/24
                print ' >> Number of unique four-center geometric derivative integrals: ',num_gderv4
                print ''
                # dumps the four-center geometric derivatives
                print ' >> Do you want to show the unique four-center geometric derivatives? (y/n) '
                if getch() == 'y':
                    four_cent_gderiv(tot_derv)
                print ''
# invalid number
else:
    print ''
    print ' !! Invalid total order of geometric derivatives: ',tot_derv

# total number of derivative integrals
num_unique = num_gderv1+num_gderv2+num_gderv3+num_gderv4
print ' >> Total number of unique derivative integrals: ',num_unique
print ' >> Used memory in MB (16 bytes for complex numbers): ',num_unique*0.0000152587890625
num_redundant = pow((3*num_atom),tot_derv)
print ' >> Total number of redundant derivative integrals (3N)^{g}: ',int(num_redundant)
print ' >> Used memory in MB (16 bytes for complex numbers): ',num_redundant*0.0000152587890625
print ' >> Sparsity (%): ',(num_redundant-num_unique)/num_redundant*100
num_non_redund = binomial(3*num_atom+tot_derv-1,tot_derv)
print ' >> Total number of non-redundant derivative integrals binomial(3N+g-1,g): ',num_non_redund
print ' >> Used memory in MB (16 bytes for complex numbers): ',num_non_redund*0.0000152587890625
print ' >> Sparsity (%): ',(num_non_redund-num_unique)/float(num_non_redund)*100
