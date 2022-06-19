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
#  This file generates the referenced results of testing geom_part_zero.F90.
#
#  2011-08-04, Bin Gao:
#  * first version

# prints creating date
from datetime import date

print ">> Generates the referenced results of testing geom_part_zero.F90 ..."

# opens geom_part_zero.h
header_file = open("geom_part_zero.h", "w")
header_file.write("  ! created by GeomPartZero.py on "+date.today().isoformat()+"\n")

# possible orders of partial geometric derivatives on bra and ket centers at first
order_geo_bra = [0,3,0,3]
order_geo_ket = [0,0,3,3]
num_poss_geo = len(order_geo_bra)

# one-center geometric derivatives
idx_one_cent = [2]
order_one_cent = 3
# initializes the indices of bra and ket centers
idx_bra_ket = [7,8]
list_idx_bk = [idx_bra_ket]
# generates the possible indices of bra and ket centers
for iatom in xrange(len(idx_bra_ket)):
    list_idx_bk.append(idx_bra_ket[0:iatom]+idx_one_cent+idx_bra_ket[iatom+1:])
# bra and ket centers are identical
list_idx_bk.append(idx_one_cent+idx_one_cent)
# number of one-center geometric derivatives
num_one_geo = len(list_idx_bk)

# two-center geometric derivatives
idx_two_cent = [1,2]
order_two_cent = [3,4]
# generates the possible indices of bra and ket centers
list_idx_bk.append([idx_two_cent[1],idx_bra_ket[0]])
list_idx_bk.append(idx_two_cent)
list_idx_bk.append([idx_two_cent[1],idx_two_cent[0]])
# number of two-center geometric derivatives
num_tot_geo = len(list_idx_bk)
num_two_geo = num_tot_geo-num_one_geo

# number of all possible terms
num_one_terms = num_one_geo*num_poss_geo
num_two_terms = num_two_geo*num_poss_geo
num_tot_terms = num_tot_geo*num_poss_geo

# number of total geometric derivatives
header_file.write("  ! number of total geometric derivatives, according to\n")
header_file.write("  ! (1) numbers of differentiated centers\n")
header_file.write("  ! (2) indices of differentiated centers\n")
header_file.write("  ! (3) different indices of bra and ket centers\n")
header_file.write("  ! (4) orders of partial geoemtric derivatives on bra and ket centers\n")
header_file.write("  integer, parameter :: NUM_TOT_GEO = "+repr(num_tot_terms)+"\n")

num_term = 8                         #how many terms per line
left_terms = num_tot_terms%num_term  #left terms in the last line
if (left_terms==0):
    left_terms = num_term

# writes the numbers of differentiated centers
header_file.write("  ! numbers of differentiated centers\n")
header_file.write("  integer, parameter :: NUM_CENTS(NUM_TOT_GEO) = (/                      &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        if (jterm<num_one_terms):
            dump_str = dump_str+repr(1).rjust(5)+",   "
        else:
            dump_str = dump_str+repr(2).rjust(5)+",   "
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    if (jterm<num_one_terms):
        dump_str = dump_str+repr(1).rjust(5)+",   "
    else:
        dump_str = dump_str+repr(2).rjust(5)+",   "
if (num_one_terms==num_tot_terms):
  dump_str = dump_str+repr(1).rjust(5)+"/)\n"
else:
  dump_str = dump_str+repr(2).rjust(5)+"/)\n"
header_file.write(dump_str)

# writes indices of differentiated centers
header_file.write("  ! indices of differentiated centers\n")
header_file.write("  integer, parameter :: IDX_GEO_CENT(NUM_TOT_GEO*2) = (/                 &\n")
str_idx_one = repr(idx_one_cent[0]).rjust(5)+","+repr(0).rjust(2)
str_idx_two = repr(idx_two_cent[0]).rjust(5)+","+repr(idx_two_cent[1]).rjust(2)
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        if (jterm<num_one_terms):
            dump_str = dump_str+str_idx_one+","
        else:
            dump_str = dump_str+str_idx_two+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    if (jterm<num_one_terms):
        dump_str = dump_str+str_idx_one+","
    else:
        dump_str = dump_str+str_idx_two+","
if (num_one_terms==num_tot_terms):
  dump_str = dump_str+str_idx_one+"/)\n"
else:
  dump_str = dump_str+str_idx_two+"/)\n"
header_file.write(dump_str)

# writes the orders of differentiated centers
header_file.write("  ! orders of differentiated centers\n")
header_file.write("  integer, parameter :: ORDER_GEO_CENT(NUM_TOT_GEO*2) = (/               &\n")
str_order_one = repr(order_one_cent).rjust(5)+","+repr(0).rjust(2)
str_order_two = repr(order_two_cent[0]).rjust(5)+","+repr(order_two_cent[1]).rjust(2)
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        if (jterm<num_one_terms):
            dump_str = dump_str+str_order_one+","
        else:
            dump_str = dump_str+str_order_two+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    if (jterm<num_one_terms):
        dump_str = dump_str+str_order_one+","
    else:
        dump_str = dump_str+str_order_two+","
if (num_one_terms==num_tot_terms):
  dump_str = dump_str+str_order_one+"/)\n"
else:
  dump_str = dump_str+str_order_two+"/)\n"
header_file.write(dump_str)

# writes the indices of bra and ket centers
header_file.write("  ! indices of bra and ket centers\n")
header_file.write("  integer, parameter :: IDX_BRA_KET(NUM_TOT_GEO*2) = (/                  &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        dump_str = dump_str+repr(list_idx_bk[jterm/num_poss_geo][0]).rjust(5)+"," \
                           +repr(list_idx_bk[jterm/num_poss_geo][1]).rjust(2)+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+repr(list_idx_bk[jterm/num_poss_geo][0]).rjust(5)+"," \
                       +repr(list_idx_bk[jterm/num_poss_geo][1]).rjust(2)+","
dump_str = dump_str+repr(list_idx_bk[num_tot_geo-1][0]).rjust(5)+"," \
                   +repr(list_idx_bk[num_tot_geo-1][1]).rjust(2)+"/)\n"
header_file.write(dump_str)

# writes the orders of partial geoemtric derivatives on bra center
header_file.write("  ! orders of partial geoemtric derivatives on bra center\n")
header_file.write("  integer, parameter :: ORDER_GEO_BRA(NUM_TOT_GEO) = (/                  &\n")
str_order_bra = []
for igeo in xrange(num_poss_geo):
    str_order_bra.append(repr(order_geo_bra[igeo]).rjust(5))
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        dump_str = dump_str+str_order_bra[jterm%num_poss_geo]+",   "
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+str_order_bra[jterm%num_poss_geo]+",   "
dump_str = dump_str+str_order_bra[num_poss_geo-1]+"/)\n"
header_file.write(dump_str)

# writes the orders of partial geoemtric derivatives on ket center
header_file.write("  ! orders of partial geoemtric derivatives on ket center\n")
header_file.write("  integer, parameter :: ORDER_GEO_KET(NUM_TOT_GEO) = (/                  &\n")
str_order_ket = []
for igeo in xrange(num_poss_geo):
    str_order_ket.append(repr(order_geo_ket[igeo]).rjust(5))
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        dump_str = dump_str+str_order_ket[jterm%num_poss_geo]+",   "
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+str_order_ket[jterm%num_poss_geo]+",   "
dump_str = dump_str+str_order_ket[num_poss_geo-1]+"/)\n"
header_file.write(dump_str)

# determines the partial derivative terms
zero_ints = []
order_part_cent = []
scatter_deriv = []
seq_part_geo = []
# one-center geometric derivatives
for iterm in xrange(num_one_geo):
    for jterm in xrange(num_poss_geo):
        num_geo_cent = list_idx_bk[iterm].count(idx_one_cent[0])
        # differentitated center is neither bra nor ket center;
        # or bra and ket centers are identical
        if (num_geo_cent==0 or num_geo_cent==2):
            zero_ints.append(True)
            order_part_cent.append("    0, 0")
            scatter_deriv.append(False)
            seq_part_geo.append(repr(0).rjust(5)+","+repr(0).rjust(2))
        # differentitated center is either bra or ket center
        else:
            zero_ints.append(False)
            place_geo_cent = list_idx_bk[iterm].index(idx_one_cent[0])
            # differentitated center is bra center
            if (place_geo_cent==0):
                order_part_cent.append(repr(order_geo_bra[jterm]+order_one_cent).rjust(5)+"," \
                                       +repr(order_geo_ket[jterm]).rjust(2))
                scatter_deriv.append(order_geo_bra[jterm]>0 or order_geo_ket[jterm]>0)
            # differentitated center is ket center
            else:
                order_part_cent.append(repr(order_geo_bra[jterm]).rjust(5)+"," \
                                       +repr(order_geo_ket[jterm]+order_one_cent).rjust(2))
                scatter_deriv.append(order_geo_ket[jterm]>0)
            seq_part_geo.append(repr(place_geo_cent+1).rjust(5)+","+repr(0).rjust(2))
# two-center geometric derivatives
for iterm in xrange(num_one_geo,num_tot_geo):
    for jterm in xrange(num_poss_geo):
        # differentitated centers are bra and ket centers
        if (list_idx_bk[iterm].count(idx_two_cent[0]) and \
            list_idx_bk[iterm].count(idx_two_cent[1])):
            zero_ints.append(False)
            place_first = list_idx_bk[iterm].index(idx_two_cent[0])
            place_second = list_idx_bk[iterm].index(idx_two_cent[1])
            # first differentitated center is bra center
            if (place_first==0):
                order_part_cent.append(repr(order_geo_bra[jterm]+order_two_cent[0]).rjust(5)+"," \
                                       +repr(order_geo_ket[jterm]+order_two_cent[1]).rjust(2))
                scatter_deriv.append(order_geo_bra[jterm]>0 or order_geo_ket[jterm]>0)
            # first differentitated center is ket center
            else:
                order_part_cent.append(repr(order_geo_bra[jterm]+order_two_cent[1]).rjust(5)+"," \
                                       +repr(order_geo_ket[jterm]+order_two_cent[0]).rjust(2))
                scatter_deriv.append(True)
            seq_part_geo.append(repr(place_first+1).rjust(5)+","+repr(place_second+1).rjust(2))
        # at least, one differentitated center is neither bra nor ket center
        else:
            zero_ints.append(True)
            order_part_cent.append("    0, 0")
            scatter_deriv.append(False)
            seq_part_geo.append(repr(0).rjust(5)+","+repr(0).rjust(2))

# writes the left information
header_file.write("  ! referenced results indicating if the total geometric derivatives are zero\n")
header_file.write("  logical, parameter :: REF_ZERO_INTS(NUM_TOT_GEO) = (/                     &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = "    "
    for jterm in xrange(iterm,iterm+num_term):
        if (zero_ints[jterm]):
            dump_str = dump_str+".true.,  "
        else:
            dump_str = dump_str+".false., "
    header_file.write(dump_str+"&\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = "    "
for jterm in xrange(iterm,num_tot_terms-1):
    if (zero_ints[jterm]):
        dump_str = dump_str+".true.,  "
    else:
        dump_str = dump_str+".false., "
if (zero_ints[num_tot_terms-1]):
    dump_str = dump_str+".true./)\n"
else:
    dump_str = dump_str+".false./)\n"
header_file.write(dump_str)
#
header_file.write("  ! referenced final orders of partial geometric derivatives by adding\n")
header_file.write("  ! those from total geometric derivatives\n")
header_file.write("  integer, parameter :: REF_ORDER_PART(NUM_TOT_GEO*2) = (/               &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        dump_str = dump_str+order_part_cent[jterm]+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+order_part_cent[jterm]+","
dump_str = dump_str+order_part_cent[num_tot_terms-1]+"/)\n"
header_file.write(dump_str)
#
header_file.write("  ! referenced results indicating if scattering the geometric derivatives later on\n")
header_file.write("  logical, parameter :: REF_SCATTER_DER(NUM_TOT_GEO) = (/                   &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = "    "
    for jterm in xrange(iterm,iterm+num_term):
        if (scatter_deriv[jterm]):
            dump_str = dump_str+".true.,  "
        else:
            dump_str = dump_str+".false., "
    header_file.write(dump_str+"&\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = "    "
for jterm in xrange(iterm,num_tot_terms-1):
    if (scatter_deriv[jterm]):
        dump_str = dump_str+".true.,  "
    else:
        dump_str = dump_str+".false., "
if (scatter_deriv[num_tot_terms-1]):
    dump_str = dump_str+".true./)\n"
else:
    dump_str = dump_str+".false./)\n"
header_file.write(dump_str)
#
header_file.write("  ! referenced sequences of bra and ket centers for partial derivative terms\n")
header_file.write("  integer, parameter :: REF_SEQ_PART(NUM_TOT_GEO*2) = (/                 &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        dump_str = dump_str+seq_part_geo[jterm]+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+seq_part_geo[jterm]+","
dump_str = dump_str+seq_part_geo[num_tot_terms-1]+"/)\n"
header_file.write(dump_str)

# closes geom_part_zero.h
header_file.close()
print ">> Please replace test_f90/geom_part_zero.h with this generated one ..."
