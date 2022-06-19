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
#  This file generates the referenced results of testing geom_part_one.F90.
#
#  2011-08-04, Bin Gao:
#  * first version

# prints creating date
from datetime import date

print ">> Generates the referenced results of testing geom_part_one.F90 ..."

# opens geom_part_one.h
header_file = open("geom_part_one.h", "w")
header_file.write("  ! created by GeomPartOne.py on "+date.today().isoformat()+"\n")

# possible orders of partial geometric derivatives on bra, ket and operator centers at first
order_geo_bra = [0,2,0,0,2,2,0,2]
order_geo_ket = [0,0,2,0,2,0,2,2]
order_geo_opt = [0,0,0,2,0,2,2,2]
num_poss_geo = len(order_geo_bra)

# one-center geometric derivatives
idx_one_cent = [1]
order_one_cent = 2
# initializes the indices of bra, ket and atomic operator centers
idx_atom_opt = [6,7,8]
num_int_cent = len(idx_atom_opt)
list_idx_int = [idx_atom_opt]
# generates the possible indices of bra, ket and operator centers
for iatom in xrange(num_int_cent):
    list_idx_int.append(idx_atom_opt[0:iatom]+idx_one_cent+idx_atom_opt[iatom+1:])
#
for iatom in xrange(num_int_cent):
    for jatom in xrange(iatom+1,num_int_cent):
        list_idx_int.append(idx_atom_opt[0:iatom]+idx_one_cent+idx_atom_opt[iatom+1:jatom] \
                            +idx_one_cent+idx_atom_opt[jatom+1:])
# bra, ket and operator centers are identical
list_idx_int.append(3*idx_one_cent)
# initializes the indices of bra, ket and non-atomic operator centers
idx_non_opt = idx_atom_opt[0:num_int_cent-1]+[-1]
list_idx_int.append(idx_non_opt)
# generates the possible indices of bra, ket and operator centers
for iatom in xrange(num_int_cent-1):
    list_idx_int.append(idx_non_opt[0:iatom]+idx_one_cent+idx_non_opt[iatom+1:])
#
for iatom in xrange(num_int_cent-1):
    for jatom in xrange(iatom+1,num_int_cent-1):
        list_idx_int.append(idx_non_opt[0:iatom]+idx_one_cent+idx_non_opt[iatom+1:jatom] \
                            +idx_one_cent+idx_non_opt[jatom+1:])
# number of one-center geometric derivatives
num_one_geo = len(list_idx_int)

# two-center geometric derivatives
idx_two_cent = [[1],[2]]
order_two_cent = [2,3]
# generates the possible indices of bra, ket and operator centers
list_idx_int.append(idx_atom_opt)
for iatom in xrange(num_int_cent):
    list_idx_int.append(idx_atom_opt[0:iatom]+idx_two_cent[0]+idx_atom_opt[iatom+1:])
#
for iatom in xrange(num_int_cent):
    for jatom in xrange(iatom):
        list_idx_int.append(idx_atom_opt[0:jatom]+idx_two_cent[1]+idx_atom_opt[jatom+1:iatom] \
                            +idx_two_cent[0]+idx_atom_opt[iatom+1:])
    for jatom in xrange(iatom+1,num_int_cent):
        list_idx_int.append(idx_atom_opt[0:iatom]+idx_two_cent[0]+idx_atom_opt[iatom+1:jatom] \
                            +idx_two_cent[1]+idx_atom_opt[jatom+1:])
#
list_idx_int.append(idx_two_cent[0]+idx_two_cent[0]+idx_two_cent[1])
list_idx_int.append(idx_two_cent[0]+idx_two_cent[1]+idx_two_cent[0])
list_idx_int.append(idx_two_cent[0]+idx_two_cent[1]+idx_two_cent[1])
list_idx_int.append(idx_two_cent[1]+idx_two_cent[0]+idx_two_cent[0])
list_idx_int.append(idx_two_cent[1]+idx_two_cent[0]+idx_two_cent[1])
list_idx_int.append(idx_two_cent[1]+idx_two_cent[1]+idx_two_cent[0])
# non-atomic operator
list_idx_int.append(idx_non_opt)
# generates the possible indices of bra, ket and operator centers
for iatom in xrange(num_int_cent-1):
    list_idx_int.append(idx_non_opt[0:iatom]+idx_two_cent[0]+idx_non_opt[iatom+1:])
#
for iatom in xrange(num_int_cent-1):
    for jatom in xrange(iatom):
        list_idx_int.append(idx_non_opt[0:jatom]+idx_two_cent[1]+idx_non_opt[jatom+1:iatom] \
                            +idx_two_cent[0]+idx_non_opt[iatom+1:])
    for jatom in xrange(iatom+1,num_int_cent-1):
        list_idx_int.append(idx_non_opt[0:iatom]+idx_two_cent[0]+idx_non_opt[iatom+1:jatom] \
                            +idx_two_cent[1]+idx_non_opt[jatom+1:])
# number of two-center geometric derivatives
num_tot_geo = len(list_idx_int)
num_two_geo = num_tot_geo-num_one_geo

# three-center geometric derivatives
idx_tri_cent = [[1],[2],[3]]
order_tri_cent = [1,2,2]
# generates the possible indices of bra, ket and operator centers
list_idx_int.append(idx_atom_opt)
for iatom in xrange(num_int_cent):
    list_idx_int.append(idx_atom_opt[0:iatom]+idx_tri_cent[0]+idx_atom_opt[iatom+1:])
#
for iatom in xrange(num_int_cent):
    for jatom in xrange(iatom):
        list_idx_int.append(idx_atom_opt[0:jatom]+idx_tri_cent[1]+idx_atom_opt[jatom+1:iatom] \
                            +idx_tri_cent[0]+idx_atom_opt[iatom+1:])
    for jatom in xrange(iatom+1,num_int_cent):
        list_idx_int.append(idx_atom_opt[0:iatom]+idx_tri_cent[0]+idx_atom_opt[iatom+1:jatom] \
                            +idx_tri_cent[1]+idx_atom_opt[jatom+1:])
#
list_idx_int.append(idx_tri_cent[0]+idx_tri_cent[0]+idx_tri_cent[1])
list_idx_int.append(idx_tri_cent[0]+idx_tri_cent[1]+idx_tri_cent[0])
list_idx_int.append(idx_tri_cent[0]+idx_tri_cent[1]+idx_tri_cent[1])
list_idx_int.append(idx_tri_cent[1]+idx_tri_cent[0]+idx_tri_cent[0])
list_idx_int.append(idx_tri_cent[1]+idx_tri_cent[0]+idx_tri_cent[1])
list_idx_int.append(idx_tri_cent[1]+idx_tri_cent[1]+idx_tri_cent[0])
#
list_idx_int.append(idx_tri_cent[0]+idx_tri_cent[1]+idx_tri_cent[2])
list_idx_int.append(idx_tri_cent[0]+idx_tri_cent[2]+idx_tri_cent[1])
list_idx_int.append(idx_tri_cent[1]+idx_tri_cent[0]+idx_tri_cent[2])
list_idx_int.append(idx_tri_cent[1]+idx_tri_cent[2]+idx_tri_cent[0])
list_idx_int.append(idx_tri_cent[2]+idx_tri_cent[0]+idx_tri_cent[1])
list_idx_int.append(idx_tri_cent[2]+idx_tri_cent[1]+idx_tri_cent[0])
# non-atomic operator
list_idx_int.append(idx_non_opt)
# generates the possible indices of bra, ket and operator centers
for iatom in xrange(num_int_cent-1):
    list_idx_int.append(idx_non_opt[0:iatom]+idx_tri_cent[0]+idx_non_opt[iatom+1:])
#
for iatom in xrange(num_int_cent-1):
    for jatom in xrange(iatom):
        list_idx_int.append(idx_non_opt[0:jatom]+idx_tri_cent[1]+idx_non_opt[jatom+1:iatom] \
                            +idx_tri_cent[0]+idx_non_opt[iatom+1:])
    for jatom in xrange(iatom+1,num_int_cent-1):
        list_idx_int.append(idx_non_opt[0:iatom]+idx_tri_cent[0]+idx_non_opt[iatom+1:jatom] \
                            +idx_tri_cent[1]+idx_non_opt[jatom+1:])
# number of three-center geometric derivatives
num_tot_geo = len(list_idx_int)
num_tri_geo = num_tot_geo-(num_one_geo+num_two_geo)

# number of all possible terms
num_one_terms = num_one_geo*num_poss_geo
num_two_terms = num_one_terms+num_two_geo*num_poss_geo
num_tot_terms = num_tot_geo*num_poss_geo

# number of total geometric derivatives
header_file.write("  ! number of total geometric derivatives, according to\n")
header_file.write("  ! (1) numbers of differentiated centers\n")
header_file.write("  ! (2) indices of differentiated centers\n")
header_file.write("  ! (3) different indices of bra, ket and operator centers\n")
header_file.write("  ! (4) orders of partial geoemtric derivatives on bra, ket and operator centers\n")
header_file.write("  integer, parameter :: NUM_TOT_GEO = "+repr(num_tot_terms)+"\n")

num_term = 6                         #how many terms per line
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
            dump_str = dump_str+repr(1).rjust(5)+",      "
        elif (jterm<num_two_terms):
            dump_str = dump_str+repr(2).rjust(5)+",      "
        else:
            dump_str = dump_str+repr(3).rjust(5)+",      "
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    if (jterm<num_one_terms):
        dump_str = dump_str+repr(1).rjust(5)+",      "
    elif (jterm<num_two_terms):
        dump_str = dump_str+repr(2).rjust(5)+",      "
    else:
        dump_str = dump_str+repr(3).rjust(5)+",      "
if (num_one_terms==num_tot_terms):
  dump_str = dump_str+repr(1).rjust(5)+"/)\n"
elif (num_two_terms==num_tot_terms):
  dump_str = dump_str+repr(2).rjust(5)+"/)\n"
else:
  dump_str = dump_str+repr(3).rjust(5)+"/)\n"
header_file.write(dump_str)

# writes indices of differentiated centers
header_file.write("  ! indices of differentiated centers\n")
header_file.write("  integer, parameter :: IDX_GEO_CENT(NUM_TOT_GEO*3) = (/                 &\n")
str_idx_one = repr(idx_one_cent[0]).rjust(5)+","+repr(0).rjust(2)+","+repr(0).rjust(2)
str_idx_two = repr(idx_two_cent[0][0]).rjust(5)+","+repr(idx_two_cent[1][0]).rjust(2) \
              +","+repr(0).rjust(2)
str_idx_tri = repr(idx_tri_cent[0][0]).rjust(5)+","+repr(idx_tri_cent[1][0]).rjust(2) \
              +","+repr(idx_tri_cent[2][0]).rjust(2)
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        if (jterm<num_one_terms):
            dump_str = dump_str+str_idx_one+","
        elif (jterm<num_two_terms):
            dump_str = dump_str+str_idx_two+","
        else:
            dump_str = dump_str+str_idx_tri+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    if (jterm<num_one_terms):
        dump_str = dump_str+str_idx_one+","
    elif (jterm<num_two_terms):
        dump_str = dump_str+str_idx_two+","
    else:
        dump_str = dump_str+str_idx_tri+","
if (num_one_terms==num_tot_terms):
  dump_str = dump_str+str_idx_one+"/)\n"
elif (num_two_terms==num_tot_terms):
  dump_str = dump_str+str_idx_two+"/)\n"
else:
  dump_str = dump_str+str_idx_tri+"/)\n"
header_file.write(dump_str)

# writes the orders of differentiated centers
header_file.write("  ! orders of differentiated centers\n")
header_file.write("  integer, parameter :: ORDER_GEO_CENT(NUM_TOT_GEO*3) = (/               &\n")
str_order_one = repr(order_one_cent).rjust(5)+","+repr(0).rjust(2)+","+repr(0).rjust(2)
str_order_two = repr(order_two_cent[0]).rjust(5)+","+repr(order_two_cent[1]).rjust(2) \
                +","+repr(0).rjust(2)
str_order_tri = repr(order_tri_cent[0]).rjust(5)+","+repr(order_tri_cent[1]).rjust(2) \
                +","+repr(order_tri_cent[2]).rjust(2)
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        if (jterm<num_one_terms):
            dump_str = dump_str+str_order_one+","
        elif (jterm<num_two_terms):
            dump_str = dump_str+str_order_two+","
        else:
            dump_str = dump_str+str_order_tri+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    if (jterm<num_one_terms):
        dump_str = dump_str+str_order_one+","
    elif (jterm<num_two_terms):
        dump_str = dump_str+str_order_two+","
    else:
        dump_str = dump_str+str_order_tri+","
if (num_one_terms==num_tot_terms):
  dump_str = dump_str+str_order_one+"/)\n"
elif (num_two_terms==num_tot_terms):
  dump_str = dump_str+str_order_two+"/)\n"
else:
  dump_str = dump_str+str_order_tri+"/)\n"
header_file.write(dump_str)

# writes the indices of bra and ket centers
header_file.write("  ! indices of bra and ket centers\n")
header_file.write("  integer, parameter :: IDX_BK_OPT(NUM_TOT_GEO*3) = (/                   &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        dump_str = dump_str+repr(list_idx_int[jterm/num_poss_geo][0]).rjust(5)+"," \
                           +repr(list_idx_int[jterm/num_poss_geo][1]).rjust(2)+"," \
                           +repr(list_idx_int[jterm/num_poss_geo][2]).rjust(2)+","
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+repr(list_idx_int[jterm/num_poss_geo][0]).rjust(5)+"," \
                       +repr(list_idx_int[jterm/num_poss_geo][1]).rjust(2)+"," \
                       +repr(list_idx_int[jterm/num_poss_geo][2]).rjust(2)+","
dump_str = dump_str+repr(list_idx_int[num_tot_geo-1][0]).rjust(5)+"," \
                   +repr(list_idx_int[num_tot_geo-1][1]).rjust(2)+"," \
                   +repr(list_idx_int[num_tot_geo-1][2]).rjust(2)+"/)\n"
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
        dump_str = dump_str+str_order_bra[jterm%num_poss_geo]+",      "
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+str_order_bra[jterm%num_poss_geo]+",      "
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
        dump_str = dump_str+str_order_ket[jterm%num_poss_geo]+",      "
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+str_order_ket[jterm%num_poss_geo]+",      "
dump_str = dump_str+str_order_ket[num_poss_geo-1]+"/)\n"
header_file.write(dump_str)

# writes the orders of partial geoemtric derivatives on operator center
header_file.write("  ! orders of partial geoemtric derivatives on operator center\n")
header_file.write("  integer, parameter :: ORDER_GEO_OPT(NUM_TOT_GEO) = (/                  &\n")
str_order_opt = []
for igeo in xrange(num_poss_geo):
    str_order_opt.append(repr(order_geo_opt[igeo]).rjust(5))
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = ""
    for jterm in xrange(iterm,iterm+num_term):
        dump_str = dump_str+str_order_opt[jterm%num_poss_geo]+",      "
    header_file.write(dump_str+" &\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = ""
for jterm in xrange(iterm,num_tot_terms-1):
    dump_str = dump_str+str_order_opt[jterm%num_poss_geo]+",      "
dump_str = dump_str+str_order_opt[num_poss_geo-1]+"/)\n"
header_file.write(dump_str)

# determines the partial derivative terms
zero_ints = []
neg_one = []
order_part_cent = []
scatter_deriv = []
seq_part_geo = []
# one-center geometric derivatives
for iterm in xrange(num_one_geo):
    for jterm in xrange(num_poss_geo):
        num_geo_cent = list_idx_int[iterm].count(idx_one_cent[0])
        # differentitated center is neither bra, ket nor operator center;
        # or bra, ket and operator centers are identical
        if (num_geo_cent==0 or num_geo_cent==3):
            zero_ints.append(True)
            neg_one.append(False)
            order_part_cent.append("    0, 0, 0")
            scatter_deriv.append(False)
            seq_part_geo.append(repr(0).rjust(5)+","+repr(0).rjust(2)+","+repr(0).rjust(2))
        else:
            zero_ints.append(False)
            # either bra, ket or operator center is the differentiated center
            if (num_geo_cent==1):
                # finds the differentiated center
                place_geo_cent = list_idx_int[iterm].index(idx_one_cent[0])
                neg_one.append(False)
            # two identical centers are the differentiated center, using translational invariance
            else:
                # finds the center which is not the differentiated one
                for iatom in xrange(num_int_cent):
                    if (list_idx_int[iterm][iatom]!=idx_one_cent[0]):
                        place_geo_cent = iatom
                        break
                neg_one.append(order_one_cent%2==1)
            # (non-)differentitated center is bra center
            if (place_geo_cent==0):
                order_part_cent.append(repr(order_geo_bra[jterm]+order_one_cent).rjust(5)+"," \
                                       +repr(order_geo_ket[jterm]).rjust(2)+","               \
                                       +repr(order_geo_opt[jterm]).rjust(2))
                scatter_deriv.append(order_geo_bra[jterm]>0 or order_geo_ket[jterm]>0 or \
                                     order_geo_opt[jterm]>0)
            # (non-)differentitated center is ket center
            elif (place_geo_cent==1):
                order_part_cent.append(repr(order_geo_bra[jterm]).rjust(5)+","                 \
                                       +repr(order_geo_ket[jterm]+order_one_cent).rjust(2)+"," \
                                       +repr(order_geo_opt[jterm]).rjust(2))
                scatter_deriv.append(order_geo_ket[jterm]>0 or order_geo_opt[jterm]>0)
            # (non-)differentitated center is operator center
            else:
                order_part_cent.append(repr(order_geo_bra[jterm]).rjust(5)+","                 \
                                       +repr(order_geo_ket[jterm]).rjust(2)+"," \
                                       +repr(order_geo_opt[jterm]+order_one_cent).rjust(2))
                scatter_deriv.append(order_geo_opt[jterm]>0)
            seq_part_geo.append(repr(place_geo_cent+1).rjust(5)+","+repr(0).rjust(2) \
                                +","+repr(0).rjust(2))
# two-center geometric derivatives
for iterm in xrange(num_one_geo,num_one_geo+num_two_geo):
    for jterm in xrange(num_poss_geo):
        num_geo_first = list_idx_int[iterm].count(idx_two_cent[0][0])
        num_geo_second = list_idx_int[iterm].count(idx_two_cent[1][0])
        # at least, one differentitated center is neither bra, ket nor operator center
        if (num_geo_first==0 or num_geo_second==0):
            zero_ints.append(True)
            neg_one.append(False)
            order_part_cent.append("    0, 0, 0")
            scatter_deriv.append(False)
            seq_part_geo.append(repr(0).rjust(5)+","+repr(0).rjust(2)+","+repr(0).rjust(2))
        else:
            zero_ints.append(False)
            # there is no identical centers
            if (num_geo_first==1 and num_geo_second==1):
                neg_one.append(False)
                place_first = list_idx_int[iterm].index(idx_two_cent[0][0])
                place_second = list_idx_int[iterm].index(idx_two_cent[1][0])
                # first differentitated center is bra center
                if (place_first==0):
                    # second differentitated center is ket center
                    if (place_second==1):
                        order_part_cent.append(repr(order_geo_bra[jterm]+order_two_cent[0]).rjust(5)+","  \
                                               +repr(order_geo_ket[jterm]+order_two_cent[1]).rjust(2)+"," \
                                               +repr(order_geo_opt[jterm]).rjust(2))
                    # second differentitated center is operator center
                    else:
                        order_part_cent.append(repr(order_geo_bra[jterm]+order_two_cent[0]).rjust(5)+"," \
                                               +repr(order_geo_ket[jterm]).rjust(2)+","                  \
                                               +repr(order_geo_opt[jterm]+order_two_cent[1]).rjust(2))
                    scatter_deriv.append(order_geo_bra[jterm]>0 or order_geo_ket[jterm]>0 or \
                                         order_geo_opt[jterm]>0)
                # first differentitated center is ket center
                elif (place_first==1):
                    # second differentitated center is bra center
                    if (place_second==0):
                        order_part_cent.append(repr(order_geo_bra[jterm]+order_two_cent[1]).rjust(5)+","  \
                                               +repr(order_geo_ket[jterm]+order_two_cent[0]).rjust(2)+"," \
                                               +repr(order_geo_opt[jterm]).rjust(2))
                        scatter_deriv.append(True)
                    # second differentitated center is operator center
                    else:
                        order_part_cent.append(repr(order_geo_bra[jterm]).rjust(5)+","                    \
                                               +repr(order_geo_ket[jterm]+order_two_cent[0]).rjust(2)+"," \
                                               +repr(order_geo_opt[jterm]+order_two_cent[1]).rjust(2))
                        scatter_deriv.append(order_geo_ket[jterm]>0 or order_geo_opt[jterm]>0)
                # first differentitated center is operator center
                else:
                    # second differentitated center is bra center
                    if (place_second==0):
                        order_part_cent.append(repr(order_geo_bra[jterm]+order_two_cent[1]).rjust(5)+"," \
                                               +repr(order_geo_ket[jterm]).rjust(2)+","                  \
                                               +repr(order_geo_opt[jterm]+order_two_cent[0]).rjust(2))
                    # second differentitated center is ket center
                    else:
                        order_part_cent.append(repr(order_geo_bra[jterm]).rjust(5)+","                    \
                                               +repr(order_geo_ket[jterm]+order_two_cent[1]).rjust(2)+"," \
                                               +repr(order_geo_opt[jterm]+order_two_cent[0]).rjust(2))
                    scatter_deriv.append(True)
                seq_part_geo.append(repr(place_first+1).rjust(5)+","+repr(place_second+1).rjust(2) \
                                    +","+repr(0).rjust(2))
            # there are identical centers which are the first or second differentiated center
            else:
                # there are identical centers which are the first differentiated center
                # performs geometric derivatives on the second differentiated center
                # by translational invariance
                if (num_geo_first==2):
                    place_geo_cent = list_idx_int[iterm].index(idx_two_cent[1][0])
                    neg_one.append(order_two_cent[0]%2==1)
                # there are identical centers which are the second differentiated center
                # performs geometric derivatives on the first differentiated center
                # by translational invariance
                else:
                    place_geo_cent = list_idx_int[iterm].index(idx_two_cent[0][0])
                    neg_one.append(order_two_cent[1]%2==1)
                # differentiated center is bra center
                if (place_geo_cent==0):
                    order_part_cent.append(repr(order_geo_bra[jterm]+sum(order_two_cent)).rjust(5)+"," \
                                           +repr(order_geo_ket[jterm]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]).rjust(2))
                # differentiated center is ket center
                elif (place_geo_cent==1):
                    order_part_cent.append(repr(order_geo_bra[jterm]).rjust(5)+","                      \
                                           +repr(order_geo_ket[jterm]+sum(order_two_cent)).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]).rjust(2))
                # differentiated center is operator center
                else:
                    order_part_cent.append(repr(order_geo_bra[jterm]).rjust(5)+","  \
                                           +repr(order_geo_ket[jterm]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]+sum(order_two_cent)).rjust(2))
                scatter_deriv.append(True)
                seq_part_geo.append(repr(place_geo_cent+1).rjust(5)+","+repr(place_geo_cent+1).rjust(2) \
                                    +","+repr(0).rjust(2))
# three-center geometric derivatives
for iterm in xrange(num_one_geo+num_two_geo,num_tot_geo):
    for jterm in xrange(num_poss_geo):
        num_geo_first = list_idx_int[iterm].count(idx_tri_cent[0][0])
        num_geo_second = list_idx_int[iterm].count(idx_tri_cent[1][0])
        num_geo_third = list_idx_int[iterm].count(idx_tri_cent[2][0])
        # at least, one differentitated center is neither bra, ket nor operator center
        if (num_geo_first==0 or num_geo_second==0 or num_geo_third==0):
            zero_ints.append(True)
            neg_one.append(False)
            order_part_cent.append("    0, 0, 0")
            scatter_deriv.append(False)
            seq_part_geo.append(repr(0).rjust(5)+","+repr(0).rjust(2)+","+repr(0).rjust(2))
        else:
            zero_ints.append(False)
            neg_one.append(False)
            place_first = list_idx_int[iterm].index(idx_tri_cent[0][0])
            place_second = list_idx_int[iterm].index(idx_tri_cent[1][0])
            place_third = list_idx_int[iterm].index(idx_tri_cent[2][0])
            # first differentitated center is bra center
            if (place_first==0):
                # second differentitated center is ket center
                if (place_second==1):
                    order_part_cent.append(repr(order_geo_bra[jterm]+order_tri_cent[0]).rjust(5)+","  \
                                           +repr(order_geo_ket[jterm]+order_tri_cent[1]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]+order_tri_cent[2]).rjust(2))
                    scatter_deriv.append(order_geo_bra[jterm]>0 or order_geo_ket[jterm]>0 or \
                                         order_geo_opt[jterm]>0)
                # second differentitated center is operator center
                else:
                    order_part_cent.append(repr(order_geo_bra[jterm]+order_tri_cent[0]).rjust(5)+","  \
                                           +repr(order_geo_ket[jterm]+order_tri_cent[2]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]+order_tri_cent[1]).rjust(2))
                    scatter_deriv.append(True)
            # first differentitated center is ket center
            elif (place_first==1):
                # second differentitated center is bra center
                if (place_second==0):
                    order_part_cent.append(repr(order_geo_bra[jterm]+order_tri_cent[1]).rjust(5)+","  \
                                           +repr(order_geo_ket[jterm]+order_tri_cent[0]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]+order_tri_cent[2]).rjust(2))
                # second differentitated center is operator center
                else:
                    order_part_cent.append(repr(order_geo_bra[jterm]+order_tri_cent[2]).rjust(5)+","  \
                                           +repr(order_geo_ket[jterm]+order_tri_cent[0]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]+order_tri_cent[1]).rjust(2))
                scatter_deriv.append(True)
            # first differentitated center is operator center
            else:
                # second differentitated center is bra center
                if (place_second==0):
                    order_part_cent.append(repr(order_geo_bra[jterm]+order_tri_cent[1]).rjust(5)+","  \
                                           +repr(order_geo_ket[jterm]+order_tri_cent[2]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]+order_tri_cent[0]).rjust(2))
                # second differentitated center is ket center
                else:
                    order_part_cent.append(repr(order_geo_bra[jterm]+order_tri_cent[2]).rjust(5)+","  \
                                           +repr(order_geo_ket[jterm]+order_tri_cent[1]).rjust(2)+"," \
                                           +repr(order_geo_opt[jterm]+order_tri_cent[0]).rjust(2))
                scatter_deriv.append(True)
            seq_part_geo.append(repr(place_first+1).rjust(5)+","+repr(place_second+1).rjust(2) \
                                +","+repr(place_third+1).rjust(2))

# writes the left information
header_file.write("  ! referenced results indicating if the total geometric derivatives are zero\n")
header_file.write("  logical, parameter :: REF_ZERO_INTS(NUM_TOT_GEO) = (/                     &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = "    "
    for jterm in xrange(iterm,iterm+num_term):
        if (zero_ints[jterm]):
            dump_str = dump_str+".true.,     "
        else:
            dump_str = dump_str+".false.,    "
    header_file.write(dump_str+"&\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = "    "
for jterm in xrange(iterm,num_tot_terms-1):
    if (zero_ints[jterm]):
        dump_str = dump_str+".true.,     "
    else:
        dump_str = dump_str+".false.,    "
if (zero_ints[num_tot_terms-1]):
    dump_str = dump_str+".true./)\n"
else:
    dump_str = dump_str+".false./)\n"
header_file.write(dump_str)
#
header_file.write("  ! referenced results indicating if the integrals will be multiplied by -1\n")
header_file.write("  logical, parameter :: REF_NEG_ONE(NUM_TOT_GEO) = (/                       &\n")
for iterm in xrange(0,num_tot_terms-left_terms,num_term):
    dump_str = "    "
    for jterm in xrange(iterm,iterm+num_term):
        if (neg_one[jterm]):
            dump_str = dump_str+".true.,     "
        else:
            dump_str = dump_str+".false.,    "
    header_file.write(dump_str+"&\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = "    "
for jterm in xrange(iterm,num_tot_terms-1):
    if (neg_one[jterm]):
        dump_str = dump_str+".true.,     "
    else:
        dump_str = dump_str+".false.,    "
if (neg_one[num_tot_terms-1]):
    dump_str = dump_str+".true./)\n"
else:
    dump_str = dump_str+".false./)\n"
header_file.write(dump_str)
#
header_file.write("  ! referenced final orders of partial geometric derivatives by adding\n")
header_file.write("  ! those from total geometric derivatives\n")
header_file.write("  integer, parameter :: REF_ORDER_PART(NUM_TOT_GEO*3) = (/               &\n")
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
            dump_str = dump_str+".true.,     "
        else:
            dump_str = dump_str+".false.,    "
    header_file.write(dump_str+"&\n")
# left terms in the last line
iterm = num_tot_terms-left_terms
dump_str = "    "
for jterm in xrange(iterm,num_tot_terms-1):
    if (scatter_deriv[jterm]):
        dump_str = dump_str+".true.,     "
    else:
        dump_str = dump_str+".false.,    "
if (scatter_deriv[num_tot_terms-1]):
    dump_str = dump_str+".true./)\n"
else:
    dump_str = dump_str+".false./)\n"
header_file.write(dump_str)
#
header_file.write("  ! referenced sequences of bra and ket centers for partial derivative terms\n")
header_file.write("  integer, parameter :: REF_SEQ_PART(NUM_TOT_GEO*3) = (/                 &\n")
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

# closes geom_part_one.h
header_file.close()
print ">> Please replace test_f90/geom_part_one.h with this generated one ..."
