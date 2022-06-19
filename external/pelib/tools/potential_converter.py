#!/usr/bin/env python
#
#   Polarizable Embedding (PE) library
#   Copyright (C) 2013, 2014 The PE library developers. See the CONTRIBUTORS file
#                      in the top-level directory of this distribution.
#
#   This tool is part of the PE library.
#
#   The PE library is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Lesser General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   The PE library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public License
#   along with the PE library. If not, see <http://www.gnu.org/licenses/>.
#
#   Contact information:
#
#   Jogvan Magnus Haugaard Olsen
#   E-mail: foeroyingur@gmail.com
#

import os
import sys
import argparse as ap


parser = ap.ArgumentParser(description='Potential Converter: convert potential files to new format',
                           usage='%(prog)s [options]',
                           fromfile_prefix_chars='@')

parser.add_argument('--version', action='version', version='Potential Converter 1.0')

parser.add_argument('-i', dest='inputfile', metavar='INPUT_FILE',
                    help='''Specify the name of the input file.''')

parser.add_argument('-o', dest='outputfile', metavar='OUTPUT_FILE',
                    help='''Specify the name of the output file.''')

parser.add_argument('--input-type', dest='inputtype', metavar='TYPE',
                    choices=['QMMM'], default='QMMM',
                    help='''Specify the format of the input file. Valid choices
                            are %(choices)s. [default: %(default)s]''')

args = parser.parse_args()

pyver = sys.version_info
if pyver[0] < 2 or pyver[0] == 2 and pyver[1] < 7:
    exit('ERROR: Python >= 2.7 required.')

if not args.inputfile:
    exit('ERROR: no input file specified.')
elif not os.path.isfile(args.inputfile):
    exit('ERROR: input file not found.')
elif not args.outputfile:
    exit('ERROR: no outputfile specified.')
elif os.path.isfile(args.outputfile):
    exit('ERROR: output file already exists.')

inputdir, inputfile = os.path.split(args.inputfile)
outputdir, outputfile = os.path.split(args.outputfile)

charge2elem = { 0:  'X',  1:  'H',  2: 'He',  3: 'Li',  4: 'Be',  5:  'B',
                6:  'C',  7:  'N',  8:  'O',  9:  'F', 10: 'Ne', 11: 'Na',
               12: 'Mg', 13: 'Al', 14: 'Si', 15:  'P', 16:  'S', 17: 'Cl',
               18: 'Ar', 19:  'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23:  'V',
               24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
               30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br',
               36: 'Kr', 37: 'Rb', 38: 'Sr', 39:  'Y', 40: 'Zr', 41: 'Nb',
               42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag',
               48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53:  'I',
               54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr',
               60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
               66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
               72: 'Hf', 73: 'Ta', 74:  'W', 75: 'Re', 76: 'Os', 77: 'Ir',
               78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi',
               84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra'}

au2aa = 0.5291772108

if args.inputtype == 'QMMM':
    fin = open('{}'.format(args.inputfile), 'r')
    line = fin.readline().split()
    unit = line[0]
    line = fin.readline().split()
    nsites = int(line[0])
    mulorder = int(line[1])
    polorder = int(line[2])
    lexlst = int(line[3])
    if len(line) == 5:
        lelems = 1
    elif len(line) == 4:
        lelems = 0
    elems = []
    coords = []
    M0s = []
    M1s = []
    M2s = []
    M3s = []
    M4s = []
    M5s = []
    alphas = []
    exclists = []
    if mulorder == 5:
        pad = 56
    elif mulorder == 4:
        pad = 35
    elif mulorder == 3:
        pad = 20
    elif mulorder == 2:
        pad = 10
    elif mulorder == 1:
        pad = 4
    elif mulorder == 0:
        pad = 1
    for line in fin:
        line = line.split()
        if not line:
            continue
        exclists.append(map(int, line[0:lexlst]))
        coords.append(map(float, line[lexlst+lelems:lexlst+lelems+3]))
        if lelems == 1:
            elems.append(int(line[lexlst]))
        elif lelems == 0:
            elems.append(0)
        if mulorder >= 0:
            M0s.append(float(line[lexlst+lelems+3]))
        if mulorder >= 1:
            M1s.append(map(float, line[lexlst+lelems+4:lexlst+lelems+7]))
        if mulorder >= 2:
            M2s.append(map(float, line[lexlst+lelems+7:lexlst+lelems+13]))
        if mulorder >= 3:
            M3s.append(map(float, line[lexlst+lelems+13:lexlst+lelems+23]))
        if mulorder >= 4:
            M4s.append(map(float, line[lexlst+lelems+23:lexlst+lelems+38]))
        if mulorder >= 5:
            M5s.append(map(float, line[lexlst+lelems+38:lexlst+lelems+59]))
        if polorder == 1:
            isoalpha = float(line[lexlst+lelems+3+pad])
            alphas.append([isoalpha, 0.0, 0.0, isoalpha, 0.0, isoalpha])
        elif polorder == 2:
            start = lexlst+lelems+pad+3
            end = lexlst+lelems+pad+9
            alphas.append(map(float, line[start:end]))
    fin.close()

    if exclists[0][0] == exclists[1][0]:
        exit('ERROR: conversion does not work for this kind of exclusion list')

    for i, excl in enumerate(exclists):
        if excl[0] != i+1:
            exit('ERROR: conversion does not work for this kind of exclusion list')

    if unit == 'AU':
        coords = [[elem * au2aa for elem in coord] for coord in coords]
    ndec = len(str(nsites))
    body = '@COORDINATES\n'
    body += '{}\n'.format(nsites)
    body += 'AA\n'
    for i, coord in enumerate(coords):
        body += '{0} {1[0]:14.8f} {1[1]:14.8f} {1[2]:14.8f}\n'.format(charge2elem[elems[i]], coord)
    if mulorder >= 0:
        body += '@MULTIPOLES\n'
    if mulorder >= 0:
        body += 'ORDER 0\n'
        body += '{0}\n'.format(nsites)
        for i, M0 in enumerate(M0s):
            body += '{0:{1}d} {2:14.8f}\n'.format(i+1, ndec, M0)
    if mulorder >= 1:
        body += 'ORDER 1\n'
        body += '{0}\n'.format(nsites)
        for i, M1 in enumerate(M1s):
            body += '{0:{1}d}'.format(i+1, ndec)
            for j in range(3):
                body += ' {0:14.8f}'.format(M1[j])
            body += '\n'
    if mulorder >= 2:
        body += 'ORDER 2\n'
        body += '{0}\n'.format(nsites)
        for i, M2 in enumerate(M2s):
            body += '{0:{1}d}'.format(i+1, ndec)
            for j in range(6):
                body += ' {0:14.8f}'.format(M2[j])
            body += '\n'
    if mulorder >= 3:
        body += 'ORDER 3\n'
        body += '{0}\n'.format(nsites)
        for i, M3 in enumerate(M3s):
            body += '{0:{1}d}'.format(i+1, ndec)
            for j in range(10):
                body += ' {0:14.8f}'.format(M3[j])
            body += '\n'
    if mulorder >= 4:
        body += 'ORDER 4\n'
        body += '{0}\n'.format(nsites)
        for i, M4 in enumerate(M4s):
            body += '{0:{1}d}'.format(i+1, ndec)
            for j in range(15):
                body += ' {0:14.8f}'.format(M4[j])
            body += '\n'
    if mulorder >= 5:
        body += 'ORDER 5\n'
        body += '{0}\n'.format(nsites)
        for i, M5 in enumerate(M5s):
            body += '{0:{1}d}'.format(i+1, ndec)
            for j in range(21):
                body += ' {0:14.8f}'.format(M5[j])
            body += '\n'
    if polorder >= 1:
        body += '@POLARIZABILITIES\n'
        body += 'ORDER 1 1\n'
        body += '{0}\n'.format(nsites)
        for i, alpha in enumerate(alphas):
            body += '{0:{1}d}'.format(i+1, ndec)
            for j in range(6):
                body += ' {0:14.8f}'.format(alpha[j])
            body += '\n'
        if exclists:
            body += 'EXCLISTS\n'
            body += '{0} {1}\n'.format(nsites, lexlst)
            for i, exclist in enumerate(exclists):
                for ex in exclist:
                    body += ' {0:{1}}'.format(ex, ndec)
                body += '\n'
        body += '\n'
    fout = open('{}'.format(args.outputfile), 'w')
    fout.write(body)
    fout.close()

else:
    exit('ERROR: potential is not supported')
    fin = open('{}'.format(args.inputfile), 'r')
    oldpot = fin.read()
    fin.close()
