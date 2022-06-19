#!/usr/bin/env python
#
#   Polarizable Embedding (PE) library
#   Copyright (C) 2013 , 2014 The PE library developers. See the CONTRIBUTORS file
#                             in the top-level directory of this distribution.
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
import math
import numpy as np
import time
import multiprocessing as mp

parser = ap.ArgumentParser(description='Cube analyzer 0.1',
                           epilog='Have a nice day :-)',
                           usage='%(prog)s [options]',
                           fromfile_prefix_chars='@')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('-i', dest='inputfiles', nargs='+', metavar='INPUT_FILE',
                    default=[],
                    help='''Specify the names of the cube files''')
parser.add_argument('-o', dest='outputfile', metavar='OUTPUT_FILE',
                    default='',
                    help='''Specify the name of the output file
                            [default: %(default)s]''')
parser.add_argument('-ref', dest='refidx', default=0, type=int,
                    metavar='REFCUBE',
                    help='''Specify the index of the reference cube to use in
                            different analysis. [default: %(default)s]''')
parser.add_argument('-add', dest='addlist', nargs='+', default=[], type=int,
                    metavar=('CUBE1', 'CUBE2'),
                    help='''Specify which cubes to add. The cubes are numbered
                            according to the input order starting from 0''')
parser.add_argument('-sub', dest='sublist', nargs='+', default=[], type=int,
                    metavar=('CUBE1', 'CUBE2'),
                    help='''Specify which cubes to subtract. The cubes are
                            numbered according to the input order starting
                            from 0''')
parser.add_argument('-mae', dest='mae', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a MAE analysis. MAE is calculated relative to
                            the reference cube in volumes between MIN times
                            vdw radius and MAX times vdw radius in STEP steps.''')
parser.add_argument('-rmsd', dest='rmsd', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a RMSD analysis. RMSD is calculated relative to
                            the reference cube in volumes between MIN times
                            vdw radius and MAX times vdw radius in STEP steps''')
parser.add_argument('-field-rmsd', dest='field_rmsd', nargs=3, default=[], type=float,
                    metavar=('MIN','MAX','STEP'),
                    help='''Do a electric field RMSD analysis. RMSD is calculated
                            relative to the reference cubes in volumes between MIN times
                            vdw radius and MAX times vdw radius in STEP steps''')
parser.add_argument('-ps', '--print-shells', dest='printshells', nargs=3, default=[],
                    type=float, metavar=('MIN','MAX','STEP'),
                    help='''Print values of all points within volumes
                            volumes between MIN times vdw radius
                            and MAX times vdw radius in STEP steps.''')
parser.add_argument('-log-prefix', dest='prefix', metavar='PREFIX',
                    default='',
                    help='''Specify prefix on log files [default: %(default)s]''')

parser.add_argument('-ncores', dest='ncores', default=1, type=int,
                    metavar='NCORES',
                    help='''Number of cores [default: %(default)s]''')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()
args = parser.parse_args()
if not args.inputfiles:
    parser.print_help()
    sys.exit()

charge2radius = { 1.0: 1.20,  2.0: 1.40,  3.0: 2.20,  4.0: 1.90,  5.0: 1.80,
                  6.0: 1.70,  7.0: 1.60,  8.0: 1.55,  9.0: 1.50, 10.0: 1.54,
                 11.0: 2.40, 12.0: 2.20, 13.0: 2.10, 14.0: 2.10, 15.0: 1.95,
                 16.0: 1.80, 17.0: 1.80, 18.0: 1.88, 19.0: 1.90,  0.0: 1.00}

aa2au = 1.8897261249935897

class Cube(object):

    def __init__(self):
        pass

    def read(self, filename):
        self.filename = filename
        self.readheader()
        self.readgrid()

    def readheader(self):
        fcube = open(self.filename, 'r')
        line = fcube.readline()
        line = fcube.readline()
        line = fcube.readline().split()
        self.natoms = int(line[0])
        self.origo = np.array([float(coord) for coord in line[1:4]])
        line = fcube.readline().split()
        self.xpoints = int(line[0])
        self.xstep = float(line[1])
        line = fcube.readline().split()
        self.ypoints = int(line[0])
        self.ystep = float(line[2])
        line = fcube.readline().split()
        self.zpoints = int(line[0])
        self.zstep = float(line[3])
        self.charges = []
        self.coords = []
        for i in xrange(self.natoms):
            line = fcube.readline().split()
            self.charges.append(float(line[1]))
            self.coords.append([float(coord) for coord in line[2:5]])
        self.grid = np.zeros((self.xpoints, self.ypoints, self.zpoints))
        fcube.close()

    def copyheader(self, cube):
        self.natoms = cube.natoms
        self.origo = cube.origo
        self.xpoints = cube.xpoints
        self.xstep = cube.xstep
        self.ypoints = cube.ypoints
        self.ystep = cube.ystep
        self.zpoints = cube.zpoints
        self.zstep = cube.zstep
        self.charges = cube.charges
        self.coords = cube.coords
        self.grid = np.zeros((self.xpoints, self.ypoints, self.zpoints))

    def readgrid(self):
        fcube = open(self.filename, 'r')
        line = fcube.readline()
        line = fcube.readline()
        line = fcube.readline().split()
        line = fcube.readline().split()
        line = fcube.readline().split()
        line = fcube.readline().split()
        for i in xrange(self.natoms):
            line = fcube.readline().split()
        data = []
        while True:
            line = [float(x) for x in fcube.readline().split()]
            if len(line) == 0:
                break
            data += line
        self.grid = np.array(data).reshape([self.xpoints, self.ypoints, self.zpoints])
        # z = 0
        # for x in xrange(self.xpoints):
        #     for y in xrange(self.ypoints):
        #         for i in xrange(int(math.ceil(float(self.zpoints) / 6.0))):
        #             line = fcube.readline().split()
        #             for value in line:
        #                 self.grid[x][y][z] = float(value)
        #                 z += 1
        #         z = 0
        fcube.close()

    def writecube(self, filename):
        self.filename = filename
        fcube = open(self.filename, 'w')
        header = ''
        header += 'CUBE file\n'
        header += 'Generated by the CUBE Analyzer 0.1\n'
        header += '{0:5d}'.format(self.natoms)
        header += '{0[0]:12.6f}{0[1]:12.6f}{0[2]:12.6f}\n'.format(self.origo)
        header += '{0:5d}'.format(self.xpoints)
        header += '{0:12.6f}{1:12.6f}{1:12.6f}\n'.format(self.xstep, 0.0)
        header += '{0:5d}'.format(self.ypoints)
        header += '{1:12.6f}{0:12.6f}{1:12.6f}\n'.format(self.ystep, 0.0)
        header += '{0:5d}'.format(self.zpoints)
        header += '{1:12.6f}{1:12.6f}{0:12.6f}\n'.format(self.zstep, 0.0)
        for charge, coord in zip(self.charges, self.coords):
            header += '{0:5d}'.format(int(charge))
            header += '{0:12.6f}'.format(charge)
            header += '{0[0]:12.6f}{0[1]:12.6f}{0[2]:12.6f}\n'.format(coord)
        fcube.write(header)
        grid = ''
        for x in xrange(self.xpoints):
            for y in xrange(self.ypoints):
                i = 0
                for z in xrange(self.zpoints):
                    i += 1
                    grid += '{0:13.5e}'.format(self.grid[x][y][z])
                    if i == 6:
                        grid += '\n'
                        i = 0
                if grid[-1] != '\n':
                    grid += '\n'
            fcube.write(grid)
            grid = ''
        fcube.close()

def rmsd_analysis(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    rmsds = []
    for i in range(len(cubelist)):
        rmsds.append([])
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            rmsds[i].append(0.0)
        grdpts.append(0)
        shells.append([round(inner, 4), round(outer, 4)])
        inner += step
        outer += step
    points = []
    point = [0.0, 0.0, 0.0]
    for x in xrange(refcub.xpoints):
        point[0] = refcub.origo[0] + x * refcub.xstep
        for y in xrange(refcub.ypoints):
            point[1] = refcub.origo[1] + y * refcub.ystep
            for z in xrange(refcub.zpoints):
                point[2] = refcub.origo[2] + z * refcub.zstep
                points.append([list(point), x, y, z])
    out_queue = mp.Queue()
    nprocs = args.ncores
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=rmsd_worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, rmsds, grdpts, refcub, cubelist, out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for rmsd, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                rmsds[j][i] += rmsd[j][i]
    for ic, cube in enumerate(cubelist):
        if args.prefix:
            name = '{0}_{1}'.format(args.prefix, cube.filename[:-5])
        else:
            name = '{0}'.format(cube.filename[:-5])
        frmsd = open('{0}_rmsd.log'.format(name), 'w')
        rmsd = 'Reference: {}\n'.format(refcub.filename)
        rmsd += '{}\n'.format(cube.filename)
        rmsd += ' Points  vdW Volume  Midpoint  RMSD / a.u.  RMSD (kJ/mol)\n'
        for ish, shell in enumerate(shells):
            inner = shell[0]
            outer = shell[1]
            rmsd += '{0:8d} '.format(grdpts[ish])
            rmsd += '{0:5.2f}-{1:<5.2f}  '.format(inner, outer)
            rmsd += '{0:5.2f}  '.format(round(inner + 0.5 * step, 4))
            rmsd += '{0:12.4e}    '.format(math.sqrt(rmsds[ic][ish] / grdpts[ish]))
            rmsd += '{0:8.4f}\n'.format(math.sqrt(rmsds[ic][ish] / grdpts[ish]) * 2625.5)
        frmsd.write(rmsd)
        frmsd.close()

def rmsd_worker(points, shells, rmsds, grdpts, refcub, cubelist, out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    rmsds[ic][ish] += (cube.grid[x][y][z] -
                                       refcub.grid[x][y][z])**2
                grdpts[ish] += 1
                break
    out_queue.put((rmsds, grdpts))

def field_rmsd_analysis(cubes, refcubs, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    norm_rmsds = []
    angle_rmsds = []
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        norm_rmsds.append(0.0)
        angle_rmsds.append(0.0)
        grdpts.append(0)
        shells.append([round(inner, 4), round(outer, 4)])
        inner += step
        outer += step
    points = []
    point = [0.0, 0.0, 0.0]
    for x in xrange(refcubs[0].xpoints):
        point[0] = refcubs[0].origo[0] + x * refcubs[0].xstep
        for y in xrange(refcubs[0].ypoints):
            point[1] = refcubs[0].origo[1] + y * refcubs[0].ystep
            for z in xrange(refcubs[0].zpoints):
                point[2] = refcubs[0].origo[2] + z * refcubs[0].zstep
                points.append([list(point), x, y, z])
    out_queue = mp.Queue()
    nprocs = args.ncores
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=field_rmsd_worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, norm_rmsds, angle_rmsds, grdpts, refcubs, cubes, out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for norm_rmsd, angle_rmsd, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            norm_rmsds[i] += norm_rmsd[i]
            angle_rmsds[i] += angle_rmsd[i]
    if args.prefix:
        name = '{0}_{1}'.format(args.prefix, cube.filename[:-5])
    else:
        name = '{0}'.format(cube.filename[:-5])
    frmsd = open('mef_rmsds.log', 'w')
    text = 'References: {0}, {1} and {2}\n'.format(refcubs[0].filename,
                                                   refcubs[1].filename,
                                                   refcubs[2].filename)
    text += '{0}, {1} and {2}\n'.format(cubes[0].filename,
                                        cubes[1].filename,
                                        cubes[2].filename)
    text += ' Points  vdW Volume  Midpoint  Norm RMSD / a.u.  Norm RMSD / V/Aa  Angle RMSD / rad.  Angle RMSD / deg.\n'
    for ish, shell in enumerate(shells):
        inner = shell[0]
        outer = shell[1]
        text += '{0:8d} '.format(grdpts[ish])
        text += '{0:5.2f}-{1:<5.2f}  '.format(inner, outer)
        text += '{0:5.2f}    '.format(round(inner + 0.5 * step, 4))
        text += '{0:12.4e}         '.format(math.sqrt(norm_rmsds[ish] / grdpts[ish]))
        text += '{0:8.4f}        '.format(math.sqrt(norm_rmsds[ish] / grdpts[ish]) * 51.4220652)
        text += '{0:12.4e}        '.format(math.sqrt(angle_rmsds[ish] / grdpts[ish]))
        text += '{0:8.4f}\n'.format(math.sqrt(angle_rmsds[ish] / grdpts[ish]) * (180.0 / np.pi))
    frmsd.write(text)
    frmsd.close()

def field_rmsd_worker(points, shells, norm_rmsds, angle_rmsds, grdpts, refcubs, cubes, out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcubs[0].coords, refcubs[0].charges)
            if include:
                cubvec = np.array([cubes[0].grid[x][y][z],
                                 cubes[1].grid[x][y][z],
                                 cubes[2].grid[x][y][z]])
                refvec = np.array([refcubs[0].grid[x][y][z],
                                 refcubs[1].grid[x][y][z],
                                 refcubs[2].grid[x][y][z]])
                nrmcub = np.linalg.norm(cubvec)
                nrmref = np.linalg.norm(refvec)
                unit_cubvec = cubvec / nrmcub
                unit_refvec = refvec / nrmref
                angle = np.arccos(np.dot(unit_cubvec, unit_refvec))
                if np.isnan(angle):
                    if (unit_cubvec == unit_refvec).all():
                        angle = 0.0
                    else:
                        angle = np.pi
                angle_rmsds[ish] += angle**2
                norm_rmsds[ish] += (nrmcub - nrmref)**2
                grdpts[ish] += 1
                break
    out_queue.put((norm_rmsds, angle_rmsds, grdpts))

def inorout(point, shell, coords, charges):
    radii = [charge2radius[charge] * aa2au for charge in charges]
    for center, radius in zip(coords, radii):
        r2 = ((point[0] - center[0])**2 +
              (point[1] - center[1])**2 +
              (point[2] - center[2])**2)
        if r2 < (shell[1] * radius)**2 and r2 >= (shell[0] * radius)**2:
            for coord, rad in zip(coords, radii):
                r2 = ((point[0] - coord[0])**2 +
                      (point[1] - coord[1])**2 +
                      (point[2] - coord[2])**2)
                if r2 > (shell[0] * rad)**2:
                    continue
                else:
                    return False
            return True
        else:
            continue
    return False

def mae_analysis(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    maes = []
    for i in range(len(cubelist)):
        maes.append([])
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            maes[i].append(0.0)
        grdpts.append(0)
        shells.append([round(inner, 4), round(outer, 4)])
        inner += step
        outer += step
    points = []
    point = [0.0, 0.0, 0.0]
    for x in xrange(refcub.xpoints):
        point[0] = refcub.origo[0] + x * refcub.xstep
        for y in xrange(refcub.ypoints):
            point[1] = refcub.origo[1] + y * refcub.ystep
            for z in xrange(refcub.zpoints):
                point[2] = refcub.origo[2] + z * refcub.zstep
                points.append([list(point), x, y, z])
    out_queue = mp.Queue()
    nprocs = args.ncores
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=mae_worker,
                          args=(points[chunksize * i:chunksize * (i + 1)],
                                shells, maes, grdpts, refcub, cubelist, out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    for mae, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                maes[j][i] += mae[j][i]
    for ic, cube in enumerate(cubelist):
        if args.prefix:
            name = '{0}_{1}'.format(args.prefix, cube.filename[:-5])
        else:
            name = '{0}'.format(cube.filename[:-5])
        fmae = open('{}_mae.log'.format(name), 'w')
        mae = 'Reference: {}\n'.format(refcub.filename)
        mae += '{}\n'.format(cube.filename)
        mae += ' Points   vdW Volume   Midpoint    MAE\n'
        for ish, shell in enumerate(shells):
            inner = shell[0]
            outer = shell[1]
            mae += '{0:8d} '.format(grdpts[ish])
            mae += '{0:5.2f}-{1:<5.2f} '.format(inner, outer)
            mae += '{0:5.2f} '.format(round(inner + 0.5 * step, 4))
            mae += '{0:12.4e}\n'.format(maes[ic][ish] / grdpts[ish])
        fmae.write(mae)
        fmae.close()

def mae_worker(points, shells, maes, grdpts, refcub, cubelist, out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    maes[ic][ish] += np.abs(cube.grid[x][y][z] -
                                       refcub.grid[x][y][z])
                grdpts[ish] += 1
                break
    out_queue.put((maes, grdpts))

def print_shells(cubelist, refcub, mini, maxi, step):
    if maxi < mini:
        exit('ERROR: maximum vdw radius is less than the specified minimum')
    vals = []
    refs = []
    for i in range(len(cubelist)):
        vals.append([])
        refs.append([])
    grdpts = []
    shells = []
    inner = mini
    outer = mini + step
    while round(outer, 4) <= round(maxi, 4):
        for i in range(len(cubelist)):
            vals[i].append([])
            refs[i].append([])
        grdpts.append(0)
        shells.append([round(inner, 4), round(outer, 4)])
        inner += step
        outer += step
    points = []
    point = [0.0, 0.0, 0.0]
    for x in xrange(refcub.xpoints):
        point[0] = refcub.origo[0] + x * refcub.xstep
        for y in xrange(refcub.ypoints):
            point[1] = refcub.origo[1] + y * refcub.ystep
            for z in xrange(refcub.zpoints):
                point[2] = refcub.origo[2] + z * refcub.zstep
                points.append([list(point), x, y, z])
    out_queue = mp.Queue()
    nprocs = args.ncores
    chunksize = int(math.ceil(len(points) / float(nprocs)))
    procs = []
    for i in xrange(nprocs):
        proc = mp.Process(target=print_shell_worker,
                          args=(i, points[chunksize * i:chunksize * (i + 1)],
                                shells, vals, refs, grdpts, refcub, cubelist,
                                out_queue))
        procs.append(proc)
        proc.start()
    results = []
    for i in xrange(nprocs):
        results.append(out_queue.get())
    for proc in procs:
        proc.join()
    results = sorted(results, key=lambda result: result[0])
    for index, val, ref, pts in results:
        for i in xrange(len(shells)):
            grdpts[i] += pts[i]
            for j in xrange(len(cubelist)):
                for k in xrange(pts[i]):
                    vals[j][i].append(val[j][i][k])
                    refs[j][i].append(ref[j][i][k])
    for ic, cube in enumerate(cubelist):
        for ish, shell in enumerate(shells):
            if args.prefix:
                name = '{0}_{1}'.format(args.prefix, cube.filename[:-5])
            else:
                name = '{0}'.format(cube.filename[:-5])
            fprt = open('{0}_shell_{1}.log'.format(name, ish), 'w')
            prt = 'Reference: {0}\n'.format(refcub.filename)
            prt += '{0}\n'.format(cube.filename)
            prt += 'Inner-outer radii: {0[0]:5.2f}-{0[1]:<5.2f}\n'.format(shell)
            prt += 'Shell midpoint: {0:5.2f}\n'.format(round(shell[0] + 0.5 * step, 4))
            prt += 'Number of grid points: {0:9d}\n'.format(grdpts[ish])
            fprt.write(prt)
            prt = ''
            for i in xrange(grdpts[ish]):
                prt += '{0:12.4e} {1:12.4e}\n'.format(vals[ic][ish][i],
                                                      refs[ic][ish][i])
            fprt.write(prt)
            fprt.close()

def print_shell_worker(index, points, shells, vals, refs, grdpts, refcub, cubelist,
                       out_queue):
    for point, x, y, z in points:
        for ish, shell in enumerate(shells):
            include = inorout(point, shell, refcub.coords, refcub.charges)
            if include:
                for ic, cube in enumerate(cubelist):
                    vals[ic][ish].append(cube.grid[x][y][z])
                    refs[ic][ish].append(refcub.grid[x][y][z])
                grdpts[ish] += 1
                break
    out_queue.put((index, vals, refs, grdpts))

if __name__ == "__main__":

    if (args.addlist or args.sublist) and (args.mae or args.rmsd or args.printshells):
        exit('ERROR: The options -add and/or -sub cannot be combined with'
             ' -mae, -rmsd or -ps.')

    cubelist = []
    for cubefile in args.inputfiles:
        cube = Cube()
        cube.read(cubefile)
        cubelist.append(cube)

    if args.addlist or args.sublist:
        newcube = Cube()
        newcube.copyheader(cubelist[0])

    if args.addlist:
        for index in args.addlist:
            newcube.grid += cubelist[index].grid

    if args.sublist:
        for index in args.sublist:
            newcube.grid -= cubelist[index].grid

    if args.mae:
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        mae_analysis(cubana, refcub, *args.mae)

    if args.rmsd:
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        rmsd_analysis(cubana, refcub, *args.rmsd)

    if args.field_rmsd:
        refcubs = []
        for i in range(3):
            refcubs.append(cubelist[i])
        cubanas = []
        for i in range(3, 6):
            cubanas.append(cubelist[i])
        field_rmsd_analysis(cubanas, refcubs, *args.field_rmsd)

    if args.printshells:
        refcub = cubelist[args.refidx]
        cubana = []
        for cube in cubelist:
            if cube == refcub:
                continue
            cubana.append(cube)
        print_shells(cubana, refcub, *args.printshells)

    if args.outputfile:
        newcube.writecube(args.outputfile)

