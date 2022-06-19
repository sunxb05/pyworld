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
#  This file contains the unit tests of Gen1Int.Tools.
#
#  2010-10-10, Bin Gao: first version

# Logging facility
from logging import getLogger
from logging import basicConfig
from logging import DEBUG as Logger_DEBUG
# Unit testing framework
from unittest import TestCase,TestLoader,TextTestRunner
#from unittest import main as TestMain

__author__ = "Bin Gao, and Andreas J. Thorvaldsen"
__copyright__ = "Copyright 2009-2012"
__credits__ = ["Radovan Bast", "Kenneth Ruud"]
__license__ = "LGPLv3"
__version__ = "0.2.1"
__maintainer__ = "Bin Gao"
__email__ = "bin.gao@uit.no"
__status__ = "Development"

# Set up a specific logger with our desired output level
logger = getLogger("test_tools")

class UnitTests(TestCase):
    """Unit tests"""
    def setUp(self):
        """Initializes unit tests"""
        from numpy import asmatrix
        from numpy.random import randn
        from numpy import trace
        # Sets the minimum number of correct digits
        self.MIN_DIGITS = 8 #12
        self.ang_bra = 10
        self.ang_ket = 12
        self.num_contr_bra = 10
        self.num_contr_ket = 12
        self.num_opt = 6
        # We set the magnetic numbers in the reverse sequence of that in Gen1Int
        self.mag_bra = xrange(self.ang_bra,-self.ang_bra-1,-1)
        self.mag_ket = xrange(self.ang_ket,-self.ang_ket-1,-1)
        num_sgto_bra = 2*self.ang_bra+1
        num_sgto_ket = 2*self.ang_ket+1
        size_sgto_int = num_sgto_bra*self.num_contr_bra*num_sgto_ket*self.num_contr_ket
        self.gint_sgto = asmatrix(xrange(size_sgto_int)).reshape([num_sgto_bra,num_sgto_ket])
        self.ref_sgto = asmatrix(xrange(size_sgto_int-1,-1,-1)).reshape([num_sgto_bra,num_sgto_ket])
        # A random matrix chosen as atomic density matrix
        self.dens_sgto = randn(num_sgto_ket,num_sgto_bra)
        self.ref_tr_sgto = trace(self.gint_sgto*self.dens_sgto)
        # We use Dalton sequence for the powers on the bra
        self.power_bra = []
        for x in xrange(self.ang_bra,-1,-1):
            for y in xrange(self.ang_bra-x,-1,-1):
                self.power_bra.append([x,y,self.ang_bra-(x+y)])
        # We use Gen1Int sequence for the powers on the ket
        self.power_ket = []
        ang_ket1 = self.ang_ket+1
        for z in xrange(ang_ket1):
            for y in xrange(ang_ket1-z):
                self.power_ket.append([self.ang_ket-(y+z),y,z])
        num_cgto_bra = (self.ang_bra+1)*(self.ang_bra+2)/2
        num_cgto_ket = (self.ang_ket+1)*(self.ang_ket+2)/2
        size_cgto_int = num_cgto_bra*num_cgto_ket
        self.gint_cgto = asmatrix(xrange(size_cgto_int)).reshape([num_cgto_bra,num_cgto_ket])
        # A random matrix chosen as atomic density matrix
        self.dens_cgto = randn(num_cgto_ket,num_cgto_bra)
        self.ref_tr_cgto = trace(self.gint_cgto*self.dens_cgto)
        return

    def testReorder(self):
        """Tests reordering integrals"""
        from Gen1Int.Tools import reorder_sgto,reorder_cgto
        # Tests reordering real solid-harmonic Gaussians
        ro_ints = reorder_sgto(self.ang_bra,self.mag_bra,
                               self.ang_ket,self.mag_ket,
                               self.gint_sgto)
        num_sgto_bra = 2*self.ang_bra+1
        num_sgto_ket = 2*self.ang_ket+1
        for i in xrange(num_sgto_bra):
            for j in xrange(num_sgto_ket):
                self.assertAlmostEqual(self.ref_sgto[i,j], \
                                       ro_ints[i,j],self.MIN_DIGITS)
        # Tests reordering Cartesian Gaussians
        ro_ints = reorder_cgto(self.ang_bra,self.power_bra,
                               self.ang_ket,self.power_ket,
                               self.gint_cgto)
        num_cgto_bra = (self.ang_bra+1)*(self.ang_bra+2)/2
        num_cgto_ket = (self.ang_ket+1)*(self.ang_ket+2)/2
        # Index of x^{l}y^{m}z^{n} with l+m+n=angm is 1+m+(2*angm+3-n)*n/2
        ang2_bra3 = 2*self.ang_bra+3
        ang2_ket3 = 2*self.ang_ket+3
        for i in xrange(num_cgto_bra):
            gint_row = self.power_bra[i][1] \
                     + (ang2_bra3-self.power_bra[i][2])*self.power_bra[i][2]/2
            for j in xrange(num_cgto_ket):
                gint_col = self.power_ket[j][1] \
                         + (ang2_ket3-self.power_ket[j][2])*self.power_ket[j][2]/2
                self.assertAlmostEqual(self.gint_cgto[gint_row,gint_col], \
                                       ro_ints[i,j],self.MIN_DIGITS)
        return

    def testTrace(self):
        """Tests the computation of trace of integrals with atomic-orbital density matrix"""
        from Gen1Int.Tools import trace_reorder_sgto,trace_reorder_cgto,trace_gto
        # Tests the trace of real solid-harmonic Gaussians
        tr_ints = trace_reorder_sgto(self.ang_bra,self.mag_bra,
                                     self.ang_ket,self.mag_ket,
                                     self.gint_sgto,self.dens_sgto)
        self.assertAlmostEqual(self.ref_tr_sgto,tr_ints,self.MIN_DIGITS)
        # Tests the trace of Cartesian Gaussians
        tr_ints = trace_reorder_cgto(self.ang_bra,self.power_bra,
                                     self.ang_ket,self.power_ket,
                                     self.gint_cgto,self.dens_cgto)
        self.assertAlmostEqual(self.ref_tr_cgto,tr_ints,self.MIN_DIGITS)
        # Tests the trace of real solid-harmonic Gaussians without reordering
        tr_ints = trace_gto(self.gint_sgto,self.dens_sgto)
        self.assertAlmostEqual(self.ref_tr_sgto,tr_ints,self.MIN_DIGITS)
        return


def test():
    """Driver of unit tests."""
    basicConfig(level=Logger_DEBUG,format="%(message)s")
    suite = TestLoader().loadTestsFromTestCase(UnitTests)
    TextTestRunner(verbosity=2).run(suite)
    #TestMain()
    return

if __name__ == '__main__':
    test()
