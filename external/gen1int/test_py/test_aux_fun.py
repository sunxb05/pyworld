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
#  This file contains the unit tests of auxiliary functions in Gen1Int.
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
logger = getLogger("test_aux_fun")

class UnitTests(TestCase):
    """Unit tests"""
    def setUp(self):
        """Initializes unit tests"""
        # Sets the minimum number of correct digits
        self.MIN_DIGITS = 9
        # Sets the threshold of ratio to the referenced result
        self.RATIO_THRSH = [1.0-pow(10.0,-self.MIN_DIGITS),1.0+pow(10.0,-self.MIN_DIGITS)]
        return

    def testBoys(self):
        """Check the accuracy of Boys function subroutine"""
        from Gen1Int.ContrInt import aux_boys_vec
        from scipy.special import gamma,gammainc
        from math import pow
        min_order_boys = 0   #minimum order of Boys functions
        max_order_boys = 30  #maximum order of Boys functions
        arg_step = 0.05      #step of argument of Boys functions
        num_steps = 400      #number of steps
        print '\n F_n(T)  n   T      Gen1Int       Reference'
        # Argument from \var(arg_step) to \var(num_steps)*\var(arg_step)
        for iarg in xrange(1,num_steps+1):
            arg_boys = arg_step*iarg
            # Gen1Int subroutine
            val_boys = aux_boys_vec(min_order_boys,max_order_boys,arg_boys)
            for iord in xrange(min_order_boys,max_order_boys+1):
                # Computes Boys function using incomplete Gamma function
                jord = iord+0.5
                ref_boys = 0.5*gamma(jord)*gammainc(jord,arg_boys)/pow(arg_boys,jord)
                ratio_to_ref = val_boys[iord]/ref_boys
                #-if ratio_to_ref<self.RATIO_THRSH[0] or ratio_to_ref>self.RATIO_THRSH[1]:
                #-    print '!!NOTE!!', iord, arg_boys, val_boys[iord], ref_boys
                #-    break  #other higher orders are inaccurate either
                self.assertAlmostEqual(ref_boys,val_boys[iord],self.MIN_DIGITS)
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
