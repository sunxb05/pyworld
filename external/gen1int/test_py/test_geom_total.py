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
#  This file contains the unit tests of geometric derivatives.
#
#  2011-03-04, Bin Gao: first version

# Logging facility
from logging import getLogger
from logging import basicConfig
from logging import DEBUG as Logger_DEBUG
# Unit testing framework
from unittest import TestCase,TestLoader,TextTestRunner
#from unittest import main as TestMain

__author__ = "Bin Gao, and Andreas J. Thorvaldsen"
__copyright__ = "Copyright 2009, 2010"
__credits__ = ["Radovan Bast", "Kenneth Ruud"]
__license__ = "LGPLv3"
__version__ = "0.2.1"
__maintainer__ = "Bin Gao"
__email__ = "bin.gao@uit.no"
__status__ = "Development"

# Set up a specific logger with our desired output level
logger = getLogger("test_geom_deriv")

class UnitTests(TestCase):
    """Unit tests"""
    def setUp(self):
        """Initializes unit tests"""
        # Sets the number of atoms
        self.n_atom = 5
        # Sets the order of geometric derivatives
        self.n_order = 4
        # Sets the number of differentiated centers
        self.n_center = 5
        return

    def testTree(self):
        """Tests recursive subroutine to generate the centers of geometric derivatives"""
        from Gen1Int.ContrInt import geom_deriv_cent_tree,geom_deriv_cent_npath
        from numpy import zeros,ones
        # Initializes
        visit_depth = self.n_order                       #depth of atom to visit
        node_list = ones(self.n_order,dtype='int32')     #selected atom nodes
        weight_node = ones(self.n_order,dtype='int32')   #weights of the selected atom nodes
        atom_cent = zeros(self.n_center,dtype='int32')   #indices of atomic centers
        atom_cent[0] = 1
        order_cent = zeros(self.n_center,dtype='int32')  #order of derivatives of the corresponding centers
        order_cent[0] = self.n_order
        num_list = 1                                     #incremental recorder of generated atomic centers
        # Computes the number of different compositions of differentiated centers
        num_paths = geom_deriv_cent_npath(self.n_atom,self.n_order,self.n_center)
        print
        print ">> Number of different paths:", num_paths
        # Generate the centers of geometric derivatives in a recursive fasion
        print num_list, '>> ', node_list, '>> ', weight_node[self.n_order-1]
        #print num_list, '>> ', ( atom_cent(ia), '(', order_cent(ia), ')', ia = 1, weight_node(n_order) )
        geom_deriv_cent_tree(self.n_atom,self.n_order,self.n_center, \
                             visit_depth,node_list,weight_node,      \
                             atom_cent,order_cent)
        print num_list, '>> ', node_list, '>> ', weight_node[self.n_order-1]
        return


def test():
    """Driver of unit tests."""
    basicConfig(level=Logger_DEBUG,format="%(message)s")
    suite = TestLoader().loadTestsFromTestCase(UnitTests)
    TextTestRunner(verbosity=2).run(suite)
    #TestMain()
    return

if __name__ == '__main__': test()
