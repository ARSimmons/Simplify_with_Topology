__author__ = 'asimmons'

from geomsimplify import *
from nose.tools import *
import unittest

## polygons = [Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)]), Polygon([(1,0),(1,1),(2,1),(2,0),(1,0)]), Polygon([(2,1),(2,2),(3,2),(3,1),(2,1)])]
##  polygon 0 and 1 shares two points (an entire border). Polygon 1 and 2 only shares one point.

class test_GeomSimplify(unittest.TestCase):

    def test_create_ring_from_arcs(self):
        todo

    def test_simplify_polygon(self):
        big one

    def test_simplify_multipolygon(self):
        todo

    def test_find_junctions_polygon(self):
        assert_false()

    def test_cut_ring_by_junctions(self):
        assert_false()

    def test_cut_polygon_by_junctions(self):
        assert_false()

    def test_cut_mpolygon_by_junctions(self):
        assert_false()