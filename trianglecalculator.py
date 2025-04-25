#! /usr/bin/env python
# encoding: utf-8

__author__ = "asimmons"


class TriangleCalculator(object):
    """
    TriangleCalculator() - Calculates the area of a triangle using the cross-product.

    """

    def __init__(self, point, index):
        # Save instance variables
        self.point = point
        self.ringIndex = index
        self.prevTriangle = None
        self.nextTriangle = None

    # enables the instantiation of 'TriangleCalculator' to be compared
    # by the calcArea().
    def __cmp__(self, other):
        return cmp(self.calcArea(), other.calcArea())

    ## calculate the effective area of a triangle given
    ## its vertices -- using the cross product
    def calcArea(self):
        # Add validation
        if not self.prevTriangle or not self.nextTriangle:
            print("ERROR:")

        p1 = self.point
        p2 = self.prevTriangle.point
        p3 = self.nextTriangle.point
        area = (
            abs(
                p1[0] * (p2[1] - p3[1])
                + p2[0] * (p3[1] - p1[1])
                + p3[0] * (p1[1] - p2[1])
            )
            / 2.0
        )
        # print "area = " + str(area) + ", point = " + str(self.point)
        return area

    def __lt__(self, other):
        # For example, if you want to compare based on the area of the triangles:
        return self.calcArea() < other.calcArea()
