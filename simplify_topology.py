__author__ = 'asimmons'


import fiona
from shapely.geometry import shape, mapping, LineString, Polygon, MultiLineString, MultiPolygon
from shapely.geometry.polygon import LinearRing
import sys
import csv
from geomsimplify import *

#Set up simplify_topology.py to be a dynamic simplification script.

#Expectations are:

#1) a user should be able to preserve boundaries between features that have a different 'iso3' code.

#2) a 'two-pass' simplification effort which first simplifies objects (line or polygon) at a different re

# solution scales,
#while preserving shared border elements, then simplifies AGAIN the shared border elements.

#Ex: 'E:\It_3\simplify\test_poly\polygon_test_iso3.shp','E:\It_3\simplify\test_poly\polygon_test_iso3_simplified.shp', CAN 1000000, USA 50000, True)

# poly for test

# poly = Polygon([(0, 0), (0, 1), (1, 1), (0, 0)])

##############################################################################################################################
#
# Note: to change the quantitization from the default, you will have to do it before running process file
#
# Linestring test:
# g = PreserveTopology();g.process_file(r'I:\It_23\simplify\test_lines_simplify_topology\sample_lines.shp',r'I:\It_23\simplify\test_lines_simplify_topology\simplified_sample.shp', 1000000, True)
#
# Multilinestring test:
#  g = PreserveTopology();g.process_file(r'I:\It_23\simplify\test_multi\four_lines_dis.shp', r'I:\It_23\simplify\test_multi\four_lines_test_simp.shp', 1000000, True)
#
###################################################################################################################################

class PreserveTopology(object):

    def process_file(self, inFile, outFile, threshold, Topology = False, DynamicThresholdFile = None):
        """
        Takes an 'inFile' of an ESRI shapefile, converts it into a Shapely geometry - simplifies.
        Returns an 'outFile' of a simplified ESRI shapefile.

        IF Topology = True
        The object to be simplified is cut into junctions.

        IF Topology = False
        The object is simplified as is.

        Note:
        - A point is considered a junction if it shares the same point
        with another shape AND has different neighbors.

        - Identical features on top of one another DO NOT have different
        neighbors, and therefor DO NOT have any junctions.

        """


        # Open input file
        # loop over each
        with fiona.open(inFile, 'r') as input:

            meta = input.meta

            if Topology:
                # declare dictJunctions as a global variable
                # key = quantitized junction points, value = 1
                dictJunctions = {}

                #create instace of Junction
                junction = Junction()

                # create dictionary of all junctions in all shapes
                junction.find_all_junctions(inFile, dictJunctions)

                dictArcThresholds = None
                if DynamicThresholdFile:
                    # create a dictionary of all arcs between junctions, keyed by the junction pair
                    # with value set to the average threshold of adjacent polygons
                    with open(DynamicThresholdFile, 'r') as iso_thresholds_file:
                        csvreader = csv.reader(iso_thresholds_file)
                        dictIsoThresholds = {}
                        for line in iso_thresholds_file:
                            if validate and len(line) != 2:
                                raise ValueError('Unexpected number of columns in iso_thresholds_file for line: ' + repr(line))
                            iso3 = line[0]
                            threshold = line[1]
                            dictIsoThresholds[iso3] = threshold

                        dictArcThresholds = junction.find_all_arc_thresholds(inFile, dictJunctions, dictIsoThresholds)

                simplify = GeomSimplify(dictJunctions, dictArcThresholds) # if you need topology
            else:
                simplify = GeomSimplify()


            # create an outFile has the same crs, schema as inFile
            with fiona.open(outFile, 'w', **meta) as output:
            # Read shapely geometries from file
            # Loop through all shapely objects
                for myGeom in input:

                    myShape = shape(myGeom['geometry'])
                    simplifiedShapes =[]
                    if isinstance(myShape, LineString):
                        line = myShape
                        simplifiedShapes = [simplify.simplify_line_topology(line, threshold)]

                    elif isinstance(myShape, MultiLineString):
                        mline = myShape
                        simplifiedShapes = [simplify.simplify_multiline_topology(mline, threshold)]

                    elif isinstance(myShape, Polygon):
                        polygon = myShape
                        simplifiedShapes = [simplify.simplify_polygon_topology(polygon, threshold)]

                    elif isinstance(myShape, MultiPolygon):
                        mpolygon = myShape
                        simplifiedShapes = [simplify.simplify_multipolygon_topology(mpolygon, threshold)]

                    else:
                        raise ValueError('Unhandled geometry type: ' + repr(myShape.type))


                    # write to outfile
                    for simpleShape in simplifiedShapes:
                        if simpleShape is not None:
                            output.write({'geometry':mapping(simpleShape), 'properties': myGeom['properties']})



def str2bool(v):
    """
    Converts strings (which all command line passed argument are) to booleans.

    """
    return v.lower() in ("yes", "true", "t", "1")


def main():
    print "number of arguments (incl. py file name): " + str(len(sys.argv))
    if len(sys.argv) != 5:
        print "Wrong amount of arguments!"
        usage()
        exit()

    geomSimplifyObject = PreserveTopology()

    inputFile = sys.argv[1]
    outputFile = sys.argv[2]
    threshold = sys.argv[3]
    topology = sys.argv[4]
    topology_convert = str2bool(topology)

    if topology_convert is False:
        geomSimplifyObject.process_file(inputFile, outputFile, float(threshold), topology_convert)
        print "Finished simplifying file (with topology NOT preserved)!"
    elif topology_convert is True:
        geomSimplifyObject.process_file(inputFile, outputFile, float(threshold), topology_convert)
        print "Finished simplifying file (topology was preserved)!"

def usage():
    print "python simplify_topology.py <input file path> <output file path> <threshold> Topology= <True> or <False>"


if __name__ == "__main__":
    main()
