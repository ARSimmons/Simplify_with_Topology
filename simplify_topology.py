__author__ = "asimmons"

import csv
import fiona
from geomsimplify import GeomSimplify
from optparse import OptionParser
from shapely.geometry import (
    shape,
    mapping,
    LineString,
    Polygon,
    MultiLineString,
    MultiPolygon,
)

################################################################################################################################################
# 1) This script is to created to simplify lines, multilines, polygons or multipolygon shapefiles.                                             #
#                                                                                                                                              #
# 2) The simplification algorithm is the Visvalingam algorithm found here  :                                                                   #
#     https://hydra.hull.ac.uk/resources/hull:8338                                                                                             #
#                                                                                                                                              #
# 3) Threshold is the area of the largest allowed triangle.                                                                                    #
#                                                                                                                                              #
# 4) Users have the ability to process different threshold levels for different countries by using a .csv file designating differing           #
#    thresholds by 'iso3' code. If you do this it is expected that the shapefile has a populated attribute field called 'iso3' .               #
#                                                                                                                                              #
# Note: to change the quantitization from the default, you will have to do it before running 'process_file'                                    #
#                                                                                                                                              #                                                                                                                                             #
################################################################################################################################################

debug = True

self_intersections_fixed = 0
validate = True


class SimplifyProcess:
    def process_file(
        self, inFile, outFile, threshold, Topology=False, DynamicThresholdFile=None
    ):
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

        # Convert threshold from str to float
        threshold = float(threshold)

        # Open input file
        # loop over each
        with fiona.open(inFile, "r") as input:

            meta = input.meta

            if Topology:
                # declare dictJunctions as a global variable
                # key = quantitized junction points, value = 1
                dictJunctions = {}

                # create instace of Junction
                simplifyObj = GeomSimplify()

                # create dictionary of all junctions in all shapes
                simplifyObj.find_all_junctions(inFile, dictJunctions)

                dictArcThresholds = None
                if DynamicThresholdFile:
                    # create a dictionary of all arcs between junctions, keyed by the junction pair
                    # with value set to the average threshold of adjacent polygons
                    with open(DynamicThresholdFile, "rb") as iso_thresholds_file:
                        csvreader = csv.reader(iso_thresholds_file)
                        dictIsoThresholds = {}
                        for line in iso_thresholds_file:
                            line_array = (
                                line.replace("\r", "").replace("\n", "").split(",")
                            )
                            if validate and len(line_array) != 2:
                                raise ValueError(
                                    "Unexpected number of columns in iso_thresholds_file for line: "
                                    + repr(line_array)
                                )
                            iso3 = line_array[0]
                            threshold = float(line_array[1])
                            dictIsoThresholds[iso3] = threshold

                        dictArcThresholds = simplifyObj.find_all_arc_thresholds(
                            inFile, dictJunctions, dictIsoThresholds
                        )

                simplify = GeomSimplify(
                    dictJunctions, dictArcThresholds
                )  # if you need topology
            else:
                simplify = GeomSimplify()

            invalid_geoms_count = 0

            # create an outFile has the same crs, schema as inFile
            with fiona.open(outFile, "w", **meta) as output:
                # Read shapely geometries from file
                # Loop through all shapely objects
                for myGeom in input:

                    myShape = shape(myGeom["geometry"])
                    simplifiedShapes = []
                    if isinstance(myShape, LineString):
                        line = myShape
                        if Topology:
                            simplifiedShapes = [
                                simplify.simplify_line_topology(line, threshold)
                            ]
                        else:
                            simplifiedShapes = [simplify.simplify_line(line, threshold)]

                    elif isinstance(myShape, MultiLineString):
                        mline = myShape
                        if Topology:
                            simplifiedShapes = [
                                simplify.simplify_multiline_topology(mline, threshold)
                            ]
                        else:
                            simplifiedShapes = [
                                simplify.simplify_multiline(mline, threshold)
                            ]

                    elif isinstance(myShape, Polygon):
                        polygon = myShape
                        if Topology:
                            simplifiedShapes = [
                                simplify.simplify_polygon_topology(polygon, threshold)
                            ]
                        else:
                            simplifiedShapes = [
                                simplify.simplify_polygon(polygon, threshold)
                            ]

                    elif isinstance(myShape, MultiPolygon):
                        mpolygon = myShape
                        if Topology:
                            simplifiedShapes = [
                                simplify.simplify_multipolygon_topology(
                                    mpolygon, threshold
                                )
                            ]
                        else:
                            simplifiedShapes = [
                                simplify.simplify_multipolygon(mpolygon, threshold)
                            ]

                    else:
                        raise ValueError(
                            "Unhandled geometry type: " + repr(myShape.type)
                        )

                    # Check for invalid geometries in shape list
                    check_invalid_geometry(simplifiedShapes)

                    # write to outfile
                    for simpleShape in simplifiedShapes:
                        if simpleShape is not None:
                            output.write(
                                {
                                    "geometry": mapping(simpleShape),
                                    "properties": myGeom["properties"],
                                }
                            )

        print(
            "Self-intersecting rings found and fixed: " + str(self_intersections_fixed)
        )

        #
        if debug and Topology:
            with open("debug_junctions.txt", "w") as output:
                for key in dictJunctions:
                    output.write(str(key))


def str2bool(v):
    """
    Converts strings (which all command line passed argument are) to booleans.

    """
    return v.lower() in ("yes", "true", "t", "1")


def main():
    parser = OptionParser()
    parser.add_option(
        "-i",
        "--input_file",
        dest="inputFile",
        help="Input file with geometries to be simplified",
        metavar="FILE",
    )
    parser.add_option(
        "-o",
        "--output_file",
        dest="outputFile",
        help="Output file with simplified geometries",
        metavar="FILE",
    )
    parser.add_option(
        "-t",
        "--threshold",
        dest="threshold",
        default=0,
        help="Threshold for simplification. Exclusive with -d",
    )
    parser.add_option(
        "-j",
        "--topology",
        action="store_true",
        dest="topology",
        default=False,
        help="Flag indicating that we should preserve topology of shapes",
    )
    parser.add_option(
        "-d",
        "--dynamic_thresholds",
        dest="dynamicThresholds",
        help="CSV file containing iso3 and corresponding threshold value. Exclusive with -t",
        metavar="FILE",
    )

    (options, args) = parser.parse_args()

    inputFile = options.inputFile
    if not inputFile:
        print("Must specify input file")
        usage()
        exit()

    outputFile = options.outputFile
    if not outputFile:
        print("Must specify output file")
        usage()
        exit()

    topology = options.topology

    threshold = options.threshold
    dynamic_thresholds = options.dynamicThresholds
    if (not threshold and not dynamic_thresholds) or (threshold and dynamic_thresholds):
        print(
            "Must set either threshold or dynamic_thresholds, but NOT both and NOT nothing."
        )
        usage()
        exit()

    # dynamic_thresholds = None
    # if len(sys.argv) == 6:
    #     dynamic_thresholds = sys.argv[5]

    geomSimplifyObject = SimplifyProcess()

    if topology is False:
        geomSimplifyObject.process_file(
            inputFile, outputFile, float(threshold), topology
        )
        print("Finished simplifying file (with topology NOT preserved)!")
    elif topology is True:
        geomSimplifyObject.process_file(
            inputFile, outputFile, float(threshold), topology, dynamic_thresholds
        )
        print("Finished simplifying file (topology was preserved)!")


def usage():
    print(
        "python simplify_topology.py <input file path> <output file path> <threshold> OR <DynamicThresholdFile=(optional) dynamic threshold csv file path>"
    )


def check_invalid_geometry(shapeList):
    # Check for self-intersecting rings, and fix
    checkedShapeList = []
    for shape in shapeList:
        if isinstance(shape, Polygon):
            shape = fix_self_intersections_polygon(shape)
        elif isinstance(shape, MultiPolygon):
            shape = fix_self_intersections_mpolygon(shape)

        # TODO
        # Check for interior rings that have points outside the exterior ring,

        # Check if more than one interior ring touches the exterior ring

        checkedShapeList.append(shape)

    return checkedShapeList


def fix_self_intersections_polygon(polygon):
    global self_intersections_fixed

    if validate:
        if not isinstance(polygon, Polygon):
            raise ValueError(
                "Non-Polygon passed to fix_self_intersections_polygon: "
                + repr(polygon.type)
            )

    shape = polygon
    if not polygon.is_valid:
        shape = polygon.buffer(0)
        self_intersections_fixed += 1

    return shape


def fix_self_intersections_mpolygon(mpolygon):
    if validate:
        if not isinstance(mpolygon, MultiPolygon):
            raise ValueError(
                "Non-MultiPolygon passed to fix_self_intersections_mpolygon: "
                + repr(mpolygon.type)
            )

    checked_polygons = []
    for polygon in mpolygon.geoms:
        if not polygon.is_valid:
            cleaned_shape = fix_self_intersections_polygon(polygon)
            for p in cleaned_shape.geoms:
                checked_polygons.append(p)
        else:
            checked_polygons.append(polygon)

    return MultiPolygon(checked_polygons)


if __name__ == "__main__":
    main()

# example usage:
# python simplify_topology.py -i input/input.shp -o output/output.shp -t 0.0001
# python simplify_topology.py -i input/input.shp -o output/output.shp -t 0.0001 -j
# python simplify_topology.py -i input/input.shp -o output/output.shp -d dynamic_thresholds.csv
