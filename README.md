# Table of Contents

[Team Members](# Team Members)

[Project Summary](#Project Summary)

## Team Members

* "Arielle Simmons" <ari.ucb.fire@gmail.com>

## Project Summary

Fiona (1.1+) and Shapely libraries simplify lines, multi-lines, polygons, and multipolygon shapefiles by an area threshold. Also provides an option to preserve topology.

Topology can be maintained when simplifying by one threshold value, or by many (i.e. 'dynamic simplification'). In the case of using 'dynamic simplification,' the user must provide a CSV (no header), containing the desired threshold values per entity (right now hardcoded to be iso3 -country- entities)

Threshold is the area of the largest allowed triangle.

The simplification algorithm is based off of M. Visvalingam and J.D. Whyatt's algorithm (1993).

More details about the Visvalingam-Whyatt algorithm can be found [here.](https://hydra.hull.ac.uk/resources/hull:8338)

### **Key points to note:**

* As of 6/16/2017 there are NO TOPOLOGY preserving rules in place for interior rings (ToDO)
* In a topology preserving process, Polygons which are smaller then the area threshold will be deleted
* UNLESS they are adjacent to another border (in which case they will delete down to a minimum of 3 vertices).
* Lines preserve their beginning and end point, thus lines CANNOT BE DELETED (regardless of the topology setting).
* Threshold units are determined by shapefile map units.
* To run from command line:

> python simplify_topology.py `<input file path>` `<output file path>` <Preserve Topology (optional) = --topology> `<threshold>` OR <DynamicThresholdFile=(optional) dynamic threshold csv file path>

**Example usage:**

> python simplify_topology.py -i input/input.shp -o output/output.shp -t 0.0001
>
> python simplify_topology.py -i input/input.shp -o output/output.shp -t 0.0001 -j
>
> python simplify_topology.py -i input/input.shp -o output/output.shp -d dynamic_thresholds.csv

![Screenshot](https://raw.github.com/ARSimmons/Simplify_with_Topology/master/dynamic_simplification.JPG)

> *THIS PROJECT IS STILL IN PROCESS 5/22/2013*
