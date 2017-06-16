
<html>
<head>
</head>
<body>
# Table of Contents
[Team Members](#team-members)

[Project Summary](#project-summary)

# <a name="team-members"></a>Team Members
* "Arielle Simmons" <ari.ucb.fire@gmail.com>
	
# <a name="project-summary"></a>Project Summary

Fiona (1.1+) and Shapely libraries used to simplify lines, multilines, polygons, and multipolygon shapefiles by an area threshold. Also provides an option to
preserve topology. Topology can be preserved when simplifying by one threshold  value, or by many (i.e. 'dynamic simplification'). In the case of using
'dynamic simplification' the user must provide a csv (no header), containing the 

The simplification algorithm is based off of M. Visvalingam and J.D. Whyatt's algorithm (1993). More details about the 
Visvalingam-Whyatt algorithm can be found here: https://hydra.hull.ac.uk/resources/hull:8338      .

Threshold is the area of the largest allowed triangle.

Key points to note:

	- As of 6/16/2017 there are NO TOPOLOGY preserving rules in place for interior rings (ToDO)
	- In a topology preserving process, Polygons which are smaller then the area threshold will be deleted
      UNLESS they are adjacent to another border (in which case they will delete down to a minimum of 3 vertices).	
	- Lines preserve their beginning and end point, thus lines CANNOT BE DELETED (regardless of the topology setting). 
	  The beginning and end points of a line feature are static throughout the 
	  simplification process.
	- threshold units are determined by shapefile map units.  
	- to run from command line: 
	
	   python  <input file path> <output file path> <IF you want to preserve Topology (optional) = --topology> <threshold> OR <DynamicThresholdFile=(optional) dynamic threshold csv file path>
	
![Screenshot](https://raw.github.com/ARSimmons/Simplify_with_Topology/master/dynamic_simplification.JPG)


*THIS PROJECT IS STILL IN PROCESS 5/22/2013*
 
</body>
</html>