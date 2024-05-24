from shapely.geometry import mapping, Polygon
import fiona

# Here's an example Shapely geometry
poly = Polygon([(0, 0), (0, 2), (1, 1), (2, 2), (2, 0), (1, 1), (0, 0)])

# Define a polygon feature geometry with one attribute
schema = {
    "geometry": "Polygon",
    "properties": {"id": "int"},
}

# Write a new Shapefile
with fiona.open(
    "E:\Tableau\2017\It_3\simplify\test_poly\bowtie.shp", "w", "ESRI Shapefile", schema
) as c:
    ## If there are multiple geometries, put the "for" loop here
    c.write(
        {
            "geometry": mapping(poly),
            "properties": {"id": 123},
        }
    )
