import shapefile


# >>> r = shapefile.Reader('shapefiles/test/polygon')
#
# >>> w = shapefile.Writer('shapefiles/test/copy')
# >>> w.fields = r.fields[1:] # skip first deletion field
# >>> # adding existing Shape objects
# >>> for shaperec in r.iterShapeRecords():
# ...     w.record(*shaperec.record)
# ...     w.shape(shaperec.shape)
# >>> # or GeoJSON dicts
# >>> for shaperec in r.iterShapeRecords():
# ...     w.record(*shaperec.record)
# ...     w.shape(shaperec.shape.__geo_interface__)
# >>> w.close()


#### cut Aus from world
# infile = r'J:\Country_shp\country.shp'
# outfile = r'J:\Country_shp\australia.shp'
#
# r = shapefile.Reader(infile)
# w = shapefile.Writer(outfile)
# w.fields = r.fields
#
# the = 12
# targetshapeRecord = r.shapeRecord(the)
# w.record(*targetshapeRecord.record)
# w.shape(targetshapeRecord.shape)
# w.record(*targetshapeRecord.record)
# w.shape(targetshapeRecord.shape.__geo_interface__)
# # w.record('crs','epsg:4326')
# r.close()
# w.close()

#### cut Aus from world
infile = r'F:\sample\tz_world_mp.shp'
outfile = r'F:\sample\simple.shp'

r = shapefile.Reader(infile)
w = shapefile.Writer(outfile)
w.fields = r.fields

the = 41
targetshapeRecord = r.shapeRecord(the)
w.record(*targetshapeRecord.record)
w.shape(targetshapeRecord.shape)

# w.record(*targetshapeRecord.record)
# w.shape(targetshapeRecord.shape.__geo_interface__)
# w.record('crs','epsg:4326')
r.close()
w.close()