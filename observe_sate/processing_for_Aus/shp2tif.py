from geopandas import *
import rasterio as rio
import rasterio.mask

shpdatafile=r'J:\Country_shp\country.shp'
rasterfile=r'J:K:\Data\Sate_ASTERGED_global_lse\MODIS\ASTER.GED.MODIS.LSE31.tif'
outfile=r'K:\Data\Sate_ASTERGED_global_lse\MODIS\LSE31_aus.tif'


src = rio.open(rasterfile)
dst_crs = src.crs


shpdata = GeoDataFrame.from_file(shpdatafile)
shpdata_crs = shpdata.to_crs(src.crs)
# geo = shpdata_crs.geometry[12]
# features = [shpdata.geometry.__geo_interface__]
# feature = [geo.__geo_interface__]
features = [shpdata_crs.geometry[12].__geo_interface__]
#

out_image, out_transform = rio.mask.mask(src, features, crop=True, nodata=src.nodata)
out_image = out_image[0,:,:]
out_meta = src.meta.copy()
out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[0],
                 "width": out_image.shape[1],
                 "transform": out_transform})
with rasterio.open(outfile, 'w', **out_meta) as dst:
    dst.write(out_image, indexes=1)


#
# affine, width, height = calcdt(src.crs, dst_crs, src.width, src.height, *src.bounds)

# kwargs = src.meta.copy()
# kwargs.update({
#     'crs': dst_crs,
#     'transform': affine,
#     'affine': affine,
#     'width': width,
#     'height': height,
#     'geotransform':(0,1,0,0,0,-1) ,
#     'driver': 'GTiff'
# })

# dst = rio.open(outfile, 'w', **kwargs)
# for i in range(1, src.count + 1):
#     reproject(
#         source = rio.band(src, i),
#         destination = rio.band(dst, i),
#         src_transform = affine,
#         src_crs = src.crs,
#         dst_transform = affine,
#         dst_crs = dst_crs,
#         dst_nodata = src.nodata,
#         resampling = Resampling.bilinear)


