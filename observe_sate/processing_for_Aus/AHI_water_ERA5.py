from myfun_image import *
from myfun_file import *
from myfun_sci import *
from pyhdf.SD import SD, SDC
from MODIS_constant import *
import time
import os
from geopandas import *
import netCDF4


targetArea = r'G:\base\type.tif'
dataset = gdal.Open(targetArea)
if dataset == None:
    print(targetArea + "文件无法打开")
ns = dataset.RasterXSize
nl = dataset.RasterYSize
im_bands = dataset.RasterCount
data = dataset.ReadAsArray(0, 0, ns, nl)
geog = dataset.GetGeoTransform()
proj = dataset.GetProjection()
nsi = np.zeros([nl, ns])
nli = np.zeros([nl, ns])
for k in range(nl):
    nli[k, :] = k
    nsi[k, :] = np.linspace(0, ns - 1, ns)
nli = np.reshape(nli,-1)
nsi = np.reshape(nsi,-1)
temp1, temp2 = imagexy2geo(dataset, nli, nsi)
lon, lat = geo2lonlat(dataset, temp1, temp2)


infile = r'e:\ahi8_sample\\world_equal_latlon.tif'
dataset = gdal.Open(infile)
if dataset == None:
    print(infile + "文件无法打开")
im_width = dataset.RasterXSize
im_height = dataset.RasterYSize
im_bands = dataset.RasterCount
im_geotrans = dataset.GetGeoTransform()
im_proj = dataset.GetProjection()
temp1,temp2=s.lonlat2geo(dataset,lat,lon)
imagey,imagex=s.geo2imagexy(dataset,temp1,temp2)
imagex = np.asarray(imagex,np.int)
imagey = np.asarray(imagey,np.int)
ind = (imagex>0)* (imagex<im_height-1) * (imagey>0) * (imagey<im_width-1)
imagex[imagex<0] = 0
imagex[imagex>im_height-1] = im_height-1
imagey[imagey<0] = 0
imagey[imagey>im_width-1] = im_width-1



indir = R'G:\ERA5\\'
outdir = R'E:\AHI8_day_01\\'

object_name = 'tcwv'
for kday in range(16,17):
    infile = indir + 'ERA5.single-level.201801%02d.nc'%kday
    dataset = netCDF4.Dataset(infile)
    data = dataset.variables[object_name][:]
    for khour in range(0,24):
        # if (khour < 21) and (khour > 9):continue
        # if (khour > 9): continue
        if (khour < 21): continue
        outfile = outdir + 'AHI_%02d'%kday+'_%02d'%khour+'00_tcwv.tif'
        print(outfile)
        tcwv = data[khour,imagex,imagey]
        tcwv = np.reshape(tcwv,[nl,ns])
        tcwv = tcwv/10.0
        write_gdalTiff(tcwv,ns,nl,1,geog,proj,outfile)



