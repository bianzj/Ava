from util.myfun import *
from geopandas import *
import rasterio as rio
from geopandas import GeoSeries
import rasterio.mask
from rasterio.warp import (reproject,Resampling, transform_bounds,calculate_default_transform as calcdt)



def read_txt2pos(infile,numskip):
    pos_ = []
    fileName_ = []
    with open(infile) as f:
        for kskip in range(numskip): lines = f.readline()
        lines = f.readline()
        while lines:
            fields = lines.split()
            fileName_.append(fields[0])
            pos_.append(fields[1:4])
            lines = f.readline()
    pos_ = np.asarray(pos_,dtype=float)
    return fileName_,pos_


sample = 'white'
band = 'nir'

wdir = r'F:\UAVpro\20200604huailai\1201\VNIR_25_result\\'
wdirout = r'F:\UAVpro\20200604huailai\1201\\'+sample+'_'+band+'_result/'

infile_loc =  r'F:\UAVpro\20200604huailai\1201\\'+sample+'_'+band+'.tif'
dataset0 = gdal.Open(infile_loc)
if dataset0 == None:
    print(infile_loc + "文件无法打开")
ns = dataset0.RasterXSize  # 栅格矩阵的列数
nl = dataset0.RasterYSize  # 栅格矩阵的行数
im_bands = dataset0.RasterCount  # 波段数
data = dataset0.ReadAsArray(0, 0, ns, nl)  # 获取数据
geog = dataset0.GetGeoTransform()  # 获取仿射矩阵信息
proj = dataset0.GetProjection()  # 获取投影信息
nsi = np.zeros([nl, ns])
nli = np.zeros([nl, ns])
for k in range(nl):
    nli[k, :] = k
    nsi[k, :] = np.linspace(0, ns - 1, ns)
nli = np.reshape(nli,-1)
nsi = np.reshape(nsi,-1)
temp1,temp2 = imagexy2geo(dataset0,nli,nsi)
lon,lat = geo2lonlat(dataset0,temp1,temp2)
lon_o = lon
lat_o = lat

fileNames = search_file(wdir,'tif')
number = len(fileNames)

minNumber = 500
maxNumber = 2000


number_threshold = ns*nl*0.75
for kimage in range(number):

    filename = fileNames[kimage][:-4]
    number = filename[18:]
    if np.int(number) < minNumber:continue
    if np.int(number) > maxNumber:continue

    infile = wdir + filename+'.tif'
    outfile = wdirout + filename+'.tif'
    print(infile)

    dataset1 = gdal.Open(infile)
    if dataset1 == None:
        print(infile + "文件无法打开")
        continue
    im_width = dataset1.RasterXSize
    im_height = dataset1.RasterYSize
    im_bands = dataset1.RasterCount
    data = dataset1.ReadAsArray(0, 0, im_width, im_height)
    im_geotrans = dataset1.GetGeoTransform()
    im_proj = dataset1.GetProjection()

    temp1, temp2 = lonlat2geo(dataset1, lat_o, lon_o)
    imagex, imagey = geo2imagexy(dataset1, temp1, temp2)
    imagex = np.asarray(imagex, np.int)
    imagey = np.asarray(imagey, np.int)
    # print(np.max(imagex),np.max(imagey))

    ind = (imagex > 0) * (imagex < im_width - 1) * (imagey > 0) * (imagey < im_height - 1)
    imagex[imagex < 0] = 0
    imagex[imagex > im_width - 1] = im_width - 1
    imagey[imagey < 0] = 0
    imagey[imagey > im_height - 1] = im_height - 1
    temp = data[imagey,imagex]

    indnew = (temp > 1) * ind
    if np.sum(indnew) <= number_threshold: continue

    data_new = np.zeros([nl, ns])
    data_new = np.reshape(data_new,[-1])
    data_new[indnew] = temp[indnew]
    data_new = np.reshape(data_new,[nl,ns])
    data_new[data_new > 10000] = 0

    write_image_gdal(data_new, ns, nl, 1, geog, proj, outfile)


