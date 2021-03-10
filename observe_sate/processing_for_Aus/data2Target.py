from processing_for_Aus.SLSTR import *
from processing_for_Aus.myfun import *

s = SLSTR()

source = r'E:\class\class_2016_Aus_1000.tif'
target = r'E:\Sample\AHI8_2000.tif'
outfile = r'E:\class\class_2016_Aus_2000.tif'

##### Accordingly, the limited data range should also be set !!!


dataset = gdal.Open(target)
if dataset == None:
    print(target + "Can not be opened")
ns = dataset.RasterXSize  # 栅格矩阵的列数
nl = dataset.RasterYSize  # 栅格矩阵的行数
im_bands = dataset.RasterCount  # 波段数
data = dataset.ReadAsArray(0, 0, ns, nl)  # 获取数据
geog = dataset.GetGeoTransform()  # 获取仿射矩阵信息
proj = dataset.GetProjection()  # 获取投影信息
nsi = np.zeros([nl, ns])
nli = np.zeros([nl, ns])
for k in range(nl):
    nli[k, :] = k
    nsi[k, :] = np.linspace(0, ns - 1, ns)
nli = np.reshape(nli,-1)
nsi = np.reshape(nsi,-1)
lon,lat= s.imagexy2geo(dataset,nli,nsi)


data = np.zeros([nl,ns])


dataset = gdal.Open(source)
if dataset == None:
    print(source + "Can not be opened")
im_width = dataset.RasterXSize  # 栅格矩阵的列数
im_height = dataset.RasterYSize  # 栅格矩阵的行数
im_bands = dataset.RasterCount  # 波段数
data_source = dataset.ReadAsArray(0, 0, im_width, im_height)  # 获取数据
im_geotrans = dataset.GetGeoTransform()  # 获取仿射矩阵信息
im_proj = dataset.GetProjection()  # 获取投影信息
temp1,temp2=s.lonlat2geo(dataset,lat,lon)
imagey,imagex=s.geo2imagexy(dataset,temp1,temp2)
imagex = np.asarray(imagex,np.int)
imagey = np.asarray(imagey,np.int)
ind = (imagex>0)* (imagex<im_height-1) * (imagey>0) * (imagey<im_width-1)
imagex[imagex<0] = 0
imagex[imagex>im_height-1] = im_height-1
imagey[imagey<0] = 0
imagey[imagey>im_width-1] = im_width-1


temp = data_source[imagex,imagey]
indnew = (temp<20)*(temp>0)*ind
data = np.reshape(data, [-1])
data[indnew] = temp[indnew]
data = np.reshape(data, [nl, ns])
data[data < 0] = 0
data[data > 20] = 0

write_gdalTiff(data,ns,nl,1,geog,proj,outfile)
