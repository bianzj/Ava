# _*_ coding: utf-8 _*_
__author__ = 'xbr'
__date__ = '2018/12/10 11:41'

#######################################
#### It failed becasue the black background can not be removed!!!
#######################################

import gdal
from processing_for_Aus.myfun import *

def get_extent(fn):
    ds = gdal.Open(fn)
    rows = ds.RasterYSize
    cols = ds.RasterXSize
    # 获取图像角点坐标
    gt = ds.GetGeoTransform()
    minx = gt[0]
    maxy = gt[3]
    maxx = gt[0] + gt[1] * rows
    miny = gt[3] + gt[5] * cols
    return (minx, maxy, maxx, miny)

wdir = r'J:\SateTDR\\'
fileNames = ope_searchfile_accept1_rej1(wdir,'S3A','ZIP')
outfile = r'J:\mosaic\Aus_20180101.tif'
band = r'\bt8_proj.tif'

in_files = fileNames[:2]

# 通过两两比较大小,将最终符合条件的四个角点坐标保存，
# 即为拼接图像的四个角点坐标
in_file = wdir + in_files[0]+band
minX, maxY, maxX, minY = get_extent(in_file)
for fn in in_files[1:]:
    in_file = wdir + fn+band
    minx, maxy, maxx, miny = get_extent(in_file)
    minX = min(minX, minx)
    maxY = max(maxY, maxy)
    maxX = max(maxX, maxx)
    minY = min(minY, miny)

# 获取输出图像的行列数
in_file =  wdir + in_files[0]+band
in_ds = gdal.Open(in_file)
gt = in_ds.GetGeoTransform()
cols = int((maxX - minX) / (gt[1]))+1000
rows = int((maxY - minY) / abs(gt[5]))+1000

# 创建输出图像
driver = gdal.GetDriverByName('gtiff')
out_ds = driver.Create(outfile,cols,rows,1)
out_ds.SetProjection(in_ds.GetProjection())
out_band = out_ds.GetRasterBand(1)

gt = list(in_ds.GetGeoTransform())
gt[0], gt[3] = minX, maxY
out_ds.SetGeoTransform(gt)

for fn in in_files:
    in_file = wdir + fn + band
    print(in_file)
    in_ds = gdal.Open(in_file)
    trans = gdal.Transformer(in_ds, out_ds, [])
    success, xyz = trans.TransformPoint(False, 0, 0)
    x, y, z = map(int, xyz)
    data = in_ds.GetRasterBand(1).ReadAsArray()
    # data = np.transpose(data)
    out_band.WriteArray(data,x,y)
out_band.FlushCache()
del in_ds, out_band, out_ds