from osgeo import gdal
from osgeo import osr, ogr
import numpy as np
import re
from scipy.linalg import lstsq
import scipy
import xlrd
import netCDF4
import cv2
import os
from geopandas import *
import rasterio as rio
from geopandas import GeoSeries
import rasterio.mask

class AHI:

    nl = 6001
    ns = 6001

    def writeTiff(self, im_data, im_width, im_height, im_bands, im_geotrans, im_proj, path):
        if 'int8' in im_data.dtype.name:
            datatype = gdal.GDT_Byte
        elif 'int16' in im_data.dtype.name:
            datatype = gdal.GDT_UInt16
        else:
            datatype = gdal.GDT_Float32

        if len(im_data.shape) == 3:
            im_bands, im_height, im_width = im_data.shape
        elif len(im_data.shape) == 2:
            im_data = np.array([im_data])
        else:
            im_bands, (im_height, im_width) = 1, im_data.shape
            # 创建文件
        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(path, im_width, im_height, im_bands, datatype)
        if (dataset != None and im_geotrans != '' and im_proj != ''):
            dataset.SetGeoTransform(im_geotrans)  # 写入仿射变换参数
            dataset.SetProjection(im_proj)  # 写入投影
        for i in range(im_bands):
            # dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
            dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
        del dataset

    def getDatafromNc(self, fileName, objectName):
        dataset = netCDF4.Dataset(fileName)
        predata = np.asarray(dataset.variables[objectName][:])
        return predata

    def getTiffData(self, fileName):
        dataset = gdal.Open(fileName)
        if dataset == None:
            print(fileName + "文件无法打开")
            return
        im_width = dataset.RasterXSize  # 栅格矩阵的列数
        im_height = dataset.RasterYSize  # 栅格矩阵的行数
        im_bands = dataset.RasterCount  # 波段数
        im_data = dataset.ReadAsArray(0, 0, im_width, im_height)  # 获取数据
        im_geotrans = dataset.GetGeoTransform()  # 获取仿射矩阵信息
        im_proj = dataset.GetProjection()  # 获取投影信息
        return im_data, im_width, im_height, im_bands, im_geotrans, im_proj

    def clipfile(self,in_file,features,out_file):

        src = rio.open(in_file)
        out_image, out_transform = rio.mask.mask(src, features, crop=True, nodata=src.nodata)
        out_image = out_image[0, :, :]
        out_meta = src.meta.copy()
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[0],
                         "width": out_image.shape[1],
                         "transform": out_transform})
        with rasterio.open(out_file, 'w', **out_meta) as dst:
            dst.write(out_image, indexes=1)


