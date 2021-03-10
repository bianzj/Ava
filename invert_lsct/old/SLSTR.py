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


class SLSTR:
    # to resize from... to ...
    adjust_xmin_n = 0
    adjust_xmax_n = 95 + 1
    adjust_xmin_o = 36
    adjust_xmax_o = 94 + 1
    # to vary from oblique to nadir
    adjust_xmin_no = 548
    adjust_xmax_no = 548 + 900
    adjust_ymin_no = 2
    adjust_ymax_no = 1200

    aodmax = 1.0
    aodmin = 0.03
    # aodvalue = 0.05 ### low
    aodvalue = 0.10 ### normal
    # aodvalue = 0.15 ### high
    aodcheck = 1 ### aodcheck=1 : 0.01-1.0 ; aodcheck = 2: aodvalue

    ns_nadir = 1500
    nl_nadir = 1200
    ns_obliq = 900
    nl_obliq = 1200
    ndvi_min = 0.054
    ndvi_max = 0.947
    wl8 = 10.8
    wl9 = 12.0

    def ifallfilesexists(self, fileDir):

        infile = fileDir + r'\S2_radiance_an.nc'
        c1 = os.path.exists(infile)
        infile = fileDir + r'\geodetic_in.nc'
        c2 = os.path.exists(infile)
        infile = fileDir + r'\geometry_to.nc'
        c3 = os.path.exists(infile)
        infile = fileDir + r'\met_tx.nc'
        c4 = os.path.exists(infile)
        infile = fileDir + r'\S2_radiance_an.nc'
        c5 = os.path.exists(infile)
        infile = fileDir + r'\S6_radiance_an.nc'
        c5 = os.path.exists(infile)
        infile = fileDir + r'\S8_BT_in.nc'
        c6 = os.path.exists(infile)
        infile = fileDir + r'\S8_BT_io.nc'
        c7 = os.path.exists(infile)
        return c1 * c2 * c3 * c4 * c5 * c6 * c7

    def getStrfromTxt(self, fileName):
        f = open(fileName, 'r')
        data = f.readlines()
        f.close()
        return data

    def getDatafromTxt(self, filename, num_pass, num_col):
        f = open(filename, 'r')
        temp = f.readlines()
        temp = np.asarray(temp)
        temp = temp[num_pass:]
        num = len(temp)
        lut = np.zeros([num, num_col])
        for k in range(num):
            if (temp[k] == ""): continue
            tempp = (re.split(r'\s+', temp[k].strip()))
            lut[k, :] = np.asarray(tempp)
        return lut

    def getDatafromNc(self, fileName, objectName):
        dataset = netCDF4.Dataset(fileName)
        predata = np.asarray(dataset.variables[objectName][:])
        return predata

    def getExceldata(self, filename):
        ExcelFile = xlrd.open_workbook(filename)
        ExcelFile.sheet_names()
        sheet = ExcelFile.sheet_by_name('Sheet1')
        return sheet

    def dataTransferlatlon(self, datafile, latlonfile, outfile):
        # infile_lat = r'F:\auxiliary\\\latlon.tif'
        dd, ns, nl, temp3, geog, proj = self.getTiffData(latlonfile)
        lon = dd[1, :, :]
        lat = dd[3, :, :]
        datanew = np.zeros([nl, ns])

        infile = datafile
        dataset = gdal.Open(infile)
        if dataset == None:
            print(infile + "文件无法打开")
        im_width = dataset.RasterXSize  # 栅格矩阵的列数
        im_height = dataset.RasterYSize  # 栅格矩阵的行数
        im_bands = dataset.RasterCount  # 波段数
        data = dataset.ReadAsArray(0, 0, im_width, im_height)  # 获取数据
        im_geotrans = dataset.GetGeoTransform()  # 获取仿射矩阵信息
        im_proj = dataset.GetProjection()  # 获取投影信息
        temp1, temp2 = self.lonlat2geo(dataset, lat, lon)
        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
        imagex = np.asarray(imagex, np.int)
        imagey = np.asarray(imagey, np.int)
        ind = (imagex > 0) * (imagex < im_height - 1) * (imagey > 0) * (imagey < im_width - 1)
        imagex[imagex < 0] = 0
        imagex[imagex > im_height - 1] = im_height - 1
        imagey[imagey < 0] = 0
        imagey[imagey > im_width - 1] = im_width - 1
        temp = data[imagex, imagey]
        data1 = np.reshape(datanew, [-1])
        data1[ind] = temp[ind]
        data1 = np.reshape(data1, [nl, ns])
        self.writeTiff(data1, ns, nl, 1, geog, proj, outfile)

    def resize(self, preArray, nl, ns):
        data = cv2.resize(preArray, (ns, nl), interpolation=cv2.INTER_NEAREST)
        return data

    def resize_file(self, infile, outfile, f=1):
        [preArray, im_width, im_height, im_bands, im_geotrans, im_proj] = self.getTiffData(infile)
        ns = im_width * f
        nl = im_height * f
        data = cv2.resize(preArray, (ns, nl), interpolation=cv2.INTER_BITS)
        self.writeTiff(data, ns, nl, 1, im_geotrans, im_proj, outfile)
        return data

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
            dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
        del dataset

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

    def resizefromSmall2Image(self, isno, preArray):
        if (isno == 1):
            newArray = cv2.resize(preArray[:, self.adjust_xmin_n:self.adjust_xmax_n],
                                  (self.ns_nadir, self.nl_nadir), cv2.INTER_NEAREST)
        else:
            newArray = cv2.resize(preArray[:, self.adjust_xmin_o:self.adjust_xmax_o],
                                  (self.ns_obliq, self.nl_obliq), cv2.INTER_NEAREST)
        return newArray

    ### get data from Nc

    def getVNIRdatafromNc_nadir(self, fileDir):

        ns = self.ns_nadir * 2
        nl = self.nl_nadir * 2
        # get red data
        fileName = fileDir + '\\S2_radiance_an.nc'
        objectName = r'S2_radiance_an'
        rad_red = self.getDatafromNc(fileName, objectName)
        fileName = fileDir + '\\S2_quality_an.nc'
        objectName = r'S2_solar_irradiance_an'
        sr_red = self.getDatafromNc(fileName, objectName)
        refl_red = rad_red / (sr_red[0] / np.pi)

        # get nir data
        fileName = fileDir + '\\S3_radiance_an.nc'
        objectName = r'S3_radiance_an'
        rad_nir = self.getDatafromNc(fileName, objectName)
        fileName = fileDir + '\\S3_quality_an.nc'
        objectName = r'S3_solar_irradiance_an'
        sr_nir = self.getDatafromNc(fileName, objectName)
        refl_nir = rad_nir / (sr_nir[0] / np.pi)

        # get nirr data
        fileName = fileDir + '\\S6_radiance_an.nc'
        objectName = r'S6_radiance_an'
        rad_ref = self.getDatafromNc(fileName, objectName)
        fileName = fileDir + '\\S6_quality_an.nc'
        objectName = r'S6_solar_irradiance_an'
        sr_ref = self.getDatafromNc(fileName, objectName)
        refl_ref = rad_ref / (sr_ref[0] / np.pi)

        outfile = fileDir + '\\refl_red.tif'
        refl_ref[refl_ref < 0] = 0
        self.writeTiff(refl_red[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\refl_ref.tif'
        self.writeTiff(refl_ref[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\refl_nir.tif'
        self.writeTiff(refl_nir[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\rad_nir.tif'
        self.writeTiff(rad_nir[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\rad_red.tif'
        self.writeTiff(rad_red[:nl, :ns], ns, nl, 1, '', '', outfile)

        return 1

    def getBLUEdatafromNc_nadir(self, fileDir):

        ns = self.ns_nadir
        nl = self.nl_nadir
        # get red data
        fileName = fileDir + '\\S1_radiance_an.nc'
        objectName = r'S1_radiance_an'
        rad_red = self.getDatafromNc(fileName, objectName)
        fileName = fileDir + '\\S1_quality_an.nc'
        objectName = r'S1_solar_irradiance_an'
        sr_red = self.getDatafromNc(fileName, objectName)
        refl_red = rad_red / (sr_red[0] / np.pi)

        outfile = fileDir + '\\refl_blue.tif'
        refl_blue=self.resize(refl_red,nl,ns)
        self.writeTiff(refl_blue[:nl, :ns], ns, nl, 1, '', '', outfile)

        return 1

    def getTIRdatafromNc_nadir(self, fileDir):
        # get red data
        fileName = fileDir + '\\S8_BT_in.nc'
        objectName = r'S8_BT_in'
        BT8 = self.getDatafromNc(fileName, objectName)
        BT8 = BT8[:self.nl_nadir, :]

        # get nir data
        fileName = fileDir + '\\S9_BT_in.nc'
        objectName = r'S9_BT_in'
        BT9 = self.getDatafromNc(fileName, objectName)
        BT9 = BT9[:self.nl_nadir, :]

        # get nirr data
        fileName = fileDir + '\\met_tx.nc'
        objectName = r'total_column_water_vapour_tx'
        tcw = self.getDatafromNc(fileName, objectName)
        tcw = tcw[0, :self.nl_nadir, :]
        tcw = self.resizefromSmall2Image(1, tcw)

        # fileName = fileDir + '\\geometry_tn.nc'
        # objectName = r'sat_zenith_tn'
        # vza = self.getDatafromNc(fileName, objectName)
        # vza =self.resizefromSmall2Image(1,vza)
        # vza = vza[:self.nl_nadir,:]
        # outfile = fileDir + '\\vza_n.tif'
        # self.writeTiff(tcw, self.ns_nadir, self.nl_nadir, 1, '', '', outfile)

        outfile = fileDir + '\\tcw.tif'
        self.writeTiff(tcw, self.ns_nadir, self.nl_nadir, 1, '', '', outfile)

        outfile = fileDir + '\\bt8.tif'
        self.writeTiff(BT8, self.ns_nadir, self.nl_nadir, 1, '', '', outfile)

        outfile = fileDir + '\\bt9.tif'
        self.writeTiff(BT9, self.ns_nadir, self.nl_nadir, 1, '', '', outfile)

        return 1

    def getVNIRdatafromNc_obliq(self, fileDir):
        ns = self.ns_obliq * 2
        nl = self.nl_obliq * 2
        # get red data
        fileName = fileDir + '\\S2_radiance_ao.nc'
        objectName = r'S2_radiance_ao'
        rad_red = self.getDatafromNc(fileName, objectName)
        fileName = fileDir + '\\S2_quality_ao.nc'
        objectName = r'S2_solar_irradiance_ao'
        sr_red = self.getDatafromNc(fileName, objectName)
        refl_red = rad_red / (sr_red[0] / np.pi)
        # rad_red = self.resize_oblique(rad_red)
        # refl_red = self.resize_oblique(refl_red)

        # get nir data
        fileName = fileDir + '\\S3_radiance_ao.nc'
        objectName = r'S3_radiance_ao'
        rad_nir = self.getDatafromNc(fileName, objectName)
        fileName = fileDir + '\\S3_quality_ao.nc'
        objectName = r'S3_solar_irradiance_ao'
        sr_nir = self.getDatafromNc(fileName, objectName)
        refl_nir = rad_nir / (sr_nir[0] / np.pi)
        # rad_nir = self.resize_oblique(rad_nir)
        # refl_nir = self.resize_oblique(refl_nir)

        # get nirr data
        fileName = fileDir + '\\S6_radiance_ao.nc'
        objectName = r'S6_radiance_ao'
        rad_ref = self.getDatafromNc(fileName, objectName)
        fileName = fileDir + '\\S6_quality_ao.nc'
        objectName = r'S6_solar_irradiance_ao'
        sr_ref = self.getDatafromNc(fileName, objectName)
        refl_ref = rad_ref / (sr_ref[0] / np.pi)
        # refl_ref = self.resize_oblique(refl_ref)

        outfile = fileDir + '\\refl_red_o.tif'
        refl_ref[refl_ref < 0] = 0
        self.writeTiff(refl_red[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\refl_ref_o.tif'
        self.writeTiff(refl_ref[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\refl_nir_o.tif'
        self.writeTiff(refl_nir[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\rad_nir_o.tif'
        self.writeTiff(rad_nir[:nl, :ns], ns, nl, 1, '', '', outfile)

        outfile = fileDir + '\\rad_red_o.tif'
        self.writeTiff(rad_red[:nl, :ns], ns, nl, 1, '', '', outfile)

        return 1

    def getTIRdatafromNc_obliq(self, fileDir):
        # get red data
        fileName = fileDir + '\\S8_BT_io.nc'
        objectName = r'S8_BT_io'
        BT8 = self.getDatafromNc(fileName, objectName)
        BT8 = BT8[:self.nl_obliq, :]

        # get nir data
        fileName = fileDir + '\\S9_BT_io.nc'
        objectName = r'S9_BT_io'
        BT9 = self.getDatafromNc(fileName, objectName)
        BT9 = BT9[:self.nl_obliq, :]

        # #get nirr data
        # fileName = fileDir + '\\met_tx.nc'
        # objectName = r'total_column_water_vapour_tx'
        # tcw = self.getDatafromNc(fileName, objectName)
        # tcw = tcw[0,:self.nl_obliq, :]
        # # tcw = np.reshape(tcw,[self.nl_nadir,-1])
        #
        # tcw = self.resizefromSmall2Image(1,tcw)
        # tcw = self.resizefromNadir2Obliq(tcw)

        # fileName = fileDir + '\\geometry_to.nc'
        # objectName = r'sat_zenith_to'
        # vza = self.getDatafromNc(fileName, objectName)
        # vza =self.resizefromSmall2Image(2,vza)
        # vza = vza[:self.nl_obliq,:]

        outfile = fileDir + '\\bt8_o.tif'
        self.writeTiff(BT8, self.ns_obliq, self.nl_obliq, 1, '', '', outfile)

        outfile = fileDir + '\\bt9_o.tif'
        self.writeTiff(BT9, self.ns_obliq, self.nl_obliq, 1, '', '', outfile)

        return 1

    def getGeodatafromNc_nadir(self, fileDir):
        # get red data

        ns = self.ns_nadir
        nl = self.nl_nadir

        fileName = fileDir + '\\geometry_tn.nc'
        objectName = r'sat_zenith_tn'
        vza = self.getDatafromNc(fileName, objectName)
        objectName = r'solar_zenith_tn'
        sza = self.getDatafromNc(fileName, objectName)
        objectName = r'sat_azimuth_tn'
        saa = self.getDatafromNc(fileName, objectName)
        objectName = r'solar_azimuth_tn'
        vaa = self.getDatafromNc(fileName, objectName)
        psi = np.abs(vaa - saa)

        vza = self.resizefromSmall2Image(1, vza)
        sza = self.resizefromSmall2Image(1, sza)
        psi = self.resizefromSmall2Image(1, psi)
        psi[psi > 180] = 360 - psi[psi > 180]

        fileName = fileDir + '\\geodetic_in.nc'
        objectName = r'latitude_in'
        lat = self.getDatafromNc(fileName, objectName)

        fileName = fileDir + '\\geodetic_in.nc'
        objectName = r'longitude_in'
        lon = self.getDatafromNc(fileName, objectName)

        fileName = fileDir + '\\flags_in.nc'
        objectName = r'cloud_in'
        cloud = self.getDatafromNc(fileName, objectName)

        outfile = fileDir + '\\lat.tif'
        self.writeTiff(lat[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\lon.tif'
        self.writeTiff(lon[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\cloud.tif'
        self.writeTiff(cloud[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\vza.tif'
        self.writeTiff(vza[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\sza.tif'
        self.writeTiff(sza[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\psi.tif'
        self.writeTiff(psi[:nl, :ns], ns, nl, 1, '', '', outfile)

        return 1

    def getGeodatafromNc_obliq(self, fileDir):
        # get red data

        ns = self.ns_obliq
        nl = self.nl_obliq

        fileName = fileDir + '\\geometry_to.nc'
        objectName = r'sat_zenith_to'
        vza = self.getDatafromNc(fileName, objectName)
        objectName = r'solar_zenith_to'
        sza = self.getDatafromNc(fileName, objectName)
        objectName = r'sat_azimuth_to'
        saa = self.getDatafromNc(fileName, objectName)
        objectName = r'solar_azimuth_to'
        vaa = self.getDatafromNc(fileName, objectName)
        psi = np.abs(vaa - saa)

        vza = self.resizefromSmall2Image(2, vza)
        sza = self.resizefromSmall2Image(2, sza)
        psi = self.resizefromSmall2Image(2, psi)
        psi[psi > 180] = 360 - psi[psi > 180]

        fileName = fileDir + '\\geodetic_io.nc'
        objectName = r'latitude_io'
        lat = self.getDatafromNc(fileName, objectName)

        fileName = fileDir + '\\geodetic_io.nc'
        objectName = r'longitude_io'
        lon = self.getDatafromNc(fileName, objectName)

        fileName = fileDir + '\\flags_io.nc'
        objectName = r'cloud_io'
        cloud = self.getDatafromNc(fileName, objectName)

        outfile = fileDir + '\\lat_o.tif'
        self.writeTiff(lat[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\lon_o.tif'
        self.writeTiff(lon[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\cloud_o.tif'
        self.writeTiff(cloud[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\vza_o.tif'
        self.writeTiff(vza[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\sza_o.tif'
        self.writeTiff(sza[:nl, :ns], ns, nl, 1, '', '', outfile)
        outfile = fileDir + '\\psi_o.tif'
        self.writeTiff(psi[:nl, :ns], ns, nl, 1, '', '', outfile)

        return 1

    def getLatLonfromPoint_nadir(self,fileDir):
        infile = fileDir + r'\lat_n_p.tif'
        [lat, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # lat = self.resizefromSmall2Image(1,lat)

        infile = fileDir + r'\lon_n_p.tif'
        [lon, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        return lat,lon,temp3,temp2,temp1

    def getNDVIofPoint_nadir(self,wdir,fileName,lat_,lon_,off=0):
        fileDir = wdir + 'data\\' + fileName.strip()

        infile = fileDir + r'\lat_v.tif'
        [lat, ns0, nl0, temp3, temp4, temp5] = self.getTiffData(infile)

        infile = fileDir + r'\lon_v.tif'
        [lon, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)

        infile = fileDir + r'\cloud.tif'
        cloud, temp1, temp2, temp3, temp4, temp5 = self.getTiffData(infile)
        infile = fileDir + r'\refl_red.tif'
        red, temp1, temp2, temp3, temp4, temp5 = self.getTiffData(infile)
        infile = fileDir + r'\refl_nir.tif'
        nir, temp1, temp2, temp3, temp4, temp5 = self.getTiffData(infile)

        numpoint = len(lat_)
        line = np.zeros(numpoint)
        for kb in range(numpoint):
            lat0 = lat_[kb]
            lon0 = lon_[kb]
            dis = (lat-lat0)*(lat-lat0)+(lon-lon0)*(lon-lon0)
            dmin = np.min(dis)
            if(dmin < 0.01):
                pos = np.where(dis == dmin)
                ns = pos[1][0]
                nl = pos[0][0]



                if((ns > ns0) or (nl>nl0) or (ns<0) or (nl<0)):
                    continue
                nsshort = np.int(ns/2)
                nlshort = np.int(nl/2)
                if(cloud[nlshort,nsshort]!=0):continue

                ndvi = (nir-red)/(nir+red)

                # nl1 = nl-off
                # nl2 = nl+off+1
                # if (nl1<0): nl1 = 0
                # if(nl2>nl0): nl2 =nl0
                # ns1 = ns-off
                # ns2 = ns+off+1
                # if (ns1<0): ns1 = 0
                # if(ns2>ns0): ns2 = ns0

                line[kb] = ndvi[nl,ns]

        return line


    def writevrt_nadir(self, infile, datafile, xfile, yfile):

        f = open(infile, 'w')
        f.write(r'<VRTDataset rasterXSize="1500" rasterYSize="1200">')
        f.write('\n')
        f.write(r'  <Metadata domain="GEOLOCATION">')
        f.write('\n')
        f.write(
            r'    <MDI key="SRS">GEOGCS["WGS 84(DD)",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AXIS["Long",EAST],AXIS["Lat",NORTH]]</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_DATASET">' + xfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_STEP">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_DATASET">' + yfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_STEP">1</MDI>')
        f.write('\n')
        f.write(r'  </Metadata>')
        f.write('\n')
        f.write(r'  <VRTRasterBand dataType = "Float32" band = "1">')
        f.write('\n')
        f.write(r'    <ColorInterp>Gray</ColorInterp >')
        f.write('\n')
        f.write(r'    <NoDataValue>0</NoDataValue >')
        f.write('\n')
        f.write(r'    <SimpleSource>')
        f.write('\n')
        f.write(r'      <SourceFilename relativeToVRT = "1" >' + datafile + r'</SourceFilename>')
        f.write('\n')
        f.write(r'      <SourceBand>1</SourceBand>')
        f.write('\n')
        f.write(r'    </SimpleSource>')
        f.write('\n')
        f.write(r'  </VRTRasterBand>')
        f.write('\n')
        f.write('</VRTDataset>')
        f.close()
        return 1

    def writevrt_obliq(self, infile, datafile, xfile, yfile):

        f = open(infile, 'w')
        f.write(r'<VRTDataset rasterXSize="900" rasterYSize="1200">')
        f.write('\n')
        f.write(r'  <Metadata domain="GEOLOCATION">')
        f.write('\n')
        f.write(
            r'    <MDI key="SRS">GEOGCS["WGS 84(DD)",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AXIS["Long",EAST],AXIS["Lat",NORTH]]</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_DATASET">' + xfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_STEP">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_DATASET">' + yfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_STEP">1</MDI>')
        f.write('\n')
        f.write(r'  </Metadata>')
        f.write('\n')
        f.write(r'  <VRTRasterBand dataType = "Float32" band = "1">')
        f.write('\n')
        f.write(r'    <ColorInterp>Gray</ColorInterp >')
        f.write('\n')
        f.write(r'    <NoDataValue>0</NoDataValue >')
        f.write('\n')
        f.write(r'    <SimpleSource>')
        f.write('\n')
        f.write(r'      <SourceFilename relativeToVRT = "1" >' + datafile + r'</SourceFilename>')
        f.write('\n')
        f.write(r'      <SourceBand>1</SourceBand>')
        f.write('\n')
        f.write(r'    </SimpleSource>')
        f.write('\n')
        f.write(r'  </VRTRasterBand>')
        f.write('\n')
        f.write('</VRTDataset>')
        f.close()
        return 1

    def writevrt_v_nadir(self, infile, datafile, xfile, yfile):

        f = open(infile, 'w')
        f.write(r'<VRTDataset rasterXSize="3000" rasterYSize="2400">')
        f.write('\n')
        f.write(r'  <Metadata domain="GEOLOCATION">')
        f.write('\n')
        f.write(
            r'    <MDI key="SRS">GEOGCS["WGS 84(DD)",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AXIS["Long",EAST],AXIS["Lat",NORTH]]</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_DATASET">' + xfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_STEP">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_DATASET">' + yfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_STEP">1</MDI>')
        f.write('\n')
        f.write(r'  </Metadata>')
        f.write('\n')
        f.write(r'  <VRTRasterBand dataType = "Float32" band = "1">')
        f.write('\n')
        f.write(r'    <ColorInterp>Gray</ColorInterp >')
        f.write('\n')
        f.write(r'    <NoDataValue>0</NoDataValue >')
        f.write('\n')
        f.write(r'    <SimpleSource>')
        f.write('\n')
        f.write(r'      <SourceFilename relativeToVRT = "1" >' + datafile + r'</SourceFilename>')
        f.write('\n')
        f.write(r'      <SourceBand>1</SourceBand>')
        f.write('\n')
        f.write(r'    </SimpleSource>')
        f.write('\n')
        f.write(r'  </VRTRasterBand>')
        f.write('\n')
        f.write('</VRTDataset>')
        f.close()
        return 1

    def writevrt_v_obliq(self, infile, datafile, xfile, yfile):

        f = open(infile, 'w')
        f.write(r'<VRTDataset rasterXSize="1800" rasterYSize="2400">')
        f.write('\n')
        f.write(r'  <Metadata domain="GEOLOCATION">')
        f.write('\n')
        f.write(
            r'    <MDI key="SRS">GEOGCS["WGS 84(DD)",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433],AXIS["Long",EAST],AXIS["Lat",NORTH]]</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_DATASET">' + xfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="X_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="PIXEL_STEP">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_DATASET">' + yfile + r'</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="Y_BAND">1</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_OFFSET">0</MDI>')
        f.write('\n')
        f.write(r'    <MDI key="LINE_STEP">1</MDI>')
        f.write('\n')
        f.write(r'  </Metadata>')
        f.write('\n')
        f.write(r'  <VRTRasterBand dataType = "Float32" band = "1">')
        f.write('\n')
        f.write(r'    <ColorInterp>Gray</ColorInterp >')
        f.write('\n')
        f.write(r'    <NoDataValue>65535</NoDataValue >')
        f.write('\n')
        f.write(r'    <SimpleSource>')
        f.write('\n')
        f.write(r'      <SourceFilename relativeToVRT = "1" >' + datafile + r'</SourceFilename>')
        f.write('\n')
        f.write(r'      <SourceBand>1</SourceBand>')
        f.write('\n')
        f.write(r'    </SimpleSource>')
        f.write('\n')
        f.write(r'  </VRTRasterBand>')
        f.write('\n')
        f.write('</VRTDataset>')
        f.close()
        return 1

    # spatial issue

    def getSRSPair(self, dataset):
        '''
        获得给定数据的投影参考系和地理参考系
        :param dataset: GDAL地理数据
        :return: 投影参考系和地理参考系
        '''
        prosrs = osr.SpatialReference()
        prosrs.ImportFromWkt(dataset.GetProjection())
        geosrs = prosrs.CloneGeogCS()
        return prosrs, geosrs

    def geo2lonlat(self, dataset, x, y):
        '''
        将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
        :param dataset: GDAL地理数据
        :param x: 投影坐标x
        :param y: 投影坐标y
        :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
        '''
        prosrs, geosrs = self.getSRSPair(dataset)
        ct = osr.CoordinateTransformation(prosrs, geosrs)
        temp = np.asarray([y, x])
        temp = np.transpose(temp)
        coords = np.asarray(ct.TransformPoints(temp))
        # coords = ct.TransformPoint(x, y)
        return coords[:, 0], coords[:, 1]

    def lonlat2geo(self, dataset, lat, lon):
        '''
            将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
            :param dataset: GDAL地理数据
            :param lon: 地理坐标lon经度
            :param lat: 地理坐标lat纬度
            :return: 经纬度坐标(lon, lat)对应的投影坐标
        '''
        # dataset = gdal.Open(fileName, gdal.GA_ReadOnly)
        prosrs, geosrs = self.getSRSPair(dataset)
        ct = osr.CoordinateTransformation(geosrs, prosrs)
        lon = np.reshape(lon, [-1])
        lat = np.reshape(lat, [-1])
        temp = np.asarray([lon, lat])
        temp = np.transpose(temp)
        # temp = np.asarray([lat[0:2],lon[0:2]])
        coords = np.asarray(ct.TransformPoints(temp))

        return coords[:, 0], coords[:, 1]

    def geo2imagexy(self, dataset, x, y):
        '''
        根据GDAL的六 参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
        :param dataset: GDAL地理数据
        :param x: 投影或地理坐标x
        :param y: 投影或地理坐标y
        :return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col)
        '''
        trans = dataset.GetGeoTransform()
        a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
        b = np.array([x - trans[0], y - trans[3]])
        c = np.linalg.solve(a, b)
        return c  # 使用numpy的linalg.solve进行二元一次方程的求解

    def imagexy2geo(self, dataset, row, col):
        '''
        根据GDAL的六参数模型将影像图上坐标（行列号）转为投影坐标或地理坐标（根据具体数据的坐标系统转换）
        :param dataset: GDAL地理数据
        :param row: 像素的行号
        :param col: 像素的列号
        :return: 行列号(row, col)对应的投影坐标或地理坐标(x, y)
        '''
        trans = dataset.GetGeoTransform()
        px = trans[0] + col * trans[1] + row * trans[2]
        py = trans[3] + col * trans[4] + row * trans[5]
        return px, py

    ### Inversion

    def inversionAOD_nadir(self, wdir,dataTime,medName,resultName):

        # fileDir = fileDir = wdir+'data\\'+fileName.strip()

        # lookup table
        infile = wdir + r'auxiliary\SLSTR_BLUE_AOD_LUT.txt'
        refllut = self.getDatafromTxt(infile, 1, 6)
        refllut = refllut[:, 5]
        refl_ = np.linspace(0, 0.15, 16)
        sza_ = np.asarray([0, 6, 12, 24, 36, 48, 54, 60, 66])
        vza_ = np.asarray([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72])
        psi_ = np.asarray([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180])
        aod_ = np.asarray([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0])
        szanum = len(sza_)
        vzanum = len(vza_)
        psinum = len(psi_)
        reflnum = len(refl_)
        aodnum = len(aod_)

        # infile = wdir+r'heihe_'+dataTime+'_rad_red.tif'
        # [rad_red,temp1,temp2,temp3,temp4] = self.getTiffData(infile)
        # infile = wdir+r'heihe_'+dataTime+'_rad_nir.tif'
        # [rad_nir,temp1,temp2,temp3,temp4] = self.getTiffData(infile)


        infile = wdir + medName + dataTime + '_refl_red_n.tif'
        [refl_red, ns, nl, nb, geog, proj] = self.getTiffData(infile)



        # ns = np.int(ns*0.5)
        # nl = np.int(nl*0.5)
        # refl_red = self.resize(refl_red, nl, ns)
        ccc = np.sum(refl_red)
        if (ccc == 0):
            return 0
        # infile = wdir+r'heihe_'+dataTime+'_refl_nir.tif'
        # [refl_nir,temp1,temp2,temp3,temp4] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_refl_ref_n.tif'
        [refl_ref, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_refl_blue_n.tif'
        [refl_blue, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_tcw_n.tif'
        [tcw, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        tcw = tcw/10.0
        # refl_ref = self.resize(refl_ref, nl, ns)
        infile = wdir + medName + dataTime + '_sza_n.tif'
        [sza, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        rd = np.pi / 180.0
        refl_red = refl_red / np.cos(sza * rd)
        refl_ref = refl_ref/ np.cos(sza*rd)
        refl_blue = refl_blue / np.cos(sza * rd)
        infile = wdir + medName + dataTime + '_vza_n.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_psi_n.tif'
        [psi, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)

        szamask = np.zeros([nl, ns], dtype=np.int)
        vzamask = np.zeros([nl, ns], dtype=np.int)
        psimask = np.zeros([nl, ns], dtype=np.int)
        reflmask = np.zeros([nl, ns], dtype=np.int)
        aod = np.zeros([nl, ns], dtype=np.float)

        # agl = np.cos(np.cos())
        tcw = tcw/(np.cos(sza*rd)*np.cos(vza*rd))
        refl_ref = refl_ref/(0.000002848*tcw*tcw-0.006428*tcw+0.898)
        refl_refd2 = refl_ref*0.5


        for krefl in range(reflnum - 1):
            ind = (refl_refd2 >= refl_[krefl]) * (refl_refd2 < refl_[krefl + 1])
            if (np.sum(ind) > 0): reflmask[ind] = krefl

        for ksza in range(szanum - 1):
            ind = (sza >= sza_[ksza]) * (sza < sza_[ksza + 1])
            if (np.sum(ind) > 0): szamask[ind] = ksza
        # this part is for pixels with sza > 48. It works but it is not very suitable
        ind = sza >= sza_[szanum - 1]
        szamask[ind] = szanum - 2

        for kvza in range(vzanum - 1):
            ind = (vza >= vza_[kvza]) * (vza < vza_[kvza + 1])
            if (np.sum(ind) > 0): vzamask[ind] = kvza

        for kpsi in range(psinum - 1):
            ind = (psi >= psi_[kpsi]) * (psi < psi_[kpsi + 1])
            if (np.sum(ind) > 0): psimask[ind] = kpsi

        num1 = vzanum * psinum * aodnum * reflnum
        num2 = psinum * aodnum * reflnum
        num3 = aodnum * reflnum

        # valid pixels
        indd = (vzamask != -1) * (refl_red > 0)
        if (np.sum(indd) <= 0):
            print('no suitable pixel can be used for vod inversion')
            return 0

        # sza1, vza1, psi1
        ind1 = szamask * num1 + vzamask * num2 + psimask * num3 + reflmask
        # sza2, vza1, psi1
        ind2 = (szamask + 1) * num1 + vzamask * num2 + psimask * num3 + reflmask
        # sza1, vza2, psi1
        ind3 = szamask * num1 + (vzamask + 1) * num2 + psimask * num3 + reflmask
        # sza2, vza2, psi1
        ind4 = (szamask + 1) * num1 + (vzamask + 1) * num2 + psimask * num3 + reflmask
        # sza1, vza1, psi2
        ind5 = szamask * num1 + (vzamask) * num2 + (psimask + 1) * num3 + reflmask
        # sza2, vza1, psi2
        ind6 = (szamask + 1) * num1 + (vzamask) * num2 + (psimask + 1) * num3 + reflmask
        # sza1, vza2, psi2
        ind7 = szamask * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + reflmask
        # sza2, vza2, psi2
        ind8 = (szamask + 1) * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + reflmask

        ind1 = ind1[indd]
        ind2 = ind2[indd]
        ind3 = ind3[indd]
        ind4 = ind4[indd]
        ind5 = ind5[indd]
        ind6 = ind6[indd]
        ind7 = ind7[indd]
        ind8 = ind8[indd]

        reflref = np.zeros([np.sum(indd), aodnum])
        aodarr = np.linspace(0, aodnum - 1, aodnum, dtype=np.int) * reflnum
        for k in range(aodnum):
            prop = (sza[indd] - sza_[szamask[indd]]) / (sza_[szamask[indd] + 1] - sza_[szamask[indd]])

            ind = ind1 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind2 + aodarr[k]
            temp2 = refllut[ind]
            refl1 = (1 - prop) * temp1 + prop * temp2

            ind = ind3 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind4 + aodarr[k]
            temp2 = refllut[ind]
            refl2 = (1 - prop) * temp1 + prop * temp2

            ind = ind5 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind6 + aodarr[k]
            temp2 = refllut[ind]
            refl3 = (1 - prop) * temp1 + prop * temp2

            ind = ind7 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind8 + aodarr[k]
            temp2 = refllut[ind]
            refl4 = (1 - prop) * temp1 + prop * temp2

            prop = (vza[indd] - vza_[vzamask[indd]]) / (vza_[vzamask[indd] + 1] - vza_[vzamask[indd]])
            refl5 = (1 - prop) * refl1 + prop * refl2
            refl6 = (1 - prop) * refl3 + prop * refl4

            prop = 1.0 * (psi[indd] - psi_[psimask[indd]]) / (psi_[psimask[indd] + 1] - psi_[psimask[indd]])
            refl7 = (1 - prop) * refl5 + prop * refl6
            reflref[:, k] = refl7

        refltemp = refl_red[indd]
        aodtemp = np.zeros(np.sum(indd))
        for k in range(np.sum(indd)):
            temp1 = reflref[k, :]
            temp2 = refltemp[k]

            if (temp2 > temp1[aodnum - 1]): aodtemp[k] = (aod_[aodnum - 1])
            if (temp2 < temp1[0]): aodtemp[k] = aod_[0]
            for kk in range(aodnum - 1):
                if temp2 >= temp1[kk] and temp2 < temp1[kk + 1]:
                    down = (temp1[kk + 1] - temp1[kk])
                    prop = (temp2 - temp1[kk]) / down
                    aodtemp[k] = (1 - prop) * aod_[kk] + prop * aod_[kk + 1]
                    break

        aod[indd] = aodtemp

        outfile = wdir + resultName + dataTime + '_aod_n.tif'
        self.writeTiff(aod, ns, nl, 1, geog, proj, outfile)

        # outfile = fileDir+'\\vza_n.tif'
        # self.writeTiff(vza,self.ns_nadir,self.nl_nadir,1,'','',outfile)
        #
        # outfile = fileDir+'\\sza_n.tif'
        # self.writeTiff(sza,self.ns_nadir,self.nl_nadir,1,'','',outfile)
        # outfile = fileDir+'\\refl6_n.tif'
        # temp = np.zeros([aodnum,self.nl_nadir,self.ns_nadir])
        # for k in range(aodnum): temp[k,indd] = reflref[:,k]
        # temp[temp<0] = 0
        # self.writeTiff(temp,self.ns_nadir,self.nl_nadir,aodnum,'','',outfile)

        return 1

    def inversionAOD_obliq(self, wdir, dataTime,medName,resultName):

        # fileDir = fileDir = wdir+'data\\'+fileName.strip()

        # lookup table
        infile = wdir + r'auxiliary\SLSTR_Red_AOD_LUT.txt'
        refllut = self.getDatafromTxt(infile, 1, 6)
        refllut = refllut[:, 5]
        refl_ = np.linspace(0, 0.15, 16)
        sza_ = np.asarray([0, 6, 12, 24, 36, 48, 54, 60, 66])
        vza_ = np.asarray([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72])
        psi_ = np.asarray([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180])
        aod_ = np.asarray([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0])
        szanum = len(sza_)
        vzanum = len(vza_)
        psinum = len(psi_)
        reflnum = len(refl_)
        aodnum = len(aod_)

        # infile = wdir+r'heihe_'+dataTime+'_rad_red.tif'
        # [rad_red,temp1,temp2,temp3,temp4] = self.getTiffData(infile)
        # infile = wdir+r'heihe_'+dataTime+'_rad_nir.tif'
        # [rad_nir,temp1,temp2,temp3,temp4] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_refl_red_n.tif'
        [refl_red, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # ns = np.int(ns*0.5)
        # nl = np.int(nl*0.5)
        # refl_red = self.resize(refl_red, nl, ns)
        ccc = np.sum(refl_red)
        if (ccc == 0):
            return 0
        # infile = wdir+r'heihe_'+dataTime+'_refl_nir.tif'
        # [refl_nir,temp1,temp2,temp3,temp4] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_refl_ref_o.tif'
        [refl_ref, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # refl_ref = self.resize(refl_ref, nl, ns)
        infile = wdir + medName + dataTime + '_sza_o.tif'
        [sza, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        rd = np.pi / 180.0
        refl_red = refl_red / np.cos(sza * rd)
        refl_ref = refl_ref / np.cos(sza * rd)
        infile = wdir + medName + dataTime + '_vza_o.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_psi_o.tif'
        [psi, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)

        szamask = np.zeros([nl, ns], dtype=np.int)
        vzamask = np.zeros([nl, ns], dtype=np.int)
        psimask = np.zeros([nl, ns], dtype=np.int)
        reflmask = np.zeros([nl, ns], dtype=np.int)
        aod = np.zeros([nl, ns], dtype=np.float)

        refl_refd2 = refl_ref / 2.0
        for krefl in range(reflnum - 1):
            ind = (refl_refd2 >= refl_[krefl]) * (refl_refd2 < refl_[krefl + 1])
            if (np.sum(ind) > 0): reflmask[ind] = krefl

        for ksza in range(szanum - 1):
            ind = (sza >= sza_[ksza]) * (sza < sza_[ksza + 1])
            if (np.sum(ind) > 0): szamask[ind] = ksza
        # this part is for pixels with sza > 48. It works but it is not very suitable
        ind = sza >= sza_[szanum - 1]
        szamask[ind] = szanum - 2

        for kvza in range(vzanum - 1):
            ind = (vza >= vza_[kvza]) * (vza < vza_[kvza + 1])
            if (np.sum(ind) > 0): vzamask[ind] = kvza

        for kpsi in range(psinum - 1):
            ind = (psi >= psi_[kpsi]) * (psi < psi_[kpsi + 1])
            if (np.sum(ind) > 0): psimask[ind] = kpsi

        num1 = vzanum * psinum * aodnum * reflnum
        num2 = psinum * aodnum * reflnum
        num3 = aodnum * reflnum

        # valid pixels
        indd = (vzamask != -1) * (refl_red > 0)
        if (np.sum(indd) <= 0):
            print('no suitable pixel can be used for vod inversion')
            return 0

        # sza1, vza1, psi1
        ind1 = szamask * num1 + vzamask * num2 + psimask * num3 + reflmask
        # sza2, vza1, psi1
        ind2 = (szamask + 1) * num1 + vzamask * num2 + psimask * num3 + reflmask
        # sza1, vza2, psi1
        ind3 = szamask * num1 + (vzamask + 1) * num2 + psimask * num3 + reflmask
        # sza2, vza2, psi1
        ind4 = (szamask + 1) * num1 + (vzamask + 1) * num2 + psimask * num3 + reflmask
        # sza1, vza1, psi2
        ind5 = szamask * num1 + (vzamask) * num2 + (psimask + 1) * num3 + reflmask
        # sza2, vza1, psi2
        ind6 = (szamask + 1) * num1 + (vzamask) * num2 + (psimask + 1) * num3 + reflmask
        # sza1, vza2, psi2
        ind7 = szamask * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + reflmask
        # sza2, vza2, psi2
        ind8 = (szamask + 1) * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + reflmask

        ind1 = ind1[indd]
        ind2 = ind2[indd]
        ind3 = ind3[indd]
        ind4 = ind4[indd]
        ind5 = ind5[indd]
        ind6 = ind6[indd]
        ind7 = ind7[indd]
        ind8 = ind8[indd]

        reflref = np.zeros([np.sum(indd), aodnum])
        aodarr = np.linspace(0, aodnum - 1, aodnum, dtype=np.int) * reflnum
        for k in range(aodnum):
            prop = (sza[indd] - sza_[szamask[indd]]) / (sza_[szamask[indd] + 1] - sza_[szamask[indd]])

            ind = ind1 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind2 + aodarr[k]
            temp2 = refllut[ind]
            refl1 = (1 - prop) * temp1 + prop * temp2

            ind = ind3 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind4 + aodarr[k]
            temp2 = refllut[ind]
            refl2 = (1 - prop) * temp1 + prop * temp2

            ind = ind5 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind6 + aodarr[k]
            temp2 = refllut[ind]
            refl3 = (1 - prop) * temp1 + prop * temp2

            ind = ind7 + aodarr[k]
            temp1 = refllut[ind]
            ind = ind8 + aodarr[k]
            temp2 = refllut[ind]
            refl4 = (1 - prop) * temp1 + prop * temp2

            prop = (vza[indd] - vza_[vzamask[indd]]) / (vza_[vzamask[indd] + 1] - vza_[vzamask[indd]])
            refl5 = (1 - prop) * refl1 + prop * refl2
            refl6 = (1 - prop) * refl3 + prop * refl4

            prop = 1.0 * (psi[indd] - psi_[psimask[indd]]) / (psi_[psimask[indd] + 1] - psi_[psimask[indd]])
            refl7 = (1 - prop) * refl5 + prop * refl6
            reflref[:, k] = refl7

        refltemp = refl_red[indd]
        aodtemp = np.zeros(np.sum(indd))
        for k in range(np.sum(indd)):
            temp1 = reflref[k, :]
            temp2 = refltemp[k]

            if (temp2 > temp1[aodnum - 1]): aodtemp[k] = (aod_[aodnum - 1])
            if (temp2 < temp1[0]): aodtemp[k] = aod_[0]
            for kk in range(aodnum - 1):
                if temp2 >= temp1[kk] and temp2 < temp1[kk + 1]:
                    down = (temp1[kk + 1] - temp1[kk])
                    prop = (temp2 - temp1[kk]) / down
                    aodtemp[k] = (1 - prop) * aod_[kk] + prop * aod_[kk + 1]
                    break

        aod[indd] = aodtemp

        outfile = wdir + resultName + dataTime + '_aod_o.tif'
        self.writeTiff(aod, ns, nl, 1, geog, proj, outfile)

        # outfile = fileDir+'\\vza_n.tif'
        # self.writeTiff(vza,self.ns_nadir,self.nl_nadir,1,'','',outfile)
        #
        # outfile = fileDir+'\\sza_n.tif'
        # self.writeTiff(sza,self.ns_nadir,self.nl_nadir,1,'','',outfile)
        # outfile = fileDir+'\\refl6_n.tif'
        # temp = np.zeros([aodnum,self.nl_nadir,self.ns_nadir])
        # for k in range(aodnum): temp[k,indd] = reflref[:,k]
        # temp[temp<0] = 0
        # self.writeTiff(temp,self.ns_nadir,self.nl_nadir,aodnum,'','',outfile)

        return 1

    def VNIRatm_nadir(self, wdir, dataTime,medName,resultName,aodcheck=2):

        # lookup table
        infile = wdir + r'auxiliary\SLSTR_Red_LUT.txt'
        temp = self.getDatafromTxt(infile, 1, 8)
        temp1 = temp[:, 5:8]
        infile = wdir + r'auxiliary\SLSTR_NIR_LUT.txt'
        temp = self.getDatafromTxt(infile, 1, 8)
        temp2 = temp[:, 5:8]
        refllut = np.hstack((temp1, temp2))

        refl_ = np.asarray([0])  # np.linspace(0,0.15,16)
        sza_ = np.asarray([0, 6, 12, 24, 36, 48, 54, 60, 66])
        vza_ = np.asarray([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72])
        psi_ = np.asarray([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180])
        aod_ = np.asarray([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0])
        szanum = len(sza_)
        vzanum = len(vza_)
        psinum = len(psi_)
        reflnum = len(refl_)
        aodnum = len(aod_)

        infile = wdir + medName + dataTime + '_refl_red_n.tif'
        [refl_red, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_refl_nir_n.tif'
        [refl_nir, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # ns = np.int(ns*0.5)
        # nl = np.int(nl*0.5)
        # refl_red = self.resize(refl_red,nl,ns)
        infile = wdir + medName + dataTime + '_rad_red_n.tif'
        [rad_red, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # rad_red = self.resize(rad_red, nl, ns)
        infile = wdir + medName + dataTime + '_rad_nir_n.tif'
        [rad_nir, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # rad_nir = self.resize(rad_nir, nl, ns)

        infile = wdir + medName + dataTime + '_sza_n.tif'
        [sza, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        rd = np.pi / 180.0
        refl_red = refl_red / np.cos(sza * rd)


        ccc = np.sum(refl_red)
        if (ccc == 0):
            return 0
        # infile = wdir+medName+dataTime+'_refl_nir.tif'
        # [refl_nir,temp1,temp2,temp3,temp4,temp5] = self.getTiffData(infile)
        # infile = wdir+medName+dataTime+'_refl_ref.tif'
        # [refl_ref,temp1,temp2,temp3,temp4,temp5] = self.getTiffData(infile)
        # infile = wdir + medName + dataTime + '_sza_n.tif'
        # [sza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_vza_n.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_psi_n.tif'
        [psi, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + '_aod_n.tif'

        c6 = os.path.exists(infile)
        if c6 ==0: return
        [aod, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)




        ###########################################################################
        # aodvalue = np.average(aod[(430 - 2):(430 + 2),(375 - 2):(375 + 2)])

        if aodcheck == 1:
            ind = (aod<=self.aodmax)*(aod>=self.aodmin)
            if (np.sum(ind)<0):
                aodvalue = self.aodvalue
            else:
                aodvalue = np.average(aod[ind])
            if aodvalue < self.aodmin: aodvalue = self.aodmin
            aod[aod != 0] = aodvalue
        if aodcheck ==2:
            aodvalue = self.aodvalue
            aod[aod != 0] = aodvalue

        # print(aodvalue)
        # aodvalue = self.aodvalue
        # aod[aod != 0] = aodvalue
        ###########################################################################
        szamask = np.zeros([nl, ns], dtype=np.int)
        vzamask = np.zeros([nl, ns], dtype=np.int)
        psimask = np.zeros([nl, ns], dtype=np.int)
        aodmask = np.zeros([nl, ns], dtype=np.int)

        for kaod in range(aodnum - 1):
            ind = (aod >= aod_[kaod]) * (aod < aod_[kaod + 1])
            if (np.sum(ind) > 0): aodmask[ind] = kaod
        ind = aod >= aod_[aodnum - 1]
        aodmask[ind] = aodnum - 2

        for ksza in range(szanum - 1):
            ind = (sza >= sza_[ksza]) * (sza < sza_[ksza + 1])
            if (np.sum(ind) > 0): szamask[ind] = ksza
        # this part is for pixels with sza > 48. It works but it is not very suitable
        ind = sza >= sza_[szanum - 1]
        szamask[ind] = szanum - 2

        for kvza in range(vzanum - 1):
            ind = (vza >= vza_[kvza]) * (vza < vza_[kvza + 1])
            if (np.sum(ind) > 0): vzamask[ind] = kvza

        for kpsi in range(psinum - 1):
            ind = (psi >= psi_[kpsi]) * (psi < psi_[kpsi + 1])
            if (np.sum(ind) > 0): psimask[ind] = kpsi

        num1 = vzanum * psinum * aodnum * reflnum
        num2 = psinum * aodnum * reflnum
        num3 = aodnum * reflnum

        # valid pixels
        indd = (vzamask != -1) * (refl_red > 0)
        if (np.sum(indd) <= 0):
            print('no suitable pixel can be used for vod inversion')
            return 0

        # sza1, vza1, psi1
        ind1 = szamask * num1 + vzamask * num2 + psimask * num3 + aodmask + 1
        # sza2, vza1, psi1
        ind2 = (szamask + 1) * num1 + vzamask * num2 + psimask * num3 + aodmask + 1
        # sza1, vza2, psi1
        ind3 = szamask * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask + 1
        # sza2, vza2, psi1
        ind4 = (szamask + 1) * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask + 1
        # sza1, vza1, psi2
        ind5 = szamask * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask + 1
        # sza2, vza1, psi2
        ind6 = (szamask + 1) * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask + 1
        # sza1, vza2, psi2
        ind7 = szamask * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask + 1
        # sza2, vza2, psi2
        ind8 = (szamask + 1) * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask + 1

        ind1 = ind1[indd]
        ind2 = ind2[indd]
        ind3 = ind3[indd]
        ind4 = ind4[indd]
        ind5 = ind5[indd]
        ind6 = ind6[indd]
        ind7 = ind7[indd]
        ind8 = ind8[indd]

        coeff1 = np.zeros([np.sum(indd), 6])
        for k in range(6):
            prop = (sza[indd] - sza_[szamask[indd]]) / (sza_[szamask[indd] + 1] - sza_[szamask[indd]])
            kk = 0
            ind = ind1 + kk
            temp1 = refllut[ind, k]
            ind = ind2 + kk
            temp2 = refllut[ind, k]
            refl1 = (1 - prop) * temp1 + prop * temp2

            ind = ind3 + kk
            temp1 = refllut[ind, k]
            ind = ind4 + kk
            temp2 = refllut[ind, k]
            refl2 = (1 - prop) * temp1 + prop * temp2

            ind = ind5 + kk
            temp1 = refllut[ind, k]
            ind = ind6 + kk
            temp2 = refllut[ind, k]
            refl3 = (1 - prop) * temp1 + prop * temp2

            ind = ind7 + kk
            temp1 = refllut[ind, k]
            ind = ind8 + kk
            temp2 = refllut[ind, k]
            refl4 = (1 - prop) * temp1 + prop * temp2

            prop = (vza[indd] - vza_[vzamask[indd]]) / (vza_[vzamask[indd] + 1] - vza_[vzamask[indd]])
            refl5 = (1 - prop) * refl1 + prop * refl2
            refl6 = (1 - prop) * refl3 + prop * refl4

            prop = 1.0 * (psi[indd] - psi_[psimask[indd]]) / (psi_[psimask[indd] + 1] - psi_[psimask[indd]])
            refl7 = (1 - prop) * refl5 + prop * refl6
            coeff1[:, k] = refl7

        #  sza1, vza1, psi1
        ind1 = szamask * num1 + vzamask * num2 + psimask * num3 + aodmask
        #  sza2, vza1, psi1
        ind2 = (szamask + 1) * num1 + vzamask * num2 + psimask * num3 + aodmask
        #  sza1, vza2, psi1
        ind3 = szamask * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask
        #  sza2, vza2, psi1
        ind4 = (szamask + 1) * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask
        # sza1, vza1, psi2
        ind5 = szamask * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask
        #  sza2, vza1, psi2
        ind6 = (szamask + 1) * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask
        #  sza1, vza2, psi2
        ind7 = szamask * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask
        #  sza2, vza2, psi2
        ind8 = (szamask + 1) * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask

        ind1 = ind1[indd]
        ind2 = ind2[indd]
        ind3 = ind3[indd]
        ind4 = ind4[indd]
        ind5 = ind5[indd]
        ind6 = ind6[indd]
        ind7 = ind7[indd]
        ind8 = ind8[indd]

        coeff2 = np.zeros([np.sum(indd), 6])
        for k in range(6):
            prop = (sza[indd] - sza_[szamask[indd]]) / (sza_[szamask[indd] + 1] - sza_[szamask[indd]])
            kk = 0
            ind = ind1 + kk
            temp1 = refllut[ind, k]
            ind = ind2 + kk
            temp2 = refllut[ind, k]
            refl1 = (1 - prop) * temp1 + prop * temp2

            ind = ind3 + kk
            temp1 = refllut[ind, k]
            ind = ind4 + kk
            temp2 = refllut[ind, k]
            refl2 = (1 - prop) * temp1 + prop * temp2

            ind = ind5 + kk
            temp1 = refllut[ind, k]
            ind = ind6 + kk
            temp2 = refllut[ind, k]
            refl3 = (1 - prop) * temp1 + prop * temp2

            ind = ind7 + kk
            temp1 = refllut[ind, k]
            ind = ind8 + kk
            temp2 = refllut[ind, k]
            refl4 = (1 - prop) * temp1 + prop * temp2

            prop = (vza[indd] - vza_[vzamask[indd]]) / (vza_[vzamask[indd] + 1] - vza_[vzamask[indd]])
            refl5 = (1 - prop) * refl1 + prop * refl2
            refl6 = (1 - prop) * refl3 + prop * refl4

            prop = 1.0 * (psi[indd] - psi_[psimask[indd]]) / (psi_[psimask[indd] + 1] - psi_[psimask[indd]])
            refl7 = (1 - prop) * refl5 + prop * refl6
            coeff2[:, k] = refl7

        prop = 1.0 * (aod[indd] - aod_[aodmask[indd]]) / (aod_[aodmask[indd] + 1] - aod_[aodmask[indd]])
        coeff = np.zeros([np.sum(indd), 6])

        for kcoeff in range(6):
            coeff[:,kcoeff] =  prop* coeff1[:, kcoeff] + (1 - prop) * coeff2[:, kcoeff]

        # coeff = coeff2
        scale = 1.0
        reflnew = coeff[:, 0] * rad_red[indd]*scale - coeff[:, 1]  #
        reflnew = reflnew / (1.0 + coeff[:, 2] * reflnew)
        # reflnew = coeff[0,:]*1.0 - coeff[1,:]  #
        # reflnew = reflnew / (1.0 + coeff[2,:] * reflnew)
        reflnewred = np.zeros([nl, ns])
        # temp =  np.zeros([self.nl_nadir,self.ns_nadir])
        reflnewred[indd] = reflnew
        # temp[indd] = coeff[0,indd]

        reflnew = coeff[:, 3] * rad_nir[indd]*scale  - coeff[:, 4]
        reflnew = reflnew / (1.0 + coeff[:, 5] * reflnew)
        reflnewnir = np.zeros([nl, ns])
        reflnewnir[indd] = reflnew

        ndvi = np.zeros([nl, ns])
        if aodcheck == 3:
            ndvi[indd] = (refl_nir[indd] - refl_red[indd]) / (refl_nir[indd] + refl_red[indd])
        else:
            ndvi[indd] = (reflnewnir[indd] - reflnewred[indd]) / (reflnewnir[indd] + reflnewred[indd])
        bf = np.zeros([nl, ns])
        bf[indd] = (reflnewnir[indd] + reflnewred[indd]) / (2.0)
        ndvi[ndvi < 0] = 0
        ndvi[ndvi > 1] = 1
        if aodcheck == 1:
            outfile = wdir + resultName + dataTime + '_ndvi_n.tif'
        else:
            outfile = wdir + resultName + dataTime + '_ndvi_n_test.tif'
        self.writeTiff(ndvi, ns, nl, 1, geog, proj, outfile)
        # outfile = wdir+resultName+dataTime+'_red_n.tif'
        # self.writeTiff(reflnewred, ns, nl, 1, geog, proj, outfile)
        # outfile = wdir+resultName+dataTime+'_nir_n.tif'
        # self.writeTiff(reflnewnir, ns, nl, 1, geog, proj, outfile)

        return 1

    def VNIRatm_obliq(self, wdir, dataTime,medName,resultName,aodcheck = 2):

        # lookup table
        infile = wdir + r'auxiliary\SLSTR_Red_LUT.txt'
        temp = self.getDatafromTxt(infile, 1, 8)
        temp1 = temp[:, 5:8]
        infile = wdir + r'auxiliary\SLSTR_NIR_LUT.txt'
        temp = self.getDatafromTxt(infile, 1, 8)
        temp2 = temp[:, 5:8]
        refllut = np.hstack((temp1, temp2))

        refl_ = np.asarray([0])  # np.linspace(0,0.15,16)
        sza_ = np.asarray([0, 6, 12, 24, 36, 48, 54, 60, 66])
        vza_ = np.asarray([0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72])
        psi_ = np.asarray([0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180])
        aod_ = np.asarray([0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0])
        szanum = len(sza_)
        vzanum = len(vza_)
        psinum = len(psi_)
        reflnum = len(refl_)
        aodnum = len(aod_)

        infile = wdir + medName + dataTime + '_refl_red_o.tif'
        [refl_red, ns, nl, nb, geog, proj] = self.getTiffData(infile)

        infile = wdir + medName + dataTime + '_refl_nir_o.tif'
        [refl_nir, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # ns = np.int(ns*0.5)
        # nl = np.int(nl*0.5)
        # refl_red = self.resize(refl_red,nl,ns)
        infile = wdir + medName + dataTime + '_rad_red_o.tif'
        [rad_red, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # rad_red = self.resize(rad_red, nl, ns)
        infile = wdir + medName + dataTime + '_rad_nir_o.tif'
        [rad_nir, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # rad_nir = self.resize(rad_nir, nl, ns)

        infile = wdir + medName + dataTime + '_sza_o.tif'
        [sza, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        rd = np.pi / 180.0
        refl_red = refl_red / np.cos(sza * rd)


        ccc = np.sum(refl_red)
        if (ccc == 0):
            return 0

        infile = wdir + medName + dataTime + '_vza_o.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_psi_o.tif'
        [psi, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + '_aod_n.tif'

        c6 = os.path.exists(infile)
        if c6 ==0: return

        [aod, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        ####################################################################
        ##THIS IS WRONG BUT IT IS A COMPROMISING
        if aodcheck == 1:
            ind = (aod<=self.aodmax)*(aod>=self.aodmin)
            if (np.sum(ind)<0):
                aodvalue = self.aodvalue
            else:
                aodvalue = np.average(aod[ind])
            if aodvalue < self.aodmin: aodvalue = self.aodmin
            aod[aod != 0] = aodvalue
        if aodcheck ==2:
            aodvalue = self.aodvalue
            aod[aod != 0] = aodvalue

        # aodvalue = self.aodvalue
        # aod[aod != 0] = aodvalue
        ####################################################################
        szamask = np.zeros([nl, ns], dtype=np.int)
        vzamask = np.zeros([nl, ns], dtype=np.int)
        psimask = np.zeros([nl, ns], dtype=np.int)
        aodmask = np.zeros([nl, ns], dtype=np.int)

        for kaod in range(aodnum - 1):
            ind = (aod >= aod_[kaod]) * (aod < aod_[kaod + 1])
            if (np.sum(ind) > 0): aodmask[ind] = kaod
        ind = aod >= aod_[aodnum - 1]
        aodmask[ind] = aodnum - 2

        for ksza in range(szanum - 1):
            ind = (sza >= sza_[ksza]) * (sza < sza_[ksza + 1])
            if (np.sum(ind) > 0): szamask[ind] = ksza
        # this part is for pixels with sza > 48. It works but it is not very suitable
        ind = sza >= sza_[szanum - 1]
        szamask[ind] = szanum - 2

        for kvza in range(vzanum - 1):
            ind = (vza >= vza_[kvza]) * (vza < vza_[kvza + 1])
            if (np.sum(ind) > 0): vzamask[ind] = kvza

        for kpsi in range(psinum - 1):
            ind = (psi >= psi_[kpsi]) * (psi < psi_[kpsi + 1])
            if (np.sum(ind) > 0): psimask[ind] = kpsi

        num1 = vzanum * psinum * aodnum * reflnum
        num2 = psinum * aodnum * reflnum
        num3 = aodnum * reflnum

        # valid pixels
        indd = (vzamask != -1) * (refl_red > 0)
        if (np.sum(indd) <= 0):
            print('no suitable pixel can be used for vod inversion')
            return 0

        # sza1, vza1, psi1
        ind1 = szamask * num1 + vzamask * num2 + psimask * num3 + aodmask + 1
        # sza2, vza1, psi1
        ind2 = (szamask + 1) * num1 + vzamask * num2 + psimask * num3 + aodmask + 1
        # sza1, vza2, psi1
        ind3 = szamask * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask + 1
        # sza2, vza2, psi1
        ind4 = (szamask + 1) * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask + 1
        # sza1, vza1, psi2
        ind5 = szamask * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask + 1
        # sza2, vza1, psi2
        ind6 = (szamask + 1) * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask + 1
        # sza1, vza2, psi2
        ind7 = szamask * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask + 1
        # sza2, vza2, psi2
        ind8 = (szamask + 1) * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask + 1

        ind1 = ind1[indd]
        ind2 = ind2[indd]
        ind3 = ind3[indd]
        ind4 = ind4[indd]
        ind5 = ind5[indd]
        ind6 = ind6[indd]
        ind7 = ind7[indd]
        ind8 = ind8[indd]

        coeff1 = np.zeros([np.sum(indd), 6])
        for k in range(6):
            prop = (sza[indd] - sza_[szamask[indd]]) / (sza_[szamask[indd] + 1] - sza_[szamask[indd]])
            kk = 0
            ind = ind1 + kk
            temp1 = refllut[ind, k]
            ind = ind2 + kk
            temp2 = refllut[ind, k]
            refl1 = (1 - prop) * temp1 + prop * temp2

            ind = ind3 + kk
            temp1 = refllut[ind, k]
            ind = ind4 + kk
            temp2 = refllut[ind, k]
            refl2 = (1 - prop) * temp1 + prop * temp2

            ind = ind5 + kk
            temp1 = refllut[ind, k]
            ind = ind6 + kk
            temp2 = refllut[ind, k]
            refl3 = (1 - prop) * temp1 + prop * temp2

            ind = ind7 + kk
            temp1 = refllut[ind, k]
            ind = ind8 + kk
            temp2 = refllut[ind, k]
            refl4 = (1 - prop) * temp1 + prop * temp2

            prop = (vza[indd] - vza_[vzamask[indd]]) / (vza_[vzamask[indd] + 1] - vza_[vzamask[indd]])
            refl5 = (1 - prop) * refl1 + prop * refl2
            refl6 = (1 - prop) * refl3 + prop * refl4

            prop = 1.0 * (psi[indd] - psi_[psimask[indd]]) / (psi_[psimask[indd] + 1] - psi_[psimask[indd]])
            refl7 = (1 - prop) * refl5 + prop * refl6
            coeff1[:, k] = refl7

        #  sza1, vza1, psi1
        ind1 = szamask * num1 + vzamask * num2 + psimask * num3 + aodmask
        #  sza2, vza1, psi1
        ind2 = (szamask + 1) * num1 + vzamask * num2 + psimask * num3 + aodmask
        #  sza1, vza2, psi1
        ind3 = szamask * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask
        #  sza2, vza2, psi1
        ind4 = (szamask + 1) * num1 + (vzamask + 1) * num2 + psimask * num3 + aodmask
        # sza1, vza1, psi2
        ind5 = szamask * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask
        #  sza2, vza1, psi2
        ind6 = (szamask + 1) * num1 + (vzamask) * num2 + (psimask + 1) * num3 + aodmask
        #  sza1, vza2, psi2
        ind7 = szamask * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask
        #  sza2, vza2, psi2
        ind8 = (szamask + 1) * num1 + (vzamask + 1) * num2 + (psimask + 1) * num3 + aodmask

        ind1 = ind1[indd]
        ind2 = ind2[indd]
        ind3 = ind3[indd]
        ind4 = ind4[indd]
        ind5 = ind5[indd]
        ind6 = ind6[indd]
        ind7 = ind7[indd]
        ind8 = ind8[indd]

        coeff2 = np.zeros([np.sum(indd), 6])
        for k in range(6):
            prop = (sza[indd] - sza_[szamask[indd]]) / (sza_[szamask[indd] + 1] - sza_[szamask[indd]])
            kk = 0
            ind = ind1 + kk
            temp1 = refllut[ind, k]
            ind = ind2 + kk
            temp2 = refllut[ind, k]
            refl1 = (1 - prop) * temp1 + prop * temp2

            ind = ind3 + kk
            temp1 = refllut[ind, k]
            ind = ind4 + kk
            temp2 = refllut[ind, k]
            refl2 = (1 - prop) * temp1 + prop * temp2

            ind = ind5 + kk
            temp1 = refllut[ind, k]
            ind = ind6 + kk
            temp2 = refllut[ind, k]
            refl3 = (1 - prop) * temp1 + prop * temp2

            ind = ind7 + kk
            temp1 = refllut[ind, k]
            ind = ind8 + kk
            temp2 = refllut[ind, k]
            refl4 = (1 - prop) * temp1 + prop * temp2

            prop = (vza[indd] - vza_[vzamask[indd]]) / (vza_[vzamask[indd] + 1] - vza_[vzamask[indd]])
            refl5 = (1 - prop) * refl1 + prop * refl2
            refl6 = (1 - prop) * refl3 + prop * refl4

            prop = 1.0 * (psi[indd] - psi_[psimask[indd]]) / (psi_[psimask[indd] + 1] - psi_[psimask[indd]])
            refl7 = (1 - prop) * refl5 + prop * refl6
            coeff2[:, k] = refl7

        prop = 1.0 * (aod[indd] - aod_[aodmask[indd]]) / (aod_[aodmask[indd] + 1] - aod_[aodmask[indd]])


        coeff = np.zeros([np.sum(indd), 6])

        for kcoeff in range(6):
            coeff[:,kcoeff] =  prop* coeff1[:, kcoeff] + (1 - prop) * coeff2[:, kcoeff]


        reflnew = coeff[:, 0] * rad_red[indd] - coeff[:, 1]
        reflnew = reflnew / (1.0 + coeff[:, 2] * reflnew)
        # reflnew = coeff[0,:]*1.0 - coeff[1,:]  #
        # reflnew = reflnew / (1.0 + coeff[2,:] * reflnew)
        reflnewred = np.zeros([nl, ns])
        # temp =  np.zeros([self.nl_nadir,self.ns_nadir])
        reflnewred[indd] = reflnew
        # temp[indd] = coeff[0,indd]

        reflnew = coeff[:, 3] * rad_nir[indd] - coeff[:, 4]
        reflnew = reflnew / (1.0 + coeff[:, 5] * reflnew)
        reflnewnir = np.zeros([nl, ns])
        reflnewnir[indd] = reflnew

        ndvi = np.zeros([nl, ns])
        if aodcheck==3:
            ndvi[indd] = (refl_nir[indd] - refl_red[indd]) / (refl_nir[indd] + refl_red[indd])
        else:
            ndvi[indd] = (reflnewnir[indd] - reflnewred[indd]) / (reflnewnir[indd] + reflnewred[indd])
        bf = np.zeros([nl, ns])
        bf[indd] = (reflnewnir[indd] + reflnewred[indd]) / (2.0)
        ndvi[ndvi < 0] = 0
        ndvi[ndvi > 1] = 1

        if aodcheck == 1:
            outfile = wdir + resultName + dataTime + '_ndvi_o.tif'
        else:
            outfile = wdir + resultName + dataTime + '_ndvi_o_test.tif'
        self.writeTiff(ndvi, ns, nl, 1, geog, proj, outfile)
        # outfile = wdir+resultName+dataTime+'_red_o.tif'
        # self.writeTiff(reflnewred, ns, nl, 1, geog, proj, outfile)
        # outfile = wdir+resultName+dataTime+'_nir_o.tif'
        # self.writeTiff(reflnewnir, ns, nl, 1, geog, proj, outfile)

        return 1

    ### LST Inversion

    def setEmis_nadir(self, wdir, dataTime,medName,resultName,aodcheck=2):

        # get soil emissivity
        infile = wdir + r'auxiliary\heihe_emis1.tif'
        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        emis1_s = emis1_s / 1000.0
        infile = wdir + r'auxiliary\heihe_emis2.tif'
        [emis2_s, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        emis2_s = emis2_s / 1000.0
        infile = wdir + r'auxiliary\heihe_classif_new.tif'
        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        classif = np.asarray(classif, np.int)
        infile = wdir + r'auxiliary\ICBP_Emi.txt'
        emis_v = self.getDatafromTxt(infile, 0, 3)

        emis1_v = emis_v[classif, 0] / 1000.0
        emis2_v = emis_v[classif, 1] / 1000.0
        F = emis_v[classif,2]

        emis1_v = np.reshape(emis1_v, [nl, ns])
        emis2_v = np.reshape(emis2_v, [nl, ns])

        if aodcheck ==1:
            infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
        else:
            infile = wdir + resultName + dataTime + r'_ndvi_n_test.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc[fvc > 1] = 1
        fvc[fvc < 0] = 0

        emis1 = emis1_s * (1 - fvc) + emis1_v * fvc+(1-emis1_s)*emis1_v*F*(1-fvc)
        emis2 = emis2_s * (1 - fvc) + emis2_v * fvc+(1-emis2_s)*emis2_v*F*(1-fvc)
        outfile = wdir + resultName + dataTime + '_emis1_n.tif'
        self.writeTiff(emis1, ns, nl, 1, geog, proj, outfile)
        outfile = wdir + resultName + dataTime + '_emis2_n.tif'
        self.writeTiff(emis2, ns, nl, 1, geog, proj, outfile)

        return 0

    def setEmis_obliq(self, wdir, dataTime,medName,resultName,aodcheck=2):

        # get soil emissivity
        infile = wdir + r'auxiliary\heihe_emis1.tif'
        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        emis1_s = emis1_s / 1000.0
        infile = wdir + r'auxiliary\heihe_emis2.tif'
        [emis2_s, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        emis2_s = emis2_s / 1000.0
        infile = wdir + r'auxiliary\heihe_classif_new.tif'
        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        classif = np.asarray(classif, np.int)
        infile = wdir + r'auxiliary\ICBP_Emi.txt'
        emis_v = self.getDatafromTxt(infile, 0, 3)

        emis1_v = emis_v[classif, 0] / 1000.0
        emis2_v = emis_v[classif, 1] / 1000.0
        F = emis_v[classif,2]
        emis1_v = np.reshape(emis1_v, [nl, ns])
        emis2_v = np.reshape(emis2_v, [nl, ns])

        if aodcheck == 1:
            infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
        else:
            infile = wdir + resultName + dataTime + r'_ndvi_o_test.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc[fvc > 1] = 1
        fvc[fvc < 0] = 0

        emis1 = emis1_s * (1 - fvc) + emis1_v * fvc + (1 - emis1_s) * emis1_v * F * (1 - fvc)
        emis2 = emis2_s * (1 - fvc) + emis2_v * fvc + (1 - emis2_s) * emis2_v * F * (1 - fvc)

        outfile = wdir + resultName + dataTime + '_emis1_o.tif'
        self.writeTiff(emis1, ns, nl, 1, geog, proj, outfile)
        outfile = wdir + resultName + dataTime + '_emis2_o.tif'
        self.writeTiff(emis2, ns, nl, 1, geog, proj, outfile)

        return 0

    def inversionLST_nadir(self, wdir, dataTime,medName,resultName,aodcheck = 2):

        infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
        coeffs = self.getDatafromTxt(infile, 0, 15)
        coeffs = coeffs[:, 2:10]

        # infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
        # coeffs = self.getDatafromTxt(infile,0,15)
        # coeffs = coeffs[:,2:10]
        infile = wdir + medName + dataTime + '_bt8_n.tif'
        [BT8, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_bt9_n.tif'
        [BT9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_tcw_n.tif'
        [tcw, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_vza_n.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + '_emis1_n.tif'
        [emis8, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + '_emis2_n.tif'
        [emis9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # emis8 = np.zeros([nl,ns])
        # emis9 = np.zeros([nl,ns])
        # emis8[:] = 0.977
        # emis9[:] = 0.980


        tcw = tcw / 10.0

        emis = (emis8 + emis9) * 0.5
        demis = emis8 - emis9
        x3 = np.zeros([nl, ns])
        x4 = np.zeros([nl, ns])
        x6 = np.zeros([nl, ns])
        x7 = np.zeros([nl, ns])
        x1 = 1
        x2 = (BT8 + BT9) / 2.0
        ind = emis != 0
        x3[ind] = (1 - emis[ind]) / emis[ind] * x2[ind]
        x4[ind] = demis[ind] / (emis[ind] * emis[ind]) * x2[ind]
        x5 = (BT8 - BT9) / 2
        x6[ind] = (1 - emis[ind]) / emis[ind] * x5[ind]
        x7[ind] = demis[ind] / (emis[ind] * emis[ind]) * x5[ind]
        x8 = (BT8 - BT9) * (BT8 - BT9)

        wvc_ = np.asarray([0, 1.25, 2.25, 3.25, 4.25, 5.25, 100])
        vza_ = np.asarray([0, 3, 14.9, 38.6, 44.5, 51.2, 58, 65, 70, 75, 80])

        LST = np.zeros([nl, ns])
        for kwvc in range(len(wvc_) - 1):
            for kvza in range(len(vza_) - 1):

                ind = (tcw >= wvc_[kwvc]) * (tcw < wvc_[kwvc + 1]) * (vza >= vza_[kvza]) * (vza < vza_[kvza + 1]) * (
                            BT8 > 0) * (emis8 > 0)
                if (np.sum(ind) <= 0): continue

                ##########################################
                ### the normal result
                ##########################################

                ind1 = kvza * 6 + kwvc
                ind2 = (kvza + 1) * 6 + kwvc
                coeff1 = coeffs[ind1, :]
                coeff2 = coeffs[ind2, :]
                vzatemp = vza[ind]
                prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                prop1 = 1 - prop2
                LSTtemp = x1 * (prop1 * coeff1[0] + coeff2[0] * prop2) + x2[ind] * (
                            prop1 * coeff1[1] + coeff2[1] * prop2) + \
                          x3[ind] * (prop1 * coeff1[2] + coeff2[2] * prop2) + x4[ind] * (
                                      prop1 * coeff1[3] + coeff2[3] * prop2) + \
                          x5[ind] * (prop1 * coeff1[4] + coeff2[4] * prop2) + x6[ind] * (
                                      prop1 * coeff1[5] + coeff2[5] * prop2) + \
                          x7[ind] * (prop1 * coeff1[6] + coeff2[6] * prop2) + x8[ind] * (
                                      prop1 * coeff1[7] + coeff2[7] * prop2)

                ##########################################
                ### intergrating with atmospheric moisture
                ##########################################
                # ind1 = kvza * 6 + kwvc
                # ind2 = (kvza + 1) * 6 + kwvc
                # coeff1 = coeffs[ind1, :]
                # coeff2 = coeffs[ind2, :]
                # vzatemp = vza[ind]
                # prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                # prop1 = 1 - prop2
                # coeff12 = np.asarray([coeff1[k] * prop1 + coeff2[k] * prop2 for k in range(8)])
                #
                # ind3 = kvza * 6 + kwvc + 1
                # ind4 = (kvza + 1) * 6 + kwvc + 1
                # coeff3 = coeffs[ind3, :]
                # coeff4 = coeffs[ind4, :]
                # prop4 = prop2
                # prop3 = prop1
                # coeff34 = np.asarray([coeff3[k] * prop3 + coeff4[k] * prop4 for k in range(8)])
                #
                # wvctemp = tcw[ind]
                # prop6 = (wvctemp - wvc_[kwvc]) / (wvc_[kwvc + 1] - wvc_[kwvc])
                # prop5 = 1 - prop6
                # coeff = np.asarray([coeff34[k, :] * prop6 + prop5 * coeff12[k, :] for k in range(8)])
                # LSTtemp = x1 * coeff[0, :] + x2[ind] * coeff[1, :] + \
                #           x3[ind] * coeff[2, :] + x4[ind] * coeff[3, :] + \
                #           x5[ind] * coeff[4, :] + x6[ind] * coeff[5, :] + \
                #           x7[ind] * coeff[6, :] + x8[ind] * coeff[7, :]

                LST[ind] = LSTtemp

        ind = LST > 335
        LST[ind] = 335
        ind = LST < 275
        LST[ind] = 0
        if aodcheck ==1:
            outfile = wdir + resultName + dataTime + '_LST_n.tif'
        else:
            outfile = wdir + resultName + dataTime + '_LST_n_test.tif'
        self.writeTiff(LST, ns, nl, 1, geog, proj, outfile)

        return 0

    def inversionLST_obliq(self, wdir, dataTime,medName,resultName,aodcheck = 2):

        infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
        coeffs = self.getDatafromTxt(infile, 0, 15)
        coeffs = coeffs[:, 2:10]

        infile = wdir + medName + dataTime + '_bt8_o.tif'
        [BT8, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_bt9_o.tif'
        [BT9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_tcw_n.tif'
        [tcw, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_vza_o.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + '_emis1_o.tif'
        [emis8, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + '_emis2_o.tif'
        [emis9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)

        # emis8 = np.zeros([nl,ns])
        # emis9 = np.zeros([nl,ns])
        # emis8[:] = 0.977
        # emis9[:] = 0.980

        tcw = tcw / 10.0

        emis = (emis8 + emis9) * 0.5
        demis = emis8 - emis9
        x3 = np.zeros([nl, ns])
        x4 = np.zeros([nl, ns])
        x6 = np.zeros([nl, ns])
        x7 = np.zeros([nl, ns])
        x1 = 1
        x2 = (BT8 + BT9) / 2.0
        ind = emis != 0
        x3[ind] = (1 - emis[ind]) / emis[ind] * x2[ind]
        x4[ind] = demis[ind] / (emis[ind] * emis[ind]) * x2[ind]
        x5 = (BT8 - BT9) / 2
        x6[ind] = (1 - emis[ind]) / emis[ind] * x5[ind]
        x7[ind] = demis[ind] / (emis[ind] * emis[ind]) * x5[ind]
        x8 = (BT8 - BT9) * (BT8 - BT9)

        wvc_ = np.asarray([0, 1.25, 2.25, 3.25, 4.25, 5.25, 100])
        vza_ = np.asarray([0, 3, 14.9, 38.6, 44.5, 51.2, 58, 65, 70, 75, 80])

        LST = np.zeros([nl, ns])
        for kwvc in range(len(wvc_) - 1):
            for kvza in range(len(vza_) - 1):

                ind = (tcw >= wvc_[kwvc]) * (tcw < wvc_[kwvc + 1]) * (vza >= vza_[kvza]) * (vza < vza_[kvza + 1]) * (
                            BT8 > 0) * (emis8 > 0)
                if (np.sum(ind) <= 0): continue


                ind1 = kvza * 6 + kwvc
                ind2 = (kvza + 1) * 6 + kwvc
                coeff1 = coeffs[ind1, :]
                coeff2 = coeffs[ind2, :]
                vzatemp = vza[ind]
                prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                prop1 = 1 - prop2
                LSTtemp = x1 * (prop1 * coeff1[0] + coeff2[0] * prop2) + x2[ind] * (
                            prop1 * coeff1[1] + coeff2[1] * prop2) + \
                          x3[ind] * (prop1 * coeff1[2] + coeff2[2] * prop2) + x4[ind] * (
                                      prop1 * coeff1[3] + coeff2[3] * prop2) + \
                          x5[ind] * (prop1 * coeff1[4] + coeff2[4] * prop2) + x6[ind] * (
                                      prop1 * coeff1[5] + coeff2[5] * prop2) + \
                          x7[ind] * (prop1 * coeff1[6] + coeff2[6] * prop2) + x8[ind] * (
                                      prop1 * coeff1[7] + coeff2[7] * prop2)

                # ### water's ...
                # ind1 = kvza * 6 + kwvc
                # ind2 = (kvza + 1) * 6 + kwvc
                # coeff1 = coeffs[ind1, :]
                # coeff2 = coeffs[ind2, :]
                # vzatemp = vza[ind]
                # prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                # prop1 = 1 - prop2
                # coeff12 = np.asarray([coeff1[k] * prop1 + coeff2[k] * prop2 for k in range(8)])
                #
                # ### water+1's
                # ind3 = kvza * 6 + kwvc + 1
                # ind4 = (kvza + 1) * 6 + kwvc + 1
                # coeff3 = coeffs[ind3, :]
                # coeff4 = coeffs[ind4, :]
                # prop4 = prop2
                # prop3 = prop1
                # coeff34 = np.asarray([coeff3[k] * prop3 + coeff4[k] * prop4 for k in range(8)])
                #
                # wvctemp = tcw[ind]
                # prop6 = (wvctemp - wvc_[kwvc]) / (wvc_[kwvc + 1] - wvc_[kwvc])
                # prop5 = 1 - prop6
                # coeff = np.asarray([coeff34[k, :] * prop6 + prop5 * coeff12[k, :] for k in range(8)])
                # LSTtemp = x1 * coeff[0, :] + x2[ind] * coeff[1, :] + \
                #           x3[ind] * coeff[2, :] + x4[ind] * coeff[3, :] + \
                #           x5[ind] * coeff[4, :] + x6[ind] * coeff[5, :] + \
                #           x7[ind] * coeff[6, :] + x8[ind] * coeff[7, :]



                LST[ind] = LSTtemp

        ind = LST > 500
        LST[ind] = 0
        if aodcheck == 1:
            outfile = wdir + resultName + dataTime + '_LST_o.tif'
        else:
            outfile = wdir + resultName + dataTime + '_LST_o_test.tif'
        self.writeTiff(LST, ns, nl, 1, geog, proj, outfile)

        return 0

    def inversionLST_nadir_EVO(self, wdir, dataTime,medName,resultName,aodcheck = 1):

        infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
        coeffs = self.getDatafromTxt(infile, 0, 15)
        coeffs = coeffs[:, 2:10]

        # infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
        # coeffs = self.getDatafromTxt(infile,0,15)
        # coeffs = coeffs[:,2:10]
        infile = wdir + medName + dataTime + '_bt8_n.tif'
        [BT8, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_bt9_n.tif'
        [BT9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_tcw_n.tif'
        [tcw, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_vza_n.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        emis8 = np.zeros([nl,ns])
        emis9 = np.zeros([nl,ns])

        if aodcheck == 1:
            infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
        else:
            infile = wdir + resultName + dataTime + r'_ndvi_n_test.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc[fvc > 1] = 1
        fvc[fvc < 0] = 0

        emis8_s = 0.968
        emis9_s = 0.970
        emis8_v = 0.967
        emis9_v = 0.968
        F = 0.65

        emis8[:] = emis8_s*(1-fvc)+emis8_v*fvc+(1-emis8_s)*emis8_v*F*(1-fvc)
        emis9[:] = emis9_s*(1-fvc)+emis9_v*fvc+(1-emis9_s)*emis9_v*F*(1-fvc)


        tcw = tcw / 10.0

        emis = (emis8 + emis9) * 0.5
        demis = emis8 - emis9
        x3 = np.zeros([nl, ns])
        x4 = np.zeros([nl, ns])
        x6 = np.zeros([nl, ns])
        x7 = np.zeros([nl, ns])
        x1 = 1
        x2 = (BT8 + BT9) / 2.0
        ind = emis != 0
        x3[ind] = (1 - emis[ind]) / emis[ind] * x2[ind]
        x4[ind] = demis[ind] / (emis[ind] * emis[ind]) * x2[ind]
        x5 = (BT8 - BT9) / 2
        x6[ind] = (1 - emis[ind]) / emis[ind] * x5[ind]
        x7[ind] = demis[ind] / (emis[ind] * emis[ind]) * x5[ind]
        x8 = (BT8 - BT9) * (BT8 - BT9)

        wvc_ = np.asarray([0, 1.25, 2.25, 3.25, 4.25, 5.25, 100])
        vza_ = np.asarray([0, 3, 14.9, 38.6, 44.5, 51.2, 58, 65, 70, 75, 80])

        LST = np.zeros([nl, ns])
        for kwvc in range(len(wvc_) - 1):
            for kvza in range(len(vza_) - 1):

                ind = (tcw >= wvc_[kwvc]) * (tcw < wvc_[kwvc + 1]) * (vza >= vza_[kvza]) * (vza < vza_[kvza + 1]) * (
                            BT8 > 0) * (emis8 > 0)
                if (np.sum(ind) <= 0): continue

                ##########################################
                ### the normal result
                ##########################################

                ind1 = kvza * 6 + kwvc
                ind2 = (kvza + 1) * 6 + kwvc
                coeff1 = coeffs[ind1, :]
                coeff2 = coeffs[ind2, :]
                vzatemp = vza[ind]
                prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                prop1 = 1 - prop2
                LSTtemp = x1 * (prop1 * coeff1[0] + coeff2[0] * prop2) + x2[ind] * (
                            prop1 * coeff1[1] + coeff2[1] * prop2) + \
                          x3[ind] * (prop1 * coeff1[2] + coeff2[2] * prop2) + x4[ind] * (
                                      prop1 * coeff1[3] + coeff2[3] * prop2) + \
                          x5[ind] * (prop1 * coeff1[4] + coeff2[4] * prop2) + x6[ind] * (
                                      prop1 * coeff1[5] + coeff2[5] * prop2) + \
                          x7[ind] * (prop1 * coeff1[6] + coeff2[6] * prop2) + x8[ind] * (
                                      prop1 * coeff1[7] + coeff2[7] * prop2)

                ##########################################
                ### intergrating with atmospheric moisture
                ##########################################
                # ind1 = kvza * 6 + kwvc
                # ind2 = (kvza + 1) * 6 + kwvc
                # coeff1 = coeffs[ind1, :]
                # coeff2 = coeffs[ind2, :]
                # vzatemp = vza[ind]
                # prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                # prop1 = 1 - prop2
                # coeff12 = np.asarray([coeff1[k] * prop1 + coeff2[k] * prop2 for k in range(8)])
                #
                # ind3 = kvza * 6 + kwvc + 1
                # ind4 = (kvza + 1) * 6 + kwvc + 1
                # coeff3 = coeffs[ind3, :]
                # coeff4 = coeffs[ind4, :]
                # prop4 = prop2
                # prop3 = prop1
                # coeff34 = np.asarray([coeff3[k] * prop3 + coeff4[k] * prop4 for k in range(8)])
                #
                # wvctemp = tcw[ind]
                # prop6 = (wvctemp - wvc_[kwvc]) / (wvc_[kwvc + 1] - wvc_[kwvc])
                # prop5 = 1 - prop6
                # coeff = np.asarray([coeff34[k, :] * prop6 + prop5 * coeff12[k, :] for k in range(8)])
                # LSTtemp = x1 * coeff[0, :] + x2[ind] * coeff[1, :] + \
                #           x3[ind] * coeff[2, :] + x4[ind] * coeff[3, :] + \
                #           x5[ind] * coeff[4, :] + x6[ind] * coeff[5, :] + \
                #           x7[ind] * coeff[6, :] + x8[ind] * coeff[7, :]

                LST[ind] = LSTtemp

        ind = LST > 335
        LST[ind] = 335
        ind = LST < 275
        LST[ind] = 0
        if aodcheck ==1:
            outfile = wdir + resultName + dataTime + '_LST_n.tif'
        else:
            outfile = wdir + resultName + dataTime + '_LST_n_test.tif'
        self.writeTiff(LST, ns, nl, 1, geog, proj, outfile)

        return 0

    def inversionLST_obliq_EVO(self, wdir, dataTime,medName,resultName,aodcheck = 1):

        infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
        coeffs = self.getDatafromTxt(infile, 0, 15)
        coeffs = coeffs[:, 2:10]

        infile = wdir + medName + dataTime + '_bt8_o.tif'
        [BT8, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_bt9_o.tif'
        [BT9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_tcw_n.tif'
        [tcw, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + medName + dataTime + '_vza_o.tif'
        [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        emis8 = np.zeros([nl,ns])
        emis9 = np.zeros([nl,ns])

        if aodcheck == 1:
            infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
        else:
            infile = wdir + resultName + dataTime + r'_ndvi_o_test.tif'

        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc[fvc > 1] = 1
        fvc[fvc < 0] = 0

        emis8_s = 0.968
        emis9_s = 0.970
        emis8_v = 0.967
        emis9_v = 0.968

        F = 0.65

        emis8[:] = emis8_s*(1-fvc)+emis8_v*fvc+(1-emis8_s)*emis8_v*F*(1-fvc)
        emis9[:] = emis9_s*(1-fvc)+emis9_v*fvc+(1-emis9_s)*emis9_v*F*(1-fvc)


        tcw = tcw / 10.0

        emis = (emis8 + emis9) * 0.5
        demis = emis8 - emis9
        x3 = np.zeros([nl, ns])
        x4 = np.zeros([nl, ns])
        x6 = np.zeros([nl, ns])
        x7 = np.zeros([nl, ns])
        x1 = 1
        x2 = (BT8 + BT9) / 2.0
        ind = emis != 0
        x3[ind] = (1 - emis[ind]) / emis[ind] * x2[ind]
        x4[ind] = demis[ind] / (emis[ind] * emis[ind]) * x2[ind]
        x5 = (BT8 - BT9) / 2
        x6[ind] = (1 - emis[ind]) / emis[ind] * x5[ind]
        x7[ind] = demis[ind] / (emis[ind] * emis[ind]) * x5[ind]
        x8 = (BT8 - BT9) * (BT8 - BT9)

        wvc_ = np.asarray([0, 1.25, 2.25, 3.25, 4.25, 5.25, 100])
        vza_ = np.asarray([0, 3, 14.9, 38.6, 44.5, 51.2, 58, 65, 70, 75, 80])

        LST = np.zeros([nl, ns])
        for kwvc in range(len(wvc_) - 1):
            for kvza in range(len(vza_) - 1):

                ind = (tcw >= wvc_[kwvc]) * (tcw < wvc_[kwvc + 1]) * (vza >= vza_[kvza]) * (vza < vza_[kvza + 1]) * (
                            BT8 > 0) * (emis8 > 0)
                if (np.sum(ind) <= 0): continue


                ind1 = kvza * 6 + kwvc
                ind2 = (kvza + 1) * 6 + kwvc
                coeff1 = coeffs[ind1, :]
                coeff2 = coeffs[ind2, :]
                vzatemp = vza[ind]
                prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                prop1 = 1 - prop2
                LSTtemp = x1 * (prop1 * coeff1[0] + coeff2[0] * prop2) + x2[ind] * (
                            prop1 * coeff1[1] + coeff2[1] * prop2) + \
                          x3[ind] * (prop1 * coeff1[2] + coeff2[2] * prop2) + x4[ind] * (
                                      prop1 * coeff1[3] + coeff2[3] * prop2) + \
                          x5[ind] * (prop1 * coeff1[4] + coeff2[4] * prop2) + x6[ind] * (
                                      prop1 * coeff1[5] + coeff2[5] * prop2) + \
                          x7[ind] * (prop1 * coeff1[6] + coeff2[6] * prop2) + x8[ind] * (
                                      prop1 * coeff1[7] + coeff2[7] * prop2)

                # ### water's ...
                # ind1 = kvza * 6 + kwvc
                # ind2 = (kvza + 1) * 6 + kwvc
                # coeff1 = coeffs[ind1, :]
                # coeff2 = coeffs[ind2, :]
                # vzatemp = vza[ind]
                # prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
                # prop1 = 1 - prop2
                # coeff12 = np.asarray([coeff1[k] * prop1 + coeff2[k] * prop2 for k in range(8)])
                #
                # ### water+1's
                # ind3 = kvza * 6 + kwvc + 1
                # ind4 = (kvza + 1) * 6 + kwvc + 1
                # coeff3 = coeffs[ind3, :]
                # coeff4 = coeffs[ind4, :]
                # prop4 = prop2
                # prop3 = prop1
                # coeff34 = np.asarray([coeff3[k] * prop3 + coeff4[k] * prop4 for k in range(8)])
                #
                # wvctemp = tcw[ind]
                # prop6 = (wvctemp - wvc_[kwvc]) / (wvc_[kwvc + 1] - wvc_[kwvc])
                # prop5 = 1 - prop6
                # coeff = np.asarray([coeff34[k, :] * prop6 + prop5 * coeff12[k, :] for k in range(8)])
                # LSTtemp = x1 * coeff[0, :] + x2[ind] * coeff[1, :] + \
                #           x3[ind] * coeff[2, :] + x4[ind] * coeff[3, :] + \
                #           x5[ind] * coeff[4, :] + x6[ind] * coeff[5, :] + \
                #           x7[ind] * coeff[6, :] + x8[ind] * coeff[7, :]



                LST[ind] = LSTtemp

        ind = LST > 500
        LST[ind] = 0
        if aodcheck == 1:
            outfile = wdir + resultName + dataTime + '_LST_o.tif'
        else:
            outfile = wdir + resultName + dataTime + '_LST_o_test.tif'
        self.writeTiff(LST, ns, nl, 1, geog, proj, outfile)

        return 0

    # def inversionLST_nadir_night(self, wdir, dataTime,dataTimeold,medName,resultName):
    #
    #     infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
    #     coeffs = self.getDatafromTxt(infile, 0, 15)
    #     coeffs = coeffs[:, 2:10]
    #
    #     # infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
    #     # coeffs = self.getDatafromTxt(infile,0,15)
    #     # coeffs = coeffs[:,2:10]
    #     wdirout = r'heihe_night\heihe_'
    #     infile = wdir + wdirout + dataTime + '_bt8_n.tif'
    #     [BT8, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #     infile = wdir + wdirout + dataTime + '_bt9_n.tif'
    #     [BT9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + wdirout + dataTime + '_tcw_n.tif'
    #     [tcw, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + wdirout + dataTime + '_vza_n.tif'
    #     [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + resultName + dataTimeold + '_emis1_n.tif'
    #     if os.path.exists(infile) == 0: return 0
    #     [emis8, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + resultName + dataTimeold + '_emis2_n.tif'
    #     [emis9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #     tcw = tcw / 10.0
    #
    #     emis = (emis8 + emis9) * 0.5
    #     demis = emis8 - emis9
    #     x3 = np.zeros([nl, ns])
    #     x4 = np.zeros([nl, ns])
    #     x6 = np.zeros([nl, ns])
    #     x7 = np.zeros([nl, ns])
    #     x1 = 1
    #     x2 = (BT8 + BT9) / 2.0
    #     ind = emis != 0
    #     x3[ind] = (1 - emis[ind]) / emis[ind] * x2[ind]
    #     x4[ind] = demis[ind] / (emis[ind] * emis[ind]) * x2[ind]
    #     x5 = (BT8 - BT9) / 2
    #     x6[ind] = (1 - emis[ind]) / emis[ind] * x5[ind]
    #     x7[ind] = demis[ind] / (emis[ind] * emis[ind]) * x5[ind]
    #     x8 = (BT8 - BT9) * (BT8 - BT9)
    #
    #     wvc_ = np.asarray([0, 1.25, 2.25, 3.25, 4.25, 5.25, 100])
    #     vza_ = np.asarray([0, 3, 14.9, 38.6, 44.5, 51.2, 58, 65, 70, 75, 80])
    #
    #     LST = np.zeros([nl, ns])
    #     for kwvc in range(len(wvc_) - 1):
    #         for kvza in range(len(vza_) - 1):
    #
    #             ind = (tcw >= wvc_[kwvc]) * (tcw < wvc_[kwvc + 1]) * (vza >= vza_[kvza]) * (vza < vza_[kvza + 1]) * (
    #                         BT8 > 0) * (emis8 > 0)
    #             if (np.sum(ind) <= 0): continue
    #
    #             # ind1 = kvza * 6 + kwvc
    #             # ind2 = (kvza + 1) * 6 + kwvc
    #             # coeff1 = coeffs[ind1, :]
    #             # coeff2 = coeffs[ind2, :]
    #             # vzatemp = vza[ind]
    #             # prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
    #             # prop1 = 1 - prop2
    #             # LSTtemp = x1 * (prop1 * coeff1[0] + coeff2[0] * prop2) + x2[ind] * (
    #             #             prop1 * coeff1[1] + coeff2[1] * prop2) + \
    #             #           x3[ind] * (prop1 * coeff1[2] + coeff2[2] * prop2) + x4[ind] * (
    #             #                       prop1 * coeff1[3] + coeff2[3] * prop2) + \
    #             #           x5[ind] * (prop1 * coeff1[4] + coeff2[4] * prop2) + x6[ind] * (
    #             #                       prop1 * coeff1[5] + coeff2[5] * prop2) + \
    #             #           x7[ind] * (prop1 * coeff1[6] + coeff2[6] * prop2) + x8[ind] * (
    #             #                       prop1 * coeff1[7] + coeff2[7] * prop2)
    #
    #             ind1 = kvza * 6 + kwvc
    #             ind2 = (kvza + 1) * 6 + kwvc
    #             coeff1 = coeffs[ind1, :]
    #             coeff2 = coeffs[ind2, :]
    #             vzatemp = vza[ind]
    #             prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
    #             prop1 = 1 - prop2
    #             coeff12 = np.asarray([coeff1[k] * prop1 + coeff2[k] * prop2 for k in range(8)])
    #
    #             ind3 = kvza * 6 + kwvc + 1
    #             ind4 = (kvza + 1) * 6 + kwvc + 1
    #             coeff3 = coeffs[ind3, :]
    #             coeff4 = coeffs[ind4, :]
    #             prop4 = prop2
    #             prop3 = prop1
    #             coeff34 = np.asarray([coeff3[k] * prop3 + coeff4[k] * prop4 for k in range(8)])
    #
    #             wvctemp = tcw[ind]
    #             prop6 = (wvctemp - wvc_[kwvc]) / (wvc_[kwvc + 1] - wvc_[kwvc])
    #             prop5 = 1 - prop6
    #             coeff = np.asarray([coeff34[k, :] * prop6 + prop5 * coeff12[k, :] for k in range(8)])
    #             LSTtemp = x1 * coeff[0, :] + x2[ind] * coeff[1, :] + \
    #                       x3[ind] * coeff[2, :] + x4[ind] * coeff[3, :] + \
    #                       x5[ind] * coeff[4, :] + x6[ind] * coeff[5, :] + \
    #                       x7[ind] * coeff[6, :] + x8[ind] * coeff[7, :]
    #
    #             LST[ind] = LSTtemp
    #
    #     outfile = wdir + r'heihe_result_night\heihe_' + dataTime + '_LST_n.tif'
    #     self.writeTiff(LST, ns, nl, 1, geog, proj, outfile)
    #
    #     return 0
    #
    # def inversionLST_obliq_night(self, wdir, dataTime,dataTimeold,medName,resultName):
    #
    #     infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
    #     coeffs = self.getDatafromTxt(infile, 0, 15)
    #     coeffs = coeffs[:, 2:10]
    #
    #     # infile = wdir + r'\auxiliary\SLSTR_SW_Coeff_Day.txt'
    #     # coeffs = self.getDatafromTxt(infile,0,15)
    #     # coeffs = coeffs[:,2:10]
    #     wdirout = r'heihe_night\heihe_'
    #     infile = wdir + wdirout + dataTime + '_bt8_o.tif'
    #     [BT8, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #     infile = wdir + wdirout + dataTime + '_bt9_o.tif'
    #     [BT9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + wdirout + dataTime + '_tcw_n.tif'
    #     [tcw, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + wdirout + dataTime + '_vza_o.tif'
    #     [vza, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + resultName + dataTimeold + '_emis1_o.tif'
    #     if os.path.exists(infile) == 0:return 0
    #     [emis8, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #     infile = wdir + resultName + dataTimeold + '_emis2_o.tif'
    #     [emis9, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #     tcw = tcw / 10.0
    #
    #     emis = (emis8 + emis9) * 0.5
    #     demis = emis8 - emis9
    #     x3 = np.zeros([nl, ns])
    #     x4 = np.zeros([nl, ns])
    #     x6 = np.zeros([nl, ns])
    #     x7 = np.zeros([nl, ns])
    #     x1 = 1
    #     x2 = (BT8 + BT9) / 2.0
    #     ind = emis != 0
    #     x3[ind] = (1 - emis[ind]) / emis[ind] * x2[ind]
    #     x4[ind] = demis[ind] / (emis[ind] * emis[ind]) * x2[ind]
    #     x5 = (BT8 - BT9) / 2
    #     x6[ind] = (1 - emis[ind]) / emis[ind] * x5[ind]
    #     x7[ind] = demis[ind] / (emis[ind] * emis[ind]) * x5[ind]
    #     x8 = (BT8 - BT9) * (BT8 - BT9)
    #
    #     wvc_ = np.asarray([0, 1.25, 2.25, 3.25, 4.25, 5.25, 100])
    #     vza_ = np.asarray([0, 3, 14.9, 38.6, 44.5, 51.2, 58, 65, 70, 75, 80])
    #
    #     LST = np.zeros([nl, ns])
    #     for kwvc in range(len(wvc_) - 1):
    #         for kvza in range(len(vza_) - 1):
    #
    #             ind = (tcw >= wvc_[kwvc]) * (tcw < wvc_[kwvc + 1]) * (vza >= vza_[kvza]) * (vza < vza_[kvza + 1]) * (
    #                         BT8 > 0) * (emis8 > 0)
    #             if (np.sum(ind) <= 0): continue
    #
    #             # ind1 = kvza * 6 + kwvc
    #             # ind2 = (kvza + 1) * 6 + kwvc
    #             # coeff1 = coeffs[ind1, :]
    #             # coeff2 = coeffs[ind2, :]
    #             # vzatemp = vza[ind]
    #             # prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
    #             # prop1 = 1 - prop2
    #             # LSTtemp = x1 * (prop1 * coeff1[0] + coeff2[0] * prop2) + x2[ind] * (
    #             #             prop1 * coeff1[1] + coeff2[1] * prop2) + \
    #             #           x3[ind] * (prop1 * coeff1[2] + coeff2[2] * prop2) + x4[ind] * (
    #             #                       prop1 * coeff1[3] + coeff2[3] * prop2) + \
    #             #           x5[ind] * (prop1 * coeff1[4] + coeff2[4] * prop2) + x6[ind] * (
    #             #                       prop1 * coeff1[5] + coeff2[5] * prop2) + \
    #             #           x7[ind] * (prop1 * coeff1[6] + coeff2[6] * prop2) + x8[ind] * (
    #             #                       prop1 * coeff1[7] + coeff2[7] * prop2)
    #
    #             ind1 = kvza * 6 + kwvc
    #             ind2 = (kvza + 1) * 6 + kwvc
    #             coeff1 = coeffs[ind1, :]
    #             coeff2 = coeffs[ind2, :]
    #             vzatemp = vza[ind]
    #             prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
    #             prop1 = 1 - prop2
    #             coeff12 = np.asarray([coeff1[k] * prop1 + coeff2[k] * prop2 for k in range(8)])
    #
    #             ind3 = kvza * 6 + kwvc + 1
    #             ind4 = (kvza + 1) * 6 + kwvc + 1
    #             coeff3 = coeffs[ind3, :]
    #             coeff4 = coeffs[ind4, :]
    #             prop4 = prop2
    #             prop3 = prop1
    #             coeff34 = np.asarray([coeff3[k] * prop3 + coeff4[k] * prop4 for k in range(8)])
    #
    #             wvctemp = tcw[ind]
    #             prop6 = (wvctemp - wvc_[kwvc]) / (wvc_[kwvc + 1] - wvc_[kwvc])
    #             prop5 = 1 - prop6
    #             coeff = np.asarray([coeff34[k, :] * prop6 + prop5 * coeff12[k, :] for k in range(8)])
    #             LSTtemp = x1 * coeff[0, :] + x2[ind] * coeff[1, :] + \
    #                       x3[ind] * coeff[2, :] + x4[ind] * coeff[3, :] + \
    #                       x5[ind] * coeff[4, :] + x6[ind] * coeff[5, :] + \
    #                       x7[ind] * coeff[6, :] + x8[ind] * coeff[7, :]
    #
    #             LST[ind] = LSTtemp
    #
    #     outfile = wdir + r'heihe_result_night\heihe_' + dataTime + '_LST_o.tif'
    #     self.writeTiff(LST, ns, nl, 1, geog, proj, outfile)
    #
    #     return 0




    #### LSCT inversion


    ### This is for evaluation


    # LSCT inversion
    def inversionLSCT_npl_agl_point_raw(self, wdir, dataTime, lat_, lon_, medName, resultName, off=2, off_=0, aodcheck=2):

        # get soil emissivity
        infile = wdir + r'auxiliary\heihe_emis1.tif'
        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        emis1_s = emis1_s / 1000.0

        # infile = wdir + r'auxiliary\heihe_emis2.tif'
        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # emis2_s = emis2_s / 1000.0

        infile = wdir + r'auxiliary\heihe_classif_new.tif'
        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        classif = np.asarray(classif, np.int)
        infile = wdir + r'auxiliary\ICBP_Emi.txt'
        emis_v = self.getDatafromTxt(infile, 0, 3)
        emis1_v = emis_v[classif, 0] / 1000.0
        # emis2_v = emis_v[classif, 1] / 1000.0
        F = emis_v[classif, 2]
        if aodcheck == 1:
            infile1 = wdir + resultName + dataTime + r'_lst_n.tif'
            infile2 = wdir + resultName + dataTime + r'_lst_o.tif'
        else:
            infile1 = wdir + resultName + dataTime + r'_lst_n_test.tif'
            infile2 = wdir + resultName + dataTime + r'_lst_o_test.tif'

        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile1)
        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile2)

        if aodcheck == 1:
            infile1 = wdir + resultName + dataTime + r'_ndvi_n.tif'
            infile2 = wdir + resultName + dataTime + r'_ndvi_o.tif'
        else:
            infile1 = wdir + resultName + dataTime + r'_ndvi_n_test.tif'
            infile2 = wdir + resultName + dataTime + r'_ndvi_o_test.tif'

        dataset = gdal.Open(infile1, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_n[fvc_n > 1] = 1
        fvc_n[fvc_n < 0] = 0

        dataset = gdal.Open(infile2, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_o[fvc_o > 1] = 1
        fvc_o[fvc_o < 0] = 0



        rad_n = self.planck(self.wl8, lst_n)
        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_n))
        rad_o = self.planck(self.wl8, lst_o)
        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_o))

        # component temperature inversion
        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n) + (1 - emis1_s) * emis1_v * F * (1 - fvc_n)])
        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o) + (1 - emis1_s) * emis1_v * F * (1 - fvc_o)])
        result = np.zeros([2, nl, ns])

        nsi = np.zeros([nl, ns])
        nli = np.zeros([nl, ns])
        for k in range(nl):
            nli[k, :] = k
            nsi[k, :] = np.linspace(0, ns - 1, ns)

        Cdvalue = 1.0
        Cdvaluep = 0.30
        puritypixe = 0.05
        # off_ = 0

        mask = np.zeros(np.shape(lst_n))
        ind = (fvc_n > 0.05) * (fvc_n < 0.95) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
                rad_o < 17)
        mask[ind] = 1
        part = 6
        ######
        ## pints
        ######
        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
        imagex = np.asarray(imagex, np.int)
        imagey = np.asarray(imagey, np.int)

        for k1 in range(off, nl - off):

            for k2 in range(off, ns - off):

                if (mask[k1, k2] != 1): continue
                temp = (imagex - k1) * (imagex - k1) + (imagey - k2) * (imagey - k2)
                ind = temp < 25
                if (np.sum(ind) < 1): continue

                nl1 = k1 - off
                nl2 = k1 + off + 1
                ns1 = k2 - off
                ns2 = k2 + off + 1

                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
                # fd[fd != 0] = 1 / fd[fd != 0]
                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
                # fd[halfpoint] = 5
                # fd[:] = 1
                fd = np.asarray(
                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
                fd = np.reshape(fd, -1) * 1.0

                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                w0 = W1[:, k1, k2]

                n0 = np.ones(len(ns_temp))
                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])

                b1 = (np.reshape(rad_n[k1, k2], -1))
                b2 = (np.reshape(rad_o[k1, k2], -1))
                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))

                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
                num_point = np.sum(ind)
                if (num_point < 12): continue

                ww = np.hstack((w11[:, ind], w22[:, ind]))
                bb = np.hstack((bb1[ind], bb2[ind]))
                #
                # ww = w11[:, ind]
                # bb = bb1[ind]

                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))

                ###############################################
                ### the least-square method
                ###############################################
                ww = np.transpose(ww)
                coeffprior = np.asarray(lstsq(ww, bb))[0]
                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ###############################################
                ### the least-square method
                ###############################################
                ### the multi-pixel method with a prior knowledge
                # bbb = sorted(b)
                # halfpoint = np.int(num_point / 2)
                # rleaf = np.average(bbb[0])
                # rsoil = np.average(bbb[1])
                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
                # bb = np.transpose(np.matrix(bb))
                # ww = np.transpose(ww)
                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
                # Cm = np.matrix(np.diag(np.ones(2 * part)))
                # Cm[0, 0] = 15.0
                # Cm[part, part] = 10.0
                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)

                #################################################
                ### result without averaging information
                #################################################
                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ################################################
                ### reult with averaing information
                ################################################
                coeffprior = np.asarray(coeffprior)
                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
                indnew = ind * (ct2 < 16.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 16.5)
                apoint = np.sum(indnew)
                if (apoint <= 1):
                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                    ctprior = np.transpose(np.matrix([cd1, cd2]))
                else:
                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
                    bb0 = bb1[offcenter]
                    fdnew = np.zeros(25)
                    # inver-distance weight
                    dismax = max(abs(bb1[indnew] - bb0))
                    if dismax != 0:
                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
                    else:
                        fdnew[:] = 1
                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])

                    # fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
                    cd1new = np.sum(ct1[indnew] * fd[indnew])
                    cd2new = np.sum(ct2[indnew] * fd[indnew])
                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))

                #################################################
                ### multi-angle
                #################################################

                nl1 = k1 - off_
                nl2 = k1 + off_ + 1
                ns1 = k2 - off_
                ns2 = k2 + off_ + 1
                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))
                w = np.transpose(w)
                # coeffprior = np.asarray(lstsq(w, b))[0]
                w = np.matrix(w)
                b = np.matrix(b)
                b = np.transpose(b)

                # if(b[0]<b[1] and conti1==1):continue

                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
                Cm = np.diag([2.0, 1.0])
                # ctprior = np.transpose(np.matrix(ctprior))
                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
                ct = np.asarray(coeff)

                if (np.min(ct) < 3): continue
                ct = self.invplanck0(self.wl8, ct)

                result[:, k1, k2] = [ct[0], ct[1]]

        result[result < 0] = 0
        result[result > 350] = 0

        if aodcheck == 1:
            outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw.tif'
        else:
            outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_test.tif'

        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)

        return 0

    def inversionLSCT_agl_point_raw(self, wdir, dataTime, lat_, lon_, medName, resultName, off=2, off_=0, aodcheck=1):

        # get soil emissivity
        infile = wdir + r'auxiliary\heihe_emis1.tif'
        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        emis1_s = emis1_s / 1000.0

        # infile = wdir + r'auxiliary\heihe_emis2.tif'
        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # emis2_s = emis2_s / 1000.0

        infile = wdir + r'auxiliary\heihe_classif_new.tif'
        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        classif = np.asarray(classif, np.int)
        infile = wdir + r'auxiliary\ICBP_Emi.txt'
        emis_v = self.getDatafromTxt(infile, 0, 3)
        emis1_v = emis_v[classif, 0] / 1000.0
        # emis2_v = emis_v[classif, 1] / 1000.0
        F = emis_v[classif, 2]

        infile = wdir + resultName + dataTime + r'_lst_n.tif'
        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + r'_lst_o.tif'
        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)

        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_n[fvc_n > 1] = 1
        fvc_n[fvc_n < 0] = 0
        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_o[fvc_o > 1] = 1
        fvc_o[fvc_o < 0] = 0

        rad_n = self.planck(self.wl8, lst_n)
        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_n))
        rad_o = self.planck(self.wl8, lst_o)
        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_o))

        # component temperature inversion
        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n) + (1 - emis1_s) * emis1_v * F * (1 - fvc_n)])
        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o) + (1 - emis1_s) * emis1_v * F * (1 - fvc_o)])
        result = np.zeros([2, nl, ns])

        nsi = np.zeros([nl, ns])
        nli = np.zeros([nl, ns])
        for k in range(nl):
            nli[k, :] = k
            nsi[k, :] = np.linspace(0, ns - 1, ns)

        Cdvalue = 0.5
        Cdvaluep = 0.30
        puritypixe = 0.05
        # off_ = 0

        mask = np.zeros(np.shape(lst_n))
        ind = (fvc_n > 0.05) * (fvc_n < 0.95) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
                rad_o < 17)
        mask[ind] = 1
        part = 6
        ######
        ## pints
        ######
        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
        imagex = np.asarray(imagex, np.int)
        imagey = np.asarray(imagey, np.int)

        for k1 in range(off, nl - off):

            for k2 in range(off, ns - off):

                if (mask[k1, k2] != 1): continue
                temp = (imagex - k1) * (imagex - k1) + (imagey - k2) * (imagey - k2)
                ind = temp < 25
                if (np.sum(ind) < 1): continue

                nl1 = k1 - off
                nl2 = k1 + off + 1
                ns1 = k2 - off
                ns2 = k2 + off + 1

                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
                # fd[fd != 0] = 1 / fd[fd != 0]
                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
                # fd[halfpoint] = 5
                # fd[:] = 1
                fd = np.asarray(
                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
                fd = np.reshape(fd, -1) * 1.0

                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                w0 = W1[:, k1, k2]

                n0 = np.ones(len(ns_temp))
                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])

                b1 = (np.reshape(rad_n[k1, k2], -1))
                b2 = (np.reshape(rad_o[k1, k2], -1))
                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))

                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
                num_point = np.sum(ind)
                if (num_point < 12): continue

                ww = np.hstack((w11[:, ind], w22[:, ind]))
                bb = np.hstack((bb1[ind], bb2[ind]))
                #
                # ww = w11[:, ind]
                # bb = bb1[ind]

                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))

                ###############################################
                ### the least-square method
                ###############################################
                ww = np.transpose(ww)
                coeffprior = np.asarray(lstsq(ww, bb))[0]
                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ###############################################
                ### the least-square method
                ###############################################
                ### the multi-pixel method with a prior knowledge
                # bbb = sorted(b)
                # halfpoint = np.int(num_point / 2)
                # rleaf = np.average(bbb[0])
                # rsoil = np.average(bbb[1])
                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
                # bb = np.transpose(np.matrix(bb))
                # ww = np.transpose(ww)
                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
                # Cm = np.matrix(np.diag(np.ones(2 * part)))
                # Cm[0, 0] = 15.0
                # Cm[part, part] = 10.0
                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)

                #################################################
                ### result without averaging information
                #################################################
                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ################################################
                ### reult with averaing information
                ################################################
                coeffprior = np.asarray(coeffprior)
                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
                apoint = np.sum(indnew)
                if (apoint <= 1):
                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                    ctprior = np.transpose(np.matrix([cd1, cd2]))
                else:
                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
                    bb0 = bb1[offcenter]
                    fdnew = np.zeros(25)
                    dismax = max(abs(bb1[indnew] - bb0))
                    if dismax != 0:
                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
                    else:
                        fdnew[:] = 1
                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])

                    # fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
                    cd1new = np.sum(ct1[indnew] * fd[indnew])
                    cd2new = np.sum(ct2[indnew] * fd[indnew])
                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))

                #################################################
                ### multi-angle
                #################################################

                # nl1 = k1 - off_
                # nl2 = k1 + off_ + 1
                # ns1 = k2 - off_
                # ns2 = k2 + off_ + 1
                # w1 = W1[:, nl1:nl2, ns1:ns2]
                # w1 = np.asarray(np.reshape(w1, (2, -1)))
                # w2 = W2[:, nl1:nl2, ns1:ns2]
                # w2 = np.asarray(np.reshape(w2, (2, -1)))
                # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
                # w = np.hstack((w1, w2))
                # b = np.hstack((b1, b2))
                # w = np.transpose(w)
                # # coeffprior = np.asarray(lstsq(w, b))[0]
                # w = np.matrix(w)
                # b = np.matrix(b)
                # b = np.transpose(b)
                #
                # # if(b[0]<b[1] and conti1==1):continue
                #
                # Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
                # Cm = np.diag([5.0, 2.0])
                # # ctprior = np.transpose(np.matrix(ctprior))
                # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
                # coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
                # # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
                # ct = np.asarray(coeff)

                nl1 = k1 - off_
                nl2 = k1 + off_ + 1
                ns1 = k2 - off_
                ns2 = k2 + off_ + 1
                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))
                w = np.transpose(w)
                ct = np.asarray(lstsq(w, b))[0]

                if (np.min(ct) < 3): continue
                ct = self.invplanck0(self.wl8, ct)

                result[:, k1, k2] = [ct[0], ct[1]]

        result[result < 0] = 0
        result[result > 350] = 0

        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_agl.tif'
        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)

        return 0

    def inversionLSCT_npl_point_raw(self, wdir, dataTime, lat_, lon_, medName, resultName, off=2, off_=0):      
        # get soil emissivity
        infile = wdir + r'auxiliary\heihe_emis1.tif'
        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        emis1_s = emis1_s / 1000.0
        
        # infile = wdir + r'auxiliary\heihe_emis2.tif'
        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # emis2_s = emis2_s / 1000.0
        
        infile = wdir + r'auxiliary\heihe_classif_new.tif'
        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        classif = np.asarray(classif, np.int)
        infile = wdir + r'auxiliary\ICBP_Emi.txt'
        emis_v = self.getDatafromTxt(infile, 0, 3)
        emis1_v = emis_v[classif, 0] / 1000.0
        F = emis_v[classif,2]
        # emis2_v = emis_v[classif, 1] / 1000.0
        
        infile = wdir + resultName + dataTime + r'_lst_n.tif'
        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + r'_lst_o.tif'
        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        
        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
        ccc = os.path.exists(infile)
        if ccc ==0: return

        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_n[fvc_n > 1] = 1
        fvc_n[fvc_n < 0] = 0
        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_o[fvc_o > 1] = 1
        fvc_o[fvc_o < 0] = 0

        rad_n = self.planck(self.wl8, lst_n)
        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_n))
        rad_o = self.planck(self.wl8, lst_o)
        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_o))

        # component temperature inversion
        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n) + (1 - emis1_s) * emis1_v * F * (1 - fvc_n)])
        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o) + (1 - emis1_s) * emis1_v * F * (1 - fvc_o)])
        result = np.zeros([2, nl, ns])
        
        nsi = np.zeros([nl, ns])
        nli = np.zeros([nl, ns])
        for k in range(nl):
           nli[k, :] = k
           nsi[k, :] = np.linspace(0, ns - 1, ns)
        
        Cdvalue = 0.5
        Cdvaluep = 0.30
        puritypixe = 0.05
        # off_ = 0
        
        mask = np.zeros(np.shape(lst_n))
        ind = (fvc_n > 0.05) * (fvc_n < 0.95) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
               rad_o < 17)
        mask[ind] = 1
        part = 6
        ######
        ## pints
        ######
        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
        imagex = np.asarray(imagex, np.int)
        imagey = np.asarray(imagey, np.int)
        
        for k1 in range(off, nl - off):
        
           for k2 in range(off, ns - off):
        
               if (mask[k1, k2] != 1): continue
               temp = (imagex - k1) * (imagex - k1) + (imagey - k2) * (imagey - k2)
               ind = temp < 25
               if (np.sum(ind) < 1): continue
        
               nl1 = k1 - off
               nl2 = k1 + off + 1
               ns1 = k2 - off
               ns2 = k2 + off + 1
        
               offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
               ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
               nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
               # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
               # fd[fd != 0] = 1 / fd[fd != 0]
               # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
               # fd[halfpoint] = 5
               # fd[:] = 1
               fd = np.asarray(
                   [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
               fd = np.reshape(fd, -1) * 1.0
        
               w1 = W1[:, nl1:nl2, ns1:ns2]
               w1 = np.asarray(np.reshape(w1, (2, -1)))
               w2 = W2[:, nl1:nl2, ns1:ns2]
               w2 = np.asarray(np.reshape(w2, (2, -1)))
               w0 = W1[:, k1, k2]
        
               n0 = np.ones(len(ns_temp))
               x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
               w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
               w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
        
               b1 = (np.reshape(rad_n[k1, k2], -1))
               b2 = (np.reshape(rad_o[k1, k2], -1))
               bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
        
               temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
               temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
               ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
               num_point = np.sum(ind)
               if (num_point < 12): continue
        
               ww = np.hstack((w11[:, ind], w22[:, ind]))
               bb = np.hstack((bb1[ind], bb2[ind]))
               #
               # ww = w11[:, ind]
               # bb = bb1[ind]
        
               w = np.hstack((w1, w2))
               b = np.hstack((b1, b2))
        
               ###############################################
               ### the least-square method
               ###############################################
               ww = np.transpose(ww)
               coeffprior = np.asarray(lstsq(ww, bb))[0]
               cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
               cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
               ct = np.transpose(np.matrix([cd1, cd2]))
        
               ###############################################
               ### the least-square method
               ###############################################
               ### the multi-pixel method with a prior knowledge
               # bbb = sorted(b)
               # halfpoint = np.int(num_point / 2)
               # rleaf = np.average(bbb[0])
               # rsoil = np.average(bbb[1])
               # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
               # bb = np.transpose(np.matrix(bb))
               # ww = np.transpose(ww)
               # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
               # Cm = np.matrix(np.diag(np.ones(2 * part)))
               # Cm[0, 0] = 15.0
               # Cm[part, part] = 10.0
               # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
               # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
        
               #################################################
               ### result without averaging information
               #################################################
               # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
               # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
               # ctprior = np.transpose(np.matrix([cd1, cd2]))
        
               ################################################
               ### reult with averaing information
               ################################################
               # coeffprior = np.asarray(coeffprior)
               # ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
               # ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
               # indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
               # apoint = np.sum(indnew)
               # if (apoint <= 1):
               #     cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
               #     cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
               #     ctprior = np.transpose(np.matrix([cd1, cd2]))
               # else:
               #     fd[indnew] = fd[indnew] / np.sum(fd[indnew])
               #     bb0 = bb1[offcenter]
               #     fdnew = np.zeros(25)
               #     dismax = max(abs(bb1[indnew] - bb0))
               #     if dismax != 0:
               #         fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
               #     else:
               #         fdnew[:] = 1
               #     fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
               #
               #     # fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
               #     cd1new = np.sum(ct1[indnew] * fd[indnew])
               #     cd2new = np.sum(ct2[indnew] * fd[indnew])
               #     ctprior = np.transpose(np.matrix([cd1new, cd2new]))
               #
               # #################################################
               # ### multi-angle
               # #################################################
               #
               # # nl1 = k1 - off_
               # # nl2 = k1 + off_ + 1
               # # ns1 = k2 - off_
               # # ns2 = k2 + off_ + 1
               # # w1 = W1[:, nl1:nl2, ns1:ns2]
               # # w1 = np.asarray(np.reshape(w1, (2, -1)))
               # # w2 = W2[:, nl1:nl2, ns1:ns2]
               # # w2 = np.asarray(np.reshape(w2, (2, -1)))
               # # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               # # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
               # # w = np.hstack((w1, w2))
               # # b = np.hstack((b1, b2))
               # # w = np.transpose(w)
               # # # coeffprior = np.asarray(lstsq(w, b))[0]
               # # w = np.matrix(w)
               # # b = np.matrix(b)
               # # b = np.transpose(b)
               # #
               # # # if(b[0]<b[1] and conti1==1):continue
               # #
               # # Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
               # # Cm = np.diag([5.0, 2.0])
               # # # ctprior = np.transpose(np.matrix(ctprior))
               # # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
               # # coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
               # # # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
               # # ct = np.asarray(coeff)
               #
               # nl1 = k1 - off_
               # nl2 = k1 + off_ + 1
               # ns1 = k2 - off_
               # ns2 = k2 + off_ + 1
               # w1 = W1[:, nl1:nl2, ns1:ns2]
               # w1 = np.asarray(np.reshape(w1, (2, -1)))
               # w2 = W2[:, nl1:nl2, ns1:ns2]
               # w2 = np.asarray(np.reshape(w2, (2, -1)))
               # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
               # w = np.hstack((w1, w2))
               # b = np.hstack((b1, b2))
               # w = np.transpose(w)
               # ct = np.asarray(lstsq(w, b))[0]
        
               if (np.min(ct) < 3): continue
               ct = self.invplanck0(self.wl8, ct)
        
               result[:, k1, k2] = [ct[0], ct[1]]
        
        result[result < 0] = 0
        result[result > 350] = 0
        
        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_npl.tif'
        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
        
        return 0

    def inversionLSCT_npl_agl_raw_EVO(self, wdir, dataTime, medName, resultName, off=2, off_=0, aodcheck=2,singlecheck = 1):

        # get soil emissivity
        # infile = wdir + r'auxiliary\heihe_emis1.tif'
        # [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # emis1_s = emis1_s / 1000.0

        # infile = wdir + r'auxiliary\heihe_emis2.tif'
        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # emis2_s = emis2_s / 1000.0

        # infile = wdir + r'auxiliary\heihe_classif_new.tif'
        # [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # classif = np.asarray(classif, np.int)
        # infile = wdir + r'auxiliary\ICBP_Emi.txt'
        # emis_v = self.getDatafromTxt(infile, 0, 2)
        # emis1_v = emis_v[classif, 0] / 1000.0
        # emis2_v = emis_v[classif, 1] / 1000.0

        if aodcheck == 1:
            infile1 = wdir + resultName + dataTime + r'_lst_n.tif'
            infile2 = wdir + resultName + dataTime + r'_lst_o.tif'
        else:
            infile1 = wdir + resultName + dataTime + r'_lst_n_test.tif'
            infile2 = wdir + resultName + dataTime + r'_lst_o_test.tif'

        # infile = wdir + resultName + dataTime + r'_lst_n.tif'
        [lst_n, ns, nl, nb, geog, proj] = self.getTiffData(infile1)
        # infile = wdir + resultName + dataTime + r'_lst_o.tif'
        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile2)

        emis1_v = np.zeros([nl, ns])
        emis1_s = np.zeros([nl, ns])
        emis1_v[:] = 0.967
        emis1_s[:] = 0.968
        F = 0.65

        if aodcheck == 1:
            infile1 = wdir + resultName + dataTime + r'_ndvi_n.tif'
            infile2 = wdir + resultName + dataTime + r'_ndvi_o.tif'
        else:
            infile1 = wdir + resultName + dataTime + r'_ndvi_n_test.tif'
            infile2 = wdir + resultName + dataTime + r'_ndvi_o_test.tif'

        # infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
        dataset = gdal.Open(infile1, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_n[fvc_n > 1] = 1
        fvc_n[fvc_n < 0] = 0
        # infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
        dataset = gdal.Open(infile2, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_o[fvc_o > 1] = 1
        fvc_o[fvc_o < 0] = 0

        rad_n = self.planck(self.wl8, lst_n)
        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_n))
        rad_o = self.planck(self.wl8, lst_o)
        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_o))

        # component temperature inversion
        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n) + (1 - emis1_s) * emis1_v * F * (1 - fvc_n)])
        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o) + (1 - emis1_s) * emis1_v * F * (1 - fvc_o)])
        result = np.zeros([2, nl, ns])

        nsi = np.zeros([nl, ns])
        nli = np.zeros([nl, ns])
        for k in range(nl):
            nli[k, :] = k
            nsi[k, :] = np.linspace(0, ns - 1, ns)

        Cdvalue = 1.0
        Cdvaluep = 0.30
        puritypixe = 0.05
        # off_ = 0

        mask = np.zeros(np.shape(lst_n))
        ind = (fvc_n > 0.05) * (fvc_n < 0.95) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
                rad_o < 17)
        mask[ind] = 1
        part = 6
        for k1 in range(off, nl - off):

            for k2 in range(off, ns - off):

                if (mask[k1, k2] != 1): continue

                nl1 = k1 - off
                nl2 = k1 + off + 1
                ns1 = k2 - off
                ns2 = k2 + off + 1

                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
                # fd[fd != 0] = 1 / fd[fd != 0]
                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
                # fd[halfpoint] = 5
                # fd[:] = 1
                fd = np.asarray(
                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
                fd = np.reshape(fd, -1) * 1.0

                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                w0 = W1[:, k1, k2]

                n0 = np.ones(len(ns_temp))
                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])

                b1 = (np.reshape(rad_n[k1, k2], -1))
                b2 = (np.reshape(rad_o[k1, k2], -1))
                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))

                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
                num_point = np.sum(ind)
                if (num_point < 6): continue

                if singlecheck == 0:
                    ww = np.hstack((w11[:, ind], w22[:, ind]))
                    bb = np.hstack((bb1[ind], bb2[ind]))
                else:
                    ww = np.hstack((w11[:, ind]))
                    bb = np.hstack((bb1[ind]))

                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))

                ###############################################
                ### the least-square method
                ###############################################
                ww = np.transpose(ww)
                coeffprior = np.asarray(lstsq(ww, bb))[0]
                # cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ###############################################
                ### the least-square method
                ###############################################
                ### the multi-pixel method with a prior knowledge
                # bbb = sorted(b)
                # halfpoint = np.int(num_point / 2)
                # rleaf = np.average(bbb[0])
                # rsoil = np.average(bbb[1])
                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
                # bb = np.transpose(np.matrix(bb))
                # ww = np.transpose(ww)
                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
                # Cm = np.matrix(np.diag(np.ones(2 * part)))
                # Cm[0, 0] = 15.0
                # Cm[part, part] = 10.0
                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)

                #################################################
                ### result without averaging information
                #################################################
                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ################################################
                ### reult with averaing information
                ################################################
                coeffprior = np.asarray(coeffprior)
                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
                apoint = np.sum(indnew)
                if (apoint <= 1):
                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                    ctprior = np.transpose(np.matrix([cd1, cd2]))
                else:
                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
                    bb0 = bb1[offcenter]
                    fdnew = np.zeros(25)
                    dismax = max(abs(bb1[indnew] - bb0))
                    if dismax != 0:
                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
                    else:
                        fdnew[:] = 1
                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])

                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
                    cd1new = np.sum(ct1[indnew] * fd[indnew])
                    cd2new = np.sum(ct2[indnew] * fd[indnew])
                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))

                #################################################
                ### multi-angle
                #################################################

                nl1 = k1 - off_
                nl2 = k1 + off_ + 1
                ns1 = k2 - off_
                ns2 = k2 + off_ + 1
                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))
                w = np.transpose(w)
                # coeffprior = np.asarray(lstsq(w, b))[0]
                w = np.matrix(w)
                b = np.matrix(b)
                b = np.transpose(b)

                # if(b[0]<b[1] and conti1==1):continue

                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
                Cm = np.diag([2.0, 1.0])
                # ctprior = np.transpose(np.matrix(ctprior))
                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
                ct = np.asarray(coeff)

                if (np.min(ct) < 3): continue
                ct = self.invplanck0(self.wl8, ct)

                result[:, k1, k2] = [ct[0], ct[1]]

        result[result < 0] = 0
        result[result > 350] = 0

        if singlecheck == 0:
            if aodcheck == 1:
                outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw.tif'
            else:
                outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_test.tif'
        else:
            if aodcheck == 1:
                outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_single.tif'
            else:
                outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_test_single.tif'
        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)

        return 0

    def inversionLSCT_agl_raw_EVO(self, wdir, dataTime, medName, resultName, off=2, off_=0):

        # get soil emissivity
        # infile = wdir + r'auxiliary\heihe_emis1.tif'
        # [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # emis1_s = emis1_s / 1000.0

        # infile = wdir + r'auxiliary\heihe_emis2.tif'
        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        # emis2_s = emis2_s / 1000.0

        # infile = wdir + r'auxiliary\heihe_classif_new.tif'
        # [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
        # classif = np.asarray(classif, np.int)
        # infile = wdir + r'auxiliary\ICBP_Emi.txt'
        # emis_v = self.getDatafromTxt(infile, 0, 2)
        # emis1_v = emis_v[classif, 0] / 1000.0
        # emis2_v = emis_v[classif, 1] / 1000.0

        infile = wdir + resultName + dataTime + r'_lst_n.tif'
        [lst_n, ns, nl, nb, geog, proj] = self.getTiffData(infile)
        infile = wdir + resultName + dataTime + r'_lst_o.tif'
        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)

        emis1_v = np.zeros([nl, ns])
        emis1_s = np.zeros([nl, ns])
        emis1_v[:] = 0.967
        emis1_s[:] = 0.968
        F = 0.65

        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_n[fvc_n > 1] = 1
        fvc_n[fvc_n < 0] = 0
        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
        data = dataset.GetRasterBand(1)
        ndvi = data.ReadAsArray()
        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
        fvc_o[fvc_o > 1] = 1
        fvc_o[fvc_o < 0] = 0

        rad_n = self.planck(self.wl8, lst_n)
        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_n))
        rad_o = self.planck(self.wl8, lst_o)
        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_o))

        # component temperature inversion
        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n) + (1 - emis1_s) * emis1_v * F * (1 - fvc_n)])
        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o) + (1 - emis1_s) * emis1_v * F * (1 - fvc_o)])
        result = np.zeros([2, nl, ns])

        nsi = np.zeros([nl, ns])
        nli = np.zeros([nl, ns])
        for k in range(nl):
            nli[k, :] = k
            nsi[k, :] = np.linspace(0, ns - 1, ns)

        Cdvalue = 0.50
        Cdvaluep = 0.30
        puritypixe = 0.05
        # off_ = 0

        mask = np.zeros(np.shape(lst_n))
        ind = (fvc_n > 0.05) * (fvc_n < 0.95) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
                rad_o < 17)
        mask[ind] = 1
        part = 6
        for k1 in range(off, nl - off):

            for k2 in range(off, ns - off):

                if (mask[k1, k2] != 1): continue

                nl1 = k1 - off
                nl2 = k1 + off + 1
                ns1 = k2 - off
                ns2 = k2 + off + 1

                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
                # fd[fd != 0] = 1 / fd[fd != 0]
                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
                # fd[halfpoint] = 5
                # fd[:] = 1
                fd = np.asarray(
                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
                fd = np.reshape(fd, -1) * 1.0

                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                w0 = W1[:, k1, k2]

                n0 = np.ones(len(ns_temp))
                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])

                b1 = (np.reshape(rad_n[k1, k2], -1))
                b2 = (np.reshape(rad_o[k1, k2], -1))
                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))

                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
                num_point = np.sum(ind)
                if (num_point < 6): continue

                ww = np.hstack((w11[:, ind], w22[:, ind]))
                bb = np.hstack((bb1[ind], bb2[ind]))

                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))

                ###############################################
                ### the least-square method
                ###############################################
                ww = np.transpose(ww)
                coeffprior = np.asarray(lstsq(ww, bb))[0]
                # cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ###############################################
                ### the least-square method
                ###############################################
                ### the multi-pixel method with a prior knowledge
                # bbb = sorted(b)
                # halfpoint = np.int(num_point / 2)
                # rleaf = np.average(bbb[0])
                # rsoil = np.average(bbb[1])
                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
                # bb = np.transpose(np.matrix(bb))
                # ww = np.transpose(ww)
                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
                # Cm = np.matrix(np.diag(np.ones(2 * part)))
                # Cm[0, 0] = 15.0
                # Cm[part, part] = 10.0
                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)

                #################################################
                ### result without averaging information
                #################################################
                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                # ctprior = np.transpose(np.matrix([cd1, cd2]))

                ################################################
                ### reult with averaing information
                ################################################
                coeffprior = np.asarray(coeffprior)
                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
                apoint = np.sum(indnew)
                if (apoint <= 1):
                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                    ctprior = np.transpose(np.matrix([cd1, cd2]))
                else:
                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
                    bb0 = bb1[offcenter]
                    fdnew = np.zeros(25)
                    dismax = max(abs(bb1[indnew] - bb0))
                    if dismax != 0:
                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
                    else:
                        fdnew[:] = 1
                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])

                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
                    cd1new = np.sum(ct1[indnew] * fd[indnew])
                    cd2new = np.sum(ct2[indnew] * fd[indnew])
                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))

                #################################################
                ### multi-angle
                #################################################

                # nl1 = k1 - off_
                # nl2 = k1 + off_ + 1
                # ns1 = k2 - off_
                # ns2 = k2 + off_ + 1
                # w1 = W1[:, nl1:nl2, ns1:ns2]
                # w1 = np.asarray(np.reshape(w1, (2, -1)))
                # w2 = W2[:, nl1:nl2, ns1:ns2]
                # w2 = np.asarray(np.reshape(w2, (2, -1)))
                # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
                # w = np.hstack((w1, w2))
                # b = np.hstack((b1, b2))
                # w = np.transpose(w)
                # # coeffprior = np.asarray(lstsq(w, b))[0]
                # w = np.matrix(w)
                # b = np.matrix(b)
                # b = np.transpose(b)
                #
                # # if(b[0]<b[1] and conti1==1):continue
                #
                # Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
                # Cm = np.diag([10.0, 10.0])
                # # ctprior = np.transpose(np.matrix(ctprior))
                # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
                # coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
                # # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
                # ct = np.asarray(coeff)

                nl1 = k1 - off_
                nl2 = k1 + off_ + 1
                ns1 = k2 - off_
                ns2 = k2 + off_ + 1
                w1 = W1[:, nl1:nl2, ns1:ns2]
                w1 = np.asarray(np.reshape(w1, (2, -1)))
                w2 = W2[:, nl1:nl2, ns1:ns2]
                w2 = np.asarray(np.reshape(w2, (2, -1)))
                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
                w = np.hstack((w1, w2))
                b = np.hstack((b1, b2))
                w = np.transpose(w)
                ct = np.asarray(lstsq(w, b))[0]

                if (np.min(ct) < 3): continue
                ct = self.invplanck0(self.wl8, ct)

                result[:, k1, k2] = [ct[0], ct[1]]

        result[result < 0] = 0
        result[result > 350] = 0

        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_agl.tif'
        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)

        return 0

    ### This is for evaluation
    def inversionLSCT_npl_raw_EVO(self, wdir, dataTime, medName,resultName,off=2, off_=0):
    
       # get soil emissivity
       # infile = wdir + r'auxiliary\heihe_emis1.tif'
       # [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
       # emis1_s = emis1_s / 1000.0
    
       # infile = wdir + r'auxiliary\heihe_emis2.tif'
       # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
       # emis2_s = emis2_s / 1000.0
    
       # infile = wdir + r'auxiliary\heihe_classif_new.tif'
       # [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
       # classif = np.asarray(classif, np.int)
       # infile = wdir + r'auxiliary\ICBP_Emi.txt'
       # emis_v = self.getDatafromTxt(infile, 0, 2)
       # emis1_v = emis_v[classif, 0] / 1000.0
       # emis2_v = emis_v[classif, 1] / 1000.0
    
       infile = wdir + resultName + dataTime + r'_lst_n.tif'
       [lst_n, ns, nl, nb, geog, proj] = self.getTiffData(infile)
       infile = wdir + resultName + dataTime + r'_lst_o.tif'
       [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    
       emis1_v = np.zeros([nl,ns])
       emis1_s = np.zeros([nl,ns])
       emis1_v[:] = 0.967
       emis1_s[:] = 0.968
       F = 0.65
    
       infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
       dataset = gdal.Open(infile, gdal.GA_ReadOnly)
       data = dataset.GetRasterBand(1)
       ndvi = data.ReadAsArray()
       fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
       fvc_n[fvc_n > 1] = 1
       fvc_n[fvc_n < 0] = 0
       infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
       dataset = gdal.Open(infile, gdal.GA_ReadOnly)
       data = dataset.GetRasterBand(1)
       ndvi = data.ReadAsArray()
       fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
       fvc_o[fvc_o > 1] = 1
       fvc_o[fvc_o < 0] = 0

       rad_n = self.planck(self.wl8, lst_n)
       rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_n))
       rad_o = self.planck(self.wl8, lst_o)
       rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s + (1 - emis1_s) * emis1_v * F * (1 - fvc_o))

       # component temperature inversion
       W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n) + (1 - emis1_s) * emis1_v * F * (1 - fvc_n)])
       W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o) + (1 - emis1_s) * emis1_v * F * (1 - fvc_o)])
       result = np.zeros([2, nl, ns])
    
       nsi = np.zeros([nl, ns])
       nli = np.zeros([nl, ns])
       for k in range(nl):
           nli[k, :] = k
           nsi[k, :] = np.linspace(0, ns - 1, ns)
    
       Cdvalue = 0.60
       Cdvaluep = 0.30
       puritypixe = 0.05
       # off_ = 0
    
       mask = np.zeros(np.shape(lst_n))
       ind = (fvc_n > 0.05) * (fvc_n < 0.95) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
               rad_o < 17)
       mask[ind] = 1
       part = 6
       for k1 in range(off, nl - off):
    
           for k2 in range(off, ns - off):
    
               if (mask[k1, k2] != 1): continue
    
               nl1 = k1 - off
               nl2 = k1 + off + 1
               ns1 = k2 - off
               ns2 = k2 + off + 1
    
               offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
               ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
               nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
               # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
               # fd[fd != 0] = 1 / fd[fd != 0]
               # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
               # fd[halfpoint] = 5
               # fd[:] = 1
               fd = np.asarray(
                   [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
               fd = np.reshape(fd, -1) * 1.0
    
               w1 = W1[:, nl1:nl2, ns1:ns2]
               w1 = np.asarray(np.reshape(w1, (2, -1)))
               w2 = W2[:, nl1:nl2, ns1:ns2]
               w2 = np.asarray(np.reshape(w2, (2, -1)))
               w0 = W1[:, k1, k2]
    
               n0 = np.ones(len(ns_temp))
               x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
               w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
               w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    
               b1 = (np.reshape(rad_n[k1, k2], -1))
               b2 = (np.reshape(rad_o[k1, k2], -1))
               bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    
               temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
               temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
               ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
               num_point = np.sum(ind)
               if (num_point < 6): continue
    
               ww = np.hstack((w11[:, ind], w22[:, ind]))
               bb = np.hstack((bb1[ind], bb2[ind]))
    
               w = np.hstack((w1, w2))
               b = np.hstack((b1, b2))
    
               ###############################################
               ### the least-square method
               ###############################################
               ww = np.transpose(ww)
               coeffprior = np.asarray(lstsq(ww, bb))[0]
               cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
               cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
               ct = np.transpose(np.matrix([cd1, cd2]))
               ### or
               # ctprior = np.transpose(np.matrix([cd1, cd2]))
    
               ###############################################
               ### the least-square method
               ###############################################
               ### the multi-pixel method with a prior knowledge
               # bbb = sorted(b)
               # halfpoint = np.int(num_point / 2)
               # rleaf = np.average(bbb[0])
               # rsoil = np.average(bbb[1])
               # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
               # bb = np.transpose(np.matrix(bb))
               # ww = np.transpose(ww)
               # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
               # Cm = np.matrix(np.diag(np.ones(2 * part)))
               # Cm[0, 0] = 15.0
               # Cm[part, part] = 10.0
               # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
               # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    
               ################################################
               ### reult with averaing information
               ################################################
               # ww = np.transpose(ww)
               # coeffprior = np.asarray(lstsq(ww, bb))[0]
               # coeffprior = np.asarray(coeffprior)
               # ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
               # ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
               # indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
               # apoint = np.sum(indnew)
               # if (apoint <= 1):
               #     cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
               #     cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
               #     ct = np.transpose(np.matrix([cd1, cd2]))
               #     ctprior = ct
               # else:
               #     fd[indnew] = fd[indnew] / np.sum(fd[indnew])
               #     bb0 = bb1[offcenter]
               #     fdnew = np.zeros(25)
               #     dismax = max(abs(bb1[indnew] - bb0))
               #     if dismax != 0:
               #         fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
               #     else:
               #         fdnew[:] = 1
               #     fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
               #
               #     # fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
               #     cd1new = np.sum(ct1[indnew] * fd[indnew])
               #     cd2new = np.sum(ct2[indnew] * fd[indnew])
               #     ct = np.transpose(np.matrix([cd1new, cd2new]))
               #     ctprior = ct
               #
               # #################################################
               # ### multi-angle
               # #################################################
               #
               # # nl1 = k1 - off_
               # # nl2 = k1 + off_ + 1
               # # ns1 = k2 - off_
               # # ns2 = k2 + off_ + 1
               # # w1 = W1[:, nl1:nl2, ns1:ns2]
               # # w1 = np.asarray(np.reshape(w1, (2, -1)))
               # # w2 = W2[:, nl1:nl2, ns1:ns2]
               # # w2 = np.asarray(np.reshape(w2, (2, -1)))
               # # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               # # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
               # # w = np.hstack((w1, w2))
               # # b = np.hstack((b1, b2))
               # # w = np.transpose(w)
               # # # coeffprior = np.asarray(lstsq(w, b))[0]
               # # w = np.matrix(w)
               # # b = np.matrix(b)
               # # b = np.transpose(b)
               # #
               # # # if(b[0]<b[1] and conti1==1):continue
               # #
               # # Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
               # # Cm = np.diag([10.0, 10.0])
               # # # ctprior = np.transpose(np.matrix(ctprior))
               # # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
               # # coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
               # # # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
               # # ct = np.asarray(coeff)
               #
               # nl1 = k1 - off_
               # nl2 = k1 + off_ + 1
               # ns1 = k2 - off_
               # ns2 = k2 + off_ + 1
               # w1 = W1[:, nl1:nl2, ns1:ns2]
               # w1 = np.asarray(np.reshape(w1, (2, -1)))
               # w2 = W2[:, nl1:nl2, ns1:ns2]
               # w2 =   np.asarray(np.reshape(w2, (2, -1)))
               # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
               # w = np.hstack((w1, w2))
               # b = np.hstack((b1, b2))
               # w = np.transpose(w)
               # ct = np.asarray(lstsq(w, b))[0]
    
               if (np.min(ct) < 3): continue
               ct = self.invplanck0(self.wl8, ct)
    
               result[:, k1, k2] = [ct[0], ct[1]]
    
       result[result < 0] = 0
       result[result > 350] = 0
    
       outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_npl.tif'
       self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    
       return 0
    
    
    
       ##################################
      #### the very normal multi-pixel method and multi-angle method
       #################################
    def inversionLSCT(self, wdir, dataTime, medName,resultName,off=2):

       # get soil emissivity
       infile = wdir + r'auxiliary\heihe_emis1.tif'
       [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
       emis1_s = emis1_s / 1000.0

       # infile = wdir + r'auxiliary\heihe_emis2.tif'
       # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
       # emis2_s = emis2_s / 1000.0

       infile = wdir + r'auxiliary\heihe_classif_new.tif'
       [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
       classif = np.asarray(classif, np.int)
       infile = wdir + r'auxiliary\ICBP_Emi.txt'
       emis_v = self.getDatafromTxt(infile, 0, 3)
       emis1_v = emis_v[classif, 0] / 1000.0
       # emis2_v = emis_v[classif, 1] / 1000.0

       infile = wdir + resultName + dataTime + r'_lst_n.tif'
       [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
       infile = wdir + resultName + dataTime + r'_lst_o.tif'
       [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)

       infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
       dataset = gdal.Open(infile, gdal.GA_ReadOnly)
       data = dataset.GetRasterBand(1)
       ndvi = data.ReadAsArray()
       fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
       fvc_n[fvc_n > 1] = 1
       fvc_n[fvc_n < 0] = 0
       infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
       dataset = gdal.Open(infile, gdal.GA_ReadOnly)
       data = dataset.GetRasterBand(1)
       ndvi = data.ReadAsArray()
       fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
       fvc_o[fvc_o > 1] = 1
       fvc_o[fvc_o < 0] = 0

       rad_n = self.planck(self.wl8, lst_n)
       rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
       rad_o = self.planck(self.wl8, lst_o)
       rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)

       # component temperature inversion
       W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
       W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
       result = np.zeros([2, nl, ns])

       nsi = np.zeros([nl, ns])
       nli = np.zeros([nl, ns])
       for k in range(nl):
           nli[k, :] = k
           nsi[k, :] = np.linspace(0, ns - 1, ns)
       Cm = np.matrix(np.diag([10.0, 5.0]))
       Cdvalue = 1.0
       Cdvaluep = 0.30
       puritypixe = 0.05
       off_ = 0

       mask = np.zeros(np.shape(lst_n))
       ind = (fvc_n > 0.05) * (fvc_n < 0.95) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
                   rad_o < 17)
       mask[ind] = 1
       part = 6
       for k1 in range(off, nl - off):

           for k2 in range(off, ns - off):


               if(mask[k1,k2]!=1):continue

               nl1 = k1 - off
               nl2 = k1 + off + 1
               ns1 = k2 - off
               ns2 = k2 + off + 1

               offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
               ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
               nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
               # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
               # fd[fd != 0] = 1 / fd[fd != 0]
               # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
               # fd[halfpoint] = 5
               # fd[:] = 1
               fd = np.asarray(
                   [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
               fd = np.reshape(fd, -1) * 1.0

               w1 = W1[:, nl1:nl2, ns1:ns2]
               w1 = np.asarray(np.reshape(w1, (2, -1)))
               w2 = W2[:, nl1:nl2, ns1:ns2]
               w2 = np.asarray(np.reshape(w2, (2, -1)))
               w0 = W1[:,k1,k2]


               n0 = np.ones(len(ns_temp))
               x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
               w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
               w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])

               b1 = (np.reshape(rad_n[k1, k2], -1))
               b2 = (np.reshape(rad_o[k1, k2], -1))
               bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))

               temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
               temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
               ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
               num_point = np.sum(ind)
               if (num_point < 12): continue

               ww = np.hstack((w11[:, ind], w22[:, ind]))
               bb = np.hstack((bb1[ind], bb2[ind]))
               #
               # ww = w11[:, ind]
               # bb = bb1[ind]

               w = np.hstack((w1, w2))
               b = np.hstack((b1, b2))

               ###############################################
               ### the least-square method
               ###############################################
               ww = np.transpose(ww)
               coeffprior = np.asarray(lstsq(ww, bb))[0]
               # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
               # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
               # ctprior = np.transpose(np.matrix([cd1, cd2]))

               ###############################################
               ### the least-square method
               ###############################################
               ### the multi-pixel method with a prior knowledge
               # bbb = sorted(b)
               # halfpoint = np.int(num_point / 2)
               # rleaf = np.average(bbb[0])
               # rsoil = np.average(bbb[1])
               # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
               # bb = np.transpose(np.matrix(bb))
               # ww = np.transpose(ww)
               # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
               # Cm = np.matrix(np.diag(np.ones(2 * part)))
               # Cm[0, 0] = 15.0
               # Cm[part, part] = 10.0
               # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
               # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)

               #################################################
               ### result without averaging information
               #################################################
               # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
               # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
               # ctprior = np.transpose(np.matrix([cd1, cd2]))

               ################################################
               ### reult with averaing information
               ################################################
               coeffprior = np.asarray(coeffprior)
               ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
               ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
               indnew = ind * (ct2 < 16.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 16.5)
               apoint = np.sum(indnew)
               if (apoint <= 1):
                   cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
                   cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
                   ctprior = np.transpose(np.matrix([cd1, cd2]))
               else:
                   fd[indnew] = fd[indnew] / np.sum(fd[indnew])
                   bb0 = bb1[offcenter]
                   fdnew = np.zeros(25)
                   # inver-distance weight
                   dismax = max(abs(bb1[indnew] - bb0))
                   if dismax != 0:
                       fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
                   else:
                       fdnew[:] = 1
                   fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])

                   # fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
                   cd1new = np.sum(ct1[indnew] * fd[indnew])
                   cd2new = np.sum(ct2[indnew] * fd[indnew])
                   ctprior = np.transpose(np.matrix([cd1new, cd2new]))

               #################################################
               ### multi-angle
               #################################################

               nl1 = k1 - off_
               nl2 = k1 + off_ + 1
               ns1 = k2 - off_
               ns2 = k2 + off_ + 1
               w1 = W1[:, nl1:nl2, ns1:ns2]
               w1 = np.asarray(np.reshape(w1, (2, -1)))
               w2 = W2[:, nl1:nl2, ns1:ns2]
               w2 = np.asarray(np.reshape(w2, (2, -1)))
               b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
               b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
               w = np.hstack((w1, w2))
               b = np.hstack((b1, b2))
               w = np.transpose(w)
               # coeffprior = np.asarray(lstsq(w, b))[0]
               w = np.matrix(w)
               b = np.matrix(b)
               b = np.transpose(b)

               # if(b[0]<b[1] and conti1==1):continue

               Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
               Cm = np.diag([2.0, 1.0])
               # ctprior = np.transpose(np.matrix(ctprior))
               pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
               coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
               # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
               ct = np.asarray(coeff)

               if (np.min(ct) < 3): continue
               ct = self.invplanck0(self.wl8, ct)

               result[:, k1, k2] = [ct[0], ct[1]]

       result[result < 0] = 0
       result[result > 350] = 0

       outfile = wdir + resultName + dataTime + '_LSCT_'+str(off)+'_comb.tif'
       self.writeTiff(result, ns, nl, 2, geog, proj, outfile)

       return 0
    #
    #
    #    ##################################
    #    #### the zhan multi-pixel method and multi-angle method
    #    #################################
    #    def inversionLSCT_npl_agl(self, wdir, dataTime, medName,resultName,off=2,off_ = 0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n_new.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o_new.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                    rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if(mask[k1,k2]!=1):continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:,k1,k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 6): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                # ww = np.transpose(ww)
    #                # coeffprior = np.asarray(lstsq(ww, bb))[0]
    #
    #
    #                ### the multi-pixel method with a prior knowledge
    #                bbb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.average(bbb[0])
    #                rsoil = np.average(bbb[1])
    #                coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                bb = np.transpose(np.matrix(bb))
    #                ww = np.transpose(ww)
    #                Cd = np.matrix(np.diag(np.ones(num_point*2)) * Cdvaluep)
    #                Cm = np.matrix(np.diag(np.ones(2*part)))
    #                Cm[0,0] = 15.0
    #                Cm[part,part] = 10.0
    #                pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                coeffprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if(apoint<=1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew]-bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0-abs(1.0*bb1[indnew]-bb0)/dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew]/np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew]*fd[indnew])/np.sum(fdnew[indnew]*fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0,12.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_'+str(off)+'_new.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_npl_agl_EVO(self, wdir, dataTime, medName, resultName, off=2, off_=0):
    #
    #        # get soil emissivity
    #        # infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        # [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        # [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        # classif = np.asarray(classif, np.int)
    #        # infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        # emis_v = self.getDatafromTxt(infile, 0, 2)
    #        # emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        emis1_s = np.zeros([nl,ns])
    #        emis1_v = np.zeros([nl,ns])
    #        emis1_s[:] = 0.968
    #        emis1_v[:] = 0.993
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.50
    #        Cdvaluep = 0.55
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 6): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                # ww = np.transpose(ww)
    #                # coeffprior = np.asarray(lstsq(ww, bb))[0]
    #
    #                ### the multi-pixel method with a prior knowledge
    #                bbb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.average(bbb[0])
    #                rsoil = np.average(bbb[1])
    #                coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                bb = np.transpose(np.matrix(bb))
    #                ww = np.transpose(ww)
    #                Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                Cm[0, 0] = 15.0
    #                Cm[part, part] = 10.0
    #                pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                coeffprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0, 12.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    ##################################
    #    #### the zhan multi-pixel method and multi-angle method
    #    #################################
    #    def inversionLSCT_npl_agl_point(self, wdir, dataTime, lat_,lon_,medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #
    #        ######
    #        ## pints
    #        ######
    #        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
    #        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
    #        imagex = np.asarray(imagex, np.int)
    #        imagey = np.asarray(imagey, np.int)
    #
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #                temp =  (imagex-k1)*(imagex-k1)+(imagey-k2)*(imagey-k2)
    #                ind = temp < 25
    #                if (np.sum(ind)<1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 6): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                # ww = np.transpose(ww)
    #                # coeffprior = np.asarray(lstsq(ww, bb))[0]
    #
    #                ### the multi-pixel method with a prior knowledge
    #                bbb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = (bbb[0])
    #                rsoil = (bbb[1])
    #                coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                bb = np.transpose(np.matrix(bb))
    #                ww = np.transpose(ww)
    #                Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                Cm[0, 0] = 15.0
    #                Cm[part, part] = 10.0
    #                pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                coeffprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([5.0, 2.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_new.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #
    #    ##################################
    #    #### the zhan multi-pixel method and multi-angle method
    #    #################################
    #    def inversionLSCT_npl_agl_point_sgl(self, wdir, dataTime, lat_,lon_,medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #
    #        ######
    #        ## pints
    #        ######
    #        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
    #        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
    #        imagex = np.asarray(imagex, np.int)
    #        imagey = np.asarray(imagey, np.int)
    #
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #                temp =  (imagex-k1)*(imagex-k1)+(imagey-k2)*(imagey-k2)
    #                ind = temp < 25
    #                if (np.sum(ind)<1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 12): continue
    #
    #                # ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                # bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #
    #                ww = ((w11[:, ind]))
    #                bb = ((bb1[ind]))
    #
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                # ww = np.transpose(ww)
    #                # coeffprior = np.asarray(lstsq(ww, bb))[0]
    #
    #                ### the multi-pixel method with a prior knowledge
    #                bbb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.min(bb)
    #                rsoil = np.max(bb)
    #                coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                bb = np.transpose(np.matrix(bb))
    #                ww = np.transpose(ww)
    #                Cd = np.matrix(np.diag(np.ones(num_point)) * Cdvaluep)
    #                Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                Cm[0, 0] = 15.0
    #                Cm[part, part] = 10.0
    #                pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                coeffprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0, 12.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_sgl.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    ##################################
    #    #### the zhan multi-pixel method and multi-angle method
    #    #################################
    #    def inversionLSCT_npl_agl_test(self, wdir, dataTime, medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n_new.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o_new.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 6): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                # ww = np.transpose(ww)
    #                # coeffprior = np.asarray(lstsq(ww, bb))[0]
    #
    #                ### the multi-pixel method with a prior knowledge
    #                bbb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.average(bbb[0])
    #                rsoil = np.average(bbb[1])
    #                coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                bb = np.transpose(np.matrix(bb))
    #                ww = np.transpose(ww)
    #                Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                Cm[0, 0] = 15.0
    #                Cm[part, part] = 10.0
    #                pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                coeffprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                # coeffprior = np.asarray(coeffprior)
    #                # ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                # ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                # indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                # apoint = np.sum(indnew)
    #                # if (apoint <= 1):
    #                #     cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                #     cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                #     ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                # else:
    #                #     fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                #     bb0 = bb1[offcenter]
    #                #     fdnew = np.zeros(25)
    #                #     dismax = max(abs(bb1[indnew] - bb0))
    #                #     if dismax != 0:
    #                #         fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                #     else:
    #                #         fdnew[:] = 1
    #                #     fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #                #
    #                #     fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                #     cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                #     cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                #     ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0, 12.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_test.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    ##################################
    #    #### the zhan multi-pixel method and multi-angle method
    #    #################################
    #    def inversionLSCT_npl_agl_point_test(self, wdir, dataTime, lat_,lon_, medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n_new.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o_new.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        ######
    #        ## pints
    #        ######
    #        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
    #        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
    #        imagex = np.asarray(imagex, np.int)
    #        imagey = np.asarray(imagey, np.int)
    #
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #                temp = (imagex - k1) * (imagex - k1) + (imagey - k2) * (imagey - k2)
    #                ind = temp < 25
    #                if (np.sum(ind) < 1): continue
    #
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 12): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                # ww = ((w11[:, ind]))
    #                # bb = ((bb1[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                # ww = np.transpose(ww)
    #                # coeffprior = np.asarray(lstsq(ww, bb))[0]
    #
    #                ### the multi-pixel method with a prior knowledge
    #                bbb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.average(bbb[0])
    #                rsoil = np.average(bbb[1])
    #                coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                bb = np.transpose(np.matrix(bb))
    #                ww = np.transpose(ww)
    #                Cd = np.matrix(np.diag(np.ones(num_point*2)) * Cdvaluep)
    #                Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                Cm[0, 0] = 15.0
    #                Cm[part, part] = 10.0
    #                pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                coeffprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                # coeffprior = np.asarray(coeffprior)
    #                # ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                # ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                # indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                # apoint = np.sum(indnew)
    #                # if (apoint <= 1):
    #                #     cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                #     cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                #     ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                # else:
    #                #     fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                #     bb0 = bb1[offcenter]
    #                #     fdnew = np.zeros(25)
    #                #     dismax = max(abs(bb1[indnew] - bb0))
    #                #     if dismax != 0:
    #                #         fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                #     else:
    #                #         fdnew[:] = 1
    #                #     fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #                #
    #                #     fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                #     cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                #     cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                #     ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0, 12.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_test.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #
    #
    #
    #
    #    def inversionLSCT_npl_agl_raw_EVO_bt9(self, wdir, dataTime, medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        # infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        # [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        # [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        # classif = np.asarray(classif, np.int)
    #        # infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        # emis_v = self.getDatafromTxt(infile, 0, 2)
    #        # emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        emis1_v = np.zeros([nl,ns])
    #        emis1_s = np.zeros([nl,ns])
    #        emis1_v[:] = 0.996
    #        emis1_s[:] = 0.977
    #
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #
    #
    #        rad_n = self.planck(self.wl9, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl9, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.50
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 6): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ww = np.transpose(ww)
    #                coeffprior = np.asarray(lstsq(ww, bb))[0]
    #                # cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ### the multi-pixel method with a prior knowledge
    #                # bbb = sorted(b)
    #                # halfpoint = np.int(num_point / 2)
    #                # rleaf = np.average(bbb[0])
    #                # rsoil = np.average(bbb[1])
    #                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                # bb = np.transpose(np.matrix(bb))
    #                # ww = np.transpose(ww)
    #                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                # Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                # Cm[0, 0] = 15.0
    #                # Cm[part, part] = 10.0
    #                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0, 15.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl9, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_test.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #
    #    def inversionLSCT_npl_agl_raw(self, wdir, dataTime, medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, ns, nl, temp3, geog, proj] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        # emis1_v = np.zeros([nl,ns])
    #        # emis1_s = np.zeros([nl,ns])
    #        # emis1_v[:] = 0.996
    #        # emis1_s[:] = 0.977
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 6): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ww = np.transpose(ww)
    #                coeffprior = np.asarray(lstsq(ww, bb))[0]
    #                # cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ### the multi-pixel method with a prior knowledge
    #                # bbb = sorted(b)
    #                # halfpoint = np.int(num_point / 2)
    #                # rleaf = np.average(bbb[0])
    #                # rsoil = np.average(bbb[1])
    #                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                # bb = np.transpose(np.matrix(bb))
    #                # ww = np.transpose(ww)
    #                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                # Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                # Cm[0, 0] = 15.0
    #                # Cm[part, part] = 10.0
    #                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0, 12.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_agl_raw(self, wdir, dataTime, medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, ns, nl, temp3, geog, proj] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        # emis1_v = np.zeros([nl,ns])
    #        # emis1_s = np.zeros([nl,ns])
    #        # emis1_v[:] = 0.996
    #        # emis1_s[:] = 0.977
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 6): continue
    #
    #                ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                bb = np.hstack((bb1[ind], bb2[ind]))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ww = np.transpose(ww)
    #                coeffprior = np.asarray(lstsq(ww, bb))[0]
    #                # cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ### the multi-pixel method with a prior knowledge
    #                # bbb = sorted(b)
    #                # halfpoint = np.int(num_point / 2)
    #                # rleaf = np.average(bbb[0])
    #                # rsoil = np.average(bbb[1])
    #                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                # bb = np.transpose(np.matrix(bb))
    #                # ww = np.transpose(ww)
    #                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                # Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                # Cm[0, 0] = 15.0
    #                # Cm[part, part] = 10.0
    #                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                # nl1 = k1 - off_
    #                # nl2 = k1 + off_ + 1
    #                # ns1 = k2 - off_
    #                # ns2 = k2 + off_ + 1
    #                # w1 = W1[:, nl1:nl2, ns1:ns2]
    #                # w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                # w2 = W2[:, nl1:nl2, ns1:ns2]
    #                # w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                # w = np.hstack((w1, w2))
    #                # b = np.hstack((b1, b2))
    #                # w = np.transpose(w)
    #                # # coeffprior = np.asarray(lstsq(w, b))[0]
    #                # w = np.matrix(w)
    #                # b = np.matrix(b)
    #                # b = np.transpose(b)
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                ct = np.asarray(lstsq(w, b))[0]
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                # Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                # Cm = np.diag([15.0, 12.0])
    #                # # ctprior = np.transpose(np.matrix(ctprior))
    #                # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                # coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                # ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_raw_agl.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    # ### This is for evaluation
    #
    #
    #    def inversionLSCT_npl_agl_point_pur(self, wdir, dataTime, lat_,lon_, medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        ######
    #        ## pints
    #        ######
    #        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
    #        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
    #        imagex = np.asarray(imagex, np.int)
    #        imagey = np.asarray(imagey, np.int)
    #
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #                temp = (imagex - k1) * (imagex - k1) + (imagey - k2) * (imagey - k2)
    #                ind = temp < 25
    #                if (np.sum(ind) < 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 12): continue
    #
    #                # ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                # bb = np.hstack((bb1[ind], bb2[ind]))
    #                #
    #                ww = w11[:, ind]
    #                bb = bb1[ind]
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ww = np.transpose(ww)
    #                coeffprior = np.asarray(lstsq(ww, bb))[0]
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ### the multi-pixel method with a prior knowledge
    #                # bbb = sorted(b)
    #                # halfpoint = np.int(num_point / 2)
    #                # rleaf = np.average(bbb[0])
    #                # rsoil = np.average(bbb[1])
    #                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                # bb = np.transpose(np.matrix(bb))
    #                # ww = np.transpose(ww)
    #                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                # Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                # Cm[0, 0] = 15.0
    #                # Cm[part, part] = 10.0
    #                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                # coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                # if(b[0]<b[1] and conti1==1):continue
    #
    #                Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                Cm = np.diag([15.0, 12.0])
    #                # ctprior = np.transpose(np.matrix(ctprior))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_pur.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_npl_point(self, wdir, dataTime, lat_, lon_, medName,resultName,off=2, off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n_new.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o_new.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        Cdvalue = 0.15
    #        Cdvaluep = 0.30
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.15) * (fvc_n < 0.85) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                rad_o < 17)
    #        mask[ind] = 1
    #        part = 6
    #        ######
    #        ## pints
    #        ######
    #        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
    #        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
    #        imagex = np.asarray(imagex, np.int)
    #        imagey = np.asarray(imagey, np.int)
    #
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #                temp = (imagex - k1) * (imagex - k1) + (imagey - k2) * (imagey - k2)
    #                ind = temp < 25
    #                if (np.sum(ind) < 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #                # fd = np.sqrt(ns_temp * ns_temp + nl_temp * nl_temp)
    #                # fd[fd != 0] = 1 / fd[fd != 0]
    #                # halfpoint = np.int(((off * 2 + 1) * (2 * off + 1) - 1) / 2)
    #                # fd[halfpoint] = 5
    #                # fd[:] = 1
    #                fd = np.asarray(
    #                    [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                fd = np.reshape(fd, -1) * 1.0
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                w0 = W1[:, k1, k2]
    #
    #                n0 = np.ones(len(ns_temp))
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])
    #
    #                b1 = (np.reshape(rad_n[k1, k2], -1))
    #                b2 = (np.reshape(rad_o[k1, k2], -1))
    #                bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                ind = (bb1 < 17) * (bb1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                if (num_point < 12): continue
    #
    #                # ww = np.hstack((w11[:, ind], w22[:, ind]))
    #                # bb = np.hstack((bb1[ind], bb2[ind]))
    #                #
    #                ww = w11[:, ind]
    #                bb = bb1[ind]
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ww = np.transpose(ww)
    #                coeffprior = np.asarray(lstsq(ww, bb))[0]
    #
    #                ###############################################
    #                ### the least-square method
    #                ###############################################
    #                ### the multi-pixel method with a prior knowledge
    #                # bbb = sorted(b)
    #                # halfpoint = np.int(num_point / 2)
    #                # rleaf = np.average(bbb[0])
    #                # rsoil = np.average(bbb[1])
    #                # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    #                # bb = np.transpose(np.matrix(bb))
    #                # ww = np.transpose(ww)
    #                # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    #                # Cm = np.matrix(np.diag(np.ones(2 * part)))
    #                # Cm[0, 0] = 15.0
    #                # Cm[part, part] = 10.0
    #                # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    #                # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)
    #
    #                #################################################
    #                ### result without averaging information
    #                #################################################
    #                # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                # ctprior = np.transpose(np.matrix([cd1, cd2]))
    #
    #                ################################################
    #                ### reult with averaing information
    #                ################################################
    #                coeffprior = np.asarray(coeffprior)
    #                ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    #                ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
    #                indnew = ind * (ct2 < 14.5) * (ct2 > 7.5) * (ct1 > 7.5) * (ct1 < 14.5)
    #                apoint = np.sum(indnew)
    #                if (apoint <= 1):
    #                    cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    #                    cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    #                    ctprior = np.transpose(np.matrix([cd1, cd2]))
    #                else:
    #                    fd[indnew] = fd[indnew] / np.sum(fd[indnew])
    #                    bb0 = bb1[offcenter]
    #                    fdnew = np.zeros(25)
    #                    dismax = max(abs(bb1[indnew] - bb0))
    #                    if dismax != 0:
    #                        fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
    #                    else:
    #                        fdnew[:] = 1
    #                    fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
    #
    #                    fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])
    #                    cd1new = np.sum(ct1[indnew] * fd[indnew])
    #                    cd2new = np.sum(ct2[indnew] * fd[indnew])
    #                    ctprior = np.transpose(np.matrix([cd1new, cd2new]))
    #
    #                ct = np.asarray(ctprior)
    #                #################################################
    #                ### multi-angle
    #                #################################################
    #
    #                # nl1 = k1 - off_
    #                # nl2 = k1 + off_ + 1
    #                # ns1 = k2 - off_
    #                # ns2 = k2 + off_ + 1
    #                # w1 = W1[:, nl1:nl2, ns1:ns2]
    #                # w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                # w2 = W2[:, nl1:nl2, ns1:ns2]
    #                # w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                # w = np.hstack((w1, w2))
    #                # b = np.hstack((b1, b2))
    #                # w = np.transpose(w)
    #                # # coeffprior = np.asarray(lstsq(w, b))[0]
    #                # w = np.matrix(w)
    #                # b = np.matrix(b)
    #                # b = np.transpose(b)
    #                #
    #                # # if(b[0]<b[1] and conti1==1):continue
    #                #
    #                # Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                # Cm = np.diag([15.0, 12.0])
    #                # # ctprior = np.transpose(np.matrix(ctprior))
    #                # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                # coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                # ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_' + str(off) + '_npl.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_raw(self, wdir, dataTime, medName,resultName,off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #        Cm = np.matrix(np.diag([10.0, 5.0]))
    #        Cdvalue = 0.15
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.1) * (fvc_n < 0.9) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                    rad_o < 17)
    #        mask[ind] = 1
    #
    #        for k1 in range(off_, nl - off_):
    #
    #            for k2 in range(off_, ns - off_):
    #
    #
    #                if(mask[k1,k2]!=1):continue
    #
    #                # nl1 = k1 - off
    #                # nl2 = k1 + off + 1
    #                # ns1 = k2 - off
    #                # ns2 = k2 + off + 1
    #
    #                # fd = np.asarray(
    #                #     [[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    #                # fd = np.reshape(fd, -1) * 1.0
    #                #
    #                # w1 = W1[:, nl1:nl2, ns1:ns2]
    #                # w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                # w2 = W2[:, nl1:nl2, ns1:ns2]
    #                # w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                # w0 = W1[:, k1, k2]
    #                #
    #                # b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                # b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                #
    #                # temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
    #                # temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
    #                # ind = (b1 < 15) * (b1 > 7.5) * (temp1 > puritypixe) * (temp1 < (1 - puritypixe)) * (temp2 > temp1) * (
    #                #             b1 > b2)
    #                # num_point = np.sum(ind)
    #                # if (num_point < 1): continue
    #
    #                # fd[ind] = fd[ind] / np.sum(fd[ind])
    #                # # fdplus =np.cos((w1[0,ind]-w0[0])*np.pi/2)
    #                # # fd[ind] = fd[ind]*fdplus/np.sum(fd[ind]*fdplus)
    #                #
    #                # w = np.stack((np.sum(w1[:, ind] * fd[ind], 1), np.sum(w2[:, ind] * fd[ind], 1)))
    #                # b = np.hstack((np.sum(b1[ind] * fd[ind]), np.sum(b2[ind] * fd[ind])))
    #                # # conti1 = 0
    #                # # if(b[0]<b[1]):conti1=1
    #                # coeffprior = np.asarray(lstsq(w, b))[0]
    #                # ctprior = np.transpose(np.matrix(coeffprior))
    #
    #                nl1 = k1 - off_
    #                nl2 = k1 + off_ + 1
    #                ns1 = k2 - off_
    #                ns2 = k2 + off_ + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                ct = np.asarray(lstsq(w, b))[0]
    #                # w = np.matrix(w)
    #                # b = np.matrix(b)
    #                # b = np.transpose(b)
    #                #
    #                # # if(b[0]<b[1] and conti1==1):continue
    #                #
    #                # Cd = np.matrix(np.diag(np.ones(2)) * Cdvalue)
    #                # # ctprior = np.transpose(np.matrix(ctprior))
    #                # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    #                # coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    #                # # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    #                # ct = np.asarray(coeff)
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_'+str(off_)+'_raw.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_agl(self, wdir, dataTime, medName,resultName,off=2,off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n_new.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o_new.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #        Cm = np.matrix(np.diag([10.0, 5.0]))
    #        Cdvalue = 0.15
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.1) * (fvc_n < 0.9) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                    rad_o < 17)
    #        mask[ind] = 1
    #
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #
    #                if(mask[k1,k2]!=1):continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                # ct = np.asarray(lstsq(w, b))[0]
    #
    #                w = np.matrix(w)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #                Cmm = np.matrix(np.diag([5.0, 5.0]))
    #                Cdd = np.matrix(np.diag(np.ones(2)) * 0.30)
    #                bb = sorted(b)
    #                halfpoint = np.int(np.size(bb) / 2)
    #                rleaf = np.average(bb[:halfpoint])
    #                rsoil = np.average(bb[halfpoint:])
    #                # rleaf = rad_o[k1,k2]
    #                # rsoil = rad_n[k1,k2]
    #                coeffprior = np.transpose(np.matrix([rsoil,rleaf]))
    #                pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cdd) * w + np.linalg.inv(Cmm))
    #                ct = pre * (np.transpose(w) * np.linalg.inv(Cdd) * b + np.linalg.inv(Cmm) * coeffprior)
    #
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_'+str(off)+'_agl.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_agl_point(self, wdir, dataTime,lat_,lon_, medName,resultName,off=0,off_=0):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n_new.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o_new.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o_new.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #        Cm = np.matrix(np.diag([10.0, 5.0]))
    #        Cdvalue = 0.15
    #        puritypixe = 0.10
    #        # off_ = 0
    #
    #        mask = np.zeros(np.shape(lst_n))
    #        ind = (fvc_n > 0.1) * (fvc_n < 0.9) * (rad_n > 7.5) * (rad_n < 17.0) * (rad_o > 7.5) * (
    #                    rad_o < 17)
    #        mask[ind] = 1
    #
    #        [temp1, temp2] = self.lonlat2geo(dataset, lat_, lon_)
    #        imagey, imagex = self.geo2imagexy(dataset, temp1, temp2)
    #        imagex = np.asarray(imagex, np.int)
    #        imagey = np.asarray(imagey, np.int)
    #
    #        for k1 in range(off, nl - off):
    #
    #            for k2 in range(off, ns - off):
    #
    #                if (mask[k1, k2] != 1): continue
    #                temp = (imagex - k1) * (imagex - k1) + (imagey - k2) * (imagey - k2)
    #                ind = temp < 25
    #                if (np.sum(ind) < 1): continue
    #
    #                nl1 = k1 - off
    #                nl2 = k1 + off + 1
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.asarray(np.reshape(w1, (2, -1)))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.asarray(np.reshape(w2, (2, -1)))
    #                b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
    #                b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
    #
    #                w = np.hstack((w1, w2))
    #                b = np.hstack((b1, b2))
    #                w = np.transpose(w)
    #                ct = np.asarray(lstsq(w, b))[0]
    #
    #                # w = np.matrix(w)
    #                # b = np.matrix(b)
    #                # b = np.transpose(b)
    #                # Cmm = np.matrix(np.diag([5.0, 5.0]))
    #                # Cdd = np.matrix(np.diag(np.ones(2)) * 0.30)
    #                # bb = sorted(b)
    #                # halfpoint = np.int(np.size(bb) / 2)
    #                # rleaf = np.average(bb[:halfpoint])
    #                # rsoil = np.average(bb[halfpoint:])
    #                # # rleaf = rad_o[k1,k2]
    #                # # rsoil = rad_n[k1,k2]
    #                # coeffprior = np.transpose(np.matrix([rsoil,rleaf]))
    #                # pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cdd) * w + np.linalg.inv(Cmm))
    #                # ct = pre * (np.transpose(w) * np.linalg.inv(Cdd) * b + np.linalg.inv(Cmm) * coeffprior)
    #
    #
    #                if (np.min(ct) < 3): continue
    #                ct = self.invplanck0(self.wl8, ct)
    #
    #                result[:, k1, k2] = [ct[0], ct[1]]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_'+str(off)+'_agl.tif'
    #        self.writeTiff(result, ns, nl, 2, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_backup(self, wdir, dataTime, medName,resultName,off=2):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 1] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl9, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl9, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        puritypixe = 0.01
    #        for k1 in range(off, nl - off):
    #            nl1 = k1 - off
    #            nl2 = k1 + off + 1
    #            for k2 in range(off, ns - off):
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #
    #                temp1 = fvc_n[k1, k2]
    #                temp2 = fvc_o[k1, k2]
    #                temp3 = rad_n[k1, k2]
    #                temp4 = rad_o[k1, k2]
    #                if (temp3 < 6.0 or temp4 < 6.0): continue
    #                if (temp1 < puritypixe or temp2 < puritypixe or temp1 > (1 - puritypixe) or temp2 > (
    #                        1 - puritypixe)): continue
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.reshape(w1, (2, -1))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.reshape(w2, (2, -1))
    #
    #                b1 = np.reshape(rad_n[nl1:nl2, ns1:ns2], -1)
    #                b2 = np.reshape(rad_o[nl1:nl2, ns1:ns2], -1)
    #                n0 = np.ones(len(ns_temp))
    #
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * ns_temp, nl_temp * nl_temp, ns_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w = np.hstack((w11, w22))
    #                b = np.hstack((b1, b2))
    #                temp1 = np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1)
    #                temp2 = np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1)
    #                temp = np.hstack((temp1, temp2))
    #                ind = (b > 6.0) * (temp > puritypixe) * (temp < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                b = b[ind]
    #                w = w[:, ind]
    #                if (num_point < 15): continue
    #
    #                w = np.transpose(w)
    #                b = np.transpose(b)
    #                # w = np.matrix(w)
    #                # b = np.matrix(b)
    #                # ww = np.transpose(w) * w
    #                # coeff = (np.linalg.inv(ww)) * np.transpose(w) * np.transpose(b)
    #
    #                u, s, vt = scipy.linalg.svd(w, 0)
    #                v = np.transpose(vt)
    #                I = np.diag(np.ones(12))
    #                ss = s * s
    #                s = np.diag(s)
    #                ss = np.diag(ss)
    #
    #                s = np.matrix(s)
    #                v = np.matrix(v)
    #                vt = np.matrix(vt)
    #                u = np.matrix(u)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                bb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.average(bb[:halfpoint])
    #                rsoil = np.average(bb[halfpoint:])
    #                coeffprior = np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0])
    #                coeffprior = np.transpose(coeffprior)
    #
    #                temp = v * (np.linalg.inv(ss + I))
    #                coeff = temp * (s * np.transpose(u) * b + vt * coeffprior)
    #
    #                coeff = np.asarray(coeff)
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ct1 = np.sum(coeff[:6, 0] * x[:, offcenter])
    #                ct2 = np.sum(coeff[6:, 0] * x[:, offcenter])
    #                if ct1 <= 0 or ct2 <= 0: continue
    #                ct1 = self.invplanck0(self.wl9, ct1)
    #                ct2 = self.invplanck0(self.wl9, ct2)
    #
    #                result[:, k1, k2] = [ct1, ct2]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT.tif'
    #        self.writeTiff(result, ns, nl, 1, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_2(self, wdir, dataTime,medName,resultName, off=2):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + resultName + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTime + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad1_n = self.planck(self.wl8, lst_n)
    #        rad1_n = np.asarray(rad1_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s))
    #        rad1_o = self.planck(self.wl8, lst_o)
    #        rad1_o = np.asarray(rad1_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s))
    #
    #        rad2_n = self.planck(self.wl9, lst_n)
    #        rad2_n = np.asarray(rad2_n * (fvc_n * emis2_v + (1 - fvc_n) * emis2_s))
    #        rad2_o = self.planck(self.wl9, lst_o)
    #        rad2_o = np.asarray(rad2_o * (fvc_o * emis2_v + (1 - fvc_o) * emis2_s))
    #
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * fvc_n])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * fvc_o])
    #        W3 = np.asarray([emis2_s * (1 - fvc_n), emis2_v * fvc_n])
    #        W4 = np.asarray([emis2_s * (1 - fvc_o), emis2_v * fvc_o])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        puritypixe = 0.01
    #        for k1 in range(off, nl - off):
    #            nl1 = k1 - off
    #            nl2 = k1 + off + 1
    #            for k2 in range(off, ns - off):
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #
    #                temp1 = fvc_n[k1, k2]
    #                temp2 = fvc_o[k1, k2]
    #                temp3 = rad1_n[k1, k2]
    #                temp4 = rad1_o[k1, k2]
    #                if (temp3 < 6.0 or temp4 < 6.0): continue
    #                if (temp1 < puritypixe or temp2 < puritypixe or temp1 > (1 - puritypixe) or temp2 > (
    #                        1 - puritypixe)): continue
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.reshape(w1, (2, -1))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.reshape(w2, (2, -1))
    #                w3 = W3[:, nl1:nl2, ns1:ns2]
    #                w3 = np.reshape(w3, (2, -1))
    #                w4 = W4[:, nl1:nl2, ns1:ns2]
    #                w4 = np.reshape(w4, (2, -1))
    #
    #                b1 = np.reshape(rad1_n[nl1:nl2, ns1:ns2], -1)
    #                b2 = np.reshape(rad1_o[nl1:nl2, ns1:ns2], -1)
    #                b3 = np.reshape(rad2_n[nl1:nl2, ns1:ns2], -1)
    #                b4 = np.reshape(rad2_o[nl1:nl2, ns1:ns2], -1)
    #
    #                n0 = np.ones(len(ns_temp))
    #
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * ns_temp, nl_temp * nl_temp, ns_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w33 = np.asarray([w3[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w44 = np.asarray([w4[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w = np.hstack((w11, w22,w33,w44))
    #                b = np.hstack((b1, b2,b3,b4))
    #                temp1 = np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1)
    #                temp2 = np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1)
    #                temp = np.hstack((temp1, temp2,temp1,temp2))
    #                ind = (b > 6.0) * (temp > puritypixe) * (temp < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                b = b[ind]
    #                w = w[:, ind]
    #                if (num_point < 15): continue
    #
    #                w = np.transpose(w)
    #                b = np.transpose(b)
    #                # w = np.matrix(w)
    #                # b = np.matrix(b)
    #                # ww = np.transpose(w) * w
    #                # coeff = (np.linalg.inv(ww)) * np.transpose(w) * np.transpose(b)
    #
    #                u, s, vt = scipy.linalg.svd(w, 0)
    #                v = np.transpose(vt)
    #                I = np.diag(np.ones(12))
    #                ss = s * s
    #                s = np.diag(s)
    #                ss = np.diag(ss)
    #
    #                s = np.matrix(s)
    #                v = np.matrix(v)
    #                vt = np.matrix(vt)
    #                u = np.matrix(u)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                bb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.average(bb[:halfpoint])
    #                rsoil = np.average(bb[halfpoint:])
    #                coeffprior = np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0])
    #                coeffprior = np.transpose(coeffprior)
    #
    #                temp = v * (np.linalg.inv(ss + I))
    #                coeff = temp * (s * np.transpose(u) * b + vt * coeffprior)
    #
    #                coeff = np.asarray(coeff)
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ct1 = np.sum(coeff[:6, 0] * x[:, offcenter])
    #                ct2 = np.sum(coeff[6:, 0] * x[:, offcenter])
    #                if ct1 <= 0 or ct2 <= 0: continue
    #                ct1 = self.invplanck0(self.wl8, ct1)
    #                ct2 = self.invplanck0(self.wl8, ct2)
    #
    #                result[:, k1, k2] = [ct1, ct2]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_2.tif'
    #        self.writeTiff(result, ns, nl, 1, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_npl(self, wdir, dataTime, medName,resultName,off=2):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + resultName + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTime + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_v)
    #
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        puritypixe = 0.01
    #        for k1 in range(off, nl - off):
    #            nl1 = k1 - off
    #            nl2 = k1 + off + 1
    #            for k2 in range(off, ns - off):
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #
    #                temp1 = fvc_n[k1, k2]
    #
    #                temp3 = rad_n[k1, k2]
    #
    #                if (temp3 < 6.0): continue
    #                if (temp1 < puritypixe): continue
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.reshape(w1, (2, -1))
    #
    #
    #                b = np.reshape(rad_n[nl1:nl2, ns1:ns2], -1)
    #
    #                n0 = np.ones(len(ns_temp))
    #
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * ns_temp, nl_temp * nl_temp, ns_temp * nl_temp])
    #                w = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #
    #                # w = np.hstack((w11))
    #                # b = np.hstack((b1))
    #                temp = np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1)
    #
    #                # temp = np.hstack((temp1))
    #                ind = (b > 6.0) * (temp > puritypixe) * (temp < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                b = b[ind]
    #                w = w[:, ind]
    #                if (num_point < 15): continue
    #
    #                # w = np.transpose(w)
    #                # coeff = np.asarray(lstsq(w, b))[0]
    #
    #                w = np.transpose(w)
    #                b = np.transpose(b)
    #                u, s, vt = scipy.linalg.svd(w, 0)
    #                v = np.transpose(vt)
    #                I = np.diag(np.ones(12))
    #                ss = s * s
    #                s = np.diag(s)
    #                ss = np.diag(ss)
    #
    #                s = np.matrix(s)
    #                v = np.matrix(v)
    #                vt = np.matrix(vt)
    #                u = np.matrix(u)
    #                b = np.matrix(b)
    #                b = np.transpose(b)
    #
    #                bb = sorted(b)
    #                halfpoint = np.int(num_point / 2)
    #                rleaf = np.average(bb[:halfpoint])
    #                rsoil = np.average(bb[halfpoint:])
    #                coeffprior = np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0])
    #                coeffprior = np.transpose(coeffprior)
    #
    #                temp = v * (np.linalg.inv(ss + I))
    #                coeff = temp * (s * np.transpose(u) * b + vt * coeffprior)
    #
    #                coeff = np.asarray(coeff)
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ct1 = np.sum(coeff[:6,0] * x[:, offcenter])
    #                ct2 = np.sum(coeff[6:,0] * x[:, offcenter])
    #
    #                if ct1 <= 0 or ct2 <= 0: continue
    #                ct1 = self.invplanck0(self.wl8, ct1)
    #                ct2 = self.invplanck0(self.wl8, ct2)
    #                result[:, k1, k2] = [ct1, ct2]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + resultName + dataTime + '_LSCT_npl.tif'
    #        self.writeTiff(result, ns, nl, 1, geog, proj, outfile)
    #
    #        return 0
    #
    #    def inversionLSCT_night(self, wdir, dataTime,dataTimeold, medName,resultName,off=2):
    #
    #        # get soil emissivity
    #        infile = wdir + r'auxiliary\heihe_emis1.tif'
    #        [emis1_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        emis1_s = emis1_s / 1000.0
    #
    #        # infile = wdir + r'auxiliary\heihe_emis2.tif'
    #        # [emis2_s, ns, nl, nb, geog, proj] = self.getTiffData(infile)
    #        # emis2_s = emis2_s / 1000.0
    #
    #        infile = wdir + r'auxiliary\heihe_classif_new.tif'
    #        [classif, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        classif = np.asarray(classif, np.int)
    #        infile = wdir + r'auxiliary\ICBP_Emi.txt'
    #        emis_v = self.getDatafromTxt(infile, 0, 2)
    #        emis1_v = emis_v[classif, 0] / 1000.0
    #        # emis2_v = emis_v[classif, 1] / 1000.0
    #
    #        infile = wdir + 'heihe_result_night\heihe_' + dataTime + r'_lst_n.tif'
    #        [lst_n, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #        infile = wdir + 'heihe_result_night\heihe_' + dataTime + r'_lst_o.tif'
    #        [lst_o, temp1, temp2, temp3, temp4, temp5] = self.getTiffData(infile)
    #
    #        infile = wdir + resultName + dataTimeold + r'_ndvi_n.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_n = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_n[fvc_n > 1] = 1
    #        fvc_n[fvc_n < 0] = 0
    #        infile = wdir + resultName + dataTimeold + r'_ndvi_o.tif'
    #        dataset = gdal.Open(infile, gdal.GA_ReadOnly)
    #        data = dataset.GetRasterBand(1)
    #        ndvi = data.ReadAsArray()
    #        fvc_o = (ndvi - self.ndvi_min) / (self.ndvi_max - self.ndvi_min)
    #        fvc_o[fvc_o > 1] = 1
    #        fvc_o[fvc_o < 0] = 0
    #
    #        rad_n = self.planck(self.wl8, lst_n)
    #        rad_n = rad_n * (fvc_n * emis1_v + (1 - fvc_n) * emis1_s)
    #        rad_o = self.planck(self.wl8, lst_o)
    #        rad_o = rad_o * (fvc_o * emis1_v + (1 - fvc_o) * emis1_s)
    #
    #        # component temperature inversion
    #        W1 = np.asarray([emis1_s * (1 - fvc_n), emis1_v * (fvc_n)])
    #        W2 = np.asarray([emis1_s * (1 - fvc_o), emis1_v * (fvc_o)])
    #        result = np.zeros([2, nl, ns])
    #
    #        nsi = np.zeros([nl, ns])
    #        nli = np.zeros([nl, ns])
    #        for k in range(nl):
    #            nli[k, :] = k
    #            nsi[k, :] = np.linspace(0, ns - 1, ns)
    #
    #        puritypixe = 0.01
    #        for k1 in range(off, nl - off):
    #            nl1 = k1 - off
    #            nl2 = k1 + off + 1
    #            for k2 in range(off, ns - off):
    #                ns1 = k2 - off
    #                ns2 = k2 + off + 1
    #
    #                ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
    #                nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
    #
    #                temp1 = fvc_n[k1, k2]
    #                temp2 = fvc_o[k1, k2]
    #                temp3 = rad_n[k1, k2]
    #                temp4 = rad_o[k1, k2]
    #                if (temp3 < 6.0 or temp4 < 6.0): continue
    #                if (temp1 < puritypixe or temp2 < puritypixe or temp1 > (1 - puritypixe) or temp2 > (
    #                        1 - puritypixe)): continue
    #
    #                w1 = W1[:, nl1:nl2, ns1:ns2]
    #                w1 = np.reshape(w1, (2, -1))
    #                w2 = W2[:, nl1:nl2, ns1:ns2]
    #                w2 = np.reshape(w2, (2, -1))
    #
    #                b1 = np.reshape(rad_n[nl1:nl2, ns1:ns2], -1)
    #                b2 = np.reshape(rad_o[nl1:nl2, ns1:ns2], -1)
    #                n0 = np.ones(len(ns_temp))
    #
    #                x = np.asarray([n0, ns_temp, nl_temp, ns_temp * ns_temp, nl_temp * nl_temp, ns_temp * nl_temp])
    #                w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(6)])
    #                w = np.hstack((w11, w22))
    #                b = np.hstack((b1, b2))
    #                temp1 = np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1)
    #                temp2 = np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1)
    #                temp = np.hstack((temp1, temp2))
    #                ind = (b > 6.0) * (temp > puritypixe) * (temp < (1 - puritypixe))
    #                num_point = np.sum(ind)
    #                b = b[ind]
    #                w = w[:, ind]
    #                if (num_point < 15): continue
    #
    #                w = np.transpose(w)
    #                coeff = np.asarray(lstsq(w, b))[0]
    #                coeff = np.asarray(coeff)
    #                offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                ct1 = np.sum(coeff[:6] * x[:, offcenter])
    #                ct2 = np.sum(coeff[6:] * x[:, offcenter])
    #
    #                # w = np.transpose(w)
    #                # b = np.transpose(b)
    #                # u, s, vt = scipy.linalg.svd(w, 0)
    #                # v = np.transpose(vt)
    #                # I = np.diag(np.ones(12))
    #                # ss = s * s
    #                # s = np.diag(s)
    #                # ss = np.diag(ss)
    #                # s = np.matrix(s)
    #                # v = np.matrix(v)
    #                # vt = np.matrix(vt)
    #                # u = np.matrix(u)
    #                # b = np.matrix(b)
    #                # b = np.transpose(b)
    #                # bb = sorted(b)
    #                # halfpoint = np.int(num_point / 2)
    #                # rleaf = np.average(bb[:halfpoint])
    #                # rsoil = np.average(bb[halfpoint:])
    #                # coeffprior = np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0])
    #                # coeffprior = np.transpose(coeffprior)
    #                # temp = v * (np.linalg.inv(ss + I))
    #                # coeff = temp * (s * np.transpose(u) * b + vt * coeffprior)
    #                # coeff = np.asarray(coeff)
    #                # offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
    #                # ct1 = np.sum(coeff[:6, 0] * x[:, offcenter])
    #                # ct2 = np.sum(coeff[6:, 0] * x[:, offcenter])
    #
    #
    #                if ct1 <= 0 or ct2 <= 0: continue
    #                ct1 = self.invplanck0(self.wl8, ct1)
    #                ct2 = self.invplanck0(self.wl8, ct2)
    #
    #                result[:, k1, k2] = [ct1, ct2]
    #
    #        result[result < 0] = 0
    #        result[result > 350] = 0
    #
    #        outfile = wdir + r'heihe_result_night\heihe_' + dataTime + '_LSCT.tif'
    #        self.writeTiff(result, ns, nl, 1, geog, proj, outfile)
    #
    #        return 0

    def planck(self, wavelength, Ts):
        c1 = 11910.439340652
        c2 = 14388.291040407

        ind = (Ts < 100) *(Ts > 10)
        Ts[ind] = Ts[ind] + 273.15

        wavelength = np.float(wavelength)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
        return rad

    def planck0(self, wavelength, Ts):
        c1 = 11910.439340652
        c2 = 14388.291040407
        

        wavelength = np.float(wavelength)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
        return rad

    def invplanck(self, wavelength, rad):
        c1 = 11910.439340652 * 10000
        c2 = 14388.291040407
        temp = c1 / (rad * np.power((wavelength), 5)) + 1
        Ts = c2 / (wavelength * np.log(temp))
        return Ts

    def invplanck0(self, wavelength, rad):
        c1 = 11910.439340652 * 10000
        c2 = 14388.291040407
        temp = c1 / (rad * np.power((wavelength), 5)) + 1
        Ts = c2 / (wavelength * np.log(temp))
        return Ts

    def fromDate2DOY(self, year, month, day):
        total = 0
        days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        for index in range(month - 1):
            total += days_of_month[index]
        temp = (year // 4 == 0 and year // 100 != 0) or (year // 400 == 0)
        if month > 2 and temp:
            total += 1
        return total + day
