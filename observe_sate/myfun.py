import numpy as np
from scipy.interpolate import griddata
from scipy.linalg import lstsq
from osgeo import gdal
from osgeo import osr,ogr
import re
import scipy
import xlrd
import os
import cv2
import struct
import netCDF4

import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
from pylab import *
########################################################################################################################
### IO and File
def search(dir):
    results = os.listdir(dir)
    return results

def search_file(dir,specstr):
    results = []
    num = len(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        specstr0 = specstr[0]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) and
                    specstr0 in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x]
    return results

def search_dir(dir,specstr):
    results = []
    num = len(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) and
                    specstr in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x]
    return results

def search_file_rej(dir,specstr,rejstr):
    results = []
    num = len(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        specstr0 = specstr[0]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x))
                    and specstr0 in x
                    and rejstr not in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and rejstr not in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isfile(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x
                    and rejstr not in x]
    return results

def search_dir_rej(dir,specstr,rejstr):
    results = []
    num = len(specstr)
    if num == 0:
        results = os.listdir(dir)
    elif num==1:
        specstr0 = specstr[0]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x))
                    and specstr0 in x
                    and rejstr not in x]
    elif num==2:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and rejstr not in x]
    elif num==3:
        specstr1 = specstr[0]
        specstr2 = specstr[1]
        specstr3 = specstr[2]
        results += [x for x in os.listdir(dir) if
                    os.path.isdir(os.path.join(dir, x)) \
                    and specstr1 in x
                    and specstr2 in x
                    and specstr3 in x
                    and rejstr not in x]
    return results

def remove_file(indir,specstr):
    for x in os.listdir(indir):
        fp = os.path.join(indir, x)
        # 如果文件存在，返回true
        if re.search(specstr, x) is not None:
            print(fp)
            os.remove(fp)

def rename_file(indir,specstr):
    for x in os.listdir(indir):
        fp = os.path.join(indir, x)
        # print(x)
        # 如果文件存在，返回true
        if re.search(specstr, x) is not None:
            [filename, hz] = os.path.splitext(x)
            outfile = indir + filename + '_test.tif'
            print(fp)
            print(outfile)
            if os.path.exists(outfile) == 1:
                os.remove(outfile)
            os.rename(fp, outfile)

def move_file(infile,outfile):
    os.rename(infile,outfile)

##########################################################################################################
### File read and write

def read_txt_float(infile):
    mydata = []
    with open(infile) as f:
        lines = f.readline()
        while lines:
            line = lines.split()
            mydata.append(line)
            lines = f.readline()
    mydata = np.asarray(mydata,dtype=float)
    return mydata

def read_txt_str(infile):
    mydata = []
    with open(infile) as f:
        lines = f.readline()
        while lines:
            line = lines.split()
            mydata.append(line)
            lines = f.readline()
    return mydata

def read_txt_array(filename, num_pass, num_col):
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

def read_excel_sheet(filename,sheetname='Sheet1'):
    ExcelFile = xlrd.open_workbook(filename)
    ExcelFile.sheet_names()
    sheet = ExcelFile.sheet_by_name(sheetname)
    return sheet

def read_excel_sheet_col(filename,col,sheetname='Sheet1'):
    ExcelFile = xlrd.open_workbook(filename)
    ExcelFile.sheet_names()
    sheet = ExcelFile.sheet_by_name(sheetname)
    colvalue = sheet.col_values(col)
    return colvalue

def read_excel_sheet_row(filename,row,sheetname='Sheet1'):
    ExcelFile = xlrd.open_workbook(filename)
    ExcelFile.sheet_names()
    sheet = ExcelFile.sheet_by_name(sheetname)
    rowvalue = sheet.row_values(row)
    return rowvalue

### binary

def read_binary(infile):
    fin = open(infile,'rb')
    temp =[]
    while True:
        fileContent = fin.read(4)
        num = len(fileContent)
        if num !=4:
            break
        tp = struct.unpack('l',fileContent)
        temp.extend(tp)
    fin.close()
    temp = np.array(temp,dtype=np.int)
    return temp

##### image
def read_image_gdal(fileName):
    dataset = gdal.Open(fileName)
    if dataset == None:
        print(fileName + "文件无法打开")
        return
    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize
    im_bands = dataset.RasterCount
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height)
    im_geotrans = dataset.GetGeoTransform()
    im_proj = dataset.GetProjection()
    return im_data,im_width,im_height,im_bands,im_geotrans,im_proj

def read_image_raw(infile,ns,nl,nb):
    fb = open(infile,'rb')
    mydata = np.zeros((ns,nl,nb))
    for kb in range(nb):
        for ks in range(ns):
            for kl in range(nl):
                arr = fb.read(2)
                elem = struct.unpack('h',arr)[0]
                mydata[ks][kl][kb] = elem
    return mydata

def read_image_Nc_group(fileName, groupName, objectName,ifscale):
    dataset = netCDF4.Dataset(fileName)
    if ifscale ==0:
        dataset.groups[groupName].variables[objectName].set_auto_maskandscale(False)
    predata = np.asarray(dataset.groups[groupName].variables[objectName][:])
    return predata

def getDatafromNc(fileName, objectName,ifscale):
    dataset = netCDF4.Dataset(fileName)
    if ifscale ==0:
        dataset.variables[objectName].set_auto_maskandscale(False)
    predata = np.asarray(dataset.variables[objectName][:])
    return predata

def write_image_gdal(im_data, im_width, im_height, im_bands, im_trans, im_proj, path):
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
    if (dataset != None and im_trans != '' and im_proj != ''):
        dataset.SetGeoTransform(im_trans)
        dataset.SetProjection(im_proj)
    for i in range(im_bands):
        dataset.GetRasterBand(i+1).WriteArray(im_data[i])
    del dataset


###########################################################################################################################
### Array and Map

rd = np.pi/180.0


def resize_data_ls(preArray, nl, ns):
    ns = np.int(ns)
    nl = np.int(nl)
    data = cv2.resize(preArray,(ns,nl),interpolation=cv2.INTER_LINEAR)
    return data

def resize_data_near(preArray, nl, ns):
    ns = np.int(ns)
    nl = np.int(nl)
    data = cv2.resize(preArray,(ns,nl),interpolation=cv2.INTER_NEAREST)
    return data

def resize_data(preArray, bs):
    [pre_nl,pre_ns] = np.shape(preArray)

    ns = np.int(pre_ns*bs)
    nl = np.int(pre_nl*bs)
    data = cv2.resize(preArray,(ns,nl),interpolation=cv2.INTER_LINEAR)
    return data

def resize_data_file(infile, ns, nl, outfile):
    # info=gdal.WarpOptions(width=ns,height=nl)
    test = gdal.Warp(outfile,infile,width=ns,height=nl)
    return 1

def getSRSPair(dataset):
    '''
    获得给定数据的投影参考系和地理参考系
    :param dataset: GDAL地理数据
    :return: 投影参考系和地理参考系
    '''
    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    return prosrs, geosrs

def geo2lonlat(dataset, x, y):
        '''
        将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
        :param dataset: GDAL地理数据
        :param x: 投影坐标x
        :param y: 投影坐标y
        :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
        '''
        prosrs, geosrs = getSRSPair(dataset)
        ct = osr.CoordinateTransformation(prosrs, geosrs)
        x = np.reshape(x, [-1])
        y = np.reshape(y, [-1])
        temp = np.asarray([x, y])
        temp = np.transpose(temp)
        coords = np.asarray(ct.TransformPoints(temp))
        return coords[:,0],coords[:,1]

def lonlat2geo(dataset,lat,lon):
    '''
        将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
        :param dataset: GDAL地理数据
        :param lon: 地理坐标lon经度
        :param lat: 地理坐标lat纬度
        :return: 经纬度坐标(lon, lat)对应的投影坐标
    '''
    # dataset = gdal.Open(fileName, gdal.GA_ReadOnly)
    prosrs, geosrs = getSRSPair(dataset)
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    lon = np.reshape(lon,[-1])
    lat = np.reshape(lat,[-1])
    temp = np.asarray([lon,lat])
    temp = np.transpose(temp)
    # temp = np.asarray([lat[0:2],lon[0:2]])
    coords = np.asarray(ct.TransformPoints(temp))

    return coords[:,0],coords[:,1]

def geo2imagexy(dataset, x, y):
    '''
    根据GDAL的六 参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
    :param dataset: GDAL地理数据
    :param x: 投影或地理坐标x
    :param y: 投影或地理坐标y
    :return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col) nl ns
    '''
    trans = dataset.GetGeoTransform()
    a = np.array([[trans[1], trans[2]], [trans[4], trans[5]]])
    b = np.array([x - trans[0], y - trans[3]])
    return np.linalg.solve(a, b)  # 使用numpy的linalg.solve进行二元一次方程的求解

def imagexy2geo(dataset, row,col):
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


def calc_azimuth(lat1, lon1, lat2, lon2):
    lat1_rad = lat1 * np.pi / 180
    lon1_rad = lon1 * np.pi / 180
    lat2_rad = lat2 * np.pi / 180
    lon2_rad = lon2 * np.pi / 180

    y = np.sin(lon2_rad - lon1_rad) * np.cos(lat2_rad)
    x = np.cos(lat1_rad) * np.sin(lat2_rad) - \
        np.sin(lat1_rad) * np.cos(lat2_rad) * np.cos(lon2_rad - lon1_rad)

    brng = np.arctan2(y, x) * 180 / np.pi

    return float((brng + 360.0) % 360.0)


### xscale = ns
### yscale = nl
def write_vrt(infile, datafile, xfile, yfile, xscale,yscale):

    f = open(infile, 'w')
    f.write(r'<VRTDataset rasterXSize="%d"'%xscale + 'rasterYSize="%d">'%yscale)
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


###########################################################################################################################
#### simple sci support
def planck(wavelength, Ts):
    c1 = 11910.439340652
    c2 = 14388.291040407
    if isinstance(Ts * 1.0, float):
        if (Ts < 100): Ts = Ts + 273.15
        wavelength = np.float(wavelength)
        Ts = np.float(Ts)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    else:
        Ts[Ts < 100] = Ts[Ts < 100] + 273.15
        wavelength = np.float(wavelength)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    return rad


def inv_planck(wavelength, rad):
    c1 = 11910.439340652 * 10000
    c2 = 14388.291040407

    temp = c1 / (rad * np.power((wavelength), 5)) + 1
    Ts = c2 / (wavelength * np.log(temp))
    return Ts

def date2DOY(year,month,day):
    days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    monthsum = np.zeros(12)

    if isinstance(day * 1.0, float):
        total = 0
        for index in range(month - 1):
            total += days_of_month[index]
        temp = (year // 4 == 0 and year // 100 != 0) or (year // 400 == 0)
        if month > 2 and temp:
            total += 1
        return total + day
    else:
        for index in range(1,12):
            monthsum[index] = monthsum[index-1]+days_of_month[index-1]
        month = np.asarray(month,dtype=np.int)
        DOY = monthsum[month-1] + day
        ind = ((year // 4 == 0) * (year // 100 != 0)) * (month>2)
        DOY[ind] = DOY[ind] + 1
        ind = ((year // 400 == 0))
        DOY[ind] = DOY[ind] + 1
    return DOY


def doy2date(year, doy):
    month_leapyear = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    month_notleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        for i in range(0, 12):
            if doy > month_leapyear[i]:
                doy -= month_leapyear[i]
            continue
            if doy <= month_leapyear[i]:
                month = i + 1
                day = doy
                break
    else:
        for i in range(0, 12):
            if doy > month_notleap[i]:
                doy -= month_notleap[i]
                continue
            if doy <= month_notleap[i]:
                month = i + 1
                day = doy
                break
    return month, day


###########################################################################################################################
### Plot


def plt_scatter(data1,data2,min1 = 275,min2 = 275,max1 = 335,max2 = 335):

    ind = (data1>min1)*(data1<max1)*(data2>min2)*(data2<max2)
    ind = np.where(ind > 0)

    dif = data2[ind]-data1[ind]
    rmse = np.sqrt(np.mean(dif*dif))
    std = np.std(dif)
    bias = np.mean(dif)
    r = np.corrcoef(data1[ind],data2[ind])
    r2 = r[0,1]*r[0,1]

    plt.plot(data1,data2,'ko',markersize=2.5)
    plt.plot([min1,max1],[min2,max2],'k-.')
    plt.text(data1,data2,)
    plt.title([rmse,r2])
    plt.xlim([min1,max1])
    plt.ylim([min2,max2])
    plt.text(min1+(max1-min1)*4/5, min2+(max2-min2)*1/5, "$RMSE$ = %2.2f$\degree$C" % rmse + "\n$Bias$ = %2.2f$\degree$C" % bias +
             "\n$r^2$ = %2.2f" % r2 + "\n$\sigma$ = %2.2f$\degree$C" % std, fontsize=14)
    plt.show()
    return 0


def plt_hist(img,img_,img_no):
    ind = (img_no > 0) * (img > 0) * (img_ > 0)
    dif0 = img[ind] * 3.1415 - img_[ind] * 3.1415
    dif = img_no[ind] * 3.1415 - img[ind] * 3.1415
    rmse = np.sqrt(np.mean(dif * dif))
    rmse0 = np.sqrt(np.mean(dif0))
    bias = np.mean(dif)
    bias0 = np.mean(dif0)
    print('nir', rmse0, rmse)
    print('nir', bias0, bias)

    kwargs = dict(histtype='stepfilled', alpha=0.5, bins=30)
    fig, axs = plt.subplots(ncols=1, figsize=(5, 4))
    plt.hist(dif, **kwargs, label='$\Delta BRF_{nir}$', color='orange')
    plt.legend()
    plt.text(-0.05, 650, 'RMSE = %4.3f\n' % rmse + r'Bias = %4.3f' % bias, fontsize=12)
    plt.xlabel('BRF Difference', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.show()

def plt_scatter_hist(data1,data2):

    dif = data2 - data1
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif * dif))
    rr = np.corrcoef(data1,data2)
    r2 = rr[1,0]*rr[1,0]

    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['figure.figsize'] = (4.0, 3.2)
    hist2d(data1, data2, bins=66, norm=LogNorm(),cmap='jet')
    plt.plot([0, 1], [0, 1], 'k-')
    plt.xlim([0,0.75])
    plt.ylim([0,0.75])
    plt.xlabel('LESS-simulated BRF ')
    plt.xlabel('model-simulated BRF ')
    plt.text(0.4,0.01,
             '$RMSE = %4.3f$\n'%rmse+
             '$Bias = %4.3f$\n'%bias+
             '$R^2 = %4.3f$'%r2,
             fontsize = 12)
    colorbar()
    show()


def plt_pie(labels,sizes,colors):

    labers = ['Raytran', 'DART', 'RGM', 'FLIGHT', 'librat', 'FliES']
    sizes1 = [131, 200, 93, 180, 84, 72]
    sizes2 = [3789, 5138, 2232, 5623, 2702, 1476]
    colors = ['red', 'blue', 'green', 'yellow', 'cyan', 'orange']
    explode = 0, 0, 0, 0, 0, 0
    patches, l_text, p_text = plt.pie(sizes1, explode=explode, labels=labers,
                                      colors=colors, autopct='%1.1f%%', shadow=False, startangle=50)
    plt.axis('equal')
    plt.show()

def plt_coutourPolar(theta,rho,z,n,maxrho):

    # ###transform data to Cartesian coordinates.
    # delta = maxrho/50
    # xx = rho*np.cos((90-theta)*(np.pi/180))
    # yy = rho*np.sin((90-theta)*(np.pi/180))
    # xi = np.linspace(-maxrho,maxrho,2*maxrho/delta)
    # yi = xi
    # [xi,yi] = np.meshgrid(xi,yi)
    # zi = griddata((xx,yy),z,(xi,yi),'cubic')
    # # fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    # plt.contourf(xi,yi,zi,n,cmap='jet')
    # plt.show()

    delta = maxrho/50
    xx = np.radians(theta)
    yy = rho
    xi = np.radians(np.arange(0,365,5))
    yi = np.arange(0,maxrho,5)
    [xi,yi] = np.meshgrid(xi,yi)
    # xi = np.radians(xi)
    zi = griddata((xx,yy),z,(xi,yi))
    fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    plt.autumn()
    cax = ax.contourf(xi, yi, zi, n,cmap='jet')
    plt.autumn()
    cb = fig.colorbar(cax)
    # cb.set_label("Pixel reflectance")
    plt.show()
