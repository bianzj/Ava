import numpy as np
import matplotlib.pylab as plt
from scipy.interpolate import griddata
from scipy.linalg import lstsq
from osgeo import gdal
from osgeo import osr,ogr
import re
import scipy
import xlrd
import os
import cv2
import netCDF4

def sci_planck( wavelength, Ts):
    c1 = 11910.439340652
    c2 = 14388.291040407
   
    if isinstance(Ts*1.0,float) :
        if (Ts < 100): Ts = Ts + 273.15
        wavelength = np.float(wavelength)
        Ts = np.float(Ts)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    else:
        Ts[Ts<100] = Ts[Ts<100] + 273.15
        wavelength = np.float(wavelength)
        rad = c1 / (np.power(wavelength, 5) * (np.exp(c2 / Ts / wavelength) - 1)) * 10000
    return rad

def sci_invplanck(wavelength,rad):
    c1 = 11910.439340652 * 10000 
    c2 = 14388.291040407

    temp = c1 / (rad * np.power((wavelength),5)) + 1
    Ts = c2 / (wavelength * np.log(temp))
    return Ts

def ope_resizeData(preArray,ns,nl):
    ns = np.int(ns)
    nl = np.int(nl)
    data = cv2.resize(preArray,(ns,nl),interpolation=cv2.INTER_LINEAR)
    return data

def ope_resizeData_(preArray,nb,ns,nl):
    ns = np.int(ns)
    nl = np.int(nl)
    result = np.zeros([nb,ns,nl])
    for k in range(nb):
        result[k,:,:] = cv2.resize(preArray[k,:,:],(ns,nl),interpolation=cv2.INTER_LINEAR)
    return result

def ope_resizeAve(preArray,ns,nl):
    [preNb,preNs,preNl] = np.shape(preArray)
    ns = np.int(ns)
    nl = np.int(nl)
    stepx = np.int(preNs/ns)
    stepy = np.int(preNl/nl)
    nb = preNb
    result = np.zeros([nb,ns,nl])
    for k in range(nb):
        for ks in range(ns):
            for kl in range(nl):
                offs = ks*stepx
                ends = ks*stepx+stepx
                offl = kl*stepy
                endl = kl*stepy+stepy
                if ends > preNs: ends = preNs
                if endl > preNl: endl = preNl
                result[k,ks,kl] = np.average(preArray[k,offs:ends,offl:endl])
    return result

def ope_resizeData_slow(preArray,ns,nl):
    ns = np.int(ns)
    nl = np.int(nl)
    data = cv2.resize(preArray,(ns,nl),interpolation=cv2.INTER_LINEAR)
    return data

def ope_resizeImage(infile,outfile,ns,nl):
    # info=gdal.WarpOptions(width=ns,height=nl)   
    test = gdal.Warp(outfile,infile,width=ns,height=nl)
    return 1

def ope_search0(dir):
    results = os.listdir(dir)
    return results

def ope_search1(dir,specstr):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x)) and specstr in x]  
    return results

def ope_search2(dir,specstr1,sepcstr2):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x)) \
        and specstr1 in x and sepcstr2 in x]  
    return results

def ope_search3(dir,specstr1,specstr2,specstr3):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x)) \
        and specstr1 in x and specstr2 in x and specstr3 in x]  
    return results

def ope_search2_rej1(dir,specstr1,specstr2,rejstr):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x)) \
        and specstr1 in x and specstr2 in x and rejstr not in x]  
    return results

def ope_search1_rej1(dir,specstr1,rejstr):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x))==0\
        and specstr1 in x and rejstr not in x]
    return results

def ope_searchfile_accept1_rej1(dir,specstr1,rejstr):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x))==0 \
        and specstr1 in x and rejstr not in x]
    return results

def ope_search2_rej2(dir,specstr1,specstr2,rejstr1,rejstr2):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x)) \
        and specstr1 in x and specstr2 in x and rejstr1 not in x and rejstr2 not in x]  
    return results

def ope_search3_rej1(dir,specstr1,specstr2,specstr3,rejstr):
    results = []
    results += [x for x in os.listdir(dir) if os.path.isfile(os.path.join(dir, x)) \
        and specstr1 in x and specstr2 in x and specstr3 in x and rejstr not in x]  
    return results

def ope_getSRSPair(dataset):
    '''
    获得给定数据的投影参考系和地理参考系
    :param dataset: GDAL地理数据
    :return: 投影参考系和地理参考系
    '''
    prosrs = osr.SpatialReference()
    prosrs.ImportFromWkt(dataset.GetProjection())
    geosrs = prosrs.CloneGeogCS()
    return prosrs, geosrs

def ope_geo2lonlat(dataset, x, y):
        '''
        将投影坐标转为经纬度坐标（具体的投影坐标系由给定数据确定）
        :param dataset: GDAL地理数据
        :param x: 投影坐标x
        :param y: 投影坐标y
        :return: 投影坐标(x, y)对应的经纬度坐标(lon, lat)
        '''
        prosrs, geosrs = ope_getSRSPair(dataset)
        ct = osr.CoordinateTransformation(prosrs, geosrs)
        coords = ct.TransformPoint(x, y)
        return coords[:2]

def ope_lonlat2geo(dataset,lat,lon):
    '''
        将经纬度坐标转为投影坐标（具体的投影坐标系由给定数据确定）
        :param dataset: GDAL地理数据
        :param lon: 地理坐标lon经度
        :param lat: 地理坐标lat纬度
        :return: 经纬度坐标(lon, lat)对应的投影坐标
    '''
    # dataset = gdal.Open(fileName, gdal.GA_ReadOnly)
    prosrs, geosrs = ope_getSRSPair(dataset)
    ct = osr.CoordinateTransformation(geosrs, prosrs)
    lon = np.reshape(lon,[-1])
    lat = np.reshape(lat,[-1])
    temp = np.asarray([lon,lat])
    temp = np.transpose(temp)
    # temp = np.asarray([lat[0:2],lon[0:2]])
    coords = np.asarray(ct.TransformPoints(temp))

    return coords[:,0],coords[:,1]

def ope_geo2imagexy(dataset, x, y):
    '''
    根据GDAL的六 参数模型将给定的投影或地理坐标转为影像图上坐标（行列号）
    :param dataset: GDAL地理数据
    :param x: 投影或地理坐标x
    :param y: 投影或地理坐标y
    :return: 影坐标或地理坐标(x, y)对应的影像图上行列号(row, col)
    '''
    trans = dataset.GetGeoTransform()
    a = np.array([[trans[2], trans[1]], [trans[5], trans[4]]])
    b = np.array([x - trans[0], y - trans[3]])
    return np.linalg.solve(a, b)  # 使用numpy的linalg.solve进行二元一次方程的求解

def ope_imagexy2geo(dataset, row,col):
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

def ope_fileRemote(indir,specstr):
    for x in os.listdir(indir):
        fp = os.path.join(indir, x)
        # 如果文件存在，返回true
        if re.search(specstr, x) is not None:
            print(fp)
            os.remove(fp)

def ope_fileRename(indir,bt):
    for x in os.listdir(indir):
        fp = os.path.join(indir, x)
        # print(x)
        # 如果文件存在，返回true
        if re.search(bt, x) is not None:
            [filename, hz] = os.path.splitext(x)
            outfile = indir + filename + '_test.tif'
            print(fp)
            print(outfile)
            if os.path.exists(outfile) == 1:
                os.remove(outfile)
            os.rename(fp, outfile)

def ope_date2DOY1(year,month,day):
    days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    monthsum = np.zeros(12)

    for index in range(1,12):
        monthsum[index] = monthsum[index-1]+days_of_month[index-1]

    month = np.asarray(month,dtype=np.int)
    DOY = monthsum[month-1] + day

    ind = ((year // 4 == 0) * (year // 100 != 0)) * (month>2)
    DOY[ind] = DOY[ind] + 1
    ind = ((year // 400 == 0))
    DOY[ind] = DOY[ind] + 1

    return DOY

def ope_date2DOY0(year, month, day):
    total = 0
    days_of_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    for index in range(month - 1):
        total += days_of_month[index]
    temp = (year // 4 == 0 and year // 100 != 0) or (year // 400 == 0)
    if month > 2 and temp:
        total += 1
    return total + day

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

def writeTiff(im_data, im_width, im_height, im_bands, im_geotrans, im_proj, path):
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

def getDatafromNc(fileName, objectName):
    dataset = netCDF4.Dataset(fileName)
    predata = np.asarray(dataset.variables[objectName][:])
    return predata

def getTiffData(fileName):
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
    return im_data,im_width,im_height,im_bands,im_geotrans,im_proj

def plt_plot(data1,data2):
    min1 = 295
    min2 = 295
    max1 = 315
    max2 = 315
    ind = (data1>min1)*(data1<max1)*(data2>min2)*(data2<max2)
    ind = np.where(ind > 0)
    dif = data2[ind]-data1[ind]
    rmse = np.sqrt(np.mean(dif*dif))
    bias = np.mean(dif)
    r = np.corrcoef(data1[ind],data2[ind])
    r2 = r[0,1]*r[0,1]
    plt.plot(data1,data2,'ko',markersize=2.5)
    plt.plot([min1,max1],[min2,max2],'k-.')
    plt.title([rmse,r2])
    plt.xlim([min1,max1])
    plt.ylim([min2,max2])
    plt.show()
    return 0 




def read_gdalTiff(fileName):
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
    return im_data,im_width,im_height,im_bands,im_geotrans,im_proj

def read_txt2array(infile):
    mydata = []
    with open(infile) as f:
        lines = f.readline()
        while lines:
            line = lines.split()
            mydata.append(line)
            lines = f.readline()
    mydata = np.asarray(mydata,dtype=float)
    return mydata

def read_txt2str(infile,spec=' '):
    mydata = []
    with open(infile) as f:
        lines = f.readline()
        while lines:
            line = lines.split(spec)
            mydata.append(line)
            lines = f.readline()
    return mydata

def read_excel2array(filename,sheetname='Sheet1'):
    ExcelFile = xlrd.open_workbook(filename)
    ExcelFile.sheet_names()
    sheet = ExcelFile.sheet_by_name(sheetname)
    return sheet

def read_excelline2array(filename,col, sheetname='Sheet1'):
    ExcelFile = xlrd.open_workbook(filename)
    ExcelFile.sheet_names()
    sheet = ExcelFile.sheet_by_name(sheetname)
    colvalue = sheet.col_values(col)
    return colvalue

def write_gdalTiff(im_data, im_width, im_height, im_bands, im_geotrans, im_proj, path):
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


from sklearn.linear_model import LinearRegression


def LSFLI(sza, vza, raa, dbt):
    ind = raa > 180
    raa[ind] = 360 - raa[ind]
    thetas = np.deg2rad(sza)
    thetaa = np.deg2rad(raa)
    thetav = np.deg2rad(abs(vza))

    kz = np.zeros(np.size(dbt))
    kz[:] = 10

    az = np.cos(thetas) * np.cos(thetav) + \
         np.sin(thetas) * np.sin(thetav) * np.cos(thetaa)

    DD = (np.tan(thetas) ^ 2 + np.tan(thetav) ^ 2 - 2 * np.tan(thetas) * np.tan(thetav) *
          np.cos(thetaa))
    D = np.sqrt(max(0, DD))

    cost = np.sqrt(D * D + (np.tan(thetav) * np.tan(thetas) * np.sin(thetaa)) ^ 2) \
           / (np.sec(thetas) + np.sec(thetav))
    cost = max(-1, min(1, cost))

    t = np.acos(cost)
    O = (1 / np.pi) * (t - np.sin(t) * np.cos(t)) * (np.sec(thetas) + np.sec(thetav))

    Klsf = ((1 + 2 * np.cos(thetav)) / ((0.96) ^ 0.5 + 2 * 0.96 * np.cos(thetav)) - 0.25 * np.cos(thetav)
            / (1 + 2 * np.cos(thetav)) + 0.15 * (1 - np.exp(-0.75 / np.cos(thetav))))
    Kli = (1 + az) * (np.sec(thetas) * (np.sec(thetav))) / (np.sec(thetav) + 1 / np.cos(thetas) - O) - 2

    A = np.asarray([kz, Klsf, Kli])

    model = LinearRegression(copy_X=True, fit_intercept=True, n_jobs=1, normalize=False)
    model.fit(A, dbt)
    result = model.predict(A)
    dif = result - dbt
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif * dif))

    return result, rmse, bias



#########################################################
def ope_reverse(preArray,nb,nl,ns,offnb):
    result = preArray*1.0
    for k in range(offnb,nb):
        matrix = preArray[k,:,:]
        # for i in range(nl //2 ):
        #     matrix[i], matrix[nl - 1 - i] = matrix[nl - 1 - i], matrix[i]
        matrix = np.flipud(matrix)
        matrix = np.fliplr(matrix)
        # for m in matrix:
        #     for j in range(ns // 2):
        #         m[j], m[ns - 1 - j] = m[ns - 1 - j], m[j]

        result[k,:,:] = matrix
    return result


def writeTiff(im_data, im_width, im_height, im_bands, im_geotrans, im_proj, path):
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
        dataset.GetRasterBand(i+1).WriteArray(im_data[i])
    del dataset


