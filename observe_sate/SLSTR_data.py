
from myfun_image import *
from myfun_sci import *
from myfun_file import *
from SLSTR_constant import *
import netCDF4
import time

#####################################
###



def ifallfilesexists(fileDir):
    infile = fileDir + r'\S2_radiance_an.nc'
    c1 = os.path.exists(infile)
    infile = fileDir + r'\geodetic_in.nc'
    c2 = os.path.exists(infile)
    infile = fileDir + r'\met_tx.nc'
    c4 = os.path.exists(infile)
    infile = fileDir + r'\S6_radiance_an.nc'
    c5 = os.path.exists(infile)
    infile = fileDir + r'\S8_BT_in.nc'
    c6 = os.path.exists(infile)
    return c1 * c2 * c4 * c5 * c6

def getDatafromNc(fileName, objectName):
    dataset = netCDF4.Dataset(fileName)
    predata = np.asarray(dataset.variables[objectName][:])
    return predata

def resize_file(infile, outfile, f=2):
    [preArray, im_width, im_height, im_bands, im_geotrans, im_proj] = read_image_gdal(infile)
    ns = im_width * f
    nl = im_height * f
    data = cv2.resize(preArray, (ns, nl), interpolation=cv2.INTER_BITS)
    write_image_gdal(data, ns, nl, 1, im_geotrans, im_proj, outfile)
    return data

def resizefromSmall2Image(ifnadir, preArray):
    if (ifnadir == 1):
        newArray = cv2.resize(preArray[:, adjust_xmin_n_g:adjust_xmax_n_g],
                              (ns_nadir_g, nl_nadir_g), cv2.INTER_NEAREST)
    else:
        newArray = cv2.resize(preArray[:, adjust_xmin_o_g:adjust_xmax_o_g],
                              (ns_obliq_g, nl_obliq_g), cv2.INTER_NEAREST)
    return newArray

def getVNIRdatafromNc(fileDir):

    ns = ns_nadir_g * 2
    nl = nl_nadir_g * 2

    # get red data
    fileName = fileDir + '\\S2_radiance_an.nc'
    objectName = r'S2_radiance_an'
    rad_red = getDatafromNc(fileName, objectName)
    fileName = fileDir + '\\S2_quality_an.nc'
    objectName = r'S2_solar_irradiance_an'
    sr_red = getDatafromNc(fileName, objectName)
    refl_red = rad_red / (sr_red[0] / np.pi)

    # get nir data
    fileName = fileDir + '\\S3_radiance_an.nc'
    objectName = r'S3_radiance_an'
    rad_nir = getDatafromNc(fileName, objectName)
    fileName = fileDir + '\\S3_quality_an.nc'
    objectName = r'S3_solar_irradiance_an'
    sr_nir = getDatafromNc(fileName, objectName)
    refl_nir = rad_nir / (sr_nir[0] / np.pi)


    outfile = fileDir + '\\red.tif'
    write_image_gdal(refl_red[:nl, :ns], ns, nl, 1, '', '', outfile)

    outfile = fileDir + '\\nir.tif'
    write_image_gdal(refl_nir[:nl, :ns], ns, nl, 1, '', '', outfile)

    return 1

def getGeodatafromNc(fileDir):
    # get red data

    ns = ns_nadir_g
    nl = nl_nadir_g

    fileName = fileDir + '\\geometry_tn.nc'
    objectName = r'sat_zenith_tn'
    vza = getDatafromNc(fileName, objectName)
    objectName = r'solar_zenith_tn'
    sza = getDatafromNc(fileName, objectName)
    objectName = r'sat_azimuth_tn'
    vaa = getDatafromNc(fileName, objectName)
    objectName = r'solar_azimuth_tn'
    saa = getDatafromNc(fileName, objectName)
    psi = np.abs(vaa - saa)

    vza = resizefromSmall2Image(1, vza)
    sza = resizefromSmall2Image(1, sza)
    psi = resizefromSmall2Image(1, psi)
    psi[psi > 180] = 360 - psi[psi > 180]

    fileName = fileDir + '\\geodetic_in.nc'
    objectName = r'latitude_in'
    lat = getDatafromNc(fileName, objectName)

    fileName = fileDir + '\\geodetic_in.nc'
    objectName = r'longitude_in'
    lon = getDatafromNc(fileName, objectName)


    outfile = fileDir + '\\lat.tif'
    write_image_gdal(lat[:nl, :ns], ns, nl, 1, '', '', outfile)

    lat2 = resize_data(lat[:nl,:ns],2)
    outfile = fileDir + '\\lat_v.tif'
    write_image_gdal(lat2, ns*2, nl*2, 1, '', '', outfile)

    outfile = fileDir + '\\lon.tif'
    write_image_gdal(lon[:nl, :ns], ns, nl, 1, '', '', outfile)

    lon2 = resize_data(lon[:nl,:ns],2)
    outfile = fileDir + '\\lon_v.tif'
    write_image_gdal(lon2, ns*2, nl*2, 1, '', '', outfile)

    outfile = fileDir + '\\vza.tif'
    write_image_gdal(vza[:nl, :ns], ns, nl, 1, '', '', outfile)
    outfile = fileDir + '\\sza.tif'
    write_image_gdal(sza[:nl, :ns], ns, nl, 1, '', '', outfile)
    outfile = fileDir + '\\psi.tif'
    write_image_gdal(psi[:nl, :ns], ns, nl, 1, '', '', outfile)

    return 1

def getTIRdatafromNc(fileDir):
    # get red data
    fileName = fileDir + '\\S8_BT_in.nc'
    objectName = r'S8_BT_in'
    BT8 = getDatafromNc(fileName, objectName)
    BT8 = BT8[:nl_nadir_g, :]

    # get nir data
    fileName = fileDir + '\\S9_BT_in.nc'
    objectName = r'S9_BT_in'
    BT9 = getDatafromNc(fileName, objectName)
    BT9 = BT9[:nl_nadir_g, :]


    outfile = fileDir + '\\bt8.tif'
    write_image_gdal(BT8, ns_nadir_g, nl_nadir_g, 1, '', '', outfile)

    outfile = fileDir + '\\bt9.tif'
    write_image_gdal(BT9, ns_nadir_g, nl_nadir_g, 1, '', '', outfile)

    return 1

def writevrt_nadir(infile, datafile, xfile, yfile):
    write_vrt(infile, datafile, xfile, yfile, 1500, 1200)
    return 1

def writevrt_obliq(infile, datafile, xfile, yfile):
    write_vrt(infile, datafile, xfile, yfile, 900, 1200)
    return 1

def writevrt_v_nadir(infile, datafile, xfile, yfile):
    write_vrt(infile, datafile, xfile, yfile, 3000, 2400)
    return 1

def writevrt_v_obliq(infile, datafile, xfile, yfile):
    write_vrt(infile,datafile,xfile,yfile,1800,2400)
    return 1

###



##############################################
######## Examples
##############################################

indir_raw = r'G:\s3b_raw\\'
indir_med = r'G:\s3b_med\\'
indir_tif = r'G:\s3b_tif\\'
indir_lst = r'G:\s3b_lst\\'
symbol = 'S3B'

# indir_raw = r'G:\s3a_raw\\'
# indir_med = r'G:\s3a_med\\'
# indir_tif = r'G:\s3a_tif\\'
# indir_lst = r'G:\s3a_lst\\'
# symbol = 'S3A'

ifunzipNC = 0
ifproj = 0
ifresample = 1
targetArea = r'G:\base\type.tif'

###############################
### from NC to tif images
###############################
if ifunzipNC ==1:
    fileNames = search_dir_rej(indir_raw,symbol,'zip')
    fileNum = len(fileNames)
    for k in range(0,fileNum):
        fileName = fileNames[k].strip()
        fileDir = indir_raw + fileName
        print(k,fileDir)
        ccc = ifallfilesexists(fileDir)
        if(ccc==0):continue
        if(os.path.exists(fileDir+'/lon_v.tif')==1):
            continue
        getVNIRdatafromNc(fileDir)
        getTIRdatafromNc(fileDir)
        getGeodatafromNc(fileDir)


###############################
### tif images projected to WGS 84 and resamppled to target study area
###############################

if ifproj ==1:
    fileNames = search_dir_rej(indir_raw, symbol, 'zip')
    fileNum = len(fileNames)
    for k in range(0,fileNum):
        fileName = fileNames[k].strip()
        fileDir = indir_raw + fileName
        print(k,fileDir)
        ccc = ifallfilesexists(fileDir)
        if(ccc==0):continue

        date_str = fileName[16:8 + 16]
        year = np.int(np.uint(date_str[:4]))
        month = np.int(np.uint(date_str[4:6]))
        day = np.int(np.uint(date_str[6:8]))
        doy = date2DOY(year,month,day)
        time_str = fileName[25:25 + 6]
        hh = np.int(np.uint(time_str[0:2]))
        mm = np.int(np.uint(time_str[2:4]))
        outfile0 = symbol+'_%4d'%year+'%03d'%doy+'_%02d'%hh+'%02d'%mm+'_'

        file ='bt8'
        outfile = indir_med + '\\' + outfile0 +file+'_proj.tif'
        if(os.path.exists(outfile)==1):
            continue


        file = 'psi'
        datafile = fileDir+'\\'+file+'.tif'
        xfile = fileDir+'\\lon.tif'
        yfile = fileDir+'\\lat.tif'
        vrtfile = fileDir+'\\'+file+'.vrt'
        writevrt_nadir(vrtfile,datafile,xfile,yfile)
        outfile = indir_med+'\\'+outfile0+file+'_proj.tif'
        dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

        file = 'sza'
        datafile = fileDir+'\\'+file+'.tif'
        xfile = fileDir+'\\lon.tif'
        yfile = fileDir+'\\lat.tif'
        vrtfile = fileDir+'\\'+file+'.vrt'
        writevrt_nadir(vrtfile,datafile,xfile,yfile)
        outfile = indir_med + '\\' + outfile0 +file+'_proj.tif'
        dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

        file = 'vza'
        datafile = fileDir+'\\'+file+'.tif'
        xfile = fileDir+'\\lon.tif'
        yfile = fileDir+'\\lat.tif'
        vrtfile = fileDir+'\\'+file+'.vrt'
        writevrt_nadir(vrtfile,datafile,xfile,yfile)
        outfile = indir_med + '\\' + outfile0 +file+'_proj.tif'
        dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

        file = 'bt8'
        datafile = fileDir+'\\'+file+'.tif'
        xfile = fileDir+'\\lon.tif'
        yfile = fileDir+'\\lat.tif'
        vrtfile = fileDir+'\\'+file+'.vrt'
        writevrt_nadir(vrtfile,datafile,xfile,yfile)
        outfile = indir_med + '\\' + outfile0 +file+'_proj.tif'
        dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

        file = 'bt9'
        datafile = fileDir+'\\'+file+'.tif'
        xfile = fileDir+'\\lon.tif'
        yfile = fileDir+'\\lat.tif'
        vrtfile = fileDir+'\\'+file+'.vrt'
        writevrt_nadir(vrtfile,datafile,xfile,yfile)
        outfile = indir_med + '\\' + outfile0 +file+'_proj.tif'
        dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)


        file = 'red'
        datafile = fileDir+'\\'+file+'.tif'
        xfile = fileDir+'\\lon_v.tif'
        yfile = fileDir+'\\lat_v.tif'
        vrtfile = fileDir+'\\'+file+'.vrt'
        writevrt_v_nadir(vrtfile,datafile,xfile,yfile)
        outfile = indir_med + '\\' + outfile0 +file+'_proj.tif'
        dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

        file = 'nir'
        datafile = fileDir+'\\'+file+'.tif'
        xfile = fileDir+'\\lon_v.tif'
        yfile = fileDir+'\\lat_v.tif'
        vrtfile = fileDir+'\\'+file+'.vrt'
        writevrt_v_nadir(vrtfile,datafile,xfile,yfile)
        outfile = indir_med + '\\' + outfile0 +file+'_proj.tif'
        dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)



if ifresample==1:


    fileNames = search_file_rej(indir_med, ['bt8'],'obliq')

    dataset = gdal.Open(targetArea)
    if dataset == None:
        print(targetArea + " ")
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
    temp1,temp2= imagexy2geo(dataset,nli,nsi)
    lon,lat = geo2lonlat(dataset,temp1,temp2)

    fileNum = np.size(fileNames)
    for k in range(fileNum):
        fileName = fileNames[k][:-9]
        fileDir = indir_med + fileName + r'\\'
        print(k, fileName)
        print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

        year = np.int(np.uint(fileName[4:8]))
        doy = np.int(np.uint(fileName[8:11]))
        time_str = fileName[12:12 + 4]
        hh = np.int(np.uint(time_str[0:2]))
        mm = np.int(np.uint(time_str[2:4]))
        infile0 = symbol+'_%4d'%year+'%03d'%doy+'_%02d'%hh+'%02d'%mm+'_'
        outfile0 = symbol + '_%4d' % year + '%03d' % doy + '_'

        ### for day data
        if hh > 22 or hh < 11:
            dayNight = 'day'
        else:
            dayNight = 'night'

        outfile1 = indir_tif + outfile0 + dayNight + '_bt1.tif'
        outfile2 = indir_tif + outfile0 + dayNight + '_bt2.tif'
        outfile3 = indir_tif + outfile0 + dayNight + '_red.tif'
        outfile4 = indir_tif + outfile0 + dayNight + '_nir.tif'
        outfile5 = indir_tif + outfile0 + dayNight + '_vza.tif'
        outfile6 = indir_tif + outfile0 + dayNight + '_sza.tif'
        outfile7 = indir_tif + outfile0 + dayNight + '_psi.tif'
        outfile8 = indir_tif + outfile0 + dayNight + '_time.tif'

        if os.path.exists(outfile1) == 1:
            [data1, temp1, temp2, temp3, geog, proj] = read_image_gdal(outfile1)
            [data2, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(outfile2)
            [data5, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(outfile5)
            [data8, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(outfile8)
            if dayNight == 'day':
                [data3, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(outfile3)
                [data4, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(outfile4)
                [data6, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(outfile6)
                [data7, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(outfile7)
        else:
            data1 = np.zeros([nl, ns])
            data2 = np.zeros([nl, ns])
            data3 = np.zeros([nl, ns])
            data4 = np.zeros([nl, ns])
            data5 = np.zeros([nl, ns])
            data6 = np.zeros([nl, ns])
            data7 = np.zeros([nl, ns])
            data8 = np.zeros([nl, ns])


        infile = indir_med + infile0+ r'bt8_proj.tif'
        dataset = gdal.Open(infile)
        if dataset == None:
            print(infile + " ")
            continue
        im_width = dataset.RasterXSize*1
        im_height = dataset.RasterYSize*1
        im_bands = dataset.RasterCount*1
        data = dataset.ReadAsArray(0, 0, im_width, im_height)
        im_geotrans = dataset.GetGeoTransform()
        im_proj = dataset.GetProjection()
        temp1,temp2=lonlat2geo(dataset,lat,lon)
        imagey,imagex=geo2imagexy(dataset,temp1,temp2)


        imagex = np.asarray(imagex,np.int)
        imagey = np.asarray(imagey,np.int)
        ind = (imagex>0)* (imagex<im_height-1) * (imagey>0) * (imagey<im_width-1)
        imagex[imagex<0] = 0
        imagex[imagex>(im_height-1)] = im_height-1
        imagey[imagey<0] = 0
        imagey[imagey>(im_width-1)] = im_width-1
        temp = data[imagex,imagey]
        indnew = (temp<500)*(temp>200)*ind
        if np.sum(indnew) == 0:continue
        data1 = np.reshape(data1,[-1])
        data1[indnew] = temp[indnew]
        data1 = np.reshape(data1,[nl,ns])
        data1[data1<0] = 0
        data1[data1 > 500] = 0

        infile = indir_med + infile0 + r'bt9_proj.tif'
        [data, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile)
        temp = data[imagex, imagey]
        indnew = (temp < 500) * (temp > 200) * ind
        data2 = np.reshape(data2, [-1])
        data2[indnew] = temp[indnew]
        data2 = np.reshape(data2, [nl, ns])
        data2[data2 < 0] = 0
        data2[data2 > 500] = 0

        data8 = np.reshape(data8, [-1])
        temp[:] = np.int(np.uint(time_str))
        data8[indnew] = temp[indnew]
        data8 = np.reshape(data8, [nl, ns])
        #
        #
        if dayNight == 'day':
            infile = indir_med + infile0 + r'red_proj.tif'
            [data, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile)
            data = resize_data_ls(data,im_height,im_width)
            temp = data[imagex, imagey]
            indnew = (temp < 500) * (temp > 0) * ind
            data3 = np.reshape(data3, [-1])
            data3[indnew] = temp[indnew]
            data3 = np.reshape(data3, [nl, ns])
            data3[data3 < 0] = 0
            data3[data3 > 500] = 0

            infile = indir_med + infile0 + r'nir_proj.tif'
            [data, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile)
            data = resize_data_ls(data,im_height,im_width)
            temp = data[imagex, imagey]
            indnew = (temp < 500) * (temp > 0) * ind
            data4 = np.reshape(data4, [-1])
            data4[indnew] = temp[indnew]
            data4 = np.reshape(data4, [nl, ns])
            data4[data4 < 0] = 0
            data4[data4 > 500] = 0

        infile = indir_med + infile0 + r'vza_proj.tif'
        [data, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile)
        temp = data[imagex, imagey]
        indnew = (temp < 500) * (temp > 0) * ind
        data5 = np.reshape(data5, [-1])
        data5[indnew] = temp[indnew]
        data5 = np.reshape(data5, [nl, ns])
        data5[data5 < 0] = 0
        data5[data5 > 500] = 0

        if dayNight == 'day':
            infile = indir_med + infile0 + r'sza_proj.tif'
            [data, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile)
            temp = data[imagex, imagey]
            indnew = (temp < 500) * (temp > 0) * ind
            data6 = np.reshape(data6, [-1])
            data6[indnew] = temp[indnew]
            data6 = np.reshape(data6, [nl, ns])
            data6[data6 < 0] = 0
            data6[data6 > 500] = 0

            infile = indir_med + infile0 + r'psi_proj.tif'
            [data, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile)
            temp = data[imagex, imagey]
            indnew = (temp < 500) * (temp > 0) * ind
            data7 = np.reshape(data7, [-1])
            data7[indnew] = temp[indnew]
            data7 = np.reshape(data7, [nl, ns])
            data7[data7 < 0] = 0
            data7[data7 > 500] = 0

        write_image_gdal(data1, ns, nl, 1, geog, proj, outfile1)
        write_image_gdal(data2, ns, nl, 1, geog, proj, outfile2)
        write_image_gdal(data5, ns, nl, 1, geog, proj, outfile5)
        write_image_gdal(data8, ns, nl, 1, geog, proj, outfile8)
        if dayNight == 'day':
            write_image_gdal(data3, ns, nl, 1, geog, proj, outfile3)
            write_image_gdal(data4, ns, nl, 1, geog, proj, outfile4)
            write_image_gdal(data6, ns, nl, 1, geog, proj, outfile6)
            write_image_gdal(data7, ns, nl, 1, geog, proj, outfile7)
