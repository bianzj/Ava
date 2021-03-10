from processing_for_Aus.SLSTR import *
import time
from processing_for_Aus.myfun import *

s = SLSTR()

wdir = r'G:\SLSTR3A\\'
wdirout = r'G:\S3A\\'


fileNames = ope_search1_rej1(wdir,'S3A','zip')
fileNum = np.size(fileNames)

infile_loc =  r'G:\base\type.tif'
dataset = gdal.Open(infile_loc)
if dataset == None:
    print(infile_loc + "文件无法打开")
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
lon,lat=s.imagexy2geo(dataset,nli,nsi)
# lon,lat = s.geo2lonlat(dataset,temp1,temp2)
lon_o = lon
lat_o = lat

# lst inversion
dataTimeold = 0
imageTimeold = 0
fileNum = np.size(fileNames)
for k in range(fileNum):
    fileName = fileNames[k]
    fileDir = wdir +  fileName + r'\\'
    print(k, fileName)
    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    ### dataTime and imageTime
    dataTime = fileName[16:8 + 16]
    dd = np.int(np.uint(dataTime[6:8]))
    # if dd > 5:continue

    imageTime = fileName[25:25 + 6]
    hh = np.int(np.uint(imageTime[0:2]))
    # if hh < 18 and hh > 6: continue  ### for day data
    # if hh >20 or hh < 5: continue  ### for night data
    if (dataTime == dataTimeold and imageTime == imageTimeold):continue
    dataTimeold = dataTime
    imageTimeold = imageTime
    outfile1 = wdirout + dataTime + '_bt8_n.tif'
    outfile2 = wdirout + dataTime + '_bt9_n.tif'
    outfile3 = wdirout + dataTime + '_tcw_n.tif'
    outfile4 = wdirout + dataTime + '_cloud_n.tif'
    outfile5 = wdirout + dataTime + '_vza_n.tif'
    outfile6 = wdirout + dataTime + '_sza_n.tif'
    outfile7 = wdirout + dataTime + '_psi_n.tif'
    outfile8 = wdirout + dataTime + '_passtime.tif'

    dataTime = np.uint(dataTime)
    imageTime = np.uint(imageTime)

    if os.path.exists(outfile1)==1:
        [data1, temp1, temp2, temp3, geog, proj] = s.getTiffData(outfile1)
        [data2, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile2)
        [data3, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile3)
        [data4, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile4)
        [data5, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile5)
        [data6, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile6)
        [data7, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile7)
        [data8, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile8)
    else:
        data1 = np.zeros([nl,ns])
        data2 = np.zeros([nl, ns])
        data3 = np.zeros([nl, ns])
        data4 = np.zeros([nl, ns])
        data5 = np.zeros([nl, ns])
        data6 = np.zeros([nl, ns])
        data7 = np.zeros([nl, ns])
        data8 = np.zeros([nl, ns])


    infile = fileDir + r'\bt8_proj.tif'
    dataset = gdal.Open(infile)
    if dataset == None:
        print(infile + "文件无法打开")
        continue
    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize
    im_bands = dataset.RasterCount
    data = dataset.ReadAsArray(0, 0, im_width, im_height)
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
    temp = data[imagex,imagey]
    indnew = (temp<500)*(temp>200)*ind
    if np.sum(indnew) == 0:continue
    data1 = np.reshape(data1,[-1])
    data1[indnew] = temp[indnew]
    data1 = np.reshape(data1,[nl,ns])
    data1[data1<0] = 0
    data1[data1 > 500] = 0

    infile = fileDir+r'bt9_proj.tif'
    [data, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile)
    temp = data[imagex, imagey]
    indnew = (temp < 500) * (temp > 200) * ind
    data2 = np.reshape(data2,[-1])
    data2[indnew] = temp[indnew]
    data2 = np.reshape(data2,[nl,ns])
    data2[data2<0] = 0
    data2[data2 > 500] = 0

    data8 = np.reshape(data8,[-1])
    temp[:] = imageTime
    data8[indnew] = temp[indnew]
    data8 = np.reshape(data8,[nl,ns])

    #
    infile = fileDir+r'tcw_proj.tif'
    [data, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile)
    temp = data[imagex, imagey]
    indnew = (temp < 500) * (temp > 0) * ind
    data3 = np.reshape(data3,[-1])
    data3[indnew] = temp[indnew]
    data3 = np.reshape(data3,[nl,ns])
    data3[data3<0] = 0
    data3[data3 > 500] = 0
    #
    infile = fileDir+r'cloud_proj.tif'
    [data, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile)
    temp = data[imagex, imagey]
    indnew = (temp < 500) * (temp > 0) * ind
    data4 = np.reshape(data4,[-1])
    data4[indnew] = temp[indnew]
    data4 = np.reshape(data4,[nl,ns])
    data4[data4<0] = 0
    data4[data4 > 500] = 0

    infile = fileDir+r'vza_proj.tif'
    [data, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile)
    temp = data[imagex, imagey]
    indnew = (temp < 500) * (temp > 0) * ind
    data5 = np.reshape(data5,[-1])
    data5[indnew] = temp[indnew]
    data5 = np.reshape(data5,[nl,ns])
    data5[data5<0] = 0
    data5[data5 > 500] = 0

    infile = fileDir+r'sza_proj.tif'
    [data, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile)
    temp = data[imagex, imagey]
    indnew = (temp < 500) * (temp > 0) * ind
    data6 = np.reshape(data6,[-1])
    data6[indnew] = temp[indnew]
    data6 = np.reshape(data6,[nl,ns])
    data6[data6<0] = 0
    data6[data6 > 500] = 0

    infile = fileDir+r'psi_proj.tif'
    [data, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile)
    temp = data[imagex, imagey]
    indnew = (temp < 500) * (temp > 0) * ind
    data7 = np.reshape(data7,[-1])
    data7[indnew] = temp[indnew]
    data7 = np.reshape(data7,[nl,ns])
    data7[data7<0] = 0
    data7[data7 > 500] = 0



    s.writeTiff(data1,ns, nl, 1, geog, proj, outfile1)
    s.writeTiff(data2,ns, nl, 1, geog, proj, outfile2)
    s.writeTiff(data3,ns, nl, 1, geog, proj, outfile3)
    s.writeTiff(data4,ns, nl, 1, geog, proj, outfile4)
    s.writeTiff(data5,ns, nl, 1, geog, proj, outfile5)
    s.writeTiff(data6,ns, nl, 1, geog, proj, outfile6)
    s.writeTiff(data7,ns, nl, 1, geog, proj, outfile7)
    s.writeTiff(data8, ns, nl, 1, geog, proj, outfile8)




