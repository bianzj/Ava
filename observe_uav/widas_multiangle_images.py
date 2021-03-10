from util.myfun import *

dir = r'D:\Code\GF3D\DATA\heihe\\'
infile_loc = dir + r'sample.tif'
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
temp1,temp2 = imagexy2geo(dataset,nli,nsi)
lon,lat = geo2lonlat(dataset,temp1,temp2)
lon_o = lon
lat_o = lat


indir1 =r'F:\WIDAS\0803\Cam01\\'
indir2 = r'F:\WIDAS\0803\Cam02\\'
outdir = r'D:\code\GF3D\DATA\heihe\\'


fileNames1 =  search_file_rej(indir1,['track11-2','TIR'],'hdr')
fileNames2 =  search_file_rej(indir2,['track11-2','TIR'],'hdr')
fileNames1 = [indir1+fileName for fileName in fileNames1]
fileNames2 = [indir2+fileName for fileName in fileNames2]
fileNames = fileNames1+fileNames2


outfile = outdir + 'sample_track11-2_tir.tif'
theBand = 0


number = len(fileNames)
data_output = np.zeros([number,nl, ns])
for kimage in range(number):
    filename = fileNames[kimage]
    infile =  filename
    print(infile)


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
    temp1, temp2 = lonlat2geo(dataset, lat_o, lon_o)
    imagex, imagey = geo2imagexy(dataset, temp1, temp2)
    imagex = np.asarray(imagex, np.int)
    imagey = np.asarray(imagey, np.int)
    # print(np.max(imagex),np.max(imagey))

    ind = (imagex > 0) * (imagex < im_width - 1) * (imagey > 0) * (imagey < im_height - 1)
    imagex[imagex < 0] = 0
    imagex[imagex > im_width - 1] = im_width - 1
    imagey[imagey < 0] = 0
    imagey[imagey > im_height - 1] = im_height - 1

    if im_bands==1:
        temp = data[imagey, imagex]
    else:
        temp = data[theBand,imagey,imagex]

    indnew = (temp > 10) * ind
    if np.sum(indnew) <= 100: continue

    data_new = np.zeros([nl,ns])
    data_new = np.reshape(data_new,[-1])
    data_new[indnew] = temp[indnew]
    data_new = np.reshape(data_new,[nl,ns])
    data_new[data_new > 20000] = np.max(temp[temp<20000])
    data_output[kimage,:,:] = data_new

write_image_gdal(data_output, ns, nl, 1, geog, proj, outfile)


