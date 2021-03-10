from processing_for_Aus.SLSTR import *
from processing_for_Aus.myfun import *

s = SLSTR()

wdir = r'J:\LANDSAT\\'
fileNames = ope_search1_rej1(wdir,'LC08','tar')
fileNum = np.size(fileNames)


infile_lse1 = 'J:\LSE\LSE1_Aus_1000.tif'
dataset1 = gdal.Open(infile_lse1)
if dataset1 == None:
    print(infile_lse1 + "文件无法打开")
im_width = dataset1.RasterXSize
im_height = dataset1.RasterYSize
im_bands = dataset1.RasterCount
lse1 = dataset1.ReadAsArray(0, 0, im_width, im_height)
im_geotrans = dataset1.GetGeoTransform()
im_proj = dataset1.GetProjection()

infile_lse2 = 'J:\LSE\LSE2_Aus_1000.tif'
[lse2, temp1, temp2, temp3, geog0, proj0] = s.getTiffData(infile_lse2)

for k in range(fileNum):

    fileName = fileNames[k]
    infile_loc = wdir+ fileName + "\\"+fileName+'_B10.TIF'

    dataset = gdal.Open(infile_loc)
    if dataset == None:
        print(infile_loc + "文件无法打开")
    ns = dataset.RasterXSize
    nl = dataset.RasterYSize
    im_bands = dataset.RasterCount
    geog = dataset.GetGeoTransform()
    proj = dataset.GetProjection()
    nsi = np.zeros([nl, ns])
    nli = np.zeros([nl, ns])
    for k in range(nl):
        nli[k, :] = k
        nsi[k, :] = np.linspace(0, ns - 1, ns)
    nli = np.reshape(nli, -1)
    nsi = np.reshape(nsi, -1)
    geox, geoy = s.imagexy2geo(dataset, nli, nsi)
    lon, lat = s.geo2lonlat(dataset, geox, geoy)



    outfile1 = wdir + fileName + "\\"+fileName+'_B10_bare_emi_M22.tif'
    outfile2 = wdir + fileName + "\\"+fileName+'_B11_bare_emi_M22.tif'


    if os.path.exists(outfile1)==1:
        [data1, temp1, temp2, temp3, geog, proj] = s.getTiffData(outfile1)
        [data2, temp1, temp2, temp3, geog, proj] = s.getTiffData(outfile2)

    else:
        data1 = np.zeros([nl,ns])
        data2 = np.zeros([nl, ns])

    temp1,temp2=s.lonlat2geo(dataset1,lat,lon)
    imagey,imagex=s.geo2imagexy(dataset1,temp1,temp2)
    imagex = np.asarray(imagex,np.int)
    imagey = np.asarray(imagey,np.int)
    ind = (imagex>0)* (imagex<im_height-1) * (imagey>0) * (imagey<im_width-1)
    imagex[imagex<0] = 0
    imagex[imagex>im_height-1] = im_height-1
    imagey[imagey<0] = 0
    imagey[imagey>im_width-1] = im_width-1


    temp = lse1[imagex,imagey]
    indnew = (temp<1000)*(temp>0)*ind
    if np.sum(indnew) == 0:continue
    data1 = np.reshape(data1,[-1])
    data1[indnew] = temp[indnew]
    data1 = np.reshape(data1,[nl,ns])
    data1[data1<0] = 0
    data1[data1 > 1000] = 0
    data1 = data1/1000.0

    temp = lse2[imagex,imagey]
    indnew = (temp<1000)*(temp>0)*ind
    if np.sum(indnew) == 0:continue
    data2 = np.reshape(data2,[-1])
    data2[indnew] = temp[indnew]
    data2 = np.reshape(data2,[nl,ns])
    data2[data2<0] = 0
    data2[data2 > 1000] = 0
    data2 = data2/1000.0




    s.writeTiff(data1,ns, nl, 1, geog, proj, outfile1)
    s.writeTiff(data2, ns, nl, 1, geog, proj, outfile2)





