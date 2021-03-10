from processing_for_Aus.myfun import *
from scipy.linalg import lstsq
import time

indir1 = r'J:\SLSTR_day_result\\'
indir2 = r'J:\AHI8_day_result\\'
indir3 = r'J:\SLSTR_day\\'
outdir = r'J:\figure\\'

fileNames = ope_search2(indir1,'LST','n')
fileNum = np.size(fileNames)



off = 13
lstmin = 250
lstmax = 350
ndvimin = 0
ndvimax = 1.0

dif_ = []
day_ = []
desp_ = []

bias1_ = []
bias2_ = []

# model = LinearRegression()
print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
for k in range(1):

    fileName = fileNames[k]
    print(fileName)
    infile1 = indir1+fileName[:off]+'LST_n.tif'
    infile2 = indir1+fileName[:off]+'LST_o.tif'
    infile3 = indir2+fileName[:off]+'LST.tif'
    infile4 = indir1+fileName[:off]+'ndvi_n.tif'
    infile5 = indir1+fileName[:off]+'ndvi_o.tif'
    infile6 = indir2+fileName[:off]+'ndvi.tif'
    infile7 = indir1+fileName[:off]+'bf_n.tif'
    infile8 = indir1+fileName[:off]+'bf_o.tif'
    infile9 = indir2+fileName[:off]+'bf.tif'
    infile10 = indir3+fileName[:off]+'cloud_n.tif'
    infile11 = indir3+fileName[:off]+'cloud_o.tif'

    [lst_slstr_n,ns2,nl2,nb,geog,proj] = read_gdalTiff(infile1)
    [lst_slstr_o, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile2)
    [lst_ahi, ns, nl, nb, geog, proj] = read_gdalTiff(infile3)
    [ndvi_slstr_n,ns2,nl2,nb,geog,proj] = read_gdalTiff(infile4)
    [ndvi_slstr_o, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile5)
    [ndvi_ahi, ns, nl, nb, geog, proj] = read_gdalTiff(infile6)
    [bf_slstr_n,ns2,nl2,nb,geog,proj] = read_gdalTiff(infile7)
    [bf_slstr_o, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile8)
    [bf_ahi, ns, nl, nb, geog, proj] = read_gdalTiff(infile9)

    [cloud_n, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile10)
    [cloud_o, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile11)


    lst_slstr_n = ope_resizeData(lst_slstr_n,ns,nl)
    lst_slstr_o = ope_resizeData(lst_slstr_o, ns, nl)
    ndvi_slstr_n = ope_resizeData(ndvi_slstr_n,ns,nl)
    ndvi_slstr_o = ope_resizeData(ndvi_slstr_o, ns, nl)
    bf_slstr_n = ope_resizeData(bf_slstr_n,ns,nl)
    bf_slstr_o = ope_resizeData(bf_slstr_o, ns, nl)
    cloud_n = ope_resizeData(cloud_n,ns,nl)
    cloud_o = ope_resizeData(cloud_o,ns,nl)

    lst_slstr_n = np.reshape(lst_slstr_n,-1)
    lst_slstr_o = np.reshape(lst_slstr_o,-1)
    lst_ahi = np.reshape(lst_ahi,-1)
    ndvi_slstr_n = np.reshape(ndvi_slstr_n,-1)
    ndvi_slstr_o = np.reshape(ndvi_slstr_o,-1)
    ndvi_ahi = np.reshape(ndvi_ahi,-1)
    bf_slstr_n = np.reshape(bf_slstr_n,-1)
    bf_slstr_o = np.reshape(bf_slstr_o,-1)
    bf_ahi = np.reshape(bf_ahi,-1)
    cloud_n = np.reshape(cloud_n,-1)
    cloud_o = np.reshape(cloud_o,-1)


    difndvi = 0.5
    ind1 = (ndvi_slstr_n < ndvimax) * (ndvi_slstr_n > ndvimin) * \
          (ndvi_slstr_o < ndvimax) * (ndvi_slstr_o > ndvimin) * \
          (ndvi_ahi < ndvimax) * (ndvi_ahi > ndvimin) * \
          (cloud_n==0) * (cloud_o==0) * (abs(ndvi_slstr_n - ndvi_slstr_o) < difndvi) * \
          (abs(ndvi_ahi - ndvi_slstr_n) < difndvi) * (abs(ndvi_ahi - ndvi_slstr_o) < difndvi)

    diflst = 15
    ind2 = (lst_slstr_n<lstmax)*(lst_slstr_n>lstmin)*\
          (lst_slstr_o<lstmax)*(lst_slstr_o>lstmin)*\
          (lst_ahi<lstmax)*(lst_ahi>lstmin)*\
          (cloud_n==0)*(cloud_o==0)*(abs(lst_slstr_n-lst_slstr_o)<diflst)*\
          (abs(lst_ahi-lst_slstr_n)<diflst)*(abs(lst_ahi-lst_slstr_o)<diflst)
    ind = ind1*ind2


    length = np.sum(ind)
    yy_ = np.zeros([length,3])
    y_ = np.transpose(np.asarray([lst_slstr_n[ind], lst_slstr_o[ind], lst_ahi[ind]]))
    for kk in range(length):
        lst1 = lst_slstr_n[ind][kk]
        lst2 = lst_slstr_o[ind][kk]
        lst3 = lst_ahi[ind][kk]
        ndvi1 = ndvi_slstr_n[ind][kk]
        ndvi2 = ndvi_slstr_o[ind][kk]
        ndvi3 = ndvi_ahi[ind][kk]
        bf1 = bf_slstr_n[ind][kk]
        bf2 = bf_slstr_o[ind][kk]
        bf3 = bf_ahi[ind][kk]

        x = np.transpose(np.asarray([[ndvi1,ndvi2,ndvi3],[bf1,bf2,bf3],[1,1,1]]))
        y = np.asarray([lst1,lst2,lst3])

        coeff = np.asarray(lstsq(x, y))[0]
        yy = np.dot(x,coeff)
        # model.fit(x,y)
        # b = model.intercept_  # 截距
        # k = model.coef_  # 回归系数
        # yy = model.predict(x)
        yy_[kk,:] = yy

    print('success!')

print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))