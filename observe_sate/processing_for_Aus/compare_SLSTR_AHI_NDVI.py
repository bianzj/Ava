
from processing_for_Aus.myfun import *

indir1 = r'E:\SLSTR_day_result\\'
indir2 = r'E:\AHI8_day_result\\'
indir3 = r'E:\SLSTR_day\\'
outdir = r'E:\figure\\'

fileNames = ope_search2(indir1,'ndvi','n')
fileNum = np.size(fileNames)

off = 13
ndvimin = 0
ndvimax = 1.0

dif_ = []
day_ = []
desp_ = []

bias1_ = []
bias2_ = []
for k in range(fileNum-1):
    fileName = fileNames[k]
    print(fileName)
    infile1 = indir1+fileName[:off]+'ndvi_n.tif'
    infile2 = indir1+fileName[:off]+'ndvi_o.tif'
    infile3 = indir2+fileName[:off]+'ndvi.tif'
    infile4 = indir3+fileName[:off]+'cloud_n.tif'
    infile5 = indir3+fileName[:off]+'cloud_o.tif'
    # outfile = outdir + fileName[:off]+'LST.tif'

    [ndvi_slstr_n, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile1)
    [ndvi_slstr_o, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile2)
    [cloud_n, ns, nl, nb, geog, proj] = read_gdalTiff(infile4)
    [cloud_o, ns, nl, nb, geog, proj] = read_gdalTiff(infile5)
    [ndvi_ahi, ns, nl, nb, geog, proj] = read_gdalTiff(infile3)

    ndvi_slstr_n = ope_resizeData(ndvi_slstr_n, ns, nl)
    ndvi_slstr_o = ope_resizeData(ndvi_slstr_o, ns, nl)
    cloud_n = ope_resizeData(cloud_n,ns,nl)
    cloud_o = ope_resizeData(cloud_o,ns,nl)

    ndvi_slstr_n = np.reshape(ndvi_slstr_n, -1)
    ndvi_slstr_o = np.reshape(ndvi_slstr_o, -1)
    ndvi_ahi = np.reshape(ndvi_ahi, -1)
    cloud_n = np.reshape(cloud_n,-1)
    cloud_o = np.reshape(cloud_o,-1)

    difndvi = 0.5
    ind = (ndvi_slstr_n < ndvimax) * (ndvi_slstr_n > ndvimin) * \
          (ndvi_slstr_o < ndvimax) * (ndvi_slstr_o > ndvimin) * \
          (ndvi_ahi < ndvimax) * (ndvi_ahi > ndvimin) * \
          (cloud_n==0) * (cloud_o==0) * (abs(ndvi_slstr_n - ndvi_slstr_o) < difndvi) * \
          (abs(ndvi_ahi - ndvi_slstr_n) < difndvi) * (abs(ndvi_ahi - ndvi_slstr_o) < difndvi)


    data1 = ndvi_ahi[ind]
    data2 = ndvi_slstr_n[ind]
    dif = data2-data1
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif*dif))
    rr = np.corrcoef(data1,data2)
    r2 = rr[1,0]*rr[1,0]
    bias1_.append(bias)
    # plt.hist2d(data1,data2,bins=100,norm=LogNorm(),cmap='jet')
    # plt.xlim([lstmin,lstmax])
    # plt.ylim([lstmin,lstmax])
    # plt.colorbar()
    # plt.text(260,320,
    #       '$RMSE = %3.2f^\circ$C\n' % rmse +
    #       '$Bias = %3.2f^\circ$C\n' % bias +
    #       '$R^2 = %3.2f$' % r2,
    #       fontsize=14)
    # plt.title('SLSTR_N-AHI')
    # plt.plot([lstmin,lstmax],[lstmin,lstmax],'k--')
    # plt.grid()
    # plt.show()
    #
    #
    data1 = ndvi_ahi[ind]
    data2 = ndvi_slstr_o[ind]
    dif = data2-data1
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif*dif))
    rr = np.corrcoef(data1,data2)
    r2 = rr[1,0]*rr[1,0]
    bias2_.append(bias)
    # plt.hist2d(data1,data2,bins=100,norm=LogNorm(),cmap='jet')
    # plt.xlim([lstmin,lstmax])
    # plt.ylim([lstmin,lstmax])
    # plt.colorbar()
    # plt.text(260,320,
    #       '$RMSE = %3.2f^\circ$C\n' % rmse +
    #       '$Bias = %3.2f^\circ$C\n' % bias +
    #       '$R^2 = %3.2f$' % r2,
    #       fontsize=14)
    # plt.title('SLSTR_N-AHI')
    # plt.plot([lstmin,lstmax],[lstmin,lstmax],'k--')
    # plt.grid()
    # plt.show()
    #
    # data1 = lst_slstr_n[ind]
    # data2 = lst_slstr_o[ind]
    # dif = data2-data1
    # bias = np.mean(dif)
    # rmse = np.sqrt(np.mean(dif*dif))
    # rr = np.corrcoef(data1,data2)
    # r2 = rr[1,0]*rr[1,0]
    # plt.hist2d(data1,data2,bins=100,norm=LogNorm(),cmap='jet')
    # plt.xlim([lstmin,lstmax])
    # plt.ylim([lstmin,lstmax])
    # plt.colorbar()
    # plt.text(260,320,
    #       '$RMSE = %3.2f^\circ$C\n' % rmse +
    #       '$Bias = %3.2f^\circ$C\n' % bias +
    #       '$R^2 = %3.2f$' % r2,
    #       fontsize=14)
    # plt.title('SLSTR_O-SLSTR_N')
    # plt.plot([lstmin,lstmax],[lstmin,lstmax],'k--')
    # plt.grid()
    # plt.show()

bias1_ = np.asarray(bias1_)
bias2_ = np.asarray(bias2_)
x = np.arange(1,32)
plt.bar(x,bias1_)
plt.bar(x,bias2_)
plt.show()

### boxplot
# data_ = ([day_,desp_,dif_])
# dff = pd.DataFrame(np.transpose(data_),
#                   columns=['day', 'desp','dif'])
# ax = sns.boxplot(x=dff["day"],
#                  y=np.asarray(dff["dif"],dtype=np.float),
#                  hue=dff["desp"],
#                  palette="Set3",)
# plt.show()