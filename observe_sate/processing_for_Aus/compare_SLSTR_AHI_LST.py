
from processing_for_Aus.myfun import *

indir1 = r'J:\SLSTR_day_result\\'
indir2 = r'J:\AHI8_day_result\\'
indir3 = r'J:\SLSTR_day\\'
outdir = r'J:\figure\\'

fileNames = ope_search2(indir1,'LST','n')
fileNum = np.size(fileNames)

off = 13
lstmin = 250
lstmax = 350

dif_ = []
day_ = []
desp_ = []

bias1_ = []
bias2_ = []
for k in range(fileNum-1):
    fileName = fileNames[k]
    print(fileName)
    infile1 = indir1+fileName[:off]+'LST_n.tif'
    infile2 = indir1+fileName[:off]+'LST_o.tif'
    infile3 = indir2+fileName[:off]+'LST.tif'
    infile4 = indir3+fileName[:off]+'cloud_n.tif'
    infile5 = indir3+fileName[:off]+'cloud_o.tif'
    outfile = outdir + fileName[:off]+'LST.tif'

    [lst_slstr_n,ns2,nl2,nb,geog,proj] = read_gdalTiff(infile1)
    [lst_slstr_o, ns2, nl2, nb, geog, proj] = read_gdalTiff(infile2)
    [cloud_n, ns, nl, nb, geog, proj] = read_gdalTiff(infile4)
    [cloud_o, ns, nl, nb, geog, proj] = read_gdalTiff(infile5)
    [lst_ahi, ns, nl, nb, geog, proj] = read_gdalTiff(infile3)

    lst_slstr_n = ope_resizeData(lst_slstr_n,ns,nl)
    lst_slstr_o = ope_resizeData(lst_slstr_o, ns, nl)
    cloud_n = ope_resizeData(cloud_n,ns,nl)
    cloud_o = ope_resizeData(cloud_o,ns,nl)

    lst_slstr_n = np.reshape(lst_slstr_n,-1)
    lst_slstr_o = np.reshape(lst_slstr_o,-1)
    lst_ahi = np.reshape(lst_ahi,-1)
    cloud_n = np.reshape(cloud_n,-1)
    cloud_o = np.reshape(cloud_o,-1)

    diflst = 15
    ind = (lst_slstr_n<lstmax)*(lst_slstr_n>lstmin)*\
          (lst_slstr_o<lstmax)*(lst_slstr_o>lstmin)*\
          (lst_ahi<lstmax)*(lst_ahi>lstmin)*\
          (cloud_n==0)*(cloud_o==0)*(abs(lst_slstr_n-lst_slstr_o)<diflst)*\
          (abs(lst_ahi-lst_slstr_n)<diflst)*(abs(lst_ahi-lst_slstr_o)<diflst)


    data1 = lst_ahi[ind]
    data2 = lst_slstr_n[ind]
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
    data1 = lst_ahi[ind]
    data2 = lst_slstr_o[ind]
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