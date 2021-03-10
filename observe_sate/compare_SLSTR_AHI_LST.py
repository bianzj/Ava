
from myfun import *

indir_slstr_lst = r'i:\S3B_lst\\'
indir_slstr_tif = r'i:\s3b_tif\\'
indir_ahi_lst = r'i:\s3b_LST_AHI\\'
indir_base = r'I:\BASE\\'
indir_ahi_tif =r'I:\S3B_TIF_AHI'


fileNames = search_file_rej(indir_slstr_lst,['lst','day'],'obliq')
fileNum = np.size(fileNames)

off = -7
lstmin = 250
lstmax = 350

dif_ = []
day_ = []
desp_ = []

infile10 = indir_base  + 'ahi_vza_china.tif'
[vza_ahi,ns2,nl2,nb,geog,proj] = read_image_gdal(infile10)
infile11 = indir_base + 'ahi_vaa_china.tif'
[vaa_ahi,ns2,nl2,nb,geog,proj] = read_image_gdal(infile11)
ns = np.int(ns2 / 2)
nl = np.int(nl2 / 2)
vza_ahi = resize_data_ls(vza_ahi,nl,ns)
vaa_ahi = resize_data_ls(vaa_ahi, nl,ns)

infile12 = indir_base + 'type.tif'
[type,ns2,nl2,nb,geog,proj] = read_image_gdal(infile12)
type = resize_data_near(type, nl,ns)


data_ahi = np.asarray([])
data_slstr = np.asarray([])
for k in range(320,fileNum-1):
    fileName = fileNames[k]
    # print(fileName)
    infile1 = indir_slstr_lst+fileName[:off]+'LST.tif'
    infile2 = indir_slstr_lst+fileName[:off]+'LST_obliq.tif'
    infile3 = indir_slstr_tif+fileName[:off]+'vza.tif'
    infile4 = indir_slstr_tif+fileName[:off]+'vaa.tif'
    infile5 = indir_slstr_tif+fileName[:off]+'vza_obliq.tif'
    infile6 = indir_slstr_tif+fileName[:off]+'vaa_obliq.tif'

    infile7 = indir_slstr_tif+fileName[:off]+'cloud.tif'
    infile8 = indir_slstr_tif+fileName[:off]+'cloud_obliq.tif'

    infile9 = indir_ahi_lst+fileName[:off]+'lst.tif'

    if os.path.exists(infile1) == 0: continue
    if os.path.exists(infile3) == 0: continue
    if os.path.exists(infile9) == 0: continue

    [lst_slstr_n,ns2,nl2,nb,geog,proj] = read_image_gdal(infile1)
    [lst_slstr_o, ns2, nl2, nb, geog, proj] = read_image_gdal(infile2)
    [vza_slstr_n, ns2, nl2, nb, geog, proj] = read_image_gdal(infile3)
    [vaa_slstr_n, ns2, nl2, nb, geog, proj] = read_image_gdal(infile4)
    [vza_slstr_o, ns2, nl2, nb, geog, proj] = read_image_gdal(infile5)
    [vaa_slstr_o, ns2, nl2, nb, geog, proj] = read_image_gdal(infile6)
    [cloud_n, ns2, nl2, nb, geog, proj] = read_image_gdal(infile7)
    [cloud_o, ns2, nl2, nb, geog, proj] = read_image_gdal(infile8)
    [lst_ahi, ns2, nl2, nb, geog, proj] = read_image_gdal(infile9)

    ns = np.int(ns2/2)
    nl = np.int(nl2/2)

    lst_slstr_n = resize_data_ls(lst_slstr_n,nl,ns)
    lst_slstr_o = resize_data_ls(lst_slstr_o,nl,ns)
    vza_slstr_n = resize_data_ls(vza_slstr_n,nl,ns)
    vza_slstr_o = resize_data_ls(vza_slstr_o,nl,ns)
    vaa_slstr_n = resize_data_ls(vaa_slstr_n,nl,ns)
    vaa_slstr_o = resize_data_ls(vaa_slstr_o,nl,ns)
    cloud_n = resize_data_ls(cloud_n,nl,ns)
    cloud_o = resize_data_ls(cloud_o,nl,ns)
    lst_ahi = resize_data_ls(lst_ahi,nl,ns)


    diflst = 5
    ind1 = (lst_slstr_n<lstmax)*(lst_slstr_n>lstmin)*\
          (lst_ahi<lstmax)*(lst_ahi>lstmin)*\
          (cloud_n==0)*\
          (abs(lst_ahi-lst_slstr_n)<diflst)*\
          (abs(vza_ahi-vza_slstr_n)<2.5) *\
          (abs(vaa_ahi - vaa_slstr_n)<2.5)*\
            (type ==  14)

    ind2 = (lst_slstr_o<lstmax)*(lst_slstr_o>lstmin)*\
          (lst_ahi<lstmax)*(lst_ahi>lstmin)*\
          (cloud_o==0)*\
          (abs(lst_ahi-lst_slstr_o)<diflst)*\
          (abs(vza_ahi-vza_slstr_o)<2.5) *\
          (abs(vaa_ahi - vaa_slstr_o)<2.5) *\
            (type ==  14)

    if np.sum(ind1) > 0:
        data_slstr = np.hstack([data_slstr,lst_slstr_n[ind1]])
        data_ahi = np.hstack([data_ahi,lst_ahi[ind1]])
    if np.sum(ind2) > 0:
        data_slstr = np.hstack([data_slstr,lst_slstr_o[ind2]])
        data_ahi = np.hstack([data_ahi,lst_ahi[ind2]])


# plt.plot(data_slstr,data_ahi,'.')
# plt.show()

zz = np.polyfit(data_slstr,data_ahi,1)
pp = np.poly1d(zz)
data_ahi_new = pp(data_slstr)
dif = data_ahi - data_ahi_new
# bias = np.average(dif)
rmse = np.sqrt(np.average(dif*dif))
# r = np.corrcoef(data_slstr, data_ahi)
# r2 = r[0, 1] * r[0, 1]
print(pp,rmse)


plt.hist2d(data_ahi,data_slstr,bins=100,norm=LogNorm(),cmap='jet')
plt.xlim([lstmin,lstmax])
plt.ylim([lstmin,lstmax])
plt.colorbar()
plt.title('SLSTR_N-AHI')
plt.plot([lstmin,lstmax],[lstmin,lstmax],'k--')
plt.grid()
plt.show()
