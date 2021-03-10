from AHI import *
from processing_for_Aus.myfun import *
from geopandas import *
import rasterio as rio

a = AHI()
wdir0 = r'F:\AHI8\L1\2018\201801\\'
outdir = r'E:\AHI8_day_01\\AHI_'
shpdatafile=r'E:\AHI8_sample\country.shp'
rasterfile = r'E:\AHI8_sample\AHI8.tif'
[data,ns,nl,nb,trans,proj]=a.getTiffData(rasterfile)
src0 = rio.open(rasterfile)
shpdata = GeoDataFrame.from_file(shpdatafile)
shpdata_crs = shpdata.to_crs(src0.crs)
features = [shpdata_crs.geometry[12].__geo_interface__]


for kd in range(16,17):

    wday = '%02d'%kd
    wdir = wdir0 + wday+'\\'
    fileNames = ope_search1(wdir,'nc')
    fileNum = np.size(fileNames)

    for k in range(0,fileNum):

        fileName = fileNames[k]
        imageTime = fileName[16:16 + 4]
        hm = np.int(imageTime[0:4])

        # if hm < 2100 and hm > 900: continue  ### for day data
        # if hm > 900 : continue
        if hm < 2100: continue
        if (hm%100)%30 != 0:continue

        # if hm <= 210 or hm > 220: continue  ### for day data

        in_file = wdir + fileName
        print(in_file)

        file_size = os.path.getsize(in_file)
        if file_size < 5000000: continue

        band = 'vza'
        out_file_temp = outdir + 'temp_' + band + '.tif'
        out_file = outdir + wday + '_' + imageTime + '_' + band + '.tif'
        data = a.getDatafromNc(in_file,'SAZ')
        a.writeTiff(data, a.ns, a.ns, 1, trans, proj, out_file_temp)
        a.clipfile(out_file_temp,features,out_file)


        band = 'sza'
        out_file_temp = outdir + 'temp_' + band + '.tif'
        out_file = outdir + wday + '_' + imageTime + '_' + band + '.tif'
        data = a.getDatafromNc(in_file,'SOZ')
        a.writeTiff(data, a.ns, a.ns, 1, trans, proj, out_file_temp)
        a.clipfile(out_file_temp,features,out_file)

        band = 'psi'
        out_file_temp = outdir + 'temp_' + band + '.tif'
        out_file = outdir + wday + '_' + imageTime + '_' + band + '.tif'
        saa = a.getDatafromNc(in_file,'SOA')
        vaa = a.getDatafromNc(in_file,'SAA')
        psi = np.abs(vaa-saa)
        a.writeTiff(psi, a.ns, a.ns, 1, trans, proj, out_file_temp)
        a.clipfile(out_file_temp,features,out_file)

        band = 'B1'
        out_file_temp = outdir + 'temp_' + band + '.tif'
        out_file = outdir + wday + '_' + imageTime + '_' + band + '.tif'
        data = a.getDatafromNc(in_file,'tbb_14')
        a.writeTiff(data, a.ns, a.ns, 1, trans, proj, out_file_temp)
        a.clipfile(out_file_temp,features,out_file)

        band = 'B2'
        out_file_temp = outdir + 'temp_' + band + '.tif'
        out_file = outdir + wday + '_' + imageTime + '_' + band + '.tif'
        data = a.getDatafromNc(in_file,'tbb_15')
        a.writeTiff(data, a.ns, a.ns, 1, trans, proj, out_file_temp)
        a.clipfile(out_file_temp,features,out_file)

        band = 'red'
        out_file_temp = outdir + 'temp_' + band + '.tif'
        out_file = outdir + wday + '_' + imageTime + '_' + band + '.tif'
        data = a.getDatafromNc(in_file,'albedo_03')
        a.writeTiff(data, a.ns, a.ns, 1, trans, proj, out_file_temp)
        a.clipfile(out_file_temp,features,out_file)

        band = 'nir'
        out_file_temp = outdir + 'temp_' + band + '.tif'
        out_file = outdir + wday + '_' + imageTime + '_' + band + '.tif'
        data = a.getDatafromNc(in_file,'albedo_04')
        a.writeTiff(data, a.ns, a.ns, 1, trans, proj, out_file_temp)
        a.clipfile(out_file_temp,features,out_file)


