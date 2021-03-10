from processing_for_Aus.myfun import *

istation = 8
lat_ = [-22.283,	-34.00206,-33.615278,-14.1592,
       -15.2588,-31.376351,	-30.1914,-12.4952,
       -13.17904178,-32.506102,	-17.11746943,-17.1507,
       -22.287,-35.6566,-43.09502,-37.4222,-34.9878]
lon_ = [133.249,140.58912,150.723611,131.3881,
        132.3706,115.71377,120.6541667,131.15005,
        130.7945459,116.966827,145.6301375,133.3502,
        133.64,148.1516,146.65452,144.0944,146.2908]
label_ = ['ASM','Calperum','CumberlandPlain','DalyUncleared',
          'DryRiver','Gingin','GWW','HowardSprints',
          'Litchfield','Ridgefield','Robson','SturtPlains',
          'TTE','Tumbarumba','Warra','WombatStateForest','Yanco']

### Tumbarumba and Warra data lost
### WombatstateForest data missing -9999

utc_darwin = 9.5
utc_adelaide = 9.5 # 10.5
utc_melbourne = 10 #11
utc_sydney = 10 #11
utc_perth = 8
utc_brisbane = 10
utc_hobart = 10 # 11


utc_t_ = [utc_darwin,utc_adelaide,utc_sydney,utc_darwin,
          utc_darwin,utc_perth,utc_perth,utc_darwin,
          utc_darwin,utc_perth,utc_brisbane,utc_darwin,
          utc_darwin,utc_sydney,utc_hobart,utc_melbourne,utc_sydney]


utc_t = utc_t_[istation]
lat = lat_[istation]
lon = lon_[istation]
label = label_[istation]

ndvi_max = 0.99
ndvi_min = 0.05
emis_s = 0.925
emis_v = 0.975
indir1 = 'J:\AHI8_day_lst\\'
indir2 = 'J:\AHI8_day\\'


# get data from measurements for validation
infile1 = 'J:\measurement\\'+label+'_fld.txt'
infile2 = 'J:\measurement\\'+label+'_flu.txt'
dl0 = read_txt2array(infile1)
ul0 = read_txt2array(infile2)

## get data from observations of AHI

infiles = ope_search1(indir1, 'LST')
infiles1 = ope_search1(indir2, 'red')
infiles2 = ope_search1(indir2, 'nir')
filenum = len(infiles)
dd_ = []
hh_ = []
mm_ = []
lst_ = []
ndvi_ = []
iffirst = 1
for k in range(0,500-1):
    filename = infiles[k]
    print(filename)
    infile = indir1 + filename
    infile1 = indir2 + infiles1[k]
    infile2 = indir2 + infiles2[k]
    dd = np.int(filename[4:6])
    hh = np.int(filename[7:9])
    mm = np.int(filename[9:11])
    if (mm != 0) and (mm!=30): continue
    dd_.append(dd)
    hh_.append(hh)
    mm_.append(mm)

    if iffirst == 1:
        dataset = gdal.Open(infile)
        if dataset == None:
            print(infile + " can not open !")
            continue
        im_width = dataset.RasterXSize
        im_height = dataset.RasterYSize
        im_bands = dataset.RasterCount
        data = dataset.ReadAsArray(0, 0, im_width, im_height)
        im_geotrans = dataset.GetGeoTransform()
        im_proj = dataset.GetProjection()
        temp1, temp2 = ope_lonlat2geo(dataset, lat, lon)
        imagey, imagex = ope_geo2imagexy(dataset, temp1, temp2)
        imagex = np.asarray(imagex, np.int)
        imagey = np.asarray(imagey, np.int)
        [data1, temp1, temp2, temp3, geog, proj] = getTiffData(infile1)
        [data2, temp1, temp2, temp3, geog, proj] = getTiffData(infile2)
        iffirst = 0
    else:
        [data, temp1, temp2, temp3, geog, proj] = getTiffData(infile)
        [data1, temp1, temp2, temp3, geog, proj] = getTiffData(infile1)
        [data2, temp1, temp2, temp3, geog, proj] = getTiffData(infile2)

    lst = data[imagey,imagex]
    nir = data2[imagey,imagex]
    red = data1[imagey, imagex]
    ndvi = 1.0*(nir-red)/(nir+red)
    lst_.append(lst)
    ndvi_.append(ndvi)
    print(dd,hh,mm,lst)


dd_ = np.asarray(dd_)
hh_ = np.asarray(hh_)
mm_ = np.asarray(mm_)
lst_ = np.reshape(np.asarray(lst_),-1)
ndvi_ = np.reshape(np.asarray(ndvi_),-1)
dl0 = np.reshape(dl0,-1)
ul0 = np.reshape(ul0,-1)


doy_=ope_date2DOY0(2018,12,dd_)
ind = np.asarray((doy_-1)*48+(hh_+utc_t)*2+mm_/30,dtype=np.int)
dl = dl0[ind]
ul = ul0[ind]

fvc = (ndvi_-ndvi_min)/(ndvi_max-ndvi_min)
emis = emis_s*(1-fvc)+(emis_v*fvc)
lst_mea = ((ul-dl*(1-emis))/(emis*5.67e-8))**0.25

### LST analysis
plt.plot(lst_mea,'ko-')
plt.plot(lst_,'ro-')
plt.ylim([260,335])
plt.show()


### LST evaluation
# lst_min = 270
# lst_max = 315
# lst_diff_limit = 15
# ind = (lst_mea > lst_min)*(lst_ > lst_min) *\
#       (abs(lst_mea-lst_)<lst_diff_limit)*(lst_<lst_max)*(lst_mea<lst_max)
# dif = lst_[ind]-lst_mea[ind]
# rmse = np.sqrt(np.mean(dif*dif))
# bias = np.mean(dif)
# print(bias,rmse)
# plt.plot(lst_mea[ind],lst_[ind],'ko')
# plt.plot([lst_min,lst_max],[lst_min,lst_max],'k-')
# plt.xlabel('Measured LST')
# plt.ylabel('Retrieved LST')
# plt.xlim([lst_min,lst_max])
# plt.ylim([lst_min,lst_max])
# plt.title(label)
# plt.grid()
# plt.show()