from myfun_image import *
from myfun_file import *
from SLSTR_constant import *
from myfun_sci import *





###################################################
####  LST number statistic
###################################################

infile_type = 'G:/base/type.tif'
symbols = ['S3A','S3B','VJ1','VNP','MOD','MYD']
indir_coeff = r'G:/coeffa/'

lon = 100.3854
lat = 38.8263
numx = 365
offs = [4,6,10]
dataset = gdal.Open(infile_type)
im_width = dataset.RasterXSize
im_height = dataset.RasterYSize
im_bands = dataset.RasterCount
data = dataset.ReadAsArray(0, 0, im_width, im_height)
im_geotrans = dataset.GetGeoTransform()
im_proj = dataset.GetProjection()
temp1, temp2 = lonlat2geo(dataset, lat, lon)
imagey, imagex = geo2imagexy(dataset, temp1, temp2)
imagex = np.int(imagex)
imagey = np.int(imagey)


numy = len(offs)
coeffa = np.zeros([numx,numy])
for ky in range(numy):
    off = offs[ky]
    for kx in range(off,numx-off):

        infile = indir_coeff+'coeffa_2019%03d'%kx+'_%02d'%off+'.tif'
        if os.path.exists(infile) == 0: continue
        print(infile)
        [data, ns, nl, nb, geog, proj] = read_image_gdal(infile)
        coeffa[kx,ky] = data[imagex,imagey]


plt.plot(coeffa[:,0],'k')
plt.plot(coeffa[:,1],'r')
plt.plot(coeffa[:,2],'b')
plt.show()

