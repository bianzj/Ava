from util.myfun import *
from geopandas import *
import rasterio as rio
from geopandas import GeoSeries
import rasterio.mask
from rasterio.warp import (reproject,Resampling, transform_bounds,calculate_default_transform as calcdt)

sample = 'sample1'
band = 'nir'

infile1 = r'F:\UAVpro\20200604huailai\1201\\'+sample+'_'+band+r'_result\track21.txt'
infile2 = r'F:\UAVpro\20200604huailai\1201\\'+sample+'_'+band+r'_result\track12.txt'
mydata1 = read_txt_array(infile1,0,6)
mydata2 = read_txt_array(infile2,0,6)

sample = mydata1[:,2]
white = mydata1[:,5]
black = mydata1[:,3]
gray = mydata1[:,4]
vza = mydata1[:,1]
vaa = mydata1[:,0]

prop = ( sample - black)/(gray-black)


# sample = mydata2[:,2]
# white = mydata2[:,5]
# black = mydata2[:,3]
# gray = mydata2[:,4]
# vza = mydata2[:,1]
# vaa = mydata2[:,0]
#
# prop = ( sample - black)/(white-black)

plt.plot(vza,prop)
plt.show()