
from processing_for_Aus.myfun import *

infile = 'J:\sample\AHI8_1000.tif'
[vza,ns2,nl2,nb,geog,proj] = read_gdalTiff(infile)
vza = vza/100.0
vza[vza<15] = None
plt.imshow(vza,cmap='jet',vmin = 15,vmax = 45)
plt.colorbar()
plt.show()
