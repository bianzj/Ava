from processing_for_Aus.SLSTR import *
from processing_for_Aus.myfun import *

s = SLSTR()

ndvi_min = 0.05
ndvi_max = 0.99

indir = r'E:\AHI8_day_01\\'
outdir = r'E:\AHI8_day_01\\'

off = 12

for kday in range(1,32):


    in_file = indir + 'AHI_%02d'%kday + '_0220_red.tif'

    if os.path.exists(in_file) ==0:continue


    [red, ns, nl, temp3, geog, proj] = s.getTiffData(in_file)
    if kday == 1:
        red_base = red*1.0
    else:
        ind = red < red_base
        red_base[ind] = red[ind]




outfile = outdir + 'AHI_red.tif'
s.writeTiff(red_base, ns, nl, 1, geog, proj, outfile)



