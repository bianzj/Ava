from AHI import *
from processing_for_Aus.myfun import *


a = AHI()
indir = r'E:\Australia_day\\'
outdir = r'E:\Australia_day_result\\'

fileNames = ope_search2_rej1(indir, 'red', 'n', 'rad')
fileNum = np.size(fileNames)

off = 13
for k in range(fileNum):

    fileName = fileNames[k]
    print(fileName)
    infile1 = indir + fileName[:off] + 'refl_red_n.tif'
    infile2 = indir + fileName[:off] + 'refl_nir_n.tif'
    infile3 = indir + fileName[:off] + 'refl_red_o.tif'
    infile4 = indir + fileName[:off] + 'refl_nir_o.tif'

    outfile1 = outdir + fileName[:off] + 'ndvi_n.tif'
    outfile2 = outdir + fileName[:off] + 'bf_n.tif'
    outfile3 = outdir + fileName[:off] + 'ndvi_o.tif'
    outfile4 = outdir + fileName[:off] + 'bf_o.tif'

    red_n,ns,nl,nb,trans,proj = a.getTiffData(infile1)
    nir_n,ns,nl,nb,trans,proj = a.getTiffData(infile2)
    red_o,ns,nl,nb,trans,proj = a.getTiffData(infile3)
    nir_o,ns,nl,nb,trans,proj = a.getTiffData(infile4)

    ndvi_n = np.zeros([nl,ns])
    ndvi_o = np.zeros([nl,ns])
    bf_n = np.zeros([nl,ns])
    bf_o = np.zeros([nl,ns])

    ind = (nir_n > 0) * (red_n > 0)
    ndvi_n[ind] = (nir_n[ind]-red_n[ind])/(nir_n[ind]+red_n[ind])
    bf_n[ind] = (nir_n[ind]+red_n[ind])/2.0
    ind = (nir_o > 0) * (red_o > 0)
    ndvi_o[ind] = (nir_o[ind]-red_o[ind])/(nir_o[ind]+red_o[ind])
    bf_o[ind] = (nir_o[ind]+red_o[ind])/2.0

    a.writeTiff(ndvi_n,ns,nl,nb,trans,proj,outfile1)
    a.writeTiff(bf_n,ns,nl,nb,trans,proj,outfile2)
    a.writeTiff(ndvi_o,ns,nl,nb,trans,proj,outfile3)
    a.writeTiff(bf_o,ns,nl,nb,trans,proj,outfile4)






