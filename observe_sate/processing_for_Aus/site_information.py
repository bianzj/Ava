from processing_for_Aus.myfun import *


indir = 'J:\measurement\\'

infiles = ope_search1(indir,'nc')

filenum = len(infiles)

for k in range(filenum):
    infile = indir + infiles[k]
    site_name=os.path.splitext(infiles[k])[0]
    dataset = netCDF4.Dataset(infile)
    print(site_name,' ',dataset.latitude,' ',dataset.longitude,' ',dataset.vegetation)