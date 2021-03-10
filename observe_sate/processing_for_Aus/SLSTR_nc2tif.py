from processing_for_Aus.SLSTR import *
from processing_for_Aus.myfun import *

wdir = r'G:\SLSTR3B\\'

s = SLSTR()

fileNames = ope_search1_rej1(wdir,'S3B','zip')
fileNum = np.size(fileNames)

for k in range(0,fileNum):

    fileName = fileNames[k].strip()
    fileDir = wdir + fileName
    print(k,fileDir)
    ccc = s.ifallfilesexists(fileDir)
    if(ccc==0):continue
    s.getBLUEdatafromNc_nadir(fileDir)

    s.getVNIRdatafromNc_nadir(fileDir)
    s.getTIRdatafromNc_nadir(fileDir)

    s.getVNIRdatafromNc_obliq(fileDir)
    s.getTIRdatafromNc_obliq(fileDir)

    s.getGeodatafromNc_obliq(fileDir)
    s.getGeodatafromNc_nadir(fileDir)

    infile = fileDir+r'\\lat.tif'
    outfile = fileDir+r'\\lat_v.tif'
    s.resize_file(infile,outfile,2)

    infile = fileDir+r'\\lon.tif'
    outfile = fileDir+r'\\lon_v.tif'
    s.resize_file(infile,outfile,2)

    infile = fileDir+r'\\lat_o.tif'
    outfile = fileDir+r'\\lat_o_v.tif'
    s.resize_file(infile,outfile,2)

    infile = fileDir+r'\\lon_o.tif'
    outfile = fileDir+r'\\lon_o_v.tif'
    s.resize_file(infile,outfile,2)
