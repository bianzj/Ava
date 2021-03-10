from processing_for_Aus.myfun import *
from processing_for_Aus.SLSTR import *

wdir = r'G:\SLSTR3A\\'

s = SLSTR()
fileNames = ope_searchfile_accept1_rej1(wdir,'S3A','ZIP')

for k in range(len(fileNames)):

    fileName = fileNames[k]
    fileDir = wdir+ fileName

    UTCtime = fileName[25:27]
    if (UTCtime > '05') & (UTCtime < '20'): continue

    print(fileDir)
    ccc = s.ifallfilesexists(fileDir)
    if(ccc==0):continue


    # file = 'refl_blue'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon.tif'
    # yfile = fileDir+'\\lat.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'psi'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon.tif'
    yfile = fileDir+'\\lat.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'sza'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon.tif'
    yfile = fileDir+'\\lat.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'vza'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon.tif'
    yfile = fileDir+'\\lat.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'bt8'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon.tif'
    yfile = fileDir+'\\lat.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'bt9'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon.tif'
    yfile = fileDir+'\\lat.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    # file = 'tcw'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon.tif'
    # yfile = fileDir+'\\lat.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'cloud'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon.tif'
    yfile = fileDir+'\\lat.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'refl_red'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon_v.tif'
    yfile = fileDir+'\\lat_v.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_v_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    file = 'refl_nir'
    datafile = fileDir+'\\'+file+'.tif'
    xfile = fileDir+'\\lon_v.tif'
    yfile = fileDir+'\\lat_v.tif'
    vrtfile = fileDir+'\\'+file+'.vrt'
    s.writevrt_v_nadir(vrtfile,datafile,xfile,yfile)
    outfile = fileDir+'\\'+file+'_proj.tif'
    dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    # file = 'refl_ref'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_v.tif'
    # yfile = fileDir+'\\lat_v.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_v_nadir(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    # file = 'rad_red'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_v.tif'
    # yfile = fileDir+'\\lat_v.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_v_nadir(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    #
    # file = 'rad_nir'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_v.tif'
    # yfile = fileDir+'\\lat_v.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_v_nadir(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    # file = 'refl_red_o'
    # datafile = fileDir + '\\' + file + '.tif'
    # xfile = fileDir + '\\lon_o_v.tif'
    # yfile = fileDir + '\\lat_o_v.tif'
    # vrtfile = fileDir + '\\' + file + '.vrt'
    # s.writevrt_v_nadir(vrtfile, datafile, xfile, yfile)
    # outfile = fileDir + '\\' + file + '_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    #
    # file = 'refl_nir_o'
    # datafile = fileDir + '\\' + file + '.tif'
    # xfile = fileDir + '\\lon_o_v.tif'
    # yfile = fileDir + '\\lat_o_v.tif'
    # vrtfile = fileDir + '\\' + file + '.vrt'
    # s.writevrt_v_nadir(vrtfile, datafile, xfile, yfile)
    # outfile = fileDir + '\\' + file + '_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    #
    # file = 'refl_ref_o'
    # datafile = fileDir + '\\' + file + '.tif'
    # xfile = fileDir + '\\lon_o_v.tif'
    # yfile = fileDir + '\\lat_o_v.tif'
    # vrtfile = fileDir + '\\' + file + '.vrt'
    # s.writevrt_v_nadir(vrtfile, datafile, xfile, yfile)
    # outfile = fileDir + '\\' + file + '_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    #
    # file = 'rad_red_o'
    # datafile = fileDir + '\\' + file + '.tif'
    # xfile = fileDir + '\\lon_o_v.tif'
    # yfile = fileDir + '\\lat_o_v.tif'
    # vrtfile = fileDir + '\\' + file + '.vrt'
    # s.writevrt_v_nadir(vrtfile, datafile, xfile, yfile)
    # outfile = fileDir + '\\' + file + '_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    #
    # file = 'rad_nir_o'
    # datafile = fileDir + '\\' + file + '.tif'
    # xfile = fileDir + '\\lon_o_v.tif'
    # yfile = fileDir + '\\lat_o_v.tif'
    # vrtfile = fileDir + '\\' + file + '.vrt'
    # s.writevrt_v_nadir(vrtfile, datafile, xfile, yfile)
    # outfile = fileDir + '\\' + file + '_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    # file = 'cloud_o'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_o.tif'
    # yfile = fileDir+'\\lat_o.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_obliq(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)

    # file = 'bt8_o'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_o.tif'
    # yfile = fileDir+'\\lat_o.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_obliq(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    #
    # file = 'bt9_o'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_o.tif'
    # yfile = fileDir+'\\lat_o.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_obliq(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    # #
    # file = 'vza_o'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_o.tif'
    # yfile = fileDir+'\\lat_o.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_obliq(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    # #
    # file = 'sza_o'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_o.tif'
    # yfile = fileDir+'\\lat_o.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_obliq(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)
    # #
    # file = 'psi_o'
    # datafile = fileDir+'\\'+file+'.tif'
    # xfile = fileDir+'\\lon_o.tif'
    # yfile = fileDir+'\\lat_o.tif'
    # vrtfile = fileDir+'\\'+file+'.vrt'
    # s.writevrt_obliq(vrtfile,datafile,xfile,yfile)
    # outfile = fileDir+'\\'+file+'_proj.tif'
    # dst_ds = gdal.Warp(outfile, vrtfile, geoloc=True, resampleAlg=gdal.GRIORA_Cubic)