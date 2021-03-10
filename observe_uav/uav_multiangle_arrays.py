from util.myfun import *



####################################################
indir = r'G:\UAVpro\20200604huailai\1201\\'
# ### 06041201 track 1
track = 1
start = np.asarray([1680,5925])
end = np.asarray([2410+1,6660+1])
# ### 06041201 track 1
# track = 2
# start = np.asarray([2525,5050])
# end = np.asarray([3275+1,5825+1])
# # # ### 06041201 track 1
# track = 3
# start = np.asarray([3375,4215])
# end = np.asarray([4125+1,4950+1])
###################################################
# indir = r'G:\UAVpro\20200604huailai\1401\\'
# # ### 06041401 track 1
# track = 1
# start = np.asarray([2010,6280])
# end = np.asarray([2710+1,6900+1])
# # ### 06041401 track 1
# track = 2
# start = np.asarray([2900,5340])
# end = np.asarray([3560+1,6020+1])
# # ### 06041401 track 1
# track = 3
# start = np.asarray([3725,4510])
# end = np.asarray([4400+1,5170+1])
###################################################
# indir = r'G:\UAVpro\20200615huailai\1001\\'
### 06151001 track 1
# track = 1
# start = np.asarray([1720,8375])
# end = np.asarray([2320+1,8975+1])
# ### 06151001 track 1
# track = 2
# start = np.asarray([2460,7610])
# end = np.asarray([3070+1,8215+1])
# ### 06151001 track 1
# track = 3
# start = np.asarray([3275,6870])
# end = np.asarray([3835+1,7445+1])
###################################################

#################### reference orth image
infile = indir + r'TIR_ref.tif'
temp,ref_width,ref_height,ref_band,ref_geotrans,ref_proj,refdataset = read_image_gdal(infile)
################## GPS DATA
infile = indir + r'cameras.txt'
fileNamePos,pos = read_txt_array(infile,2)
################ INPUT FILE
indir = indir + r'TIR_result\\sample1\\'
outdirmed = indir + r'TIR_result\\sample1\\'

numVZA = 60
numVAA = 360
outdir = indir
for kline in range(2):

    dbt = np.zeros([numVZA,numVAA])
    count = np.zeros([numVZA,numVAA], dtype=np.int)

    for kimage in range(start[kline], end[kline],5):

        ####################################################
        ##### this file
        infile = "%05d"%(kimage)
        fileName = indir + infile+'.tif'
        if 1-os.path.exists(fileName):
            continue
        print(infile)
        #######################################################
        ### camera position in this file
        kpos = fileNamePos.index(infile)
        if kpos < 0: continue
        camera_lat = pos[kpos, 1]
        camera_lon = pos[kpos, 0]
        camera_height = pos[kpos, 2]


        #######################################################
        #### position for each pixel in this file
        dataset1 = gdal.Open(fileName)
        temp_width1 = dataset1.RasterXSize
        temp_height1 = dataset1.RasterYSize
        temp_geotrans1 = dataset1.GetGeoTransform()
        temp_proj1 = dataset1.GetProjection()
        temp_data1 = dataset1.ReadAsArray(0, 0, temp_width1, temp_height1)
        temp_x1 = np.zeros([temp_height1, temp_width1], dtype=np.int)
        temp_y1 = np.zeros([temp_height1, temp_width1], dtype=np.int)
        temp = np.arange(0, temp_width1)
        for k in range(temp_height1):
            temp_x1[k, :] = k
            temp_y1[k, :] = temp
        temp_geox1, temp_geoy1 = imagexy2geo(dataset1, temp_x1, temp_y1)
        temp_lon1, temp_lat1 = geo2lonlat(dataset1, temp_geox1, temp_geoy1)
        ind = (temp_data1 < 3000)
        temp_data1[ind] = 0
        ind = (temp_data1 > 9000)
        temp_data1[ind] = 0

        #########################################################
        #### calculate viewing angle for each pixel in this file
        camera_geox1, camera_geoy1 = lonlat2geo(dataset1, camera_lon, camera_lat)
        camera_y1, camera_x1 = geo2imagexy(dataset1, camera_geox1, camera_geoy1)

        ref_geogx1, ref_geogy1 = lonlat2geo(refdataset, temp_lon1, temp_lat1)
        ref_y, ref_x = geo2imagexy(refdataset, ref_geogx1, ref_geogy1)
        distance = np.sqrt((temp_geox1 - camera_geox1) * (temp_geox1 - camera_geox1) +
                           (temp_geoy1 - camera_geoy1) * (temp_geoy1 - camera_geoy1))
        temp_zenith = np.arctan(distance / camera_height) / np.pi * 180.0
        temp_azimuth = np.arctan2((temp_lon1 - camera_lon), (temp_lat1 - camera_lat)) / np.pi * 180.0
        temp_azimuth = ((temp_azimuth + 360.0 + 180.0) % 360.0)
        temp_azimuth = np.reshape(temp_azimuth, [temp_height1, temp_width1])

        temp_zenith = np.asarray(temp_zenith,dtype=np.int)
        temp_azimuth = np.asarray(temp_azimuth,dtype= np.int)

        ind = (temp_data1 > 3000) *(temp_data1 < 9000)

        dbt[temp_zenith[ind],temp_azimuth[ind]] = dbt[temp_zenith[ind],temp_azimuth[ind]] + temp_data1[ind]
        count[temp_zenith[ind], temp_azimuth[ind]] = count[temp_zenith[ind], temp_azimuth[ind]] +1



    ind = count > 0
    dbt[ind] = dbt[ind] / count[ind]

    outfile = outdir + r'dbt%1d' % track + '_%1d.tif'%kline
    write_image_gdal(dbt, numVAA, numVZA, 1, '', '', outfile)


