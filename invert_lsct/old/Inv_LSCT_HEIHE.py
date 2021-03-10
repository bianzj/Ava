import time
from old.SLSTR import *

wdir = r'G:\Work\LSCT_SLSTR\\'

s = SLSTR()


dataTimeOld = 0
imageTimeOld = 0
x = os.listdir(wdir+r'heihe_result\\')
num_x = len(x)
xx = os.listdir(wdir+r'heihe\\')
num_xx = len(xx)
off1 = 0
off2 = np.int(num_x/5*1)
off3 = np.int(num_x/5*2)
off4 = np.int(num_x/5*3)
off5 = np.int(num_x/5*4)
off6 = num_x

lat_ = np.asarray([38.0473,38.8399,37.838,38.0142,
                   38.8555,38.7583,38.97514,
                   41.99317639,42.0012,42.15555])
lon_ = np.asarray([100.4643,98.94060,101.1160,100.2421,
                   100.3722,100.30556 , 100.4464,
                   101.1238647,101.1374,100.9833])

lat_ = np.asarray([38.0473,38.8399,37.838,38.0142,
                   38.8555,38.7659,38.97514,
                   41.99317639,42.0012,42.15555])
lon_ = np.asarray([100.4643,98.94060,101.1160,100.2421,
                   100.3722,100.3201 , 100.4464,
                   101.1238647,101.1374,100.9833])


num_points = len(lat_)

### 100.3201   38.7659

for i in range(off1,off6):
# for i in range(off4, off6):
    ix = x[i]
    dataTime = ix[6:14]
    imageTime = ix[15:21]
    timetime = dataTime+'_'+imageTime
    # if dataTime != '20170715':continue
    # if np.int(dataTime) < 20170815 : continue
    if (dataTime==dataTimeOld) and (imageTime==imageTimeOld):continue
    dataTimeOld = dataTime
    imageTimeOld = imageTime
    print(ix)

    print(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime( time.time())))

    infile = wdir + r'heihe_result\heihe_' + timetime + '_lst_n.tif'
    ccc = os.path.exists(infile)
    if (ccc == 0): continue
    infile = wdir + r'heihe_result\heihe_' + timetime + '_lst_o.tif'
    ccc = os.path.exists(infile)
    if (ccc == 0): continue

    medName = r'heihe\heihe_'
    resultName = r'heihe_result\heihe_'
    c1 = os.path.exists(wdir + medName + timetime + '_bt8_o.tif')
    c2 = os.path.exists(wdir + medName + timetime + '_bt8_n.tif')
    c3 = os.path.exists(wdir + medName + timetime + '_refl_ref_n.tif')
    c4 = os.path.exists(wdir + medName + timetime + '_refl_ref_o.tif')
    c5 = os.path.exists(wdir + medName + timetime + '_refl_blue_n.tif')
    c6 = os.path.exists(wdir + resultName + timetime + '_ndvi_o.tif')
    ccc = c1 * c2 * c3 * c4 * c5*c6
    if ccc == 0: continue


    s.inversionLSCT(wdir,timetime,medName, resultName,2)

    # s.inversionLSCT_npl_point_raw(wdir,timetime,lat_,lon_,medName, resultName,2,0)
    # s.inversionLSCT_agl_point_raw(wdir,timetime,lat_,lon_,medName, resultName,2,0,aodcheck = 1)
    # s.inversionLSCT_npl_agl_point_raw(wdir,timetime,lat_,lon_,medName, resultName,2,0,aodcheck=1)
    # s.inversionLSCT_npl_agl_point_raw(wdir,timetime,lat_,lon_,medName, resultName,2,0,aodcheck=1)
    #
    # s.inversionLSCT_npl_agl_point_pur(wdir,timetime,lat_,lon_,2,0)
    # s.inversionLSCT_npl_agl_point_sgl(wdir,timetime,lat_,lon_,2,0)
    # s.inversionLSCT_npl_agl_point(wdir, timetime, lat_, lon_, 2, 0)
    # s.inversionLSCT_npl_agl_point_test(wdir,timetime,lat_,lon_,2,0)
    # s.inversionLSCT_agl_point(wdir, timetime, lat_, lon_, 0, 0)
    # s.inversionLSCT_npl_point(wdir, timetime, lat_, lon_, 2, 0)










