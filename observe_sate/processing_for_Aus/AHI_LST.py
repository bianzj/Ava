from processing_for_Aus.SLSTR import *
from processing_for_Aus.myfun import *

s = SLSTR()

ndvi_min = 0.05
ndvi_max = 0.99

in_file_red ='E:\AHI8_day_09\AHI_red.tif'
[red_base, ns, nl, nb, geog, proj] = s.getTiffData(in_file_red)

indir1 = r'E:\AHI8_day_09\\'
indir2 = r'E:\AHI8_day_LST_09\\'
indir3 = r'G:\LSE\\'
indir4 = r'G:\class\\'
indir5 = r'E:\AHI8_sample\\'
outdir = r'E:\AHI8_day_LST_09\\'

fileNames = ope_search2(indir1, 'B1', 'AHI')
fileNum = np.size(fileNames)

infile10 = indir3 + 'LSE1_Aus_1000.tif'
infile11 = indir3 + 'LSE2_Aus_1000.tif'
infile13 = indir3 + r'ICBP_Emi.txt'
infile12 = indir4 + 'class_2016_Aus_1000.tif'

infile14 = indir5 + 'AHI_SW_Coeff_Day.txt'
coeffs = s.getDatafromTxt(infile14, 0, 15)
coeffs = coeffs[:, 2:10]

[emis1_s, ns, nl, nb, geog, proj] = s.getTiffData(infile10)
emis1_s = emis1_s / 1000.0
[emis2_s, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile11)
emis2_s = emis2_s / 1000.0
[classif, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile12)
classif = np.asarray(classif, np.int)
emis_v = s.getDatafromTxt(infile13, 0, 2)

off = 12
for k in range(0,fileNum):
    fileName = fileNames[k]
    print(fileName)

    infile1 = indir1 + fileName[:off] + 'B1.tif'
    infile2 = indir1 + fileName[:off] + 'B2.tif'
    infile5 = indir1 + fileName[:off] + 'red.tif'
    infile6 = indir1 + fileName[:off] + 'nir.tif'
    infile8 = indir1 + fileName[:off] + 'vza.tif'
    infile7 = indir1 + fileName[:off-3] + '00_tcwv.tif'

    if os.path.exists(infile7) ==0:continue

    dd = fileName[4:6]
    hh = fileName[7:9]
    mm_prop = np.int(fileName[9:11])/60.0

    if (hh =='23'):
        dd = '%02d'%(np.int(dd)+1)
        hh = '00'
    else:
        hh = '%02d'%(np.int(hh))

    infile9 = indir1 + 'AHI_'+dd +'_'+hh+ '00_tcwv.tif'

    [BT8, ns, nl, nb, geog, proj] = s.getTiffData(infile1)
    [BT9, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile2)
    [tcw1, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile7)
    [tcw2, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile9)
    [vza, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile8)
    [red, ns2, nl2, temp3, temp4, temp5] = s.getTiffData(infile5)
    [nir, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(infile6)

    tcw = tcw1 * (1 - mm_prop) + tcw2 * mm_prop
    ndvi = (nir-red)/(nir+red)
    fvc = (ndvi - ndvi_min) / (ndvi_max - ndvi_min)
    ind = np.isnan(ndvi)
    fvc[ind] = 0

    fvc[fvc > 1] = 1
    fvc[fvc < 0] = 0

    emis1_v = emis_v[classif, 0] / 1000.0
    emis2_v = emis_v[classif, 1] / 1000.0
    emis1_v = ope_resizeData(emis1_v,ns,nl)
    emis2_v = ope_resizeData(emis2_v, ns,nl)
    emis1_s = ope_resizeData(emis1_s, ns,nl)
    emis2_s = ope_resizeData(emis2_s,ns,nl)

    emis8 = emis1_s*(1-fvc) + emis1_v*fvc
    emis9 = emis2_s*(1-fvc) + emis2_v*fvc
    # emis8 = ope_resizeData(emis8,ns,nl)
    # emis9 = ope_resizeData(emis9,ns,nl)

    emis = (emis8 + emis9) * 0.5
    demis = emis8 - emis9
    x3 = np.zeros([nl, ns])
    x4 = np.zeros([nl, ns])
    x6 = np.zeros([nl, ns])
    x7 = np.zeros([nl, ns])
    x1 = 1
    x2 = (BT8 + BT9) / 2.0
    ind = emis != 0
    x3[ind] = (1 - emis[ind]) / emis[ind] * x2[ind]
    x4[ind] = demis[ind] / (emis[ind] * emis[ind]) * x2[ind]
    x5 = (BT8 - BT9) / 2
    x6[ind] = (1 - emis[ind]) / emis[ind] * x5[ind]
    x7[ind] = demis[ind] / (emis[ind] * emis[ind]) * x5[ind]
    x8 = (BT8 - BT9) * (BT8 - BT9)

    wvc_ = np.asarray([0, 1.25, 2.25, 3.25, 4.25, 5.25, 100])
    vza_ = np.asarray([0, 3, 14.9, 38.6, 44.5, 51.2, 58, 65, 70, 75, 80])

    LST = np.zeros([nl, ns])
    for kwvc in range(len(wvc_) - 1):
        for kvza in range(len(vza_) - 1):

            ind = (tcw >= wvc_[kwvc]) * (tcw < wvc_[kwvc + 1]) * (vza >= vza_[kvza]) * (vza < vza_[kvza + 1]) * (
                    BT8 > 0.5) * (emis8 > 0.5)
            if (np.sum(ind) <= 0): continue

            ##########################################
            ### the normal result
            ##########################################

            ind1 = kvza * 6 + kwvc
            ind2 = (kvza + 1) * 6 + kwvc
            coeff1 = coeffs[ind1, :]
            coeff2 = coeffs[ind2, :]
            vzatemp = vza[ind]
            prop2 = (vzatemp - vza_[kvza]) / (vza_[kvza + 1] - vza_[kvza])
            prop1 = 1 - prop2
            LSTtemp = x1      * (prop1 * coeff1[0] + coeff2[0] * prop2) + \
                      x2[ind] * (prop1 * coeff1[1] + coeff2[1] * prop2) + \
                      x3[ind] * (prop1 * coeff1[2] + coeff2[2] * prop2) + \
                      x4[ind] * (prop1 * coeff1[3] + coeff2[3] * prop2) + \
                      x5[ind] * (prop1 * coeff1[4] + coeff2[4] * prop2) + \
                      x6[ind] * (prop1 * coeff1[5] + coeff2[5] * prop2) + \
                      x7[ind] * (prop1 * coeff1[6] + coeff2[6] * prop2) + \
                      x8[ind] * (prop1 * coeff1[7] + coeff2[7] * prop2)

            LST[ind] = LSTtemp



    ind = LST > 355
    LST[ind] = 355
    ind = LST < 275
    LST[ind] = 0
    ind = ((red - red_base) > 0.050) *(LST < 285)
    LST[ind] = 0

    outfile = outdir + fileName[:off]+'LST.tif'
    s.writeTiff(LST, ns, nl, 1, geog, proj, outfile)



