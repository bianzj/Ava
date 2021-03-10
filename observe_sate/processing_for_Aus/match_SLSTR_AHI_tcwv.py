from processing_for_Aus.SLSTR import *
from processing_for_Aus.myfun import *

s = SLSTR()
indir1 = r'J:\SLSTR3A_day\\'

indir2 = r'J:\AHI8_day\\'
indir3 = r'J:\AHI8_day_lst\\'
outdir = r'J:\AHI8_day_selected\\'
month_str = '12'

fileNames = ope_search1(indir1, 'passtime')
fileNum = np.size(fileNames)

for k in range(fileNum):
    print(fileNames[k])
    fileName = fileNames[k]
    day_str = fileName[10:12]
    day = np.int(day_str)
    infile = indir1 + fileName
    [passtime, ns2, nl2, nb, geog, proj] = s.getTiffData(infile)
    ns = np.int(ns2 / 2)
    nl = np.int(nl2 / 2)
    passtime = s.resize(passtime, nl, ns)
    passtime_unique = np.unique(passtime)

    outfile1 = outdir + 'Aus_2018' + month_str + day_str + '_tcwv.tif'


    if os.path.exists(outfile1) == 1:
        [data1, nl, ns, temp3, geog, proj] = s.getTiffData(outfile1)
    else:
        data1 = np.zeros([nl, ns])


    for kk in range(len(passtime_unique)):
        passtime_one = passtime_unique[kk]
        if passtime_one == 0: continue
        hh = passtime_one // 10000
        m = ((passtime_one // 1000) % 10)
        mm = (passtime_one // 100) % 10
        mm_prop = mm / 10.0
        hh_str1 = "%02d" % hh
        mm_str1 = "%02d" % (m * 10)
        day_str1 = "%02d" % day
        hh_str2 = "%02d" % hh
        mm_str2 = "%02d" % ((m + 1) * 10)
        day_str2 = "%02d" % day

        if m == 5:
            mm_str2 = "00"
            hh_str2 = "%02d" % (hh + 1)
            if hh == 23:
                hh_str2 = "00"
                day_str2 = "%02d" % (day + 1)

        name1 = "AHI_" + day_str1 + "_" + hh_str1 + '00'
        name2 = "AHI_" + day_str2 + "_" + hh_str2 + '00'
        print(name1)
        print(name2)

        infile1 = indir2 + name1 + '_tcwv.tif'
        infile8 = indir2 + name2 + '_tcwv.tif'

        if os.path.exists(infile1) == 0: continue
        if os.path.exists(infile8) == 0: continue

        ind = (passtime == passtime_one)
        [temp1,ns,nl,nb,geog,proj] = s.getTiffData(infile1)
        [temp2,ns,nl,nb,geog,proj] = s.getTiffData(infile8)
        data1[ind] = temp1[ind]*(1-mm_prop)+temp2[ind]*mm_prop



    s.writeTiff(data1, ns, nl, nb, geog, proj, outfile1)










