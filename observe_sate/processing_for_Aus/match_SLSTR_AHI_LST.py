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
    infile = indir1+fileName
    [passtime,ns2,nl2,nb,geog,proj] = s.getTiffData(infile)
    ns = np.int(ns2/2)
    nl = np.int(nl2/2)
    passtime = s.resize(passtime,nl,ns)
    passtime_unique = np.unique(passtime)

    outfile1 = outdir + 'Aus_2018' + month_str + day_str + '_lst.tif'
    outfile2 = outdir + 'Aus_2018' + month_str + day_str + '_red.tif'
    outfile3 = outdir + 'Aus_2018' + month_str + day_str + '_nir.tif'
    outfile4 = outdir + 'Aus_2018' + month_str + day_str + '_vza.tif'
    outfile5 = outdir + 'Aus_2018' + month_str + day_str + '_psi.tif'



    if os.path.exists(outfile4) == 1:
        [data1, nl, ns, temp3, geog, proj] = s.getTiffData(outfile1)
        [data2, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile2)
        [data3, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile3)
        [data4, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile4)
        [data5, temp1, temp2, temp3, temp4, temp5] = s.getTiffData(outfile5)


    else:
        data1 = np.zeros([nl, ns])
        data2 = np.zeros([nl,ns])
        data3 = np.zeros([nl,ns])
        data4 = np.zeros([nl, ns])
        data5 = np.zeros([nl, ns])



    for kk in range(len(passtime_unique)):
        passtime_one = passtime_unique[kk]
        if passtime_one ==0: continue
        hh = passtime_one//10000
        m = ((passtime_one//1000) % 10)
        mm = (passtime_one//100) % 10
        mm_prop = mm/10.0
        hh_str1 = "%02d"%hh
        mm_str1 = "%02d"%(m*10)
        day_str1 = "%02d"%day
        hh_str2 = "%02d"%hh
        mm_str2 = "%02d"%((m+1)*10)
        day_str2 = "%02d"%day

        if m == 5:
            mm_str2 = "00"
            hh_str2 = "%02d"%(hh+1)
            if hh == 23:
                hh_str2 = "00"
                day_str2 = "%02d"%(day+1)


        name1 = "AHI_"+day_str1+"_"+hh_str1+mm_str1
        name2 = "AHI_"+day_str2+"_"+hh_str2+mm_str2
        print(name1)
        print(name2)


        infile1 = indir3+name1+'_lst.tif'
        infile3 = indir2+name1+'_red.tif'
        infile4 = indir2+name1+'_nir.tif'
        infile5 = indir2 + name1 + '_vza.tif'
        infile6 = indir2 + name1 + '_psi.tif'


        infile8 = indir3+name2+'_lst.tif'
        infile10 = indir2+name2+'_red.tif'
        infile11 = indir2+name2+'_nir.tif'
        infile12 = indir2 + name2 + '_vza.tif'
        infile13 = indir2 + name2 + '_psi.tif'


        if os.path.exists(infile1) == 0:continue
        if os.path.exists(infile8) == 0: continue
        if os.path.exists(infile4) == 0: continue
        if os.path.exists(infile11) == 0: continue

        ind = (passtime == passtime_one)
        # [temp1,ns,nl,nb,geog,proj] = s.getTiffData(infile1)
        # [temp2,ns,nl,nb,geog,proj] = s.getTiffData(infile8)
        # data1[ind] = temp1[ind]*(1-mm_prop)+temp2[ind]*mm_prop
        #
        # [temp1, ns, nl, nb, geog, proj] = s.getTiffData(infile3)
        # [temp2, ns, nl, nb, geog, proj] = s.getTiffData(infile10)
        # data2[ind] = temp1[ind] * (1 - mm_prop) + temp2[ind] * mm_prop
        # # s.writeTiff(data2,ns,nl,nb,geog,proj,outfile2)
        #
        # [temp1, ns, nl, nb, geog, proj] = s.getTiffData(infile4)
        # [temp2, ns, nl, nb, geog, proj] = s.getTiffData(infile11)
        # data3[ind] = temp1[ind] * (1 - mm_prop) + temp2[ind] * mm_prop

        [temp1, ns, nl, nb, geog, proj] = s.getTiffData(infile5)
        [temp2, ns, nl, nb, geog, proj] = s.getTiffData(infile12)
        data4[ind] = temp1[ind] * (1 - mm_prop) + temp2[ind] * mm_prop

        [temp1, ns, nl, nb, geog, proj] = s.getTiffData(infile6)
        [temp2, ns, nl, nb, geog, proj] = s.getTiffData(infile13)
        data5[ind] = temp1[ind] * (1 - mm_prop) + temp2[ind] * mm_prop



    # s.writeTiff(data1, ns, nl, nb, geog, proj, outfile1)
    # s.writeTiff(data2, ns, nl, nb, geog, proj, outfile2)
    # s.writeTiff(data3, ns, nl, nb, geog, proj, outfile3)
    s.writeTiff(data4, ns, nl, nb, geog, proj, outfile4)
    s.writeTiff(data5, ns, nl, nb, geog, proj, outfile5)









