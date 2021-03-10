from processing_for_Aus.myfun import *

indir1 = r'J:\multiple_result\\'
indir2 = r'J:\SLSTR3A_day\\'
indir3 = r'J:\AHI8_day_lst\\'
indir4 = r'J:\AHI8_day\\'
outdir = r'J:\temporal_result\\'

infile1 = indir1 + 'Aus_20180102_9VTa.tif'
infile2 = indir1 + 'Aus_20180102_9VTb.tif'
infile3 = indir1 + 'Aus_20180102_9VTc.tif'
[coeffa, nl, ns, temp3, geog, proj] = getTiffData(infile1)
[coeffb, nl, ns, temp3, geog, proj] = getTiffData(infile2)
[coeffc, nl, ns, temp3, geog, proj] = getTiffData(infile3)

infile4 = indir2 + 'Aus_20181202_passtime.tif'
[passtime,nl,ns,temp3,geog,proj] = getTiffData(infile4)
passtime = ope_resizeData(passtime,nl/2,ns/2)
passtime_unique = np.unique(passtime)

infile8 = indir2 + 'Aus_20181201_refl_red_n.tif'
infile9 = indir2 + 'Aus_20181201_refl_nir_n.tif'
[red, nl, ns, temp3, geog, proj] = getTiffData(infile8)
[nir, nl, ns, temp3, geog, proj] = getTiffData(infile9)
ndvi0 = (nir-red)/(nir+red)
bf0 = (nir+red)/2.0
ndvi0 = ope_resizeData(ndvi0,nl/2,ns/2)
bf0 = ope_resizeData(bf0,nl/2,ns/2)


iflast = 0
ifnext = 0
ifcurr = 1
dd = 2

if iflast == 1:
    kday = dd-1
    for khour in range(21,24):
        for kmin in range(0,6):

            infile5 = indir4 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_red.tif'
            infile6 = indir4 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_nir.tif'
            infile7 = indir3 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_lst.tif'
            outfile = outdir + 'AHI_%02d_' % dd + '%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_test.tif'

            if not os.path.exists(infile5): continue

            [red, nl, ns, temp3, geog, proj] = getTiffData(infile5)
            [nir, nl, ns, temp3, geog, proj] = getTiffData(infile6)
            [lst_obs, nl, ns, temp3, geog, proj] = getTiffData(infile7)
            ndvi = (nir-red)/(nir+red)
            bf = (nir+red)/2.0
            lst_sim = lst_obs-coeffa*(ndvi-ndvi0)-coeffb*(bf-bf0)

            distance_ = np.zeros([ns,nl])
            for kk in range(len(passtime_unique)):
                passtime_one = passtime_unique[kk]
                if passtime_one == 0: continue
                hh = passtime_one // 10000
                m = ((passtime_one // 1000) % 10)
                mm = (passtime_one // 100) % 10
                distance = (24 - khour - 1) * 60.0 + (6 - kmin - 1) * 10  + hh * 60 + m * 10 + mm
                ind = passtime == passtime_one
                distance_[ind] = distance


            result = np.stack([lst_sim,distance_])
            writeTiff(result, ns, nl, 2, geog, proj, outfile)

    iflast = 0

if ifcurr == 1:
    kday = dd
    for khour in range(0, 3):
        for kmin in range(0, 6):

            infile5 = indir4 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_red.tif'
            infile6 = indir4 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_nir.tif'
            infile7 = indir3 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_lst.tif'
            outfile = outdir + 'AHI_%02d_' % dd + '%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_test.tif'

            if not os.path.exists(infile5): continue
            [red, nl, ns, temp3, geog, proj] = getTiffData(infile5)
            [nir, nl, ns, temp3, geog, proj] = getTiffData(infile6)
            [lst_obs, nl, ns, temp3, geog, proj] = getTiffData(infile7)
            ndvi = (nir - red) / (nir + red)
            bf = (nir + red) / 2.0
            lst_sim = lst_obs - coeffa * (ndvi - ndvi0) - coeffb * (bf - bf0)

            distance_ = np.zeros([ns, nl])
            for kk in range(len(passtime_unique)):
                passtime_one = passtime_unique[kk]
                if passtime_one == 0: continue
                hh = passtime_one // 10000
                m = ((passtime_one // 1000) % 10)
                mm = (passtime_one // 100) % 10
                distance = abs((khour*60+kmin*10)-(hh*60+m*10+mm))
                ind = passtime == passtime_one
                distance_[ind] = distance


            result = np.stack([lst_sim, distance_])
            writeTiff(result, ns, nl, 2, geog, proj, outfile)
    ifcurr = 0

if ifnext == 1:
    kday = dd+1
    for khour in range(0, 9):
        for kmin in range(0, 6):

            infile5 = indir4 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_red.tif'
            infile6 = indir4 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_nir.tif'
            infile7 = indir3 + 'AHI_%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_lst.tif'
            outfile = outdir + 'AHI_%02d_' % dd + '%02d_' % kday + '%02d' % khour + '%1d0' % kmin + '_test.tif'
            if not os.path.exists(infile5): continue
            [red, nl, ns, temp3, geog, proj] = getTiffData(infile5)
            [nir, nl, ns, temp3, geog, proj] = getTiffData(infile6)
            [lst_obs, nl, ns, temp3, geog, proj] = getTiffData(infile7)
            ndvi = (nir - red) / (nir + red)
            bf = (nir + red) / 2.0
            lst_sim = lst_obs - coeffa * (ndvi - ndvi0) - coeffb * (bf - bf0)

            distance_ = np.zeros([ns, nl])
            for kk in range(len(passtime_unique)):
                passtime_one = passtime_unique[kk]
                if passtime_one == 0: continue
                hh = passtime_one // 10000
                m = ((passtime_one // 1000) % 10)
                mm = (passtime_one // 100) % 10

                distance = (24-hh-1)*60.0 + (6-m-1) * 10 +(10-mm) + khour * 60 + kmin *10
                ind = passtime == passtime_one
                distance_[ind] = distance


            result = np.stack([lst_sim, distance_])
            writeTiff(result, ns, nl, 2, geog, proj, outfile)
    ifnext = 0
