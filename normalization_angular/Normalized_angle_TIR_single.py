from util.myfun import *
from fitting_kernels import *
from SLSTR_constant import *
# from myfun_plot import *
os.environ['PROJ_LIB'] = r'C:\Users\zzz\anaconda3\envs\patrol\Library\share\proj'


indir_s3a = r'G:/S3A_tif/'
indir_s3b = r'G:/S3B_tif/'
indirs2 = [indir_s3a,indir_s3b]
indir_s3a = r'G:/S3A_LST/'
indir_s3b = r'G:/S3B_LST/'
indirs1 = [indir_s3a,indir_s3b]

outdira = 'G:/coeffa/'
outdirbk = 'G:/coeffbk/'


infile_type = 'G:/base/type.tif'
[type, ns, nl, nb, geog, proj] = read_image_gdal(infile_type)
symbols = ['S3A','S3B']

off = 4
num = 365
rd = np.pi/180.0

ifcoeffa = 0
ifcoeffb = 1
th = 270*270.0

if ifcoeffa == 1:
    iffirst = 1
    mark = 0
    nb = off*2+1
    lst_nadir = np.zeros([nb * 2, nl, ns])
    lst_obliq = np.zeros([nb * 2, nl, ns])
    vza_nadir = np.zeros([nb * 2, nl, ns])
    sza_nadir = np.zeros([nb * 2, nl, ns])
    psi_obliq = np.zeros([nb * 2, nl, ns])
    psi_nadir = np.zeros([nb * 2, nl, ns])
    vza_obliq = np.zeros([nb * 2, nl, ns])
    sza_obliq = np.zeros([nb * 2, nl, ns])
    for kdoy in range(off,num-off):
        print(kdoy)
        a = np.zeros([nl,ns])
        if iffirst == 1:
            for ksate in range(2):
                indir1 = indirs1[ksate]
                indir2 = indirs2[ksate]
                symbol = symbols[ksate]
                kstep = nb*ksate
                for krun in range(nb):
                    kday = kdoy-off+krun
                    infile_nadir_lst = indir1 + symbol +'_2019%03d'%(kday+1)+'_night_lst.tif'
                    infile_obliq_lst = indir1 + symbol +'_2019%03d'%(kday+1)+'_night_lst_obliq.tif'
                    infile_nadir_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_night_vza.tif'
                    infile_obliq_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_night_vza_obliq.tif'
                    infile_nadir_sza = indir2 + symbol +'_2019%03d'%(kday+1)+'_night_sza.tif'
                    infile_obliq_sza = indir2 + symbol +'_2019%03d'%(kday+1)+'_night_sza_obliq.tif'
                    infile_nadir_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_night_vza.tif'
                    infile_obliq_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_night_vza_obliq.tif'

                    infile_cloud = indir2+ symbol +'_2019%03d'%(kday+1)+'_night_cloud.tif'



                    if ((os.path.exists(infile_nadir_lst) ==1) * (os.path.exists(infile_obliq_lst)==1)):
                        [lst1,ns,nl,temp,geog,proj] = read_image_gdal(infile_nadir_lst)
                        [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
                        [vza1,ns,nl,temp,geog,proj] = read_image_gdal(infile_nadir_vza)
                        [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)
                        vza_nadir[krun+kstep,:] = vza1*1.0
                        vza_obliq[krun+kstep,:] = vza2*1.0
                        lst_nadir[krun+kstep,:] = lst1*1.0
                        lst_obliq[krun+kstep,:] = lst2*1.0
            iffirst = 0
        else:
            for ksate in range(2):
                indir1 = indirs1[ksate]
                indir2 = indirs2[ksate]
                symbol = symbols[ksate]
                kstep = nb * ksate
                kday = kdoy +off
                infile_nadir_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_night_lst.tif'
                infile_obliq_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_night_lst_obliq.tif'
                infile_nadir_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_night_vza.tif'
                infile_obliq_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_night_vza_obliq.tif'
                infile_nadir_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_night_sza.tif'
                infile_obliq_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_night_sza_obliq.tif'
                infile_nadir_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_night_vza.tif'
                infile_obliq_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_night_vza_obliq.tif'
                if ((os.path.exists(infile_nadir_lst) == 1) * (os.path.exists(infile_obliq_lst) == 1)):
                    [lst1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_lst)
                    [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
                    [vza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vza)
                    [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)
                    vza_nadir[mark+kstep, :] = vza1 * 1.0
                    vza_obliq[mark+kstep, :] = vza2 * 1.0
                    lst_nadir[mark+kstep, :] = lst1 * 1.0
                    lst_obliq[mark+kstep, :] = lst2 * 1.0
                else:
                    vza_nadir[mark+kstep, :] = 0.0
                    vza_obliq[mark+kstep, :] = 0.0
                    lst_nadir[mark+kstep, :] = 0.0
                    lst_obliq[mark+kstep, :] = 0.0
            print(kdoy, mark)
            mark = mark + 1
            if mark >= nb: mark = 0

        for kl in range(nl):
            for ks in range(ns):
                lst1 = lst_nadir[:,kl,ks]
                lst2 = lst_obliq[:,kl,ks]
                temp = lst1*lst2
                ind = (temp > 1)
                if np.sum(ind) < 3: continue
                lst1 = lst1[ind]
                lst2 = lst2[ind]
                vza1 = vza_nadir[ind,kl,ks]
                vza2 = vza_obliq[ind,kl,ks]

                # y = lst1 - lst2
                # theta1 = 1- np.cos(vza1*rd)
                # theta2 = 1- np.cos(vza2*rd)
                # x = theta1*lst2 - theta2*lst1
                # a[kl, ks] = np.mean(y / x)
                # temp = kernels_VinRL_a(lst1,lst2,vza1,vza2)

                a[kl,ks] = fitting_Vin(lst1,lst2,vza1,vza2)


        outfile = outdira +'coeffa_2019%03d'%(kdoy+1)+'_%02d'%off+'.tif'
        write_image_gdal(a,ns,nl,1,geog,proj,outfile)


if ifcoeffb == 1:
    iffirst = 1
    mark = 0
    nb = off*2+1
    lst_nadir = np.zeros([nb * 2, nl, ns])
    lst_obliq = np.zeros([nb * 2, nl, ns])
    vza_nadir = np.zeros([nb * 2, nl, ns])
    sza_nadir = np.zeros([nb * 2, nl, ns])
    psi_obliq = np.zeros([nb * 2, nl, ns])
    psi_nadir = np.zeros([nb * 2, nl, ns])
    vza_obliq = np.zeros([nb * 2, nl, ns])
    sza_obliq = np.zeros([nb * 2, nl, ns])
    sza_obliq = np.zeros([nb * 2, nl, ns])
    rads = np.zeros([nb * 2, nl, ns])
    clouds = np.zeros([nb*2,nl,ns])
    coeffas = np.zeros([nb*2,nl,ns])
    for kdoy in range(off,num-off):
        print(kdoy)
        b = np.zeros([nl,ns])
        k = np.zeros([nl,ns])
        if iffirst == 1:
            for ksate in range(2):
                indir1 = indirs1[ksate]
                indir2 = indirs2[ksate]
                symbol = symbols[ksate]
                kstep = nb*ksate
                for krun in range(nb):
                    kday = kdoy-off+krun
                    infile_nadir_lst = indir1 + symbol +'_2019%03d'%(kday+1)+'_day_lst.tif'
                    infile_obliq_lst = indir1 + symbol +'_2019%03d'%(kday+1)+'_day_lst_obliq.tif'
                    infile_nadir_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_day_vza.tif'
                    infile_obliq_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_day_vza_obliq.tif'
                    infile_nadir_sza = indir2 + symbol +'_2019%03d'%(kday+1)+'_day_sza.tif'
                    infile_obliq_sza = indir2 + symbol +'_2019%03d'%(kday+1)+'_day_sza_obliq.tif'
                    infile_nadir_psi = indir2 + symbol +'_2019%03d'%(kday+1)+'_day_psi.tif'
                    infile_obliq_psi = indir2 + symbol +'_2019%03d'%(kday+1)+'_day_psi_obliq.tif'
                    infile_cloud = indir2+ symbol +'_2019%03d'%(kday+1)+'_day_cloud.tif'
                    infile_rad = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_rad.tif'

                    infile_ca = outdira + 'coeffa_2019%03d' % (kdoy+1) + '_%02d' % off + '.tif'

                    if ((os.path.exists(infile_nadir_lst) ==1) * (os.path.exists(infile_obliq_lst)==1)):
                        [lst1,ns,nl,temp,geog,proj] = read_image_gdal(infile_nadir_lst)
                        [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
                        [vza1,ns,nl,temp,geog,proj] = read_image_gdal(infile_nadir_vza)
                        [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)

                        [sza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_sza)
                        [sza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_sza)
                        [psi1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_psi)
                        [psi2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_psi)
                        [rad, ns, nl, temp, geog, proj] = read_image_gdal(infile_rad)
                        [coeffa, ns, nl, temp, geog, proj] = read_image_gdal(infile_ca)
                        # [cloud, ns, nl, temp, geog, proj] = read_image_gdal(infile_cloud)
                        vza_nadir[krun+kstep,:] = vza1*1.0
                        vza_obliq[krun+kstep,:] = vza2*1.0
                        sza_nadir[krun+kstep,:] = sza1*1.0
                        sza_obliq[krun+kstep,:] = sza2*1.0
                        psi_nadir[krun+kstep,:] = psi1*1.0
                        psi_obliq[krun+kstep,:] = psi2*1.0
                        lst_nadir[krun+kstep,:] = lst1*1.0
                        lst_obliq[krun+kstep,:] = lst2*1.0
                        rads[krun+kstep,:] = rad*1.0
                        coeffas[krun + kstep, :] = coeffa * 1.0
                        # clouds[krun + kstep, :] = cloud * 1.0
            iffirst = 0
        else:
            for ksate in range(2):
                indir1 = indirs1[ksate]
                indir2 = indirs2[ksate]
                symbol = symbols[ksate]
                kstep = nb * ksate
                kday = kdoy +off
                infile_nadir_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_day_lst.tif'
                infile_obliq_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_day_lst_obliq.tif'
                infile_nadir_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza.tif'
                infile_obliq_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza_obliq.tif'
                infile_nadir_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_sza.tif'
                infile_obliq_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_sza_obliq.tif'
                infile_nadir_psi = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_psi.tif'
                infile_obliq_psi = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_psi_obliq.tif'
                infile_cloud = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_cloud.tif'
                infile_rad = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_rad.tif'
                infile_ca = outdira + 'coeffa_2019%03d' % (kdoy + 1) + '_%02d' % off + '.tif'

                if ((os.path.exists(infile_nadir_lst) == 1) * (os.path.exists(infile_obliq_lst) == 1)):
                    [lst1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_lst)
                    [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
                    [vza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vza)
                    [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)
                    [sza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_sza)
                    [sza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_sza)
                    [psi1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_psi)
                    [psi2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_psi)
                    # [cloud, ns, nl, temp, geog, proj] = read_image_gdal(infile_cloud)
                    [rad, ns, nl, temp, geog, proj] = read_image_gdal(infile_rad)
                    [coeffa, ns, nl, temp, geog, proj] = read_image_gdal(infile_ca)


                    vza_nadir[mark+kstep, :] = vza1 * 1.0
                    vza_obliq[mark+kstep, :] = vza2 * 1.0
                    lst_nadir[mark+kstep, :] = lst1 * 1.0
                    lst_obliq[mark+kstep, :] = lst2 * 1.0
                    sza_nadir[mark + kstep, :] = sza1 * 1.0
                    sza_obliq[mark + kstep, :] = sza2 * 1.0
                    psi_nadir[mark + kstep, :] = psi1 * 1.0
                    psi_obliq[mark + kstep, :] = psi2 * 1.0
                    # clouds[mark + kstep, :] = cloud * 1.0
                    rads[mark + kstep, :] = rad * 1.0
                    coeffas[mark + kstep, :] = coeffa * 1.0
                else:
                    vza_nadir[mark+kstep, :] = 0.0
                    vza_obliq[mark+kstep, :] = 0.0
                    lst_nadir[mark+kstep, :] = 0.0
                    lst_obliq[mark+kstep, :] = 0.0
                    sza_nadir[mark + kstep, :] = 0.0
                    sza_obliq[mark + kstep, :] = 0.0
                    psi_nadir[mark + kstep, :] = 0.0
                    psi_obliq[mark + kstep, :] = 0.0
                    rads[mark + kstep, :] = 0.0
                    coeffas[mark + kstep, :] = 0.0
            print(kdoy, mark)
            mark = mark + 1
            if mark >= nb: mark = 0

        for kl in range(nl):
            for ks in range(ns):
                lst1 = lst_nadir[:,kl,ks]
                lst2 = lst_obliq[:,kl,ks]
                temp = 1.0*lst1*lst2
                ind = (temp > th)
                if np.sum(ind) < 3: continue
                lst1 = lst1[ind]
                lst2 = lst2[ind]
                y = lst1 - lst2
                vza1 = vza_nadir[ind,kl,ks]
                vza2 = vza_obliq[ind,kl,ks]
                sza1 = sza_nadir[ind,kl,ks]
                sza2 = sza_obliq[ind,kl,ks]
                psi1 = psi_nadir[ind,kl,ks]
                psi2 = psi_obliq[ind,kl,ks]
                rad = rads[ind,kl,ks]
                coeffa = coeffas[ind, kl, ks]
                temp = fitting_RLE(lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,coeffa,rad)

                b[kl,ks] = temp[0]
                k[kl,ks] = temp[1]

        b[b> 1000] = 0
        b[b<-1000] = 0
        k[k<-1000] = 0
        k[k> 1000] = 0
        outfile = outdirbk +'coeffb_2019%03d'%(kdoy+1)+'_%02d'%off+'.tif'
        write_image_gdal(b,ns,nl,1,geog,proj,outfile)
        outfile = outdirbk +'coeffk_2019%03d'%(kdoy+1)+'_%02d'%off+'.tif'
        write_image_gdal(k,ns,nl,1,geog,proj,outfile)

