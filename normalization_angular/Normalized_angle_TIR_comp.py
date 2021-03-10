from util.myfun import *
from fitting_kernels import *
# from myfun_plot import *
# os.environ['PROJ_LIB'] = r'C:\Users\zzz\anaconda3\envs\patrol\Library\share\proj'
import time
from multiprocessing import Pool, cpu_count
import multiprocessing






def coeff_invert(kdoy,off,l1,l2,s1,s2,outdir):

    indir_s3a = r'F:/S3A_tif/'
    indir_s3b = r'H:/S3B_tif/'
    indirs2 = [indir_s3a, indir_s3b]
    indir_s3a = r'F:/S3A_LST/'
    indir_s3b = r'H:/S3B_LST/'
    indirs1 = [indir_s3a, indir_s3b]

    infile_type = r'D:\data\lst\base/type.tif'
    [type, ns, nl, nb, geog, proj] = read_image_gdal(infile_type)
    symbols = ['S3A', 'S3B']
    dayNight = 'day'
    num = 365
    rd = np.pi / 180.0

    nll = l2-l1
    nss = s2-s1


    ifcoeffa = 1
    ifcoeffb = 1
    th = 220

    if ifcoeffb == 1:
        iffirst = 1
        mark = 0
        nb = off*2+1

        lst_nadir = np.zeros([nb * 2, nll, nss])
        lst_obliq = np.zeros([nb * 2, nll, nss])
        vza_nadir = np.zeros([nb * 2, nll, nss])
        sza_nadir = np.zeros([nb * 2, nll, nss])
        psi_obliq = np.zeros([nb * 2, nll, nss])
        psi_nadir = np.zeros([nb * 2, nll, nss])
        vza_obliq = np.zeros([nb * 2, nll, nss])
        sza_obliq = np.zeros([nb * 2, nll, nss])
        clouds = np.zeros([nb*2,nll,nss])


        rads = np.zeros([nb * 2, nll, nss])
        coeffas = np.zeros([nb*2,nll,nss])

        print(kdoy)
        a = np.zeros([nll,nss])
        b = np.zeros([nll,nss])
        k = np.zeros([nll,nss])
        if iffirst == 1:
            for ksate in range(2):
                indir1 = indirs1[ksate]
                indir2 = indirs2[ksate]
                symbol = symbols[ksate]
                kstep = nb*ksate
                for krun in range(nb):
                    kday = kdoy-off+krun
                    infile_nadir_lst = indir1 + symbol +'_2019%03d'%(kday+1)+'_'+dayNight+'_lst.tif'
                    infile_obliq_lst = indir1 + symbol +'_2019%03d'%(kday+1)+'_'+dayNight+'_lst_obliq.tif'
                    infile_nadir_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_'+dayNight+'_vza.tif'
                    infile_obliq_vza = indir2 + symbol +'_2019%03d'%(kday+1)+'_'+dayNight+'_vza_obliq.tif'
                    infile_nadir_sza = indir2 + symbol +'_2019%03d'%(kday+1)+'_'+dayNight+'_sza.tif'
                    infile_obliq_sza = indir2 + symbol +'_2019%03d'%(kday+1)+'_'+dayNight+'_sza_obliq.tif'
                    infile_nadir_vaa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_'+dayNight+'_vaa.tif'
                    infile_obliq_vaa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_'+dayNight+'_vaa_obliq.tif'
                    infile_nadir_saa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_'+dayNight+'_saa.tif'
                    infile_obliq_saa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_'+dayNight+'_saa_obliq.tif'
                    infile_cloud = indir2+ symbol +'_2019%03d'%(kday+1)+'_'+dayNight+'_cloud.tif'
                    infile_rad = indir2 + symbol + '_2019%03d' % (kday + 1) + '_'+dayNight+'_rad.tif'

                    infile_ca = outdir + '2019%03d' % (kdoy+1)  + '_nighta.tif'
                    print(kdoy, infile_nadir_lst)
                    if ((os.path.exists(infile_nadir_lst) ==1) * (os.path.exists(infile_nadir_vaa)==1)):
                        [lst1,ns,nl,temp,geog,proj] = read_image_gdal(infile_nadir_lst)
                        [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
                        [vza1,ns,nl,temp,geog,proj] = read_image_gdal(infile_nadir_vza)
                        [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)

                        [sza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_sza)
                        [sza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_sza)
                        [vaa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vaa)
                        [vaa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vaa)
                        [saa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_saa)
                        [saa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_saa)
                        [rad, ns, nl, temp, geog, proj] = read_image_gdal(infile_rad)

                        [cloud, ns, nl, temp, geog, proj] = read_image_gdal(infile_cloud)


                        vza_nadir[krun+kstep,:] = vza1[l1:l2,s1:s2]*1.0
                        vza_obliq[krun+kstep,:] = vza2[l1:l2,s1:s2]*1.0
                        sza_nadir[krun+kstep,:] = sza1[l1:l2,s1:s2]*1.0
                        sza_obliq[krun+kstep,:] = sza2[l1:l2,s1:s2]*1.0
                        psi_nadir[krun + kstep, :] = np.abs(saa1[l1:l2,s1:s2] - vaa1[l1:l2,s1:s2]) * 1.0
                        psi_obliq[krun + kstep, :] = np.abs(saa2[l1:l2,s1:s2] - vaa2[l1:l2,s1:s2]) * 1.0
                        lst_nadir[krun + kstep, :] = eliminate_edge(lst1[l1:l2, s1:s2]) * 1.0
                        lst_obliq[krun + kstep, :] = eliminate_edge(lst2[l1:l2, s1:s2]) * 1.0
                        rads[krun+kstep,:] = rad[l1:l2,s1:s2]*1.0

                        clouds[krun + kstep, :] = cloud[l1:l2,s1:s2] * 1.0
            iffirst = 0

        for kl in range(nll):
            for ks in range(nss):
                lst1 = lst_nadir[:,kl,ks]
                lst2 = lst_obliq[:,kl,ks]
                cloud0 = clouds[:,kl,ks]
                ind = (lst1>0)*(lst2>0)*(cloud0<=0.05)

              #  print(np.sum(ind))
                if np.sum(ind) < 2: continue

                lst1 = lst1[ind]
                lst2 = lst2[ind]
                y = lst1 - lst2
                vza1 = vza_nadir[ind,kl,ks]
                vza2 = vza_obliq[ind,kl,ks]
                sza1 = sza_nadir[ind,kl,ks]
                sza2 = sza_obliq[ind,kl,ks]
                psi1 = psi_nadir[ind,kl,ks]
                psi2 = psi_obliq[ind,kl,ks]
                rad = rads[ind,kl,ks]/3600.0
                temp = fitting_VinRLE_rad_full(lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,rad)
                a[kl,ks] = temp[0]
                b[kl,ks] = temp[1]
                k[kl,ks] = temp[2]
        a[a> 1000] = 0
        a[a<-1000] = 0
        b[b> 1000] = 0
        b[b<-1000] = 0
        k[k<-1000] = 0
        k[k> 1000] = 0
        outfile = outdir +'2019%03d'%(kdoy+1)+'_daya.tif'
        write_image_gdal(a,nss,nll,1,geog,proj,outfile)
        outfile = outdir +'2019%03d'%(kdoy+1)+'_dayb.tif'
        write_image_gdal(b,nss,nll,1,geog,proj,outfile)
        outfile = outdir +'2019%03d'%(kdoy+1)+'_dayk.tif'
        write_image_gdal(k,nss,nll,1,geog,proj,outfile)

def doMultiple(param):
    return coeff_invert(np.int(param[0]), np.int(param[1]),np.int(param[2]),np.int(param[3]),np.int(param[4]),np.int(param[5]),param[6])

if __name__ == '__main__':
    multiprocessing.freeze_support()

    off = 15
    l1 = 550
    l2 = 750
    s1 = 900
    s2 = 1100
    outdir = r'D:\data\an/Ecoeff_comp/'

    # param = []
    # for i in range(15,16,30):
    #     param.append([i, off,l1,l2,s1,s2,outdir])
    # doMultiple(param[0])

    ####multiple
    start = time.time()
    param = []
    for i in range(15,360,30):
        param.append([i, off,l1,l2,s1,s2,outdir])
    pool = Pool(min(48, cpu_count()))
    for p in param:
        pool.apply_async(doMultiple, (p,))
    pool.close()
    pool.join()
    end = time.time()
    print('time cost: ', end - start, 's')
