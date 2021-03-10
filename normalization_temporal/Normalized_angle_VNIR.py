from util.myfun import *
from fitting_kernels import *
# from myfun_plot import *
from multiprocessing import Pool, cpu_count
import multiprocessing
import time

def coeffAndRmse(kdoy, off):
# def coeffAndRmse(kdoy, nl, ns, iffirst, mark):

    indir_s3a = r'F:/S3A_tif/'
    indir_s3b = r'H:/S3B_tif/'
    symbols = ['S3A', 'S3B']
    outdira = r'F:/S3A_coeff/'
    indir_vnir = [indir_s3a, indir_s3b]
    indir_tif = [indir_s3a,indir_s3b]



    infile_type = 'F:/LST_baseandTcw/type.tif'
    [type, ns, nl, nb, geog, proj] = read_image_gdal(infile_type)
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
    cloud = np.zeros([nb*2,nl,ns])
    print(kdoy)
    coeff0 = np.zeros([nl, ns])
    coeff_vol = np.zeros([nl, ns])
    coeff_hs = np.zeros([nl, ns])
    rmse0 = np.zeros([nl, ns])
    if iffirst == 1:
        for ksate in range(2):
            indir1 = indir_vnir[ksate]
            indir2 = indir_tif[ksate]
            symbol = symbols[ksate]
            kstep = nb * ksate
            for krun in range(nb):
                kday = kdoy - off + krun

                infile_nadir_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_day_red.tif'
                infile_obliq_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_day_red_obliq.tif'
                infile_nadir_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza.tif'
                infile_obliq_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza_obliq.tif'
                infile_nadir_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_sza.tif'
                infile_obliq_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_sza_obliq.tif'
                infile_nadir_vaa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vaa.tif'
                infile_obliq_vaa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vaa_obliq.tif'
                infile_nadir_saa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_saa.tif'
                infile_obliq_saa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_saa_obliq.tif'
                infile_cloud = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_cloud.tif'

                if ((os.path.exists(infile_nadir_lst) == 1) * (os.path.exists(infile_obliq_lst) == 1)):
                    [lst1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_lst)
                    [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
                    [vza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vza)
                    [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)

                    [sza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_sza)
                    [sza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_sza)
                    [vaa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vaa)
                    [vaa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vaa)
                    [saa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_saa)
                    [saa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_saa)
                    [cloud0, ns, nl, temp, geog, proj] = read_image_gdal(infile_cloud)

                    vza_nadir[krun + kstep, :] = vza1 * 1.0
                    vza_obliq[krun + kstep, :] = vza2 * 1.0
                    sza_nadir[krun + kstep, :] = sza1 * 1.0
                    sza_obliq[krun + kstep, :] = sza2 * 1.0
                    psi_nadir[krun + kstep, :] = np.abs(saa1 - vaa1) * 1.0
                    psi_obliq[krun + kstep, :] = np.abs(saa2 - vaa2) * 1.0
                    lst_nadir[krun + kstep, :] = lst1 * 1.0
                    lst_obliq[krun + kstep, :] = lst2 * 1.0
                    cloud[krun + kstep, :] = cloud0 * 1.0

        iffirst = 0
    else:
        coeff0 = np.zeros([nl, ns])
        coeff_vol = np.zeros([nl, ns])
        coeff_hs = np.zeros([nl, ns])
        for ksate in range(2):
            indir1 = indir_vnir[ksate]
            indir2 = indir_tif[ksate]
            symbol = symbols[ksate]
            kstep = nb * ksate
            kday = kdoy + off

            infile_nadir_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_day_red.tif'
            infile_obliq_lst = indir1 + symbol + '_2019%03d' % (kday + 1) + '_day_red_obliq.tif'
            # infile_nadir_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza.tif'
            # infile_obliq_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza_obliq.tif'
            infile_nadir_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_sza.tif'
            infile_obliq_sza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_sza_obliq.tif'
            infile_nadir_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza.tif'
            infile_obliq_vza = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vza_obliq.tif'
            infile_nadir_vaa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vaa.tif'
            infile_obliq_vaa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_vaa_obliq.tif'
            infile_nadir_saa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_saa.tif'
            infile_obliq_saa = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_saa_obliq.tif'
            infile_cloud = indir2 + symbol + '_2019%03d' % (kday + 1) + '_day_cloud.tif'

            if ((os.path.exists(infile_nadir_lst) == 1) * (os.path.exists(infile_obliq_lst) == 1)):
                [lst1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_lst)
                [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
                [vza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vza)
                [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)

                [sza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_sza)
                [sza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_sza)
                [vaa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vaa)
                [vaa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vaa)
                [saa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_saa)
                [saa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_saa)
                [cloud0, ns, nl, temp, geog, proj] = read_image_gdal(infile_cloud)

                vza_nadir[mark + kstep, :] = vza1 * 1.0
                vza_obliq[mark + kstep, :] = vza2 * 1.0
                lst_nadir[mark + kstep, :] = lst1 * 1.0
                lst_obliq[mark + kstep, :] = lst2 * 1.0
                sza_nadir[mark + kstep, :] = sza1 * 1.0
                sza_obliq[mark + kstep, :] = sza2 * 1.0
                psi_nadir[mark + kstep, :] = np.abs(vaa1 - saa1) * 1.0
                psi_obliq[mark + kstep, :] = np.abs(vaa2 - saa2) * 1.0
                cloud[mark + kstep, :] = cloud0 * 1.0
            else:
                vza_nadir[mark + kstep, :] = 0.0
                vza_obliq[mark + kstep, :] = 0.0
                lst_nadir[mark + kstep, :] = 0.0
                lst_obliq[mark + kstep, :] = 0.0
                sza_nadir[mark + kstep, :] = 0.0
                sza_obliq[mark + kstep, :] = 0.0
                psi_nadir[mark + kstep, :] = 0.0
                psi_obliq[mark + kstep, :] = 0.0
                cloud[mark + kstep, :] = 0.0

        print(kdoy, mark)
        mark = mark + 1
        if mark >= nb: mark = 0

    for kl in range(nl):
        for ks in range(ns):
            if type[kl, ks] >= 17: continue
            lst1 = lst_nadir[:, kl, ks]
            lst2 = lst_obliq[:, kl, ks]
            ind1 = (lst1 > 0.001) * (cloud[:, kl, ks] == 0) * (lst1 < 1)
            ind2 = (lst2 > 0.001) * (cloud[:, kl, ks] == 0) * (lst2 < 1)
            if np.sum(ind1) + np.sum(ind2) < 5: continue
            lst1 = lst1[ind1]
            lst2 = lst2[ind2]
            vza1 = vza_nadir[ind1, kl, ks]
            vza2 = vza_obliq[ind2, kl, ks]
            sza1 = sza_nadir[ind1, kl, ks]
            sza2 = sza_obliq[ind2, kl, ks]
            psi1 = psi_nadir[ind1, kl, ks]
            psi2 = psi_obliq[ind2, kl, ks]
            lst = np.hstack([lst1, lst2])
            vza = np.hstack([vza1, vza2])
            sza = np.hstack([sza1, sza2])
            psi = np.hstack([psi1, psi2])

            coeff, rmse = fitting_LSRROSS(sza, vza, psi, lst)

            coeff0[kl, ks] = coeff[0]
            coeff_vol[kl, ks] = coeff[1]
            coeff_hs[kl, ks] = coeff[2]
            rmse0[kl, ks] = rmse

    outfile = outdira + '2019%03d' % (kdoy + 1) + 'red0_%02d' % off + '.tif'
    write_image_gdal(coeff0, ns, nl, 1, geog, proj, outfile)
    outfile = outdira + '2019%03d' % (kdoy + 1) + 'redvol_%02d' % off + '.tif'
    write_image_gdal(coeff_vol, ns, nl, 1, geog, proj, outfile)
    outfile = outdira + '2019%03d' % (kdoy + 1) + 'redhs_%02d' % off + '.tif'
    write_image_gdal(coeff_hs, ns, nl, 1, geog, proj, outfile)
    outfile = outdira + '2019%03d' % (kdoy + 1) + 'redrmse_%02d' % off + '.tif'
    write_image_gdal(rmse0, ns, nl, 1, geog, proj, outfile)

def doMultiple(param):
    return coeffAndRmse(int(param[0]), int(param[1]))


# os.environ['PROJ_LIB'] = r'C:\Users\zzz\anaconda3\envs\patrol\Library\share\proj'

if __name__ == '__main__':
    multiprocessing.freeze_support()
    # indir_tif = r'G:/S3B_tif/'
    # indir_vnir = r'G:/S3B_tif/'
    #
    # outdira = r'G:/S3B_coeff/'
    # symbol = 'S3B'
    #
    # infile_type = 'G:/base/type.tif'
    # [type, ns, nl, nb, geog, proj] = read_image_gdal(infile_type)

    # off = 4
    # num = 365


    off = 4
    num = 31
    rd = np.pi/180.0

    ifcoeffa = 0
    ifcoeffb = 1
    th = 0.001

    ####multiple
    start = time.time()
    param = []
    for i in range(8,20):
        param.append([i, off])
    print(param[0])
    # doMultiple(param[1])
    pool = Pool(min(6, cpu_count()))
    for p in param:
        pool.apply_async(doMultiple, (p,))
    pool.close()
    pool.join()
    end = time.time()
    print('time cost: ', end - start, 's')

    ####single
    # start = time.time()
    # param = [5, off, indir_vnir, indir_tif, symbol, outdira]
    # doMultiple(param)
    # end = time.time()
    # print('time cost: ',end - start,'s')