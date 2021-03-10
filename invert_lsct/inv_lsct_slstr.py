import numpy as np
from scipy import interpolate
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
import time
from myfun import *
from multiple_lsct import *

if __name__ == '__main__':
    indir_lst = 'F:/S3A_LST/'
    indir_tif = 'F:/S3A_TIF/'
    indir_ems = 'H:/EMIS/'
    outdir_lsct = 'F:/S3A_LSCT/'

    subdir = np.asarray(['2019_01','2019_02','2019_03','2019_04','2019_05','2019_06',
              '2019_07','2019_08','2019_09','2019_10','2019_11','2019_12'])

    num_month = len(subdir)

    infile_soil_emis = indir_ems+'emis1.tif'
    infile_type = indir_ems+'type.tif'
    infile_veg_emis = indir_ems +'ICBP_Emi.txt'

    [emis_s,emis_v,geog,proj] = file2Array_static(infile_soil_emis,infile_veg_emis,infile_type)


    for km in range(0,3):
        indir_lst_sub = indir_lst+subdir[km]+'/'
        indir_tif_sub = indir_tif+subdir[km]+'/'
        outdir_lsct_sub = outdir_lsct+subdir[km]+'/'


        files = search_file_rej(indir_lst_sub,'tif','obliq')
        num_files = len(files)
        for kf in range(num_files):
            infile = files[kf][:-8]
            infile_lst = indir_lst_sub+infile+'_lst.tif'
            infile_lst_obliq = indir_lst_sub+infile+'_lst_obliq.tif'
            infile_red = indir_tif_sub+infile+'_red.tif'
            infile_red_obliq = indir_tif_sub+infile+'_red_obliq.tif'
            infile_nir = indir_tif_sub+infile+'_nir.tif'
            infile_nir_obliq = indir_tif_sub+infile+'_nir_obliq.tif'
            outfile_lsct = outdir_lsct_sub+infile+'_lsct.tif'
            print(outfile_lsct)

            c1 = os.path.exists(infile_lst)
            c2 = os.path.exists(infile_lst_obliq)
            c3 = os.path.exists(infile_nir)
            c4 = os.path.exists(infile_nir_obliq)
            c5 = os.path.exists(infile_red)
            c6 = os.path.exists(infile_red_obliq)
            ccc = c1 * c2 * c3 * c4 * c5 * c6
            if ccc == 0: continue

            [lst_n,ndvi_n,geog,prop] = file2Array_dynamic(infile_lst,infile_red,infile_nir)
            [lst_o,ndvi_o,geog,prop] = file2Array_dynamic(infile_lst_obliq,infile_red_obliq,infile_nir_obliq)
            # lsct = inv_mpl_agl_multiple(lst_n,lst_o,emis_s,emis_v,ndvi_n,ndvi_o,wl=10.8)
            lsct = inv_mpl_agl(lst_n, lst_o, emis_s, emis_v, ndvi_n, ndvi_o, wl=10.8)
            [ns,nl,nb] = np.shape(lsct)
            write_image_gdal(lsct,ns,nl,nb,geog,proj,outfile_lsct)










