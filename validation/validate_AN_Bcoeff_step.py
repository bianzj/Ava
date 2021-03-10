from util.myfun import *
from AngularNormalization.fitting_kernels import *

indir_coeff = 'J:/AN/Bcoeff_comp/'
indir_lst = 'J:/s3a_lst/'
indir_tif = 'J:/s3a_tif/'


files = search_file(indir_lst,['day','lst'])

fileNum = np.size(files)
symbol = 'S3A'
print('filenum:',fileNum)

l1 = 550
l2 = 750
s1 = 900
s2 = 1100


# for kmonth in range(12):
#     data1_old = np.asarray([])
#     data2_old = np.asarray([])
#     data1_new = np.asarray([])
#     data2_new = np.asarray([])
#     infile_coeffa = indir_coeff + '2019%03d' % (kmonth * 30 + 16) + '_nighta.tif'
#     [coeffa, ns, nl, nb, geog, proj] = read_image_gdal(infile_coeffa)
#     for kday in range(31):
#         kdoy = date2DOY(2019,kmonth + 1,kday+1)-1
#         infile_nadir_lst = indir_lst + symbol + '_2019%03d' % (kdoy + 1) + '_night_lst.tif'
#         infile_obliq_lst = indir_lst + symbol + '_2019%03d' % (kdoy + 1) + '_night_lst_obliq.tif'
#         infile_nadir_vza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_night_vza.tif'
#         infile_obliq_vza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_night_vza_obliq.tif'
#         infile_nadir_cloud = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_night_cloud.tif'
#         infile_obliq_cloud = indir_tif +  symbol + '_2019%03d' % (kdoy + 1) + '_night_cloud_obliq.tif'
#         if ((os.path.exists(infile_nadir_lst) == 1) * (os.path.exists(infile_obliq_lst) == 1)):
#             [lst1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_lst)
#             [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
#             [vza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vza)
#             [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)
#             [cloud1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_cloud)
#             [cloud2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_cloud)
#             vza1 = vza1[l1:l2, s1:s2] * 1.0
#             vza2 = vza2[l1:l2, s1:s2] * 1.0
#             lst1 = lst1[l1:l2, s1:s2] * 1.0
#             lst2 = lst2[l1:l2, s1:s2] * 1.0
#             cloud1 = cloud1[l1:l2, s1:s2] * 1.0
#             cloud2 = cloud2[l1:l2, s1:s2] * 1.0
#             nll = l2 - l1
#             nss = s2 - s1
#             LSF = model_LSF(vza1,coeffa)
#             ind = LSF!=0
#             lst1_new = np.zeros([nll,nss])
#             lst1_new[ind] = lst1[ind]+(LSF[ind])
#             LSF = model_LSF(vza2,coeffa)
#             lst2_new = np.zeros([nll,nss])
#             ind = LSF != 0
#             lst2_new[ind] = lst2[ind]+(LSF[ind])
#             ind = (cloud1 == 0.0) *(cloud2 ==0 )
#             data1_old = np.hstack([data1_old,lst1[ind]])
#             data2_old = np.hstack([data2_old,lst2[ind]])
#             data1_new = np.hstack([data1_new,lst1_new[ind]])
#             data2_new = np.hstack([data2_new,lst2_new[ind]])
#     plt_scatter(data1_old,data2_old,'%02d'%(kmonth+1)+'_old',15)
#     plt_scatter(data1_new,data2_new,'%02d'%(kmonth+1)+'_new',15)

for kmonth in range(12):
    data1_old = np.asarray([])
    data2_old = np.asarray([])
    data1_new = np.asarray([])
    data2_new = np.asarray([])
    infile_coeffa = indir_coeff + '2019%03d' % (kmonth * 30 + 16) + '_nighta.tif'
    [coeffa, ns, nl, nb, geog, proj] = read_image_gdal(infile_coeffa)
    for kday in range(31):
        kdoy = date2DOY(2019,kmonth + 1,kday+1)-1
        infile_nadir_lst = indir_lst + symbol + '_2019%03d' % (kdoy + 1) + '_night_lst.tif'
        infile_obliq_lst = indir_lst + symbol + '_2019%03d' % (kdoy + 1) + '_night_lst_obliq.tif'
        infile_nadir_vza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_night_vza.tif'
        infile_obliq_vza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_night_vza_obliq.tif'
        infile_nadir_cloud = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_night_cloud.tif'
        infile_obliq_cloud = indir_tif +  symbol + '_2019%03d' % (kdoy + 1) + '_night_cloud_obliq.tif'
        if ((os.path.exists(infile_nadir_lst) == 1) * (os.path.exists(infile_obliq_lst) == 1)):
            [lst1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_lst)
            [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
            [vza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vza)
            [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)
            [cloud1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_cloud)
            [cloud2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_cloud)
            vza1 = vza1[l1:l2, s1:s2] * 1.0
            vza2 = vza2[l1:l2, s1:s2] * 1.0
            lst1 = lst1[l1:l2, s1:s2] * 1.0
            lst2 = lst2[l1:l2, s1:s2] * 1.0
            cloud1 = cloud1[l1:l2, s1:s2] * 1.0
            cloud2 = cloud2[l1:l2, s1:s2] * 1.0
            nll = l2 - l1
            nss = s2 - s1
            LSF = model_LSF(vza1,coeffa)
            ind = LSF!=0
            lst1_new = np.zeros([nll,nss])
            lst1_new[ind] = lst1[ind]+(LSF[ind])
            LSF = model_LSF(vza2,coeffa)
            lst2_new = np.zeros([nll,nss])
            ind = LSF != 0
            lst2_new[ind] = lst2[ind]+(LSF[ind])
            ind = (cloud1 == 0.0) *(cloud2 ==0 )
            data1_old = np.hstack([data1_old,lst1[ind]])
            data2_old = np.hstack([data2_old,lst2[ind]])
            data1_new = np.hstack([data1_new,lst1_new[ind]])
            data2_new = np.hstack([data2_new,lst2_new[ind]])
    plt_scatter(data1_old,data2_old,'%02d'%(kmonth+1)+'_old',15)
    plt_scatter(data1_new,data2_new,'%02d'%(kmonth+1)+'_new',15)




