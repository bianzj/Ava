from util.myfun import *
from AngularNormalization.fitting_kernels import *

indir_coeff = 'J:/AN/Ecoeff_step/'
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

nll = l2 - l1
nss = s2 - s1
nlld2 = np.int(nll / 2)
nssd2 = np.int(nss / 2)

for kmonth in range(12):
    data1_old = np.asarray([])
    data2_old = np.asarray([])
    data1_new = np.asarray([])
    data2_new = np.asarray([])
    infile_coeffa = indir_coeff + '2019%03d' % (kmonth * 30 + 16) + '_nighta.tif'
    [coeffa, ns, nl, nb, geog, proj] = read_image_gdal(infile_coeffa)
    coeffa = resize_data_ls(coeffa,nlld2,nssd2)
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

            vza1 = resize_data_ls(vza1,nlld2,nssd2)
            vza2 = resize_data_ls(vza2,nlld2,nssd2)
            lst1 = resize_data_ls(lst1,nlld2,nssd2)
            lst2 = resize_data_ls(lst2,nlld2,nssd2)
            cloud1 = resize_data_ls(cloud1,nlld2,nssd2)
            cloud2 = resize_data_ls(cloud2,nlld2,nssd2)

            Vin = model_Vin(vza1,coeffa)
            ind = Vin!=0
            lst1_new = np.zeros([nlld2,nssd2])
            lst1_new[ind] = lst1[ind]/(1+Vin[ind])
            Vin = model_LSF(vza2,coeffa)
            lst2_new = np.zeros([nlld2,nssd2])
            ind = Vin != 0
            lst2_new[ind] = lst2[ind]/(1+Vin[ind])
            ind = (cloud1 == 0.0) *(cloud2 ==0 )
            data1_old = np.hstack([data1_old,lst1[ind]])
            data2_old = np.hstack([data2_old,lst2[ind]])
            data1_new = np.hstack([data1_new,lst1_new[ind]])
            data2_new = np.hstack([data2_new,lst2_new[ind]])
    plt_scatter(data1_old,data2_old,'%02d'%(kmonth+1)+'_old',30)
    plt_scatter(data1_new,data2_new,'%02d'%(kmonth+1)+'_new',30)

# for kmonth in range(12):
#     data1_old = np.asarray([])
#     data2_old = np.asarray([])
#     data1_new = np.asarray([])
#     data2_new = np.asarray([])
#     infile_coeffa = indir_coeff + '2019%03d' % (kmonth * 30 + 16) + '_daya.tif'
#     [coeffa, ns, nl, nb, geog, proj] = read_image_gdal(infile_coeffa)
#     infile_coeffb = indir_coeff + '2019%03d' % (kmonth * 30 + 16) + '_dayb.tif'
#     [coeffb, ns, nl, nb, geog, proj] = read_image_gdal(infile_coeffa)
#     infile_coeffk = indir_coeff + '2019%03d' % (kmonth * 30 + 16) + '_dayk.tif'
#     [coeffk, ns, nl, nb, geog, proj] = read_image_gdal(infile_coeffa)
#
#     for kday in range(31):
#
#         kdoy = date2DOY(2019,kmonth + 1,kday+1)-1
#
#
#         infile_nadir_cloud = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_cloud.tif'
#         infile_obliq_cloud = indir_tif +  symbol + '_2019%03d' % (kdoy + 1) + '_day_cloud_obliq.tif'
#         infile_nadir_lst = indir_lst + symbol + '_2019%03d' % (kdoy + 1) + '_day_lst.tif'
#         infile_obliq_lst = indir_lst + symbol + '_2019%03d' % (kdoy + 1) + '_day_lst_obliq.tif'
#         infile_nadir_vza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_vza.tif'
#         infile_obliq_vza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_vza_obliq.tif'
#         infile_nadir_sza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_sza.tif'
#         infile_obliq_sza = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_sza_obliq.tif'
#         infile_nadir_vaa = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_vaa.tif'
#         infile_obliq_vaa = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_vaa_obliq.tif'
#         infile_nadir_saa = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_saa.tif'
#         infile_obliq_saa = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_saa_obliq.tif'
#         infile_rad = indir_tif + symbol + '_2019%03d' % (kdoy + 1) + '_day_rad.tif'
#
#
#
#         if ((os.path.exists(infile_nadir_lst) == 1) * (os.path.exists(infile_nadir_vaa) == 1)):
#             [lst1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_lst)
#             [lst2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_lst)
#             [vza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vza)
#             [vza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vza)
#
#             [sza1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_sza)
#             [sza2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_sza)
#             [vaa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_vaa)
#             [vaa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_vaa)
#             [saa1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_saa)
#             [saa2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_saa)
#             [rad, ns, nl, temp, geog, proj] = read_image_gdal(infile_rad)
#             [cloud1, ns, nl, temp, geog, proj] = read_image_gdal(infile_nadir_cloud)
#             [cloud2, ns, nl, temp, geog, proj] = read_image_gdal(infile_obliq_cloud)
#
#
#             psi1 = np.abs(saa1[l1:l2, s1:s2] - vaa1[l1:l2, s1:s2]) * 1.0
#             psi2= np.abs(saa2[l1:l2, s1:s2] - vaa2[l1:l2, s1:s2]) * 1.0
#             rad = rad[l1:l2, s1:s2] * 1.0/3600
#             vza1 = vza1[l1:l2, s1:s2] * 1.0
#             vza2 = vza2[l1:l2, s1:s2] * 1.0
#             lst1 =  eliminate_edge(lst1[l1:l2, s1:s2]) * 1.0
#             lst2 =  eliminate_edge(lst2[l1:l2, s1:s2]) * 1.0
#             cloud1 = cloud1[l1:l2, s1:s2] * 1.0
#             cloud2 = cloud2[l1:l2, s1:s2] * 1.0
#             sza1 = sza1[l1:l2,s1:s2]*1.0
#             sza2 = sza2[l1:l2,s1:s2]*1.0
#
#
#             Vin,RLE = model_VinRLE(sza1,vza1,psi1,rad,coeffa,coeffb,coeffk)
#             lst1_new = (lst1-RLE)/(1+Vin)
#             Vin,RLE =  model_VinRLE(sza2,vza2,psi2,rad,coeffa,coeffb,coeffk)
#             lst2_new = (lst2-RLE)/(1+Vin)
#
#             ind = (cloud1 == 0.0) *(cloud2 ==0 )
#             data1_old = np.hstack([data1_old,lst1[ind]])
#             data2_old = np.hstack([data2_old,lst2[ind]])
#             data1_new = np.hstack([data1_new,lst1_new[ind]])
#             data2_new = np.hstack([data2_new,lst2_new[ind]])
#
#
#     plt_scatter(data1_old,data2_old,'%02d'%(kmonth+1)+'_old',15)
#     plt_scatter(data1_new,data2_new,'%02d'%(kmonth+1)+'_new',15)






