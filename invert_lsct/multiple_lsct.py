from myfun import *
from multiprocessing import Pool,cpu_count

def file2Array_dynamic(infile_lst,infile_red,infile_nir):

    [lst, temp1, temp2, temp3, geog, proj] = read_image_gdal(infile_lst)

    [red, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile_red)
    [nir, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile_nir)
    ndvi = red*1.0
    ind = (red >0) *(nir >0)
    ndvi[ind] = (nir[ind]-red[ind])/(nir[ind]+red[ind])

    return lst,ndvi,geog,proj


def file2Array_static(infile_emis_s,infile_emis_v,infile_classif):

    number_colume_emis = 2
    time_emissivity = 1000.0
    # get soil emissivity
    [emis_s, ns, nl, nb, geog, proj] = read_image_gdal(infile_emis_s)
    emis_s = emis_s / time_emissivity

    [classif, temp1, temp2, temp3, temp4, temp5] = read_image_gdal(infile_classif)
    classif = np.asarray(classif, np.int)

    emis_v = read_txt_array(infile_emis_v, 0, number_colume_emis)
    emis_v = emis_v[classif, 0] / time_emissivity
    return emis_s,emis_v,geog,proj



def inv_mpl_agl(lst_n,lst_o,emis_s,emis_v,ndvi_n,ndvi_o,wl=10.5,ndvi_min = 0.05,ndvi_max = 0.95,off=2,off_=0):

    radiance_max_threshold = 20
    radiance_min_threshold = 4.5
    temperature_max_threshold = 350
    temperature_min_threshold = 250
    default_nan_temperature = 0
    soil_radiance_variance = 2
    leaf_radiance_variance = 1
    sensor_noise = 0.5
    pure_soil_threshold = 0.05
    pure_leaf_threshold = 0.95
    fvc_max_threshold = 0.99
    fvc_min_threshold = 0.01
    fd0 =  np.asarray([[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    number_component = 2
    number_angle     = 2
    number_value_threshold = 12
    

    [nl,ns] = np.shape(lst_n)
    result = np.zeros([number_component, nl, ns])

    ######### FVC
    fvc_n = (ndvi_n - ndvi_min) / (ndvi_max - ndvi_min)
    fvc_n[fvc_n > fvc_max_threshold] = fvc_max_threshold
    fvc_n[fvc_n < fvc_min_threshold] = fvc_min_threshold

    fvc_o = (ndvi_o - ndvi_min) / (ndvi_max - ndvi_min)
    fvc_o[fvc_o > fvc_max_threshold] = fvc_max_threshold
    fvc_o[fvc_o < fvc_min_threshold] = fvc_min_threshold

    ####### The radiance
    rad_n = planck(wl, lst_n)
    rad_n = rad_n * (fvc_n * emis_v + (1 - fvc_n) * emis_s)
    rad_o = planck(wl, lst_o)
    rad_o = rad_o * (fvc_o * emis_v + (1 - fvc_o) * emis_s)

    # component effective emissivity matrix
    W1 = np.asarray([emis_s * (1 - fvc_n), emis_v * (fvc_n)])
    W2 = np.asarray([emis_s * (1 - fvc_o), emis_v * (fvc_o)])


    nsi = np.zeros([nl, ns])
    nli = np.zeros([nl, ns])
    for k in range(nl):
        nli[k, :] = k
        nsi[k, :] = np.linspace(0, ns - 1, ns)

    mask = np.zeros([nl,ns])
    ind = (fvc_n > pure_soil_threshold) * (fvc_n < pure_leaf_threshold) *\
          (rad_n > radiance_min_threshold) * (rad_n < radiance_max_threshold) *\
          (rad_o > radiance_min_threshold) * (rad_o < radiance_max_threshold)
    mask[ind] = 1

    part = 6

    for k1 in range(off, nl - off):
        for k2 in range(off, ns - off):

            if (mask[k1, k2] != 1): continue

            nl1 = k1 - off
            nl2 = k1 + off + 1
            ns1 = k2 - off
            ns2 = k2 + off + 1

            offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)
            ns_temp = np.reshape(nsi[nl1:nl2, ns1:ns2] - k2, -1)
            nl_temp = np.reshape(nli[nl1:nl2, ns1:ns2] - k1, -1)
            fd = fd0
            fd = np.reshape(fd, -1) * 1.0

            ### w1 for nadir and w2 for oblique
            w1 = W1[:, nl1:nl2, ns1:ns2]
            w1 = np.asarray(np.reshape(w1, (2, -1)))
            w2 = W2[:, nl1:nl2, ns1:ns2]
            w2 = np.asarray(np.reshape(w2, (2, -1)))

            n0 = np.ones(len(ns_temp))
            x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
            w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
            w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])

            b1 = (np.reshape(rad_n[k1, k2], -1))
            b2 = (np.reshape(rad_o[k1, k2], -1))
            bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
            bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))

            temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
            temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
            ind = (bb1 < radiance_max_threshold) * (bb1 > radiance_min_threshold) * \
                  (temp1 > pure_soil_threshold) * (temp1 < pure_leaf_threshold) *\
                  (temp2 > pure_soil_threshold) *(temp2 <pure_leaf_threshold)*\
                  (bb2 < radiance_max_threshold) * (bb2 > radiance_min_threshold)

            num_val_point = np.sum(ind)
            if (num_val_point < number_value_threshold): continue

            ww = np.hstack((w11[:, ind], w22[:, ind]))
            bb = np.hstack((bb1[ind], bb2[ind]))

            w = np.hstack((w1, w2))
            b = np.hstack((b1, b2))

            ###############################################
            ### least-square method
            ###############################################
            ww = np.transpose(ww)

            coeffprior = np.asarray(lstsq(ww, bb,rcond=None))[0]
            # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
            # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
            # ctprior = np.transpose(np.matrix([cd1, cd2]))

            ###############################################
            ### Bayes method with priori information
            ###############################################
            # bbb = sorted(b)
            # halfpoint = np.int(num_point / 2)
            # rleaf = np.average(bbb[0])
            # rsoil = np.average(bbb[1])
            # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
            # bb = np.transpose(np.matrix(bb))
            # ww = np.transpose(ww)
            # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
            # Cm = np.matrix(np.diag(np.ones(2 * part)))
            # Cm[0, 0] = 15.0
            # Cm[part, part] = 10.0
            # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
            # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)

            #################################################
            ### first result without averaging information
            #################################################
            # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
            # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
            # ctprior = np.transpose(np.matrix([cd1, cd2]))

            ################################################
            ### first reult with averaing information
            ################################################
            coeffprior = np.asarray(coeffprior)
            ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
            ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])
            indnew = ind * (ct2 < radiance_max_threshold) * (ct2 > radiance_min_threshold) *\
                     (ct1 > radiance_min_threshold) * (ct1 < radiance_max_threshold)
            apoint = np.sum(indnew)
            if (apoint <= 1):
                cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
                cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
                ctprior = np.transpose(np.matrix([cd1, cd2]))
            else:
                fd[indnew] = fd[indnew] / np.sum(fd[indnew])

                # inver-distance weight but not used
                # bb0 = bb1[offcenter]
                # fdnew = np.zeros(25)
                # dismax = max(abs(bb1[indnew] - bb0))
                # if dismax != 0:
                #     fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
                # else:
                #     fdnew[:] = 1
                # fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
                # fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])

                cd1new = np.sum(ct1[indnew] * fd[indnew])
                cd2new = np.sum(ct2[indnew] * fd[indnew])
                ctprior = np.transpose(np.matrix([cd1new, cd2new]))

            #################################################
            ### second result by combining one-point result
            #################################################

            nl1 = k1 - off_
            nl2 = k1 + off_ + 1
            ns1 = k2 - off_
            ns2 = k2 + off_ + 1
            w1 = W1[:, nl1:nl2, ns1:ns2]
            w1 = np.asarray(np.reshape(w1, (2, -1)))
            w2 = W2[:, nl1:nl2, ns1:ns2]
            w2 = np.asarray(np.reshape(w2, (2, -1)))
            b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
            b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
            w = np.hstack((w1, w2))
            b = np.hstack((b1, b2))
            w = np.transpose(w)

            #################################################
            ### least-squared method
            #################################################
            # coeffprior = np.asarray(lstsq(w, b))[0]

            #################################################
            ### Bayes method
            #################################################

            w = np.matrix(w)
            b = np.matrix(b)
            b = np.transpose(b)

            Cd = np.matrix(np.diag(np.ones(number_component)) * sensor_noise)
            Cm = np.diag([soil_radiance_variance, leaf_radiance_variance])
            pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
            coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
            ### another solution for using priori information
            #coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
            component_radiance = np.asarray(coeff)

            if (np.min(component_radiance) < radiance_min_threshold): continue
            if (np.max(component_radiance) > radiance_max_threshold): continue

            component_temperature = inv_planck(wl, component_radiance)

            result[:, k1, k2] = [component_temperature[0], component_temperature[1]]

    result[result < temperature_min_threshold] = default_nan_temperature
    result[result > temperature_max_threshold] = default_nan_temperature

    return result


def inv_process(bb,ww,b,w,x,ind,wl=10.8):
 #   return 1
    part = 6
    radiance_max_threshold = 17
    radiance_min_threshold = 6.5
    temperature_max_threshold = 350
    temperature_min_threshold = 250
    default_nan_temperature = 0
    soil_radiance_variance = 2
    leaf_radiance_variance = 1
    sensor_noise = 0.5
    pure_soil_threshold = 0.05
    pure_leaf_threshold = 0.95
    fvc_max_threshold = 0.99
    fvc_min_threshold = 0.01
    fd0 = np.asarray([[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    number_component = 2
    fd0 = np.asarray([[1, 4, 7, 4, 1], [4, 16, 26, 16, 4], [7, 26, 41, 26, 7], [4, 16, 26, 16, 4], [1, 4, 7, 4, 1]])
    number_angle = 2
    off = 2
    off_ = 1
    number_value_threshold = 12
    offcenter = np.int(((2 * off + 1) * (2 * off + 1) - 1) / 2)

    fd = fd0
    fd = np.reshape(fd, -1) * 1.0


    ###############################################
    ### least-square method
    ###############################################
    ww = np.transpose(ww)
    coeffprior = np.asarray(lstsq(ww, bb,rcond=None))[0]
    # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    # ctprior = np.transpose(np.matrix([cd1, cd2]))

    ###############################################
    ### Bayes method with priori information
    ###############################################
    # bbb = sorted(b)
    # halfpoint = np.int(num_point / 2)
    # rleaf = np.average(bbb[0])
    # rsoil = np.average(bbb[1])
    # coeffprior = np.transpose(np.matrix([rsoil, 0, 0, 0, 0, 0, rleaf, 0, 0, 0, 0, 0]))
    # bb = np.transpose(np.matrix(bb))
    # ww = np.transpose(ww)
    # Cd = np.matrix(np.diag(np.ones(num_point * 2)) * Cdvaluep)
    # Cm = np.matrix(np.diag(np.ones(2 * part)))
    # Cm[0, 0] = 15.0
    # Cm[part, part] = 10.0
    # pre = np.linalg.inv(np.transpose(ww) * np.linalg.inv(Cd) * ww + np.linalg.inv(Cm))
    # ctprior = pre * (np.transpose(ww) * np.linalg.inv(Cd) * bb + np.linalg.inv(Cm) * coeffprior)

    #################################################
    ### first result without averaging information
    #################################################
    # cd1 = np.sum(coeffprior[:6] * x[:, offcenter])
    # cd2 = np.sum(coeffprior[6:] * x[:, offcenter])
    # ctprior = np.transpose(np.matrix([cd1, cd2]))

    ################################################
    ### first reult with averaing information
    ################################################
    coeffprior = np.asarray(coeffprior)
    ct1 = np.asarray([np.sum(coeffprior[:part] * x[:part, j]) for j in range(25)])
    ct2 = np.asarray([np.sum(coeffprior[part:2 * part] * x[:part, j]) for j in range(25)])

    indnew = ind * (ct2 < radiance_max_threshold) * (ct2 > radiance_min_threshold) * \
             (ct1 > radiance_min_threshold) * (ct1 < radiance_max_threshold)
    apoint = np.sum(indnew)
    if (apoint <= 1):
        cd1 = np.sum(coeffprior[:part] * x[:, offcenter])
        cd2 = np.sum(coeffprior[part:] * x[:, offcenter])
        ctprior = np.transpose(np.matrix([cd1, cd2]))
    else:
        fd[indnew] = fd[indnew] / np.sum(fd[indnew])

        # inver-distance weight but not used
        # bb0 = bb1[offcenter]
        # fdnew = np.zeros(25)
        # dismax = max(abs(bb1[indnew] - bb0))
        # if dismax != 0:
        #     fdnew[indnew] = 1.0 - abs(1.0 * bb1[indnew] - bb0) / dismax
        # else:
        #     fdnew[:] = 1
        # fdnew[indnew] = fdnew[indnew] / np.sum(fdnew[indnew])
        # fd[indnew] = (fdnew[indnew] * fd[indnew]) / np.sum(fdnew[indnew] * fd[indnew])

        cd1new = np.sum(ct1[indnew] * fd[indnew])
        cd2new = np.sum(ct2[indnew] * fd[indnew])
        ctprior = np.transpose(np.matrix([cd1new, cd2new]))

    #################################################
    ### second result by combining one-point result
    #################################################


    #################################################
    ### least-squared method
    #################################################
    # coeffprior = np.asarray(lstsq(w, b))[0]

    #################################################
    ### Bayes method
    #################################################

    w = np.matrix(w)
    b = np.matrix(b)
    b = np.transpose(b)

    Cd = np.matrix(np.diag(np.ones(number_component)) * sensor_noise)
    Cm = np.diag([soil_radiance_variance, leaf_radiance_variance])
    pre = np.linalg.inv(np.transpose(w) * np.linalg.inv(Cd) * w + np.linalg.inv(Cm))
    coeff = pre * (np.transpose(w) * np.linalg.inv(Cd) * b + np.linalg.inv(Cm) * ctprior)
    ### another solution for using priori information
    # coeff = ctprior + np.transpose(Cd * np.transpose(w) * pre) * (b - np.transpose(w) * ctprior)
    component_radiance = np.asarray(coeff)

    if (np.min(component_radiance) < radiance_min_threshold): return [0,0]
    if (np.max(component_radiance) > radiance_max_threshold): return [0,0]

    component_temperature = inv_planck(wl, component_radiance)

    result = [component_temperature[0], component_temperature[1]]

    return result

def inv_mpl_agl_multiple(lst_n, lst_o, emis_s, emis_v, ndvi_n, ndvi_o, wl=10.5, ndvi_min=0.05, ndvi_max=0.95, off=2, off_=0):
    radiance_max_threshold = 17
    radiance_min_threshold = 6.5
    temperature_max_threshold = 350
    temperature_min_threshold = 250
    default_nan_temperature = 0
    soil_radiance_variance = 2
    leaf_radiance_variance = 1
    sensor_noise = 0.5
    pure_soil_threshold = 0.05
    pure_leaf_threshold = 0.95
    fvc_max_threshold = 0.99
    fvc_min_threshold = 0.01

    number_component = 2
    number_angle = 2
    part = 6
    number_value_threshold = 12

    [nl, ns] = np.shape(lst_n)
    result = np.zeros([number_component, nl, ns])
    result1 = np.zeros([ nl, ns])
    result2 = np.zeros([ nl, ns])
    ######### FVC
    fvc_n = (ndvi_n - ndvi_min) / (ndvi_max - ndvi_min)
    fvc_n[fvc_n > fvc_max_threshold] = fvc_max_threshold
    fvc_n[fvc_n < fvc_min_threshold] = fvc_min_threshold

    fvc_o = (ndvi_o - ndvi_min) / (ndvi_max - ndvi_min)
    fvc_o[fvc_o > fvc_max_threshold] = fvc_max_threshold
    fvc_o[fvc_o < fvc_min_threshold] = fvc_min_threshold

    ####### The radiance
    rad_n = planck(wl, lst_n)
    rad_n = rad_n * (fvc_n * emis_v + (1 - fvc_n) * emis_s)
    rad_o = planck(wl, lst_o)
    rad_o = rad_o * (fvc_o * emis_v + (1 - fvc_o) * emis_s)

    # component effective emissivity matrix
    W1 = np.asarray([emis_s * (1 - fvc_n), emis_v * (fvc_n)])
    W2 = np.asarray([emis_s * (1 - fvc_o), emis_v * (fvc_o)])

    nsi = np.zeros([nl, ns])
    nli = np.zeros([nl, ns])
    for k in range(nl):
        nli[k, :] = k
        nsi[k, :] = np.linspace(0, ns - 1, ns)

    mask = np.zeros([nl, ns])
    ind0 = (fvc_n > pure_soil_threshold) * (fvc_n < pure_leaf_threshold) * \
          (rad_n > radiance_min_threshold) * (rad_n < radiance_max_threshold) * \
          (rad_o > radiance_min_threshold) * (rad_o < radiance_max_threshold) * \
          (nli>off) *(nli<(nl-off))*\
          (nsi>off)*(nsi<(ns-off))

    ind0 = np.where(ind0 > 0)
    ind0 = np.asarray(ind0)
    #valuable_position = nli[ind0]*ns+nsi[ind0]
    ns_temp = np.asarray([-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2,-2,-1,0,1,2])
    nl_temp = np.asarray(np.repeat([-2,-1,0,1,2],5))


    tasks = []
    inds = []
    [num_col,num_position]  = np.shape(ind0)
    pool = Pool(min(18, cpu_count()))
    num_position = 10
    num = 0
    for kn in range(num_position):
        k1 = np.int(ind0[0,kn])
        k2 = np.int(ind0[1,kn])
        nl1 = k1 - off
        nl2 = k1 + off + 1
        ns1 = k2 - off
        ns2 = k2 + off + 1

        w1 = W1[:, nl1:nl2, ns1:ns2]
        w1 = np.asarray(np.reshape(w1, (2, -1)))
        w2 = W2[:, nl1:nl2, ns1:ns2]
        w2 = np.asarray(np.reshape(w2, (2, -1)))

        bb1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
        bb2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))

        temp1 = (np.reshape(fvc_n[nl1:nl2, ns1:ns2], -1))
        temp2 = (np.reshape(fvc_o[nl1:nl2, ns1:ns2], -1))
        ind = (bb1 < radiance_max_threshold) * (bb1 > radiance_min_threshold) * \
              (temp1 > pure_soil_threshold) * (temp1 < pure_leaf_threshold) * \
              (temp2 > pure_soil_threshold) * (temp2 < pure_leaf_threshold) * \
              (bb2 < radiance_max_threshold) * (bb2 > radiance_min_threshold)

        num_val_point = np.sum(ind)
        if (num_val_point < number_value_threshold): continue



        ### w1 for nadir and w2 for oblique

        n0 = np.ones(len(ns_temp))
        x = np.asarray([n0, ns_temp, nl_temp, ns_temp * nl_temp, ns_temp * ns_temp, nl_temp * nl_temp])
        w11 = np.asarray([w1[j, :] * x[i, :] for j in range(2) for i in range(part)])
        w22 = np.asarray([w2[j, :] * x[i, :] for j in range(2) for i in range(part)])

        ww = np.hstack((w11[:, ind], w22[:, ind]))
        bb = np.hstack((bb1[ind], bb2[ind]))

        nl1 = k1 - off_
        nl2 = k1 + off_ + 1
        ns1 = k2 - off_
        ns2 = k2 + off_ + 1
        w1 = W1[:, nl1:nl2, ns1:ns2]
        w1 = np.asarray(np.reshape(w1, (2, -1)))
        w2 = W2[:, nl1:nl2, ns1:ns2]
        w2 = np.asarray(np.reshape(w2, (2, -1)))
        b1 = (np.reshape(rad_n[nl1:nl2, ns1:ns2], -1))
        b2 = (np.reshape(rad_o[nl1:nl2, ns1:ns2], -1))
        w = np.hstack((w1, w2))
        b = np.hstack((b1, b2))
        w = np.transpose(w)

        # w = np.hstack((W1[:,k1,k2], W2[:,k1,k2]))
        # b = np.hstack((rad_n[k1,k2], rad_o[k1,k2]))

        #arg = (bb,ww,b,w,x,ind)
        #lsct = inv_process(bb,ww,b,w,x,ind)

        tasks.append(pool.apply_async(inv_process,(bb,ww,b,w,x,ind,)))
        inds.append(ind0[:,kn])
        #num = num + 1


    pool.close()
    pool.join()

    if (len(tasks) < 1): return result
    dataMerge = np.asarray([task.get() for task in tasks])
    inds = np.asarray(inds)


    result1[inds] = dataMerge[:,0]
    result2[inds] = dataMerge[:,1]
    result[0,:,:] = result1
    result[1,:,:] = result2

    result[result < temperature_min_threshold] = default_nan_temperature
    result[result > temperature_max_threshold] = default_nan_temperature
    ind = np.isnan(result)
    result[ind] = 0

    return result

