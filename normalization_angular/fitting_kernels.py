
from util.myfun import *
from sklearn.linear_model import LinearRegression
from scipy.optimize import root,fsolve
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.optimize import fmin


rd = np.pi/180.0

############################################
####  Kernels
############################################

def kernel_LSFLI(sza, vza, raa):

    sza
    ind = raa > 180
    raa[ind] = 360 - raa[ind]
    thetas = np.deg2rad(sza)
    thetaa = np.deg2rad(raa)
    thetav = np.deg2rad(abs(vza))

    kz = np.zeros(np.size(vza))
    kz[:] = 1.0

    az = np.cos(thetas) * np.cos(thetav) + \
         np.sin(thetas) * np.sin(thetav) * np.cos(thetaa)

    DD = (np.power(np.tan(thetas), 2) + np.power(np.tan(thetav), 2) - 2 * np.tan(thetas) * np.tan(thetav) *
          np.cos(thetaa))
    DD[DD < 0] = 0
    D = np.sqrt(DD)

    cost = np.sqrt(D * D + np.power(np.tan(thetav) * np.tan(thetas) * np.sin(thetaa), 2)) \
           / (1.0 / np.cos(thetas) + 1.0 / np.cos(thetav))
    cost[cost < -1] = -0.999
    cost[cost > 1] = 0.999

    t = np.arccos(cost)
    O = (1 / np.pi) * (t - np.sin(t) * np.cos(t)) * (1.0 / np.cos(thetas) + 1.0 / np.cos(thetav))

    # Klsf = ((1 + 2 * np.cos(thetav)) / (np.sqrt(0.96) + 2 * 0.96 * np.cos(thetav)) - 0.25 * np.cos(thetav)
    #         / (1 + 2 * np.cos(thetav)) + 0.15 * (1 - np.exp(-0.75 / np.cos(thetav))))

    Klsf = 0.069*np.cos(thetav)*np.cos(thetav)-0.215*np.cos(thetav)+1.176

    # if sza <= 75:
    Kli = (1 + az) * (1.0 / np.cos(thetas) * (1.0 / np.cos(thetav))) / (1.0 / np.cos(thetav) + 1 / np.cos(thetas) - O) - 2
    Kli[sza >= 75] = 0
    # else:
    #     Kli = 0

    return kz,Klsf,Kli

def kernel_LSF(vza):
    thetav = np.deg2rad(abs(vza))
    Klsf = 0.069*np.cos(thetav)*np.cos(thetav)-0.215*np.cos(thetav)+1.176
    return Klsf



def kernel_LSRRoss(sza, vza, raa):

    ind = raa > 180
    raa[ind] = 360 - raa[ind]
    thetas = np.deg2rad(sza)
    thetaa = np.deg2rad(raa)
    thetav = np.deg2rad(abs(vza))

    kz = np.zeros(np.size(vza))
    kz[:] = 1.0

    az = np.cos(thetas) * np.cos(thetav) + \
         np.sin(thetas) * np.sin(thetav) * np.cos(thetaa)

    DD = (np.power(np.tan(thetas), 2) + np.power(np.tan(thetav), 2) - 2 * np.tan(thetas) * np.tan(thetav) *
          np.cos(thetaa))
    DD[DD < 0] = 0
    D = np.sqrt(DD)

    cost = np.sqrt(D * D + np.power(np.tan(thetav) * np.tan(thetas) * np.sin(thetaa), 2)) \
           / (1.0 / np.cos(thetas) + 1.0 / np.cos(thetav))
    cost[cost < -1] = -0.999
    cost[cost > 1] = 0.999

    t = np.arccos(cost)
    O = (1 / np.pi) * (t - np.sin(t) * np.cos(t)) * (1.0 / np.cos(thetas) + 1.0 / np.cos(thetav))
    Klsr = O - (1.0 / np.cos(thetas))-(1.0/np.cos(thetav))+0.5*(1+az)*(1.0 / np.cos(thetas))*(1.0/np.cos(thetav))

    angle = np.arccos(az)

    Kross = ((np.pi/2.0 - angle)*az + np.sin(angle))/(np.cos(thetav)+np.cos(thetas))-np.pi/4.0
    # Klsf = ((1 + 2 * np.cos(thetav)) / (np.sqrt(0.96) + 2 * 0.96 * np.cos(thetav)) - 0.25 * np.cos(thetav)
    #         / (1 + 2 * np.cos(thetav)) + 0.15 * (1 - np.exp(-0.75 / np.cos(thetav))))
    # Kli = (1 + az) * (1.0 / np.cos(thetas) * (1.0 / np.cos(thetav))) / (1.0 / np.cos(thetav) + 1 / np.cos(thetas) - O) - 2

    return kz,Klsr,Kross

############################################3
#### forward model
###########################################3

def model_Vin(x,coeffa):

    lst1 = x[0]
    lst2 = x[1]
    vza1 = x[2]
    vza2 = x[3]
    phi1 = 1 - np.cos(vza1*rd)
    phi2 = 1 - np.cos(vza2*rd)
    x = phi1*lst2 - phi2*lst1
    return x*coeffa

def model_RLE(x,b,k):

    a = x[8]
    rad = x[9]
    lst1 = x[0]
    lst2 = x[1]
    vza1 = x[2]
    vza2 = x[3]
    sza1 = x[4]
    sza2 = x[5]
    psi1 = x[6]
    psi2 = x[7]


    thetas1 = np.deg2rad(sza1)
    thetaa1 = np.deg2rad(psi1)
    thetav1 = np.deg2rad(abs(vza1))

    thetas2 = np.deg2rad(sza2)
    thetaa2 = np.deg2rad(psi2)
    thetav2 = np.deg2rad(abs(vza2))


    f1 = np.sqrt((np.power(np.tan(thetas1), 2) + np.power(np.tan(thetav1), 2) - 2 * np.tan(thetas1) * np.tan(thetav1) *
                 np.cos(thetaa1)))
    f2 = np.sqrt((np.power(np.tan(thetas2), 2) + np.power(np.tan(thetav2), 2) - 2 * np.tan(thetas2) * np.tan(thetav2) *
                 np.cos(thetaa2)))
    fn = np.tan(thetas1)
    # phi1 = 1 - np.cos(thetav1)
    # phi2 = 1 - np.cos(thetav2)
    phi1 = 0.069*np.cos(vza1*rd)*np.cos(vza1*rd)-0.215*np.cos(vza1*rd)+1.176
    phi2 = 0.069*np.cos(vza2*rd)*np.cos(vza2*rd)-0.215*np.cos(vza2*rd)+1.176

    return a * (phi1*lst2 - phi2*lst1) + b * rad *(np.exp(-k*f1)-np.exp(-k*f2))/(1-np.exp(-k*fn))

################################333
###  fun
####################################3

def fun_Vin(x,lst1,lst2,vza1,vza2):

    phi1 = 1 - np.cos(vza1*rd)
    phi2 = 1 - np.cos(vza2*rd)
    y = lst1 - lst2
    res = (phi1*lst2 - phi2*lst1)*x - y

    return np.sum(res**2)

def fun_LSF(x,lst1,lst2,vza1,vza2):

    phi1 = 0.069*np.cos(vza1*rd)*np.cos(vza1*rd)-0.215*np.cos(vza1*rd)+1.176
    phi2 = 0.069*np.cos(vza2*rd)*np.cos(vza2*rd)-0.215*np.cos(vza2*rd)+1.176
    y = lst1 - lst2
    res = (phi1*lst2 - phi2*lst1)*x - y

    return np.sum(res**2)

def fun_VinRLE_rad_bk(x,lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,a,rad):

    b = x[0]
    k = x[1]
    thetas1 = np.deg2rad(sza1)
    thetaa1 = np.deg2rad(psi1)
    thetav1 = np.deg2rad(abs(vza1))

    thetas2 = np.deg2rad(sza2)
    thetaa2 = np.deg2rad(psi2)
    thetav2 = np.deg2rad(abs(vza2))


    f1 = np.sqrt((np.power(np.tan(thetas1), 2) + np.power(np.tan(thetav1), 2) - 2 * np.tan(thetas1) * np.tan(thetav1) *
                 np.cos(thetaa1)))
    f2 = np.sqrt((np.power(np.tan(thetas2), 2) + np.power(np.tan(thetav2), 2) - 2 * np.tan(thetas2) * np.tan(thetav2) *
                 np.cos(thetaa2)))
    fn = np.tan(thetas1)
    phi1 = 1 - np.cos(thetav1)
    phi2 = 1 - np.cos(thetav2)




    res = a * (phi1*lst2 - phi2*lst1) + b * rad *(np.exp(-k*f1)-np.exp(-k*f2))/(1-np.exp(-k*fn)) - (lst1 - lst2)
    return np.sum(res*res)

def fun_VinRLE_rad_full(x,lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,rad):

    a = x[0]
    b = x[1]
    k = x[2]
    thetas1 = np.deg2rad(sza1)
    thetaa1 = np.deg2rad(psi1)
    thetav1 = np.deg2rad(abs(vza1))

    thetas2 = np.deg2rad(sza2)
    thetaa2 = np.deg2rad(psi2)
    thetav2 = np.deg2rad(abs(vza2))


    f1 = np.sqrt((np.power(np.tan(thetas1), 2) + np.power(np.tan(thetav1), 2) - 2 * np.tan(thetas1) * np.tan(thetav1) *
                 np.cos(thetaa1)))
    f2 = np.sqrt((np.power(np.tan(thetas2), 2) + np.power(np.tan(thetav2), 2) - 2 * np.tan(thetas2) * np.tan(thetav2) *
                 np.cos(thetaa2)))
    fn = np.tan(thetas1)
    phi1 = 1 - np.cos(thetav1)
    phi2 = 1 - np.cos(thetav2)

    mea = lst1 - lst2
    mod = a * (phi1 * lst2 - phi2 * lst1) + b * rad * (np.exp(-k * f1) - np.exp(-k * f2)) / (1 - np.exp(-k * fn))
    ind = sza1 > 80
    mod_night = a * (phi1 * lst2 - phi2 * lst1)
    mod[ind] = mod_night[ind]


    res = mod - mea

    return np.sum(res*res)

def fun_VinRLE_ta_bk(x, lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, a, ta):

    b = x[0]
    k = x[1]
    thetas1 = np.deg2rad(sza1)
    thetaa1 = np.deg2rad(psi1)
    thetav1 = np.deg2rad(abs(vza1))

    thetas2 = np.deg2rad(sza2)
    thetaa2 = np.deg2rad(psi2)
    thetav2 = np.deg2rad(abs(vza2))


    f1 = np.sqrt((np.power(np.tan(thetas1), 2) + np.power(np.tan(thetav1), 2) - 2 * np.tan(thetas1) * np.tan(thetav1) *
                 np.cos(thetaa1)))
    f2 = np.sqrt((np.power(np.tan(thetas2), 2) + np.power(np.tan(thetav2), 2) - 2 * np.tan(thetas2) * np.tan(thetav2) *
                 np.cos(thetaa2)))
    fn = np.tan(thetas1)
    phi1 = 1 - np.cos(thetav1)
    phi2 = 1 - np.cos(thetav2)


    res = a * (phi1*lst2 - phi2*lst1) + b * ta * (np.exp(-k * f1) - np.exp(-k * f2)) / (1 - np.exp(-k * fn)) - (lst1 - lst2)
    return np.sum(res*res)

def fun_VinRLE_ta_full(x, lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, ta):

    a = x[0]
    b = x[1]
    k = x[2]
    thetas1 = np.deg2rad(sza1)
    thetaa1 = np.deg2rad(psi1)
    thetav1 = np.deg2rad(abs(vza1))

    thetas2 = np.deg2rad(sza2)
    thetaa2 = np.deg2rad(psi2)
    thetav2 = np.deg2rad(abs(vza2))


    f1 = np.sqrt((np.power(np.tan(thetas1), 2) + np.power(np.tan(thetav1), 2) - 2 * np.tan(thetas1) * np.tan(thetav1) *
                 np.cos(thetaa1)))
    f2 = np.sqrt((np.power(np.tan(thetas2), 2) + np.power(np.tan(thetav2), 2) - 2 * np.tan(thetas2) * np.tan(thetav2) *
                 np.cos(thetaa2)))
    fn = np.tan(thetas1)
    phi1 = 1 - np.cos(thetav1)
    phi2 = 1 - np.cos(thetav2)

    mea = lst1 - lst2
    mod = a * (phi1 * lst2 - phi2 * lst1) + b * ta * (np.exp(-k * f1) - np.exp(-k * f2)) / (1 - np.exp(-k * fn))
    ind = sza1 > 80
    mod_night = a * (phi1 * lst2 - phi2 * lst1)
    mod[ind] = mod_night[ind]

    res = mod - mea

    return np.sum(res*res)

def fun_LSFRLE_ta_bk(x, lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, a, ta):

    b = x[0]
    k = x[1]
    thetas1 = np.deg2rad(sza1)
    thetaa1 = np.deg2rad(psi1)
    thetav1 = np.deg2rad(abs(vza1))

    thetas2 = np.deg2rad(sza2)
    thetaa2 = np.deg2rad(psi2)
    thetav2 = np.deg2rad(abs(vza2))


    f1 = np.sqrt((np.power(np.tan(thetas1), 2) + np.power(np.tan(thetav1), 2) - 2 * np.tan(thetas1) * np.tan(thetav1) *
                 np.cos(thetaa1)))
    f2 = np.sqrt((np.power(np.tan(thetas2), 2) + np.power(np.tan(thetav2), 2) - 2 * np.tan(thetas2) * np.tan(thetav2) *
                 np.cos(thetaa2)))
    fn = np.tan(thetas1)
    # phi1 = 1 - np.cos(thetav1)
    # phi2 = 1 - np.cos(thetav2)

    phi1 = 0.069*np.cos(vza1*rd)*np.cos(vza1*rd)-0.215*np.cos(vza1*rd)+1.176
    phi2 = 0.069*np.cos(vza2*rd)*np.cos(vza2*rd)-0.215*np.cos(vza2*rd)+1.176


    res = a * (phi1*lst2 - phi2*lst1) + b * ta * (np.exp(-k * f1) - np.exp(-k * f2)) / (1 - np.exp(-k * fn)) - (lst1 - lst2)
    return np.sum(res*res)

def fun_LSFRLE_ta_full(x, lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, ta):

    a = x[0]
    b = x[1]
    k = x[2]
    thetas1 = np.deg2rad(sza1)
    thetaa1 = np.deg2rad(psi1)
    thetav1 = np.deg2rad(abs(vza1))

    thetas2 = np.deg2rad(sza2)
    thetaa2 = np.deg2rad(psi2)
    thetav2 = np.deg2rad(abs(vza2))


    f1 = np.sqrt((np.power(np.tan(thetas1), 2) + np.power(np.tan(thetav1), 2) - 2 * np.tan(thetas1) * np.tan(thetav1) *
                 np.cos(thetaa1)))
    f2 = np.sqrt((np.power(np.tan(thetas2), 2) + np.power(np.tan(thetav2), 2) - 2 * np.tan(thetas2) * np.tan(thetav2) *
                 np.cos(thetaa2)))
    fn = np.tan(thetas1)
    # phi1 = 1 - np.cos(thetav1)
    # phi2 = 1 - np.cos(thetav2)

    phi1 = 0.069*np.cos(vza1*rd)*np.cos(vza1*rd)-0.215*np.cos(vza1*rd)+1.176
    phi2 = 0.069*np.cos(vza2*rd)*np.cos(vza2*rd)-0.215*np.cos(vza2*rd)+1.176

    mea = lst1 - lst2
    mod = a * (phi1 * lst2 - phi2 * lst1) + b * ta * (np.exp(-k * f1) - np.exp(-k * f2)) / (1 - np.exp(-k * fn))
    ind = sza1 > 80
    mod_night = a * (phi1 * lst2 - phi2 * lst1)
    mod[ind] = mod_night[ind]

    res = mod - mea

    return np.sum(res*res)

#############################################
### fitting USING FMIN
###############################################3

def fitting_Vin(lst1, lst2, vza1,vza2):

    res = fmin(fun_Vin, x0 = [-0.05],args = (lst1,lst2,vza1,vza2),disp = 0)

    return res

def fitting_LSF(lst1, lst2, vza1,vza2):

    res = fmin(fun_LSF, x0 = [-0.05],args = (lst1,lst2,vza1,vza2),disp = 0)

    return res


def fitting_VinRLE_rad_bk(lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, a, rad):

    res = fmin(fun_VinRLE_rad_bk, x0 = [0.5,0.5],args = (lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,a,rad),disp = 0)

    return res

def fitting_VinRLE_rad_full(lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, rad):

    res = fmin(fun_VinRLE_rad_full, x0 = [0.5,0.5,0.5],args = (lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,rad),disp = 0)

    return res

def fitting_VinRLE_ta_bk(lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, a, rad):

    res = fmin(fun_VinRLE_ta_bk, x0 = [0.5,0.5],args = (lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,a,rad),disp = 0)

    return res

def fitting_VinRLE_ta_full(lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, rad):

    res = fmin(fun_VinRLE_ta_full, x0 = [0.5,0.5,0.5],args = (lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,rad),disp = 0)

    return res

def fitting_LSFRLE_ta_bk(lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, a, ta):

    res = fmin(fun_LSFRLE_ta_bk, x0 = [0.5,0.5],args = (lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,a,ta),disp = 0)

    return res

def fitting_LSFRLE_ta_full(lst1, lst2, vza1, vza2, sza1, sza2, psi1, psi2, ta):

    res = fmin(fun_LSFRLE_ta_full, x0 = [0.5,0.5,0.5],args = (lst1,lst2,vza1,vza2,sza1,sza2,psi1,psi2,ta),disp = 0)

    return res


########################################################
####  FITTING USING LINEARREGRESSION
#######################################################

def fitting_LSRROSS(sza, vza, raa, dbt):

    std = np.std(dbt)
    kz,Klsf,Kli = kernel_LSRRoss(sza, vza, raa)
    A = np.asarray([kz, Klsf, Kli])
    A = np.transpose(A)
    model = LinearRegression(copy_X=True, fit_intercept=True, n_jobs=-1, normalize=False)
    model.fit(A, dbt)
    newdbt = model.predict(A)
    coeffs = model.coef_
    coeffs[0] = model.intercept_
    dif = newdbt - dbt
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif * dif))
    # maxdif = np.max(dbt)-np.min(dbt)
    # kz,Klsf,Kli = kernel_LSRRoss(sza, np.asarray([0]), np.asarray([0]))
    # A = np.asarray([kz, Klsf, Kli])
    # A = np.transpose(A)
    # predicted = model.predict(A)

    return coeffs,rmse

def fitting_LSFLI_ta(sza,vza1,vza2,raa1,raa2,lst1,lst2,ta):

    kz,Klsf1,Kli1 = kernel_LSFLI(sza, vza1, raa1)
    kz,Klsf2,Kli2 = kernel_LSFLI(sza, vza2, raa2)

    W = np.asarray([Klsf1-Klsf2, Kli1-Kli2])
    W = np.transpose(W)
    y = (lst1 - lst2)/ta
    model = LinearRegression(copy_X=True, fit_intercept=True, n_jobs=-1, normalize=False)
    model.fit(W, y)
    newdbt = model.predict(W)
    coeffs = model.coef_
    coeffs[0] = model.intercept_
    dif = newdbt - y
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif * dif))
    # maxdif = np.max(dbt)-np.min(dbt)
    # kz,Klsf,Kli = kernel_LSRRoss(sza, np.asarray([0]), np.asarray([0]))
    # A = np.asarray([kz, Klsf, Kli])
    # A = np.transpose(A)
    # predicted = model.predict(A)

    return coeffs,rmse

def fitting_LSF_ta(vza1,vza2,lst1,lst2,ta):


    klsf1 = kernel_LSF(vza1)
    klsf2 = kernel_LSF(vza2)

    W = np.asarray([klsf1-klsf2])
    W = np.transpose(W)
    y = (lst1 - lst2)/ta
    model = LinearRegression(copy_X=True, fit_intercept=True, n_jobs=-1, normalize=False)
    model.fit(W, y)
    newdbt = model.predict(W)
    coeffs = model.coef_
    coeffs[0] = model.intercept_
    dif = newdbt - y
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif * dif))
    # maxdif = np.max(dbt)-np.min(dbt)
    # kz,Klsf,Kli = kernel_LSRRoss(sza, np.asarray([0]), np.asarray([0]))
    # A = np.asarray([kz, Klsf, Kli])
    # A = np.transpose(A)
    # predicted = model.predict(A)

    return coeffs,rmse


def fitting_LI_ta(sza,vza1,vza2,raa1,raa2,lst1,lst2,ta,coeffa):


    kz,Klsf1,Kli1 = kernel_LSFLI(sza, vza1, raa1)
    kz,Klsf2,Kli2 = kernel_LSFLI(sza, vza2, raa2)

    W = np.asarray([Kli1-Kli2])
    W = np.transpose(W)
    y = (lst1 - lst2)/ta - coeffa(Klsf1-Klsf2)
    model = LinearRegression(copy_X=True, fit_intercept=True, n_jobs=-1, normalize=False)
    model.fit(W, y)
    newdbt = model.predict(W)
    coeffs = model.coef_
    coeffs[0] = model.intercept_
    dif = newdbt - y
    bias = np.mean(dif)
    rmse = np.sqrt(np.mean(dif * dif))
    # maxdif = np.max(dbt)-np.min(dbt)
    # kz,Klsf,Kli = kernel_LSRRoss(sza, np.asarray([0]), np.asarray([0]))
    # A = np.asarray([kz, Klsf, Kli])
    # A = np.transpose(A)
    # predicted = model.predict(A)

    return coeffs,rmse