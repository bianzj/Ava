
from util.myfun import *
from sklearn.linear_model import LinearRegression
from scipy.optimize import root,fsolve
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy.optimize import fmin


############################################3
#### residual function
###########################################3

def fun_DTC(x, t, mea):

    T0 = x[0]
    Ta = x[1]
    tm = x[2]
    tsr = x[3]

    omg = 4/3.0*(tm-tsr)
    mod = T0 + Ta*np.cos(np.pi/omg*(t-tm))
    res = mod - mea

    return np.sum(res*res)

def model_DTC(x,t):
    T0 = x[0]
    Ta = x[1]
    tm = x[2]
    tsr = x[3]

    omg = 4 / 3.0 * (tm - tsr)
    mod = T0 + Ta * np.cos(np.pi / omg * (t - tm))

    return mod

#############################################
### fitting USING FMIN
###############################################3

def fitting_DTC(t, mea):

    res = fmin(fun_DTC, x0 = [6,16,14,6],args = (t, mea),disp = 0)

    return res

