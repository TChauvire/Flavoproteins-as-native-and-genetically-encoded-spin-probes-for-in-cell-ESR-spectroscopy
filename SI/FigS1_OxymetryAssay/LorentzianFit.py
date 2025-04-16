# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 12:10:59 2024

@author: tim_t
"""
from ImportMultipleFiles import ImportMultipleNameFiles, OpenDoubleIntegral2, OpenDoubleIntegralKin2
import matplotlib.pyplot as plt
import numpy as np
from os import path, makedirs, getcwd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from os import path, makedirs, getcwd

folder = getcwd() + "\\250321_LiPC\\"  # '\\2D\\'  # "\\ThirdTry\\"
ListOfFiles = ImportMultipleNameFiles(folder, Extension='.DSC')
# FirstDeriv, ZeroDeriv, Header, _, _ = OpenDoubleIntegralKin2(
#     ListOfFiles[0], Scaling=None, polyorder=[1, 1], window=50, col=[0, 1, 2, 3, 4, 5, 6, 7])
FirstDeriv, ZeroDeriv, Header, _, _ = OpenDoubleIntegral2(
    ListOfFiles, Scaling=None, polyorder=[1, 1], window_length=50)


def ComputeRMSE(y, yfit, p0):
    '''
    Compute the normalized residual sum of square of the residual of a function
    or Root Mean Square of Error (RMSE)
    See by instance :
    https://statisticsbyjim.com/regression/root-mean-square-error-rmse/
    Script written by Timothée Chauviré 10/26/2023

    Parameters
    ----------
    y : experimental data
        TYPE : Numpy data array
    yfit : experimental data
        TYPE : Numpy data array
    p0 : paremeters used for the fit
        TYPE : Numpy data array
    Returns
    -------
    RMSE : normalized residual sum of square or Root mean square of error
        TYPE : real float value
    '''
    NumMeas = y.shape[0]
    NumParams = len(p0)
    resnorm = np.sum((y-yfit)**2)
    RMSE = np.sqrt(resnorm/(NumMeas - NumParams))
    return RMSE

# Fit Slope with a gaussian derivative curve: (According to Poole 1979:


def gaussianderiv(x, a, b, c):
    fun = np.exp(0.5)*a*(x-b)/(0.5*c)*np.exp(-0.5*((x-b)/(0.5*c))**2)
    return fun

# Fit Slope with a Loretzian derivative curve: (According to Poole 1979:
# Y'(H) = [16y'm*(H-H0)/(1/2DHpp)]/[3+((H-H0)/(1/2DHpp))^2]^2 (lorentzian)


def lorentzianderiv(x, a, b, c):
    fun = 16*a*((x-b)/(0.5*c))/((3+((x-b)/(0.5*c))**2)**2)
    return fun


# According to Peric et al 1994, "the Derivative Voigt Lineshape, [...], can be
# well approximated by the variably weighted sum of the derivative of the
# gaussian and lorentzian lineshapes ref (4,8).
# def voigtderiv(x, a, b, c, d):
#     fun = d*gaussianderiv(x, a, b, c)+(1-d)*lorentzianderiv(x, a, b, c)
#     return fun
# p1 = [1, center, DeltaH, 1]
# popt3, pcov3 = curve_fit(voigtderiv, x[:,].ravel(), y.ravel(),
#                          p0=p1, sigma=None, absolute_sigma=False,
#                          check_finite=None)
# yfit3 = voigtderiv(x, popt3[0], popt3[1], popt3[2],  popt3[3]).ravel()
# l4 = plt.plot(x.ravel(), yfit3.ravel(), 'g',
#               label='Voigt Fit', linewidth=2)
FitParam = np.full((7, 2*len(ListOfFiles)), np.nan, dtype="float")
FitParamHeader = list(np.zeros((8,)))

plt.close('all')
for i in range(len(ListOfFiles)):
    fileID = path.split(ListOfFiles[i])
    # plot the first derivative
    x = FirstDeriv[:, 4*i]
    y = FirstDeriv[:, 4*i+1]

    max_x = x[np.argmax(y)]
    min_x = x[np.argmin(y)]
    center = max_x+(min_x-max_x)/2
    DeltaH = min_x-max_x
    p0 = [1, center, DeltaH]
    popt1, pcov1 = curve_fit(gaussianderiv, x[:,].ravel(), y.ravel(),
                             p0=p0, sigma=None, absolute_sigma=False,
                             check_finite=None)
    yfit1 = gaussianderiv(x, popt1[0], popt1[1], popt1[2]).ravel()
    perr1 = np.sqrt(np.diag(pcov1))
    RMSE1 = ComputeRMSE(y, yfit1[:,], popt1)
    popt2, pcov2 = curve_fit(lorentzianderiv, x[:,].ravel(), y.ravel(),
                             p0=p0, sigma=None, absolute_sigma=False,
                             check_finite=None)
    yfit2 = lorentzianderiv(x, popt2[0], popt2[1], popt2[2]).ravel()
    perr2 = np.sqrt(np.diag(pcov2))
    RMSE2 = ComputeRMSE(y, yfit2[:,], popt2)
    FitParam[0, 2*i] = popt1[0]
    FitParam[1, 2*i] = perr1[0]
    FitParam[2, 2*i] = popt1[1]
    FitParam[3, 2*i] = perr1[1]
    FitParam[4, 2*i] = popt1[2]
    FitParam[5, 2*i] = perr1[2]
    FitParam[6, 2*i] = RMSE1
    FitParam[0, 2*i+1] = popt2[0]
    FitParam[1, 2*i+1] = perr2[0]
    FitParam[2, 2*i+1] = popt2[1]
    FitParam[3, 2*i+1] = perr2[1]
    FitParam[4, 2*i+1] = popt2[2]
    FitParam[5, 2*i+1] = perr2[2]
    FitParam[6, 2*i+1] = RMSE2
    plt.figure(i)
    plt.title(fileID[1])
    l1 = plt.plot(x.ravel(), y.ravel(), 'k',
                  label='First Derivative', linewidth=2)
    l2 = plt.plot(x.ravel(), yfit1.ravel(), 'r',
                  label='Gaussian Fit', linewidth=2)
    l3 = plt.plot(x.ravel(), yfit2.ravel(), 'b',
                  label='Lorentzian Fit', linewidth=2)
    plt.xlabel("Magnetic Field [G]", fontsize=14, weight='bold')
    plt.ylabel("Intensity [a.u.]", fontsize=14, weight='bold')
    plt.legend(fontsize=16)
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.tight_layout()
    makedirs(fileID[0] + '\\Figures\\', exist_ok=True)
    plt.savefig(fileID[0] + '\\Figures\\figurefit{0}'.format(i+1))


FitParamHeader[0] = 'Name of the Sample'
FitParamHeader[1] = 'a'  # amplitude factor
FitParamHeader[2] = 'Error (a)'
FitParamHeader[3] = 'Center Field H0 (G)'
FitParamHeader[4] = 'Error (Center Field)'
FitParamHeader[5] = 'Linewidth H(P-P) (G)'
FitParamHeader[6] = 'Error Linewidth H(P-P)'
FitParamHeader[7] = 'RMSE'

del center, DeltaH, figManager, fileID, i, l1, l2, l3, max_x, min_x, p0, pcov1
del pcov2, perr1, perr2, popt1, popt2, RMSE1, RMSE2, x, y, yfit1, yfit2
