# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 12:10:59 2024

@author: tim_t
"""
from ImportMultipleFiles_FigS8 import ImportMultipleNameFiles, MaxLengthOfFiles
from ImportMultipleFiles_FigS8 import OpenDavies, eprload, automatic_phase
from ImportMultipleFiles_FigS8 import basecorr1D, datasmooth
import numpy as np
from os import getcwd

# Import Davies Data
folder = getcwd() + '\\RawData\\ENDORMeasurements\\Davies\\'
ListOfFiles = ImportMultipleNameFiles(folder, Extension='.DSC')
maxlen = MaxLengthOfFiles(ListOfFiles)
fulldata, HeaderDavies = OpenDavies(
    ListOfFiles, Scaling=None, polyorder=1, window_length=50)

# Import Mims Data with tau averaging
folder2 = getcwd() + '\\RawData\\ENDORMeasurements\\Mims\\'
ListOfFiles2 = ImportMultipleNameFiles(folder2, Extension='.DSC')
y, x, par = eprload(ListOfFiles2[0], Scaling=None)
newy = np.sum(y[:, :], axis=1)
new_data2, _ = automatic_phase(vector=newy, pivot1=int(
    newy.shape[0]/2), funcmodel='minfunc')
data_real = -1*new_data2.real
data_imag = -1*new_data2.imag
data_real_new, _, _, _ = basecorr1D(
    x=x[:, 0], y=data_real, polyorder=1, window=50)
ncol = 7
npts = x[:, 0].shape[0]
fulldata2 = np.full((npts, ncol*len(ListOfFiles2)), np.nan, dtype=float)
fulldata2[0:npts, 0] = x[0:npts, 0]
fulldata2[0:npts, 1] = data_real[0:npts]
fulldata2[0:npts, 2] = data_imag[0:npts]
fulldata2[0:npts, 3] = np.absolute(data_real+1j*data_imag)
fulldata2[0:npts, 4] = data_real_new[0:npts]
fulldata2[0:npts, 5] = datasmooth(data_real_new[0:npts],
                                  window_length=4,
                                  method='flat')
d2 = fulldata2[0:npts, 5]
fulldata2[0:npts, 6] = (d2) / np.amax(d2)
HeaderMims = list(np.zeros((7,)))
HeaderMims[0] = par['XNAM']
HeaderMims[1] = par['TITL']+str("_real")
HeaderMims[2] = par['TITL']+str("_imag")
HeaderMims[3] = par['TITL']+str("_abs")
HeaderMims[4] = par['TITL']+str("_real_bckg")
HeaderMims[5] = par['TITL']+str("_real_4ptsSmoothed")
HeaderMims[6] = par['TITL']+str("_real_Normalized")

del d2, data_imag, data_real, data_real_new, folder, folder2, maxlen, ncol
del new_data2, newy, npts, par, x, y
