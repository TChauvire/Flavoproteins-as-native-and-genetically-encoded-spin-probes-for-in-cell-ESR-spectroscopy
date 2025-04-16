# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 12:10:59 2024

@author: tim_t
"""
from ImportMultipleFiles_FigS11 import ImportMultipleNameFiles
from ImportMultipleFiles_FigS11 import OpenDoubleIntegralKin
import numpy as np
from os import path, getcwd
import matplotlib.pyplot as plt

folder = getcwd() + '\\2D\\'
plt.close('all')
ListOfFiles = ImportMultipleNameFiles(folder, Extension='.DSC')
nfile = len(ListOfFiles)
IntgValue = np.full((800, 2*nfile), np.nan)
HeaderInt = list(np.zeros(2*nfile,))
ConversionFactor = 0.434  # Value coming from a tempo scale M=1G
# See Excel file in the same folder "Calculation_ESRDoubleIntegral.xlsx"
for i in range(len(ListOfFiles)):
    nfile = i
    Filename = ListOfFiles[nfile]
    _, _, _, IntgValueKin, _ = OpenDoubleIntegralKin(
        Filename, Scaling=None, polyorder=[1, 1], window=50)
    fileID = path.split(ListOfFiles[nfile])

    # Export all the data in a single matrix:
    nslice = IntgValueKin.shape[0]
    IntgValue[:nslice, 2*i] = IntgValueKin[:nslice, 0]
    IntgValue[:nslice, 2*i+1] = IntgValueKin[:nslice, 1]
    HeaderInt[2*i] = 'Time (s)'
    HeaderInt[2*i+1] = fileID[1] + '_DoubleIntegral'
del fileID, Filename, folder, i, nfile, nslice
