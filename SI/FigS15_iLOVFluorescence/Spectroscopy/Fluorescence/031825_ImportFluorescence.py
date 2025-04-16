# -*- coding: utf-8 -*-
"""
Script to import multiple files from Fluorescence experiments.
Raw data were recorded with a Fluorimeter Fluoromax from Jobin-Yvon
(ref)
This spectrometer is using origin software so every spectrum were converted 
to ascci format .dat
First column is wavelength, second column is intensity (cps)
Created on Tue Mar 18 13:29:05 2025

@author: tim_t
"""

from ImportMultipleFiles import ImportMultipleNameFiles
import matplotlib.pyplot as plt
import numpy as np
from os import path, makedirs, getcwd
import numpy.polynomial.polynomial as pl
folder = getcwd()
ListOfFiles = ImportMultipleNameFiles(folder, Extension='.dat')
# Check for the maximum length datafiles
for file in ListOfFiles:
    list1 = []
    maxlen = 0
    list1 = np.loadtxt(file, dtype=float, delimiter='\t')
    if maxlen < len(list1):
        maxlen = len(list1)

# Initial the export variable and tjhe header
fulldata = np.full((maxlen, len(ListOfFiles)*4), np.nan, dtype=float)
header = [None] * len(ListOfFiles)*4
HeaderInt = list(np.zeros((len(ListOfFiles),)))
IntgValue = list(np.zeros((len(ListOfFiles),)))
# import and plot the data
for i in range(len(ListOfFiles)):
    fileID = path.split(ListOfFiles[i])
    # Import the data
    x, y = [], []
    data = np.loadtxt(ListOfFiles[i], dtype=float, delimiter='\t')
    x = data[:, 0]
    y = data[:, 1]
    npts = len(x)
    fulldata[0:npts, 4*i] = x
    fulldata[0:npts, 4*i+1] = y
    # subtract line
    i1 = int(np.argwhere(x == 650))
    i2 = int(np.argwhere(x == 700))
    c, stats = pl.polyfit(x[i1:i2], y[i1:i2], deg=0, full=True)
    ypoly = pl.polyval(x, c)
    ycorr = y-ypoly
    fulldata[0:npts, 4*i+2] = ycorr
    fulldata[0:npts, 4*i+3] = ycorr/np.max(ycorr)
    # Assign header of the data
    header[4*i] = 'wavelength [nm]'
    header[4*i+1] = fileID[1]
    header[4*i+2] = str(fileID[1] + '_Bckgcorrected')
    header[4*i+3] = str(fileID[1] + '_normalized')
    # Plot the data in a single figure
    plt.figure(i)
    fig, axes = plt.subplots(1, 1)
    fig.suptitle('Fluorescence_' + fileID[1])
    l1, = axes.plot(x.ravel(), y.ravel(), 'k', linewidth=2)
    axes.set_xlabel("Wavelength [nm]", fontsize=14, weight='bold')
    axes.set_ylabel("Intensity [cps]", fontsize=14, weight='bold')
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.subplots_adjust(left=0.046, right=0.976, top=0.945, bottom=0.085,
                        hspace=0.2, wspace=0.2)
    plt.tight_layout()
    makedirs(fileID[0] + '\\Figures\\', exist_ok=True)
    plt.savefig(fileID[0] + '\\Figures\\figure{0}'.format(i+1))
    # Extract the integral value of the data
    IntgValue[i] = np.trapz(y.ravel(), x.ravel())
    HeaderInt[i] = fileID[1]

plt.close('all')
