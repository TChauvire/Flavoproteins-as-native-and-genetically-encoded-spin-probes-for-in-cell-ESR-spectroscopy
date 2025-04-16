# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 17:48:00 2023

@author: tim_t
"""

import numpy as np
import matplotlib.pyplot as plt
from DeerTools import importasciiSVDData, importasciiTimeData
from DeerTools import ImportMultipleFiles
from os import path, makedirs, getcwd

folder1 = getcwd() + '\\Deer\\Pr\\'
folder2 = getcwd() + '\\Deer\\Data\\'
Listoffilename = ImportMultipleFiles(folder1, 'txt')
# Initiate Variable:
filenametest = Listoffilename[0]
fileID = path.split(filenametest)
TITL1 = str(fileID[0][:-2]+'Reconstructed\\Reconstructed_'
            + fileID[1][3:])
y, yfit, t, _ = importasciiTimeData(TITL1)
npts = y.shape[0]
ncol = 5*len(Listoffilename)
ncol2 = 5*len(Listoffilename)
ExportDistData = np.full((5*npts, ncol), np.nan, dtype=float)
ExportTimeData = np.full((5*npts, ncol2), np.nan, dtype=float)
HeaderDist = list(np.zeros(ncol,))
HeaderTime = list(np.zeros(ncol2,))
# Gaussianreport = np.full((10, ncol2), np.nan, dtype=float)

plt.close('all')
for i in range(len(Listoffilename)):
    filename = Listoffilename[i]
    fileID = path.split(filename)
    # Get Time Domain data and Time domain reconstructed Data
    TITL1 = str(fileID[0][:-2]+'Reconstructed\\Reconstructed_'
                + fileID[1][3:])
    y1, yfit1, t1, _ = importasciiTimeData(TITL1)
    t1 = (t1)*1000  # Time domain in nanoseconds
    npts = y1.shape[0]
    y1norm = y1/np.max(y1)
    yfit1_norm = yfit1/np.max(yfit1)
    # Export the time domain
    ExportTimeData[0:npts, 5*i] = t1[0:npts,]
    ExportTimeData[0:npts, 5*i+1] = y1[0:npts,]
    ExportTimeData[0:npts, 5*i+2] = y1norm[0:npts,]
    ExportTimeData[0:npts, 5*i+3] = yfit1[0:npts,]
    ExportTimeData[0:npts, 5*i+4] = yfit1_norm[0:npts,]
    HeaderTime[5*i] = 'Time (ns)'
    HeaderTime[5*i+1] = fileID[1][3:]+'_S(t)'
    HeaderTime[5*i+2] = fileID[1][3:]+'_S(t)_Norm'
    HeaderTime[5*i+3] = fileID[1][3:]+'_SVDFit'
    HeaderTime[5*i+4] = fileID[1][3:]+'_SVDFit_Norm'

    # Get SVD Distance Domain spectra:
    P1, r1, header, P1_down, P1_up = importasciiSVDData(
        Listoffilename[i], uncertainty=True)
    # Get the factor to normalize the distance domain (Total probability =1)
    intgValue = np.trapz(P1.ravel(), r1.ravel())
    # # Export the data:
    npts = P1.shape[0]
    ExportDistData[0:npts, 5*i] = r1[0:npts,]
    ExportDistData[0:npts, 5*i+1] = P1[0:npts,]/intgValue
    ExportDistData[0:npts, 5*i+2] = P1_down[0:npts,]/intgValue
    ExportDistData[0:npts, 5*i+3] = P1_up[0:npts,]/intgValue
    ExportDistData[0:npts, 5*i+4] = P1[0:npts,]/np.max(P1[0:npts,])

    HeaderDist[5*i] = 'Distance (nm)'
    HeaderDist[5*i+1] = fileID[1]+'_SVDMethod_Pr'
    HeaderDist[5*i+2] = fileID[1]+'_SVDMethod_Prdown'
    HeaderDist[5*i+3] = fileID[1]+'_SVDMethod_Prup'
    HeaderDist[5*i+4] = fileID[1]+'_SVDMethod_Norm'

    # Plot the data : Time Domain // Distance Domain // Pake Pattern
    plt.figure(i)
    fig, axes = plt.subplots(1, 2)
    font = {'family': 'tahoma',
            'weight': 'normal',
            'size': 16}
    plt.rc('grid', linestyle="--", color='grey')
    plt.rc('font', **font)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    for j, ax in enumerate(axes):
        axes[j].grid()
    fig.suptitle(fileID[1], fontsize=20, fontweight='bold')
    l1, = axes[0].plot(t1, y1, 'k', label='Signal', linewidth=1.5)
    l2, = axes[0].plot(t1, yfit1, 'r', label='SVD Reconstructed Signal',
                       linewidth=1.5)
    axes[0].set_xlabel("Time Domain [ns]", fontsize=18, fontweight='bold')
    axes[0].set_ylabel("S(t)", fontsize=18, fontweight='bold')
    axes[0].legend(loc='upper right')
    # PR plot:
    axes[1].fill_between(x=r1, y1=P1_down/intgValue, y2=P1_up/intgValue,
                         color="red", alpha=1, lw=0.1, interpolate=True)
    l4, = axes[1].plot(r1, P1/intgValue, 'k', label='SF-SVD Method',
                       linewidth=1.5)
    axes[1].set_xlabel("Distance Domain [nm]", fontsize=18, fontweight='bold')
    axes[1].set_ylabel("Normalized P(r)", fontsize=18, fontweight='bold')
    axes[1].legend(loc='upper right')
    axes[1].set_xlim(2.0, 8.0)

    #
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    fig.set_dpi(300)
    makedirs(folder2[:-6] + '\\Figures\\', exist_ok=True)
    size = [19.2,  9.51]
    # plt.tight_layout()
    plt.gcf().set_size_inches(size[0], size[1])
    plt.savefig(folder2[:-6] + '\\Figures\\figure_DistNorm{0}'.format(i))

del ax, axes, fig, figManager, fileID, filename, filenametest, folder1, folder2
del font, header, i, j, l1, l2
del l4, ncol, ncol2, npts, P1, P1_down, P1_up, r1, size, t1, TITL1, yfit
del y, y1, y1norm, yfit1, yfit1_norm, t
