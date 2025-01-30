# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 17:48:00 2023

@author: tim_t
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
from DeerTools import importasciiSVDData, importasciiTimeData
from DeerTools import ImportMultipleFiles
from os import path, makedirs, getcwd

folder1 = getcwd() + '\\Deer\\Pr\\'
folder2 = getcwd() + '\\Deer\\Data\\'
Listoffilename = ImportMultipleFiles(folder1, 'txt')

# Initiate Variable:
timestackedvalue = [0, 0.4, 0.1, 0.2]
diststackedvalue = [0, 1.8, 0.65, 1.3]
color = ['b', 'k', 'g', 'r']
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


plt.close('all')
# This first loop generate a stacked plot of the time domain of the  4pDeer
fig, axes = plt.subplots(1, 1)
for i in range(len(Listoffilename)):
    filename = Listoffilename[i]
    fileID = path.split(filename)
    # Get Time Domain data and Time domain reconstructed Data
    TITL1 = str(fileID[0][:-2]+'Reconstructed\\Reconstructed_'
                + fileID[1][3:])
    y1, yfit1, t1, _ = importasciiTimeData(TITL1)
    t1 = (t1)*1000  # -t1[0]Time domain in nanoseconds
    npts = y1.shape[0]
    y1norm = y1/np.max(y1)
    yfit1_norm = yfit1/np.max(yfit1)
    # Export Data in a single table for a simple copy and paste toward a
    # plotter. (The following 10 lines could be taking out so...)
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
    font = {'family': 'tahoma',
            'weight': 'normal',
            'size': 18}
    plt.rc('grid', linestyle="--", color='grey')
    plt.rc('font', **font)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    nameline1, = axes.plot(t1/1000, y1+timestackedvalue[i], color=color[i],
                           label='Signal', linewidth=2.0)
    nameline2, = axes.plot(t1/1000, yfit1+timestackedvalue[i], color='sienna',
                           linestyle='dashed', linewidth=3,
                           label='SVD Reconstructed Signal')
    axes.set_xlabel("Time [us]", fontsize=20, fontweight='bold')
    axes.set_ylabel("Offset Amplitude S(t)", fontsize=20, fontweight='bold')
    # axes.legend(loc='upper right')
    axes.set_xlim(-0.01, 4.85)
    fig.suptitle(fileID[1], fontsize=20, fontweight='bold')
    axes.grid(visible=True, which='major')
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
fig.set_dpi(300)
folder3 = folder2[:-6] + '\\Figures\\'
makedirs(folder3, exist_ok=True)
size = [12,  8]
# plt.tight_layout()
plt.gcf().set_size_inches(size[0], size[1])
# Save automatically the figures associated to each files for record
plt.savefig(folder3 + 'figure_Time.pdf', transparent=True)
plt.savefig(folder3 + 'figure_Time.tiff', transparent=True)
# This second loop generate a stacked plot of the distance domain of the 4pDeer
fig2, axes2 = plt.subplots(1, 1)
for i in range(len(Listoffilename)):
    filename = Listoffilename[i]
    fileID = path.split(filename)

    # # Get SVD Distance Domain spectra:
    P1, r1, header, P1_down, P1_up = importasciiSVDData(
        Listoffilename[i], uncertainty=True)
    intgValue = np.trapz(P1.ravel(), r1.ravel())
    # # Export the data:
    npts = P1.shape[0]
    # Export Data in a single table for a simple copy and paste toward a
    # plotter. (The following 10 lines could be taking out so...)
    ExportDistData[0:npts, 5*i] = r1[0:npts,]
    ExportDistData[0:npts, 5*i+1] = P1[0:npts,]
    ExportDistData[0:npts, 5*i+2] = P1_down[0:npts,]
    ExportDistData[0:npts, 5*i+3] = P1_up[0:npts,]
    ExportDistData[0:npts, 5*i+4] = P1[0:npts,]/np.max(P1[0:npts,])
    HeaderDist[5*i] = 'Distance (nm)'
    HeaderDist[5*i+1] = fileID[1]+'_SVDMethod_Pr'
    HeaderDist[5*i+2] = fileID[1]+'_SVDMethod_Prdown'
    HeaderDist[5*i+3] = fileID[1]+'_SVDMethod_Prup'
    HeaderDist[5*i+4] = fileID[1]+'_SVDMethod_Norm'

    # Plot the data : Time Domain // Distance Domain // Pake Pattern
    font = {'family': 'tahoma',
            'weight': 'normal',
            'size': 18}
    plt.rc('grid', linestyle="--", color='grey')
    plt.rc('font', **font)
    plt.rc('xtick', labelsize=18)
    plt.rc('ytick', labelsize=18)
    fig2.suptitle(fileID[1], fontsize=20, fontweight='bold')

    # PR plot:
    axes2.fill_between(x=r1, y1=(P1_down/intgValue+diststackedvalue[i]),
                       y2=(P1_up/intgValue+diststackedvalue[i]),
                       color="red", alpha=0.5, lw=0.1, interpolate=True)
    l4, = axes2.plot(r1, (P1/intgValue+diststackedvalue[i]), color=color[i],
                     label='SF-SVD Method',
                     linewidth=3)
    axes2.set_yticks([0.])
    axes2.set_xlabel("Distance [nm]", fontsize=20, fontweight='bold')
    axes2.set_ylabel("P(r)", fontsize=20, fontweight='bold')
    plt.locator_params(axis='y', nbins=6)
    axes2.set_xlim(2, 9)
    plt.locator_params(axis='x', nbins=8)
    axes2.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    axes2.grid(visible=True, which='major')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
fig.set_dpi(300)
folder3 = folder2[:-6] + '\\Figures\\'
makedirs(folder3, exist_ok=True)
size = [12, 8]
# plt.tight_layout()
plt.gcf().set_size_inches(size[0], size[1])
# Save automatically the figures associated to each files for record
plt.savefig(folder3 + 'figure_Dist.pdf', transparent=True)
plt.savefig(folder3 + 'figure_Dist.tiff', transparent=True)

del axes, fig, figManager, fileID, filename, filenametest, folder1, folder2
del font, header, i, axes2, color, diststackedvalue, fig2, folder3, l4, ncol
del ncol2, npts, P1, P1_down, P1_up, r1, size, t1, TITL1, yfit, y, y1, y1norm
del yfit1, yfit1_norm, t, intgValue, nameline1, nameline2, timestackedvalue
