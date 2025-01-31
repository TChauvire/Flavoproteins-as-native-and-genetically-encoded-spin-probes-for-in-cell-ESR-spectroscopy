# -*- coding: utf-8 -*-
"""
1) Import 4P DEER data from QBand
2) Plot the data and export the 4PDEER param table
3) import the Peter Borbat 4PDeer data and achieve bckg subtraction before
SVD treatment.

@author: tim_t
"""
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from ImportMultipleFiles import eprload
from os import path, makedirs, getcwd
from DeerTools import BckgndSubtractionOfRealPart, ExportOneAscii
from DeerTools import ExportParamTable, GetPakePattern, ImportMultipleFiles

plt.close('all')
# File location in folder
folder1 = getcwd() + '\\Deer\\Data\\'
# Importing global path to files
ListOfFiles1 = ImportMultipleFiles(
    FolderPath=folder1, Extension='DSC')

# Initialize dictionnaries
DataDict1, ParamDict1 = {}, {}

for j in range(len(ListOfFiles1)):
    d, p, d2, p2 = {}, {}, {}, {}
    filename = ListOfFiles1[j]
    # Import data and achieve baseline subtraction
    # Two parameters can be ajusted: 1) the truncation at the end of the time
    # domain is defined by the ratio Percent_tmax. Percent_tmax of the data
    # is kept, and at the end of the time domain, 1-Percent_tmax is truncated.
    # 2) truncsize defined the area the background fit is done. It's start from
    # the end of the time domain. So for a 256pts time domain, with a truncsize
    # of 1/2, the fitting range of the background is [128-256pts] in the time
    # domain. With a truncsize of 3/4, the fitting range would be [64-256].
    d, p = BckgndSubtractionOfRealPart(
        filename, d, p, Dimensionorder=1, Percent_tmax=79/80,
        mode='polyexp', truncsize=3/4)
    DataDict1.update(deepcopy(d))
    ParamDict1.update(deepcopy(p))
    ExportOneAscii(filename, d, p, mode='time')

# Comparison of the data before and after background subtraction
for j in range(len(ListOfFiles1)):
    filename = ListOfFiles1[j]
    _, x, par = eprload(filename, Scaling=None)
    npts = x.shape[0]
    fileID = path.split(filename)
    TITL = fileID[1]
    fulldata1 = DataDict1.get(fileID[1])
    Exp_parameter1 = ParamDict1.get(fileID[1])
    itmin = Exp_parameter1.get('itmin')
    itmax = Exp_parameter1.get('itmax')
    time1 = fulldata1[:itmax, 0].ravel()/1000
    x1 = x[:itmax,].ravel()/1000
    y1 = fulldata1[:itmax, 5].ravel()
    plt.figure(j)
    fig, axes = plt.subplots(2, 1)
    font = {'family': 'tahoma',
            'weight': 'normal',
            'size': 16}
    plt.rc('grid', linestyle="--", color='grey')
    plt.rc('font', **font)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    for i, ax in enumerate(axes):
        axes[i].grid()
    plt.suptitle(fileID[1], fontsize=18, fontweight='bold')
    y = DataDict1.get(TITL)[0:itmax, 3].ravel()
    bckg = DataDict1.get(TITL)[0:itmax, 4].ravel()
    l1 = axes[0].plot(x1, y, 'k', linewidth=2, label='Raw data')
    l2 = axes[0].plot(x1, bckg, 'r', linewidth=2, label='Bckg')
    axes[0].set_xlabel('Time [us]', fontsize=16, fontweight='bold')
    axes[0].set_ylabel('S(t) [a.u.]', fontsize=16, fontweight='bold')
    axes[0].legend()
    l1 = axes[1].plot(time1, y1, 'k', linewidth=2,
                      label='Bckg subtracted data')
    axes[1].set_xlabel('Time [us]', fontsize=16, fontweight='bold')
    axes[1].set_ylabel('S(t) [a.u.]', fontsize=16, fontweight='bold')
    axes[1].legend()
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.tight_layout()
    fig.set_dpi(200)
    makedirs(folder1 + '\\Figures\\', exist_ok=True)
    size = [19.2,  9.51]
    plt.gcf().set_size_inches(size[0], size[1])
    plt.savefig(folder1 + '\\Figures\\figure{0}'.format(j))
# Export 4PDEER parameter tables
ExportParamTable(FolderPath=folder1, ParamDict=ParamDict1)

# Comparison of the data and its Pakepattern
for i in range(len(ListOfFiles1)):
    filename = ListOfFiles1[i]
    fileID = path.split(filename)
    _, x, par = eprload(filename, Scaling=None)
    npts = x.shape[0]
    TITL = fileID[1]
    Exp_parameter1 = ParamDict1.get(TITL)
    itmax1 = Exp_parameter1.get('itmax')
    x1 = DataDict1.get(TITL)[0:npts, 0].ravel()
    y1 = DataDict1.get(TITL)[0:npts, 2].ravel()
    bckg1 = DataDict1.get(TITL)[0:npts, 4].ravel()
    ynew1 = DataDict1.get(TITL)[0:npts, 5].ravel()
    ModDepth = Exp_parameter1.get('ModDepth')
    tmin1 = Exp_parameter1.get('zerotime')
    tmax1 = Exp_parameter1.get('tmax')
    itmin1 = Exp_parameter1.get('itmin')
    ynew1 = np.real(ynew1[itmin1:itmax1,])
    plt.figure(j+i+2)
    PakeSpectraPhased1, PakeAbs1, freq1 = GetPakePattern(
        x[itmin1:itmax1,], ynew1)
    l1, = plt.plot(freq1, PakeAbs1, 'k', linewidth=2,
                   label='Pake Pattern Data')
    plt.title(fileID[1], fontsize=18, fontweight='bold')
    plt.grid()
    plt.xlabel('Frequency [MHz]', fontsize=16, fontweight='bold')
    plt.ylabel('F(t) [a.u]', fontsize=16, fontweight='bold')
    plt.title(fileID[1], fontsize=18, fontweight='bold')
    plt.legend()
    plt.grid()
    plt.xlabel('Frequency [MHz]', fontsize=16, fontweight='bold')
    plt.ylabel('F(t) [a.u]', fontsize=16, fontweight='bold')
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.tight_layout()
    fig.set_dpi(200)
    makedirs(folder1 + '\\Figures\\', exist_ok=True)
    size = [19.2,  9.51]
    plt.gcf().set_size_inches(size[0], size[1])
    plt.savefig(folder1 + '\\Figures\\figure_fit{0}'.format(i+j+2))
