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
from ImportMultipleFiles import eprload
from os import path, makedirs, getcwd
from DeerTools import GetPakePattern, ImportMultipleFiles
from DeerTools import ComputeRMSE, ExponentialCorr1D_DEER
from automatic_phase import automatic_phase
import csv
# from SVD_scripts import get_KUsV, process_data, rminrmaxcalculus
# import deerlab as dl

plt.close('all')
# Get the files recorded in 2D format (Tau averaging)
folder1 = getcwd() + '\\Deer\\RawData\\2D\\'
ListOfFiles1 = ImportMultipleFiles(FolderPath=folder1, Extension='DSC')


def ExportDeerParam(x, newy, yfit, Percent_tmax, p0):
    itmin = np.abs(x).argmin()
    RelBckgndDecay = 1 - yfit[-1] / yfit[itmin]
    npts = x.shape[0]
    itmax = int(np.floor(Percent_tmax*npts)-1)
    tmax = float(x[itmax,])
    tmin = float(x[itmin,])
    center = int(np.floor(npts/2))
    sigma_noise = np.std(newy.imag[itmin:itmax,])/np.max(newy.real)
    sigma_noise_half = np.std(newy.imag[center-int(np.floor(npts/4)):
                                        center+int(np.floor(npts/4))],)/np.max(newy.real)
    # NoiseLevel = (sigma_noise, sigma_noise_half)
    RMSE = ComputeRMSE(newy.real/np.max(newy.real),
                       yfit/np.max(newy.real), p0)
    FinalData = (new_data.real[:itmax,] - yfit[:itmax,]
                 )/np.max(new_data.real[:itmax,])
    Mod_Depth = FinalData[itmin-1:itmin+1].mean()
    DEER_SNR = (Mod_Depth//sigma_noise, Mod_Depth//sigma_noise_half)
    # Paramlist = ['tmax', 'tmin', 'NoiseLevel', 'ModDepth', 'DEER_SNR',
    # 'RelBckgndDecay', 'RMSE']
    Param = [tmax, tmin, sigma_noise_half, Mod_Depth, DEER_SNR[1],
             RelBckgndDecay, RMSE]
    return Param


ListOfFiles = list()
shape = (7, len(ListOfFiles1))
FinalParam = np.full(shape, np.nan, dtype=float)
# This loop plot both time domain and background subtraction impact, and Pake
# Pattern so the frequency domain after FFT transformation.
for i in range(len(ListOfFiles1)):
    filename = ListOfFiles1[i]
    fileID = path.split(filename)
    ListOfFiles += [fileID[1]]
    # Import data and achieve baseline subtraction
    y, x, par = eprload(filename, Scaling=None)

    newy = np.mean(y[:, :], axis=1)
    pivot = int(np.floor(newy.shape[0]/2))
    new_data, _ = automatic_phase(vector=newy, pivot1=pivot,
                                  funcmodel='minfunc')
    newx = x[:, 0]/1000
    npts = x.shape[0]
    Percent_tmax = 78/80
    itmax = int(np.floor(Percent_tmax*npts))
    yfit, p0, _ = ExponentialCorr1D_DEER(x=newx, y=np.real(new_data),
                                         Dimensionorder=2,
                                         Percent_tmax=Percent_tmax,
                                         mode='polyexp', truncsize=8/8)
    newx = x[:itmax, 0]/1000
    Bckg_type = str('polyexp2_1.0')
    FinalData = (new_data.real[:itmax,] - yfit[:itmax,]
                 )/np.max(new_data.real[:itmax,])
    ynew = FinalData-FinalData[-20:].mean()
    Param = ExportDeerParam(newx*1000, newy, yfit, Percent_tmax, p0)
    FinalParam[0:7, i] = Param
    plt.figure(2*i)
    fig, axes = plt.subplots(2, 1)
    font = {'family': 'tahoma',
            'weight': 'normal',
            'size': 16}
    plt.rc('grid', linestyle="--", color='grey')
    plt.rc('font', **font)
    plt.rc('xtick', labelsize=16)
    plt.rc('ytick', labelsize=16)
    for j, ax in enumerate(axes):
        axes[j].grid()
    plt.suptitle(fileID[1], fontsize=18, fontweight='bold')
    l1 = axes[0].plot(newx, new_data.real[:itmax,], 'k',
                      linewidth=2, label='Raw data')
    l2 = axes[0].plot(newx, yfit[:itmax,], 'r', linewidth=2, label='Bckg')
    axes[0].set_xlabel('Time [us]', fontsize=16, fontweight='bold')
    axes[0].set_ylabel('S(t) [a.u.]', fontsize=16, fontweight='bold')
    axes[0].legend()
    l1 = axes[1].plot(newx, ynew, 'k', linewidth=2,
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
    plt.savefig(folder1 + '\\Figures\\figure{0}'.format(2*i))
    # Let's get the Pake pattern graph now:
    plt.figure(2*i+1)
    fig, axes = plt.subplots(1, 1)
    PakeSpectraPhased1, PakeAbs1, freq1 = GetPakePattern(
        newx[6:,]*1000, ynew[6:,])
    l1, = axes.plot(freq1, PakeAbs1, 'k', linewidth=2,
                    label='Pake Pattern Data')
    plt.xlabel('Frequency [MHz]', fontsize=16, fontweight='bold')
    plt.ylabel('F(t) [a.u]', fontsize=16, fontweight='bold')
    plt.title(fileID[1], fontsize=18, fontweight='bold')
    plt.legend()
    plt.grid()
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    plt.tight_layout()
    fig.set_dpi(200)
    makedirs(folder1 + '\\Figures\\', exist_ok=True)
    size = [19.2,  9.51]
    plt.gcf().set_size_inches(size[0], size[1])
    plt.savefig(folder1 + '\\Figures\\figure_Pake{0}'.format(2*i+1))

    # Export the backgnd subtracted data:
    fullnewfilename = str(filename+'.txt')
    header = '\n'.join(['Filename: ' + fileID[1],
                        'Bckg_type: ' + Bckg_type,
                        'FirstColumn: ' + 'Time[us]',
                        'SecondColumn: ' + 'BckgSubtractedData'])
    data = np.column_stack([newx, ynew])
    makedirs(path.dirname(fullnewfilename), exist_ok=True)
    np.savetxt(fullnewfilename, data, fmt=['%.5e', '%.15e'],
               delimiter='\t', header=header)

fullfilename = path.normpath(path.join(folder1, 'DeerParam.csv'))
with open(fullfilename, 'w') as file:
    wr = csv.writer(file, dialect='excel',
                    delimiter='\t', lineterminator='\n')
    wr.writerow(ListOfFiles)
    for row in FinalParam:
        wr.writerow(row)
