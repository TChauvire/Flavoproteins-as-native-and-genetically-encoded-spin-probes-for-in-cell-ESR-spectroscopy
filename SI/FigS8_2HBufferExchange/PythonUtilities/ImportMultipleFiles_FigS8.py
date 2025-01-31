# -*- coding: utf-8 -*-
"""
Importing Multiple files
Return a list of name and foldername
Created on Mon Jul  6 11:57:42 2020

@author: TC229401
"""
from os import walk, getcwd, path
from datasmooth import datasmooth
from eprload_BrukerBES3T import eprload_BrukerBES3T
from eprload_BrukerESP import eprload_BrukerESP
from basecorr1D import basecorr1D
from automatic_phase import automatic_phase
import numpy as np


def eprload(FileName=None, Scaling=None, *args, **kwargs):
    if FileName[-4:].upper() in ['.DSC', '.DTA']:
        data, abscissa, par = eprload_BrukerBES3T(FileName, Scaling)
    elif FileName[-4:].lower() in ['.spc', '.par']:
        data, abscissa, par = eprload_BrukerESP(FileName, Scaling)
    else:
        data, abscissa, par = None, None, None
        raise ValueError("Can\'t Open the File {0} ".format(str(FileName)) +
                         "because the extension isn\'t a Bruker extension " +
                         ".DSC, .DTA, .spc, or .par!")
    return data, abscissa, par


def ImportMultipleNameFiles(FolderPath=getcwd(), Extension=None,
                            *args, **kwargs):
    ListOfFiles = []
    for root, dirs, files in walk(FolderPath):
        for file in files:
            if file.endswith(Extension):
                ListOfFiles.append(path.normpath(
                    path.join(root, file)))
    return ListOfFiles


def MaxLengthOfFiles(ListOfFiles, *args, **kwargs):
    maxlen = 0
    for file in ListOfFiles:
        data, abscissa, par = eprload(file, Scaling=None)
        if maxlen < data.shape[0]:
            maxlen = data.shape[0]
    return maxlen


def OpenMultipleFiles(ListOfFiles, Scaling=None, polyorder=0,
                      window=20, *args, **kwargs):
    maxlen = MaxLengthOfFiles(ListOfFiles, *args, **kwargs)
    ncol = 4
    fulldata = np.full((maxlen, 4*len(ListOfFiles)), np.nan)
    Header = list(np.zeros((ncol*len(ListOfFiles),)))
    for file in ListOfFiles:
        i = ListOfFiles.index(file)
        data, abscissa, par = eprload(FileName=file, Scaling=Scaling)
        if (data.shape[0] != np.ravel(data).shape[0]):
            raise ValueError(
                'The file {0} is\'t a column vector'.format(par['TITL']))
        else:
            data = np.ravel(data)
            data, _, _, _ = basecorr1D(
                x=abscissa, y=data, polyorder=polyorder, window=window)
            npts = abscissa.shape[0]
            fulldata[0:npts, 4*i] = abscissa[0:npts, 0]
            fulldata[0:npts, 4*i+1] = data[0:npts]
            newdata = datasmooth(data[0:npts], window_length=1, method='binom')
            fulldata[0:npts, 4*i+2], _, _, _ = basecorr1D(
                x=abscissa, y=newdata, polyorder=polyorder, window=window)
            fulldata[0:npts, 4*i+3] = fulldata[0:npts, 4*i+2] / \
                np.max(fulldata[0:npts, 4*i+2])
            Header[ncol*i] = par['XNAM']
            Header[ncol*i+1] = par['TITL']
            Header[ncol*i+2] = par['TITL']+str("_1ptsSmoothed_BkgdCorrected")
            Header[ncol*i+3] = par['TITL'] + \
                str("_4ptsSmoothed_BkgdCorrected_Normalized")
    return fulldata, Header


def OpenDavies(ListOfFiles, Scaling=None, polyorder=0, window=20,
               *args, **kwargs):
    maxlen = MaxLengthOfFiles(ListOfFiles, *args, **kwargs)
    ncol = 7
    fulldata = np.full((maxlen, ncol*len(ListOfFiles)),
                       np.nan, dtype=float)
    Header = list(np.zeros((ncol*len(ListOfFiles),)))
    for file in ListOfFiles:
        i = ListOfFiles.index(file)
        data, x, par = eprload(FileName=file, Scaling=None)
        if (data.shape[0] != np.ravel(data).shape[0]):
            raise ValueError(
                'The file {0} is\'t a column vector'.format(par['TITL']))
        else:
            data = np.ravel(data)
            # First Phase the data
            new_data, _ = automatic_phase(vector=data, pivot1=int(
                data.shape[0]/2), funcmodel='minfunc')
            data_real = new_data.real
            data_imag = new_data.imag
            data_real_new, _, _, _ = basecorr1D(
                x=x, y=-data_real, polyorder=polyorder, window=window)
            npts = x.shape[0]
            fulldata[0:npts, ncol*i] = x[0:npts, 0]
            fulldata[0:npts, ncol*i+1] = -data_real[0:npts]
            fulldata[0:npts, ncol*i+2] = data_imag[0:npts]
            fulldata[0:npts, ncol*i+3] = np.absolute(-data_real+1j*data_imag)
            fulldata[0:npts, ncol*i+4] = data_real_new[0:npts]
            fulldata[0:npts, ncol*i+5] = datasmooth(data_real_new[0:npts],
                                                    window_length=4,
                                                    method='flat')
            d2 = fulldata[0:npts, ncol*i+5]
            fulldata[0:npts, ncol*i+6] = (d2) / np.amax(d2)
            Header[ncol*i] = par['XNAM']
            Header[ncol*i+1] = par['TITL']+str("_real")
            Header[ncol*i+2] = par['TITL']+str("_imag")
            Header[ncol*i+3] = par['TITL']+str("_abs")
            Header[ncol*i+4] = par['TITL']+str("_real_bckg")
            Header[ncol*i+5] = par['TITL']+str("_real_4ptsSmoothed")
            Header[ncol*i+6] = par['TITL']+str("_real_Normalized")
    return fulldata, Header
