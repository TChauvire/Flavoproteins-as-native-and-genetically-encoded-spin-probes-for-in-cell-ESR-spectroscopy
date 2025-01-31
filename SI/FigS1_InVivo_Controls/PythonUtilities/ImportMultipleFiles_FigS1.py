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
import numpy as np
from scipy.integrate import cumulative_trapezoid


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


def OpenDoubleIntegral(ListOfFiles, Scaling=None, polyorder=[0, 0],
                       window=20, *args, **kwargs):
    maxlen = MaxLengthOfFiles(ListOfFiles, *args, **kwargs)
    ncol = 4
    FirstDeriv = np.full((maxlen, 4*len(ListOfFiles)), np.nan)
    ZeroDeriv = np.full((maxlen, 4*len(ListOfFiles)), np.nan)
    Header = list(np.zeros((ncol*len(ListOfFiles),)))
    HeaderInt = list(np.zeros((len(ListOfFiles),)))
    IntgValue = list(np.zeros((len(ListOfFiles),)))
    bckg, _, _ = eprload(FileName=ListOfFiles[0], Scaling=Scaling)

    for file in ListOfFiles:
        i = ListOfFiles.index(file)
        data2, abscissa2, par = eprload(FileName=file, Scaling=Scaling)
        npts = abscissa2.shape[0]
        ninit = 0  # int(np.floor(npts/4))
        data = data2[0+ninit:npts-ninit]
        abscissa = abscissa2[0+ninit:npts-ninit]
        if (data.shape[0] != np.ravel(data).shape[0]):
            raise ValueError(
                'The file {0} is\'t a column vector'.format(par['TITL']))
        else:
            # data = np.ravel(data) - np.ravel(bckg) # line to adjust if you
            # want to automatically subtract the background to the other files
            data, _, _, _ = basecorr1D(
                x=abscissa, y=data, polyorder=polyorder[0], window=window)
            npts = abscissa.shape[0]
            FirstDeriv[0:npts, 4 *
                       i] = abscissa[0:npts, 0]
            FirstDeriv[0:npts, 4*i+1] = data[0:npts,]
            newdata = datasmooth(data[0:npts], window_length=4, method='binom')
            FirstDeriv[0:npts, 4*i +
                       2] = newdata[0:npts,]
            FirstDeriv[0:npts, 4*i+3] = FirstDeriv[0:npts, 4*i+2] / \
                np.max(FirstDeriv[0:npts, 4*i+2])
            Header[ncol*i] = par['XNAM']
            Header[ncol*i+1] = par['TITL']
            Header[ncol*i+2] = par['TITL']+str("_4ptsSmoothed_BkgdCorrected")
            Header[ncol*i+3] = par['TITL'] + \
                str("_4ptsSmoothed_BkgdCorrected_Normalized")
            # achieve first integral of the data
            newdata2 = cumulative_trapezoid(data.ravel(),
                                            abscissa.ravel(), initial=None)
            npts2 = newdata2.shape[0]
            newdata2, _, _, _ = basecorr1D(
                x=abscissa[0:npts2,].ravel(), y=newdata2[:npts2,].ravel(),
                polyorder=polyorder[1], window=window)
            ZeroDeriv[0:npts2, 4 *
                      i] = abscissa[0:npts2, 0]
            ZeroDeriv[0:npts2, 4*i +
                      1] = newdata2[0:npts2,]
            newdata2 = datasmooth(
                newdata2[0:npts2,], window_length=4, method='binom')
            ZeroDeriv[0:npts2, 4*i +
                      2] = newdata2[0:npts2,]
            ZeroDeriv[0:npts2, 4*i+3] = newdata2[0:npts2,] / \
                np.max(newdata2[0:npts2,])
            # integrate the zero derivative for half of the spectrum
            init = int(np.floor(npts2/8))
            end = int(np.floor(7*npts2/8))
            IntgValue[i] = np.trapz(ZeroDeriv[init:end, 4*i+1].ravel(),
                                    abscissa[init+1:end+1].ravel())
            HeaderInt[i] = par['TITL']
    return FirstDeriv, ZeroDeriv, Header, IntgValue, HeaderInt
