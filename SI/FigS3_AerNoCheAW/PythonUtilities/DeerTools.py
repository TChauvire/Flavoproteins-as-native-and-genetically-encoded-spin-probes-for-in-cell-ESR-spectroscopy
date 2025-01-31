# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 12:37:19 2023

@author: tim_t

General toolbox scripts for 1) importing DEER Bruker data, 2) autophase the
complex data 3) generate background subtraction, 4) importing result
from the DEER analysis matlab software, 5) analyze distance domain data,
6) generate table for publications, 7) ...

Scripts to import selected information from DEER_Analysis software
(https://epr.ethz.ch/software.html) to construct the table for the publication:
    [ref to insert]
The script rely on the EPR_ESR Suite:
(https://github.com/TChauvire/EPR_ESR_Suite).

# For the data :
# 1) First Column is time
# 2) Second column is the raw data (complex))
# 3) Third column is the phase corrected data (complex))
# 4) Fourth column is the real part of the phase corrected data (real))
# 5) Fifth column are the background data (real)
# 5) Sixth column would be the background subtracted data (real)
# To call a column, you have to use the title of the data to access
TITL = path.split(Filename)[1]
# the global table : datatable = DataDictionnary.get(TITL)

For the parameters :
1)
2) Noise level estimated from the root-mean-square
amplitude of the imaginary part after phase correction.
Parameters

"""

from os import walk, getcwd, path, makedirs
import csv
from automatic_phase import automatic_phase
from basecorr1D import basecorr1D, error_vandermonde
from eprload_BrukerBES3T import eprload_BrukerBES3T
from eprload_BrukerESP import eprload_BrukerESP
from fdaxis import fdaxis
from windowing import windowing
import numpy as np
import numpy.polynomial.polynomial as pl
from scipy.optimize import curve_fit  # 1) used for gaussian analysis of the
# distance domain
# 2) used for the backgnd subtraction in the time domain
# import deerlab as dl (To use in case to determine the zerotime)
#############################################################################
# Initialization of the variables:
DataDict, ParamDict = {}, {}


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


def ImportMultipleFiles(FolderPath=getcwd(), Extension='txt'):
    ListOfFiles = []
    for root, dirs, files in walk(FolderPath):
        for file in files:
            if file.endswith(Extension):
                ListOfFiles.append(path.normpath(path.join(root, file)))
    return ListOfFiles


def ExponentialCorr1D_DEER(x=None, y=None, Dimensionorder=3, Percent_tmax=9/10,
                           mode='strexp', truncsize=3/4, *args, **kwargs):
    '''
    Function that achieve a baseline correction by fitting a function
    parameterized by a streched exponential for a supposed homogeneous
    three-dimensional solution (d=3) or a stretched exponential for other
    dimensions as described by multiple publications.
    See by example : (dx.doi.org/10.1039/c9cp06111h )

    .. math::

    B(t) = \exp\left(-\kappa \vert t\vert^{d}\right)

    k is the decay rate constant of the background and d is
    the so-called fractal dimension of the exponential

    The fitting is done on the last 3/4 points of the data.
    Script written by Timothée Chauviré
    (https://github.com/TChauvire/EPR_ESR_Suite/), 10/18/2023

    Parameters
    ----------
    x : abscissa of the data, TYPE : numpy data array, column vector
        DESCRIPTION. The default is None.
    y : data which baseline has to be corrected,
        TYPE : numpy data array, column vector
        It has to have the same size than x
        DESCRIPTION. The default is None.
    Dimensionorder : order of so-called fractal dimension of the exponential
        TYPE : Integer, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    ynew : baseline data array
        TYPE: numpy data array, same shape as input data y
    (k,d) : coefficient used for the exponential fit
        TYPE : tuple of real values, coefficient of the streched exponential
    perr : error coefficient obtained from the covariance matrix
        (perr = np.sqrt(np.diag(pcov)))
        TYPE : diagonalized 2-D array
    mode='strexp', stretched exponential ackground subtraction
    mode='poly', polynomial fit of the logarithmic data
    '''
    shape = y.shape
    if x.shape[0] != np.ravel(x).shape[0]:
        raise ValueError('x must be a column vector. ExponentialCorr1D_DEER'
                         'function does not work on 2D arrays.')
    else:
        x = np.ravel(x)

    if y.shape[0] != np.ravel(y).shape[0]:
        raise ValueError('y must be a column vector. ExponentialCorr1D_DEER'
                         'function does not work on 2D arrays.')
    else:
        y = np.ravel(y)

    if y.shape[0] != x.shape[0]:
        raise ValueError('x and y must be column vector of the same size.')
    yfit = np.full(y.shape, np.nan)
    npts = x.shape[0]
    npts_new = int(np.floor(npts*truncsize))
    itmax = int(np.floor(Percent_tmax*npts))
    xfitinit = (np.array(x[(npts-npts_new):itmax])).ravel()
    # xfitinit = np.array(x[-npts_new:tmax]).ravel()
    yfit = (y/np.max(y)).ravel().real
    yfitinit = (np.array(yfit[(npts-npts_new):itmax])).ravel()
    # yfitinit = np.array(yfit[-npts_new:tmax]).ravel()

    def strexp(x, ModDepth, decay, stretch):
        a = (1-ModDepth)*(np.exp((-1)*(np.abs(decay*x))) ** (stretch/3))
        return a

    def strexp2(x, ModDepth, decay):
        a = (1-ModDepth)*(np.exp((-1)*(np.abs(decay*x))) **
                          (Dimensionorder/3))
        return a

    # # Add parameters
    p0_1 = [0.3, 0.25, Dimensionorder]
    b_1 = ([0, 0, 0], [1, 200, 6])  # (lowerbound,upperbound)  bounds=b,
    # Add parameters
    p0_2 = [0.3, 0.25]
    b_2 = ([0, 0], [1, 200])  # (lowerbound,upperbound)  bounds=b,

    if mode == 'strexp':
        poptarray, pcov = curve_fit(strexp, xfitinit, yfitinit, p0=p0_1,
                                    sigma=None, absolute_sigma=False,
                                    check_finite=None, bounds=b_1)
        perr = np.sqrt(np.diag(pcov))
        yfit2 = (strexp(x, poptarray[0],
                 poptarray[1], poptarray[2]))*np.max(y)
        yfit2 = yfit2.reshape(shape)
        return yfit2, poptarray, perr
    if mode == 'strexp_fixed':
        poptarray, pcov = curve_fit(strexp2, xfitinit, yfitinit, p0=p0_2,
                                    sigma=None, absolute_sigma=False,
                                    check_finite=None, bounds=b_2)
        perr = np.sqrt(np.diag(pcov))
        yfit2 = (strexp2(x, poptarray[0], poptarray[1]))*np.max(y)
        poptarray = np.append(poptarray, Dimensionorder)
        yfit2 = yfit2.reshape(shape)
        return yfit2, poptarray, perr
    if mode == 'poly':
        c, stats = pl.polyfit(xfitinit, yfitinit,
                              deg=Dimensionorder, full=True)
        ypoly = pl.polyval(x, c)*np.max(y)
        error_parameters, _ = error_vandermonde(x, residuals=stats[0], rank=3)
        return ypoly, c, error_parameters
    if mode == 'polyexp':
        c, stats = pl.polyfit(xfitinit, np.log(yfitinit),
                              deg=Dimensionorder, full=True)
        ypoly = np.exp(pl.polyval(x, c))*np.max(y)
        error_parameters, _ = error_vandermonde(x, residuals=stats[0], rank=3)
        return ypoly, c, error_parameters


def BckgndSubtractionOfRealPart(Filename, DataDict, ParamDict, Scaling=None,
                                Dimensionorder=3, Percent_tmax=9/10,
                                mode='strexp', truncsize=3/4, *args, **kwargs):
    '''

    ----------
    ListOfFiles : TYPE
        DESCRIPTION.
    Scaling : TYPE, optional
        DESCRIPTION. The default is None.
    Dimensionorder : TYPE, optional
        DESCRIPTION. The default is 1.
    Percent_tmax : TYPE, optional
        DESCRIPTION. The default is 9/10.
    *args : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    BckgndSubtractedData : TYPE
        DESCRIPTION.
    Modulation_Depth : TYPE
        DESCRIPTION.
    RelBckgndDecay : TYPE
        DESCRIPTION.
    NoiseLevel : TYPE
        DESCRIPTION.

    '''
    fileID = path.split(Filename)
    data, abscissa, par = eprload(Filename, Scaling)
    if (data.shape[0] != np.ravel(data).shape[0]):
        raise ValueError(
            'The file {0} is\'t a column vector'.format(par['TITL']))
    else:
        npts = abscissa.shape[0]
        new_data = np.full(npts, np.nan, dtype="complex_")
        pivot = int(np.floor(data.shape[0]/2))
        new_data, _ = automatic_phase(vector=data, pivot1=pivot,
                                      funcmodel='minfunc')
        data_real = np.ravel(new_data.real)
        data_imag = np.ravel(new_data.imag)
        abscissa = np.ravel(abscissa)
        # Achieve background correction of the real part :
        # newdata_real = datasmooth(
        #     data_real[0:npts], window_length=10, method='binom')
        # itmin = np.argmax(newdata_real)
        # newx = dl.correctzerotime(data_real, abscissa)
        newx = abscissa
        itmin = np.abs(newx).argmin()
        # itmin = newx([newx==0])
        if itmin > 50:  # The determination of zerotime didn't work, do nothing
            itmin = 0
            newx = abscissa
        tmin = newx[itmin]
        data_bckgnd, p0, perr = ExponentialCorr1D_DEER(x=newx, y=data_real,
                                                       Dimensionorder=Dimensionorder,
                                                       Percent_tmax=Percent_tmax, mode=mode,
                                                       truncsize=truncsize)
        w = int(np.floor(npts/2))
        # Achieve automatic base line correction correction of the imaginary
        # part :
        data_imag_new, _, _, _ = basecorr1D(x=newx, y=data_imag,
                                            polyorder=1, window=w)
        if np.floor(Percent_tmax*npts)-1 <= npts:
            itmax = int(np.floor(Percent_tmax*npts)-1)
        else:
            itmax = npts
        RelBckgndDecay = 1 - data_bckgnd[itmax] / data_bckgnd[itmin]
        FinalData = (data_real - data_bckgnd)/np.max(data_real)
        FinalData = FinalData-FinalData[itmax-20:itmax].mean()
        # BckgndSubtractedData = (data_real - data_bckgnd)/(np.max(data_real - data_bckgnd))
        # Two noises-level are computed :
        # 1) With the full imaginary part "sigma_noise"
        # 2) With half of the imaginary part "sigma_noise_half"
        center = int(np.floor(npts/2))
        sigma_noise = np.std(data_imag_new[itmin:itmax,])/np.max(data_real)
        sigma_noise_half = np.std(data_imag_new[center-int(np.floor(npts/4)):
                                                center+int(np.floor(npts/4))],
                                  )/np.max(data_real)
        NoiseLevel = (sigma_noise, sigma_noise_half)
        # Calculate the Root mean square of error
        RMSE = ComputeRMSE(data_real/np.max(data_real),
                           data_bckgnd/np.max(data_real), p0)
        # Let's create a global dictionnary for storing all the data and the
        # parameters:
        # TITL = str(par['TITL'])
        TITL = fileID[1]
        fulldata = np.full((5*npts, 50), np.nan, dtype="complex_")
        fulldata[0:npts, 0] = newx.ravel()
        fulldata[0:npts, 1] = data.ravel()
        fulldata[0:npts, 2] = new_data.ravel()
        fulldata[0:npts, 3] = new_data.real.ravel()
        fulldata[0:npts, 4] = data_bckgnd.ravel()
        fulldata[0:npts, 5] = FinalData.ravel()

        DataDict.update({TITL: fulldata})
        Header = list(np.zeros((50,)))
        Header[0] = str('Time [ns]')
        Header[1] = str(TITL+"_rawData")
        Header[2] = str(TITL+"_phased")
        Header[3] = str(TITL+"_phased_real")
        Header[4] = str(TITL+"_background_real")
        Header[5] = str(TITL+"_backgroundsubtracted_real")
        HeaderTITL = str(TITL+'_Header')
        DataDict[HeaderTITL] = Header
        Exp_parameter = {'RelBckgndDecay': RelBckgndDecay, 'tmin':  tmin,
                         'NoiseLevel': NoiseLevel, 'tmax': abscissa[itmax],
                         'itmax': itmax, 'itmin': itmin, 'RMSE': RMSE}
        # Assign the modulation depth in the parameters
        if mode == 'strexp':
            Bckg_type = str(mode+'_'+str(Percent_tmax)+'_'+str(truncsize))
            DEER_SNR = (p0[0]/sigma_noise, p0[0]/sigma_noise_half)
            Exp_parameter.update({'ModDepth': p0[0], 'decay': p0[1],
                                  'stretch': p0[2], 'Bckg_type': Bckg_type,
                                  'DEER_SNR': DEER_SNR, 'polyparam': ''})

        if mode == 'strexp_fixed':
            Bckg_type = str(mode+'_'+str(Percent_tmax)+'_'+str(truncsize))
            DEER_SNR = (p0[0]/sigma_noise, p0[0]/sigma_noise_half)
            Exp_parameter.update({'ModDepth': p0[0], 'decay': p0[1],
                                  'stretch': p0[2], 'Bckg_type': Bckg_type,
                                  'DEER_SNR': DEER_SNR, 'polyparam': ''})
        elif mode == 'poly':
            Bckg_type = str(mode+str(len(p0)-1)+'_'+str(Percent_tmax)+'_' +
                            str(truncsize))
            Mod_Depth = FinalData[itmin-1:itmin+1].mean()
            DEER_SNR = (Mod_Depth/sigma_noise, Mod_Depth/sigma_noise_half)
            Exp_parameter.update({'ModDepth': Mod_Depth, 'polyparam': p0,
                                  'Bckg_type': Bckg_type,
                                  'DEER_SNR': DEER_SNR})
        elif mode == 'polyexp':
            Bckg_type = str(mode+str(len(p0)-1)+'_'+str(Percent_tmax)+'_' +
                            str(truncsize))
            Mod_Depth = FinalData[itmin-1:itmin+1].mean()
            DEER_SNR = (Mod_Depth/sigma_noise, Mod_Depth/sigma_noise_half)
            Exp_parameter.update({'ModDepth': Mod_Depth, 'polyparam': p0,
                                  'Bckg_type': Bckg_type,
                                  'DEER_SNR': DEER_SNR})
        ParamDict.update({TITL: Exp_parameter})
    return DataDict, ParamDict


def GetPakePattern(t, y):
    tmax = np.max(t)
    npts = y.shape[0]
    newt = np.linspace(0, tmax*5, 5*npts)/1000  # Time axis in us
    freq = fdaxis(TimeAxis=newt)
    win = windowing(window_type='exp+', N=npts, alpha=3)
    y2 = np.zeros((npts*5,), dtype="complex_")
    y2[0:npts,] = y[0:npts,]*win[0:npts,]
    PakeSpectra = np.fft.fftshift(np.fft.fft(y2))
    Pivot = int(np.floor(PakeSpectra.shape[0]/2))
    PakeSpectraPhased, _ = automatic_phase(vector=PakeSpectra, pivot1=Pivot,
                                           funcmodel='minfunc')
    PakeAbs, _, _, _ = basecorr1D(x=freq, y=np.absolute(PakeSpectraPhased),
                                  polyorder=1, window=200)
    return PakeSpectraPhased, PakeAbs, freq


def ExportParamTable(FolderPath=getcwd(), ParamDict=None):
    # '''
    # The script export a table of the important parameters deduced from the
    # export files created by Deer-Analysis. (Typically the files are named :
    #                                   "Filename_res.txt")
    # In the following order, it To redescribe

    #     Returns
    #     -------
    #     None.
    # To do, check if a version with https://xlsxwriter.readthedocs.io/
    # would be better
    #     '''
    ListOfFiles = list(ParamDict.keys())
    Paramlist = ['tmax', 'tmin', 'NoiseLevel', 'ModDepth', 'DEER_SNR',
                 'RelBckgndDecay', 'RMSE', 'polyparam']
    PolyOrder = [0]*len(ListOfFiles)
    for i in range(len(ListOfFiles)):
        PolyOrder[i] = len(ParamDict[ListOfFiles[i]]['polyparam'])
    maxn = max(PolyOrder)
    shape = (int(len(Paramlist)+maxn), len(ListOfFiles))
    Param = np.full(shape, np.nan, dtype=float)
    for i in range(len(ListOfFiles)):
        Value = []
        for j in range(len(Paramlist)):
            Value = ParamDict[ListOfFiles[i]][Paramlist[j]]
            if type(Value) == tuple:
                Param[j, i] = Value[1]
            elif type(Value) == np.ndarray:
                for k in range(len(Value)):
                    Param[j+k, i] = Value[k]
            else:
                try:
                    Param[j, i] = Value
                except:
                    type(Value) == str
    fullfilename = path.normpath(path.join(FolderPath, 'DeerParam.csv'))
    with open(fullfilename, 'w') as file:
        wr = csv.writer(file, dialect='excel',
                        delimiter='\t', lineterminator='\n')
        wr.writerow(ListOfFiles)
        for row in Param:
            wr.writerow(row)
    return


def ExportOneAscii(fullpath_to_filename, DataDict, ParamDict, mode='time'):
    'Export a .txt ascii datafile at the same location of the datafolder file.'
    'append the keyword "waveletname_denoised" at the end of ascii filename'
    fileID = path.split(fullpath_to_filename)
    Exp_parameter = ParamDict.get(fileID[1])
    itmax = Exp_parameter.get('itmax')
    fulldata = DataDict.get(fileID[1])
    if mode == 'time':
        time = fulldata[:, 0].ravel()
        y = fulldata[:, 5].ravel()
        cleanedy = y[~np.isnan(y)][:itmax].real
        cleanedx = (time[~np.isnan(time)][:itmax].real)/1000
        newfilename = str(fileID[1]+'.txt')
        fullnewfilename = path.join(fileID[0], newfilename)
    elif mode == 'timeden':
        time = fulldata[:, 0].ravel()
        ydenoised = fulldata[:, 5].ravel()
        cleanedy = ydenoised[~np.isnan(ydenoised)][:itmax].real
        cleanedx = (time[~np.isnan(time)][:itmax].real)/1000
        wavename = Exp_parameter.get('wavename')
        mode = Exp_parameter.get('denoisedmode')
        newfilename = str(fileID[1]+'_'+wavename+'_'+mode+'_denoised.txt')
    fullnewfilename = path.join(fileID[0], newfilename)
    Bckg_type = Exp_parameter.get('Bckg_type')
    # HeaderTITL = str(TITL+'_Header')
    # Header = DataDict.get(HeaderTITL)
    header = '\n'.join(['Filename: ' + newfilename,
                        'Bckg_type: ' + Bckg_type,
                       'FirstColumn: ' + 'Time[us]',
                        'SecondColumn: ' + 'BckgSubtractedData'])
    if cleanedy.shape == cleanedx.shape:
        data = np.column_stack([cleanedx, cleanedy])
        makedirs(path.dirname(fullnewfilename), exist_ok=True)
        np.savetxt(fullnewfilename, data, fmt=['%.5e', '%.15e'],
                   delimiter='\t', header=header)
    # , footer='', comments='# ', encoding=None)
    return fullnewfilename


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


def importasciiSVDData(fullpath_to_filename, uncertainty=False):
    Header = {}
    x, y, y_up, y_down = [], [], [], []
    if uncertainty == False:
        with open(fullpath_to_filename, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if lines[i][0] == '%':
                    pairs = lines[i].split('\t')
                    Header.update({str(i): pairs})
                else:
                    pairs = lines[i].split('\t')
                    x.append(float(pairs[0]))
                    y.append(float(pairs[1]))
        return np.asarray(y), np.asarray(x), Header
    else:
        with open(fullpath_to_filename, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if lines[i][0] == '%':
                    pairs = lines[i].split('\t')
                    Header.update({str(i): pairs})
                else:
                    pairs = lines[i].split('\t')
                    x.append(float(pairs[0]))
                    y.append(float(pairs[1]))
                    y_down.append(float(pairs[2]))
                    y_up.append(float(pairs[3]))
        return np.asarray(y), np.asarray(x), Header, np.asarray(y_down), np.asarray(y_up)


def importasciiTimeData(fullpath_to_filename):
    Header = {}
    t, y1, y2 = [], [], []
    with open(fullpath_to_filename, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i][0] == '%':
                pairs = lines[i].split('\t')
                Header.update({str(i): pairs})
            else:
                pairs = lines[i].split('\t')
                t.append(float(pairs[0]))
                y1.append(float(pairs[1]))
                y2.append(float(pairs[2]))
    return np.asarray(y1), np.asarray(y2), np.asarray(t), Header
