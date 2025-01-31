# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 17:48:00 2023

@author: tim_t
"""

import numpy as np
import matplotlib.pyplot as plt
from os import path, makedirs, getcwd, walk

folder0 = getcwd()
folder1 = getcwd() + '\\Pr\\'


# Define two useful functions for importing data
def ImportMultipleNameFiles(FolderPath=getcwd(), Extension=None,
                            *args, **kwargs):
    ListOfFiles = []
    for root, dirs, files in walk(FolderPath):
        for file in files:
            if file.endswith(Extension):
                ListOfFiles.append(path.normpath(
                    path.join(root, file)))
    return ListOfFiles


def importasciiSVDData(fullpath_to_filename, uncertainty=False):
    Header = {}
    x, y, y_up, y_down = [], [], [], []
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
    return np.asarray(y), np.asarray(x), Header, np.asarray(y_down), \
        np.asarray(y_up)


# Define a gaussian function to generate gaussian curve associated to a
# specific distance. The only parameter to choose is the width (wid).
def Gaussian(x, params):
    yfit = np.zeros_like(x)
    ctr = params[0]
    wid = params[1]
    amp = params[2]
    yfit = yfit + amp*np.sqrt(1/(2*np.pi))*(1/wid) * \
        np.exp(-0.5*((x - ctr)/wid)**2)
    return yfit


# Start to import the distance reference in the txt file N5histogram.txt
Listoffilename = ImportMultipleNameFiles(folder0, 'txt')
Header = {}
x, y = [], []
fit = np.zeros((1024,))
r = np.linspace(1.5, 10, 1024)
with open(Listoffilename[0], 'r') as f:
    lines = f.read().splitlines()  # lines()
    for i in range(len(lines)):
        pairs = lines[i].split('  ')
        x.append(float(pairs[0]))
        y.append(float(pairs[1]))

for i in range(len(x)):
    width = 0.15
    params = [x[i]/10, width, y[i]]
    yfit = Gaussian(r, params)
    fit += yfit

# normalize the distance domain of the homology model (Total sum of
# probability = 1)
intgValue = np.trapz(fit.ravel(), r.ravel())
fitnorm = fit/intgValue

# plot the figures
plt.close('all')
color = ['b', 'k', 'g', 'r']
Listoffilename2 = ImportMultipleNameFiles(folder1, 'txt')


font = {'family': 'tahoma',
        'weight': 'normal',
        'size': 24}
plt.rc('grid', linestyle="--", color='grey')
plt.rc('font', **font)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)


# Start the import of the SVD distance domain data and compare it to the
# reconstituted distance domain data of the homology domain.
for i in range(len(Listoffilename2)):
    fig, axes = plt.subplots(1, 1)
    filename = Listoffilename2[i]
    fileID = path.split(filename)
    l1, = axes.plot(r, fitnorm, color='brown', label='Homology Model',
                    linewidth=2)
    axes.set_xlabel("Distance [nm]", fontsize=20, fontweight='bold')
    axes.set_ylabel("P(r)", fontsize=20, fontweight='bold')
    axes.set_xlim(2.0, 9)
    axes.grid(visible=True, which='major')
    # # Get SVD Distance Domain spectra:
    P1, r1, header, P1_down, P1_up = importasciiSVDData(
        Listoffilename2[i], uncertainty=True)
    intgValue2 = np.trapz(P1.ravel(), r1.ravel())
    fig.suptitle(fileID[1], fontsize=20, fontweight='bold')
    # PR plot:
    axes.fill_between(x=r1, y1=(P1_down)/intgValue2, y2=(P1_up)/intgValue2,
                      color="red", alpha=0.5, lw=0.1, interpolate=True)
    l2, = axes.plot(r1, P1/intgValue2, color=color[i], label='SF-SVD Method',
                    linewidth=2)
    axes.set_yticks([0.])
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    fig.set_dpi(300)
    folder3 = folder0 + '\\Figures\\'
    makedirs(folder3, exist_ok=True)
    size = [12, 8]
    plt.gcf().set_size_inches(size[0], size[1])
    plt.savefig(folder3 + 'figure_Dist_{0}.pdf'.format(i), transparent=True)
    plt.savefig(folder3 + 'figure_Dist_{0}.tiff'.format(i), transparent=True)
