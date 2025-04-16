# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 14:37:37 2024

@author: tim_t
"""

from ImportMultipleFiles_FigS11 import ImportMultipleNameFiles
from ImportMultipleFiles_FigS11 import OpenDoubleIntegralKin2
from os import path, makedirs, getcwd
import matplotlib.pyplot as plt

folder = getcwd() + '\\2D\\'

ListOfFiles = ImportMultipleNameFiles(folder, Extension='.DSC')
nfile = 2  # Select the right file here. Check ListOfFiles variable for a ref.
col = 10*[0]  # An average on 10 spectrum was done.
for i in range(10):
    col[i] = 36+i  # starting point has to be determined manually to correspond
    # to the right time in the kinetic. To be optimized in the future!
print(col)

Filename = ListOfFiles[nfile]
FirstDeriv, ZeroDeriv, Header, IntgValue, HeaderInt = OpenDoubleIntegralKin2(
    Filename, Scaling=None, col=col, polyorder=[1, 1], window=200)  #

fileID = path.split(ListOfFiles[nfile])
# plot the first derivative
x = FirstDeriv[:, 0]/10
y = FirstDeriv[:, 1]
x2 = ZeroDeriv[:, 0]/10
y2 = ZeroDeriv[:, 1]
plt.figure()
fig, axes = plt.subplots(2, 1)
fig.suptitle(fileID[1])
l1, = axes[0].plot(x.ravel(), y.ravel(), 'k',
                   label='First Derivative', linewidth=2)
axes[0].set_xlabel("Magnetic Field [mT]", fontsize=14, weight='bold')
axes[0].set_ylabel("Intensity [a.u.]", fontsize=14, weight='bold')
# axes[0].legend(loc='upper right')
axes[0].set_title('First Derivative')
l2, = axes[1].plot(x2.ravel(), y2.ravel(), 'k',
                   label='First Derivative', linewidth=2)
axes[1].set_xlabel("Magnetic Field [mT]", fontsize=14, weight='bold')
axes[1].set_ylabel("Intensity [a.u.]", fontsize=14, weight='bold')
# axes[0].legend(loc='upper right')
axes[1].set_title('Zero Derivative')
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.subplots_adjust(left=0.046, right=0.976, top=0.945, bottom=0.085,
                    hspace=0.2, wspace=0.2)
plt.tight_layout()
makedirs(fileID[0] + '\\Figures\\', exist_ok=True)
plt.savefig(fileID[0] + '\\Figures\\figure{0}'.format(nfile))


del axes, col, Filename, fig, figManager, fileID, l1, l2, x, x2, nfile, folder
del i, y, y2
