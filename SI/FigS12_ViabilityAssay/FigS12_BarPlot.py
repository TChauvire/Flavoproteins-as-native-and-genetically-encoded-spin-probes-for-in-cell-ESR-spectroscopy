# -*- coding: utf-8 -*-
"""
Bar plot figure of the viability assay

Created on Sat Mar 22 17:12:25 2025

@author: Timothee Chauvire, tsc84@cornell.edu
"""

import matplotlib.pyplot as plt
from os import getcwd, makedirs
plt.close('all')
font = {'family': 'tahoma',
        'weight': 'normal',
        'size': 18}
# plt.rc('lines', linewidth=18)
plt.rc('font', **font)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)


name = ['No Light', '3s', '30s', '180s']
A = [8030000, 3750000, 640000, 110000]
color = ['black', 'blue', 'blue', 'blue']
edgecolor = ['dimgrey', 'dimgrey', 'dimgrey', 'dimgrey']
fig, ax = plt.subplots()
ax.bar(name, A, color=color, edgecolor=edgecolor)

ax.set_xlabel("TarCheWCheA-iLOV", fontsize=20, weight='bold')
ax.set_ylabel(r'CFU/mL/OD$_{600}$', fontsize=20, weight='bold')
ax.set_ylim(10000, 20000000)
ax.set_yscale('log')
ax.tick_params(axis='both', width=1.5)
for axis in ['bottom', 'left']:
    ax.spines[axis].set_linewidth(1.5)
# ax.spines['both'].set_linewidth(2)
# Show graph
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
# Save Graph
folder = getcwd() + '\\Figure\\'
makedirs(folder, exist_ok=True)
size = [12,  8]
# plt.tight_layout()
plt.gcf().set_size_inches(size[0], size[1])
plt.savefig(folder + 'figure_ViabilityAssay.pdf', transparent=True)
plt.savefig(folder + 'figure_ViabilityAssay.tiff', transparent=True)
