# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:41:25 2023

Import Multiple cw-ESR files
"""

from datasmooth import *
from ImportMultipleFiles import ImportMultipleNameFiles, OpenMultipleFiles, eprload
from basecorr1D import *

import matplotlib.pyplot as plt
import numpy as np
from os import getcwd
#############################################################################
folder = getcwd()+"\\250321_LiPC"
fulldata, Header = [], []
ListOfFiles = ImportMultipleNameFiles(FolderPath=folder, Extension='DSC')
fulldata, Header = OpenMultipleFiles(ListOfFiles=ListOfFiles, Scaling=None,
                                     polyorder=1, window=100)
