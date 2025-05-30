# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 12:10:59 2024

@author: tim_t
"""
from ImportMultipleFiles import ImportMultipleNameFiles, MaxLengthOfFiles, OpenDavies
import numpy as np
from os import getcwd
folder = getcwd()
ListOfFiles = ImportMultipleNameFiles(folder, Extension='.DSC')
maxlen = MaxLengthOfFiles(ListOfFiles)
fulldata, header = OpenDavies(
    ListOfFiles, Scaling=None, polyorder=1, window_length=50)
