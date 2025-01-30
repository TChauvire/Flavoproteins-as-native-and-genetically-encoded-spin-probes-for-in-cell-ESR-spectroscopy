# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 12:10:59 2024

@author: tim_t
"""
from ImportMultipleFiles_Fig1 import ImportMultipleNameFiles, OpenMultipleFiles
from os import getcwd
folder = getcwd() + '\\RawData\\'
ListOfFiles = ImportMultipleNameFiles(folder, Extension='.DSC')
fulldata, header = OpenMultipleFiles(
    ListOfFiles, Scaling=None, polyorder=2, window_length=100)
