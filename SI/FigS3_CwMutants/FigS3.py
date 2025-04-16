# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:41:25 2023

Import Multiple cw-ESR files
"""

from ImportMultipleFiles import ImportMultipleNameFiles, OpenMultipleFiles
from os import getcwd
#############################################################################
folder = getcwd()+"//RawData//"
fulldata, Header = [], []
ListOfFiles = ImportMultipleNameFiles(FolderPath=folder, Extension='DSC')
fulldata, Header = OpenMultipleFiles(ListOfFiles=ListOfFiles, Scaling=None,
                                     polyorder=1, window=100)
