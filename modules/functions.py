#Functions for HI tools

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Collection of functions that I use
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import os

#Get file dir
this_dir,this_filename = os.path.split(__file__)
hitoolsdir = this_dir[:-7]
filedir = os.path.join(aperinfodir,"files")


def plot_hi_dwarfs():
    """
    Function for plotting HI dwarfs in context
    ALFALFA as parent / blind sample
    Highlight SHIELD, HUDs
    Add LITTLE THINGS, FIGGS, VLA/ANGST for context
    """
    #get ALFALFA data
    a100 = ascii.read(os.path.join(filedir,"a100.code12.table2.190808.csv"))
    a100_sdss_xmatch = ascii.read(os.path.join(filedir,"durbala2020-table1.21-Sep-2020.fit"))
    
