#Survey Predictions

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Module for making predictions for HI 
surveys
Class-based
"""

import numpy as np
import astropy.units as u

class Survey(object):
    def __init__(self, sens, spec_res, beam_fwhm, footprint):
        """
        Initialize the survey object
        sens: Quantity giving the sensitivity; e.g., u.mJy / u.beam
        spec_res: Quantity giving the spectral resolution
        beam_fwhm: Quantity giving the beam size (FWHM)
        footprint: Quantity giving the area of the survey
        """
        self.sens = sens
        self.spec_res = spec_res
        self.beam_fwhm = beam_fwhm
        self.footprint = footprint

        #as part of initializitation, calcule 3sig N_HI (over 20 km/s)
        
