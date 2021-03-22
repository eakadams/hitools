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
    def __init__(self, sens, spec_res, beam_maj, beam_min, footprint,
                 res=3):
        """
        Initialize the survey object
        Take both major & minor beam to be able to account 
        for Apertif decl dep (roughly)
        For most other programs, maj=min, just have to enter twice
        Inputs:
        sens: Quantity giving the sensitivity; e.g., u.mJy / u.beam
        spec_res: Quantity giving the spectral resolution
        beam_maj: Quantity giving the beam major axis (FWHM)
        beam_min: Quantify giving the beam minor axis (FWHM)
        footprint: Quantity giving the area of the survey
        res (optional): Float giving number of maj beams for to resolve
        """
        self.sens = sens
        self.spec_res = spec_res
        self.beam_maj = beam_maj
        self.beam_min = beam_min
        self.footprint = footprint

        #as part of initializitation, calcule 3sig N_HI (over 20 km/s)
        sens_20 = hitools.get_sens_for_res(self.sens,
                                           self.spec_res,20*u.km/u.s)
        self.nhi_3sig = hitools.get_nhi(3,sens_20,20*u.km/u.s,
                                        self.beam_fwhm)


    def predic_himf():
        """
        Predictions based on a HIMF
        """
        print("HIMF predictions not yet enabled")

    def predic_a100():
        """
        Predictions based on a100
        https://ui.adsabs.harvard.edu/abs/2020AJ....160..271D/abstract
        Table 2
        """
        print("ALFALFA-based predictions not yet enabled")

    def predic_hunt():
        """
        Predictions based on Mstar-SFR work from L. Hunt
        This is specific to SDSS spring footprint w/in Apertif wide survey
        """
        print("Predictions from Leslie not yet enabled")
