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
import hitools as hit
from astropy.io import ascii

class Survey(object):
    def __init__(self, sens, spec_res, beam_maj, beam_min, footprint,
                 nbeam=3, snint=5, snres=3):
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
        nbeam (optional): Float giving number of maj beams for to resolve. 
                          Default is 3
        snint (opt): S/N for detection (flux); default=5
        snres (opt): S/N for res detection; default = 3
        """
        self.sens = sens
        self.spec_res = spec_res
        self.beam_maj = beam_maj
        self.beam_min = beam_min
        self.footprint = footprint
        self.nbeam = nbeam
        self.snint = snint
        self.snres = snres

        #as part of initializitation, calcule 3sig N_HI (over 20 km/s)
        sens_20 = hit.get_sens_for_res(self.sens,
                                           self.spec_res,20*u.km/u.s)
        self.nhi_sens = hit.get_nhi(self.snres,sens_20,20*u.km/u.s,
                                        self.beam_maj,self.beam_min)


    def predic_himf(self, mlow,mhigh,w50,
                    alpha=-1.25,  mstar=10**9.94, phi=4.5e-3):
        """
        Predictions based on a HIMF
        described by alpha, phi, mstar
        Default are global values from Jones+2018
        Mstar in solar masses, not log units
        Inputs:
        - mlow (Quantity): Lower bound of HI masses to consider, solar masses
        - mhigh (Quantity): Upper bound of HI masses to consider, solar masses
        - W50 (Quantity): Velocity width to consider for detection, km/s
        """

        #Calculate detection limit for survey and provided w50
        sdet = self.snint * w50 * hit.get_sens_for_res(self.sens,self.spec_res,w50) * u.beam

        #Want to provide estimates both for detections resolved
        #Thus need to calculate both ways
        #for setting distance limits, should I do that by min or max in mass?
        #maybe do both and give range?
        #also presume that user will limit mass ranges provided here

        #do point source detection first, easiest
        #integrate_himf not set up for quantities, so use values
        ##check my numbers
        #print(mlow.value, mhigh.value, sdet.to(u.Jy*u.km/u.s).value,
        #                                      self.footprint.to(u.deg*u.deg).value)
        dmax,Ndet = hit.integrate_himf_proper(phi,mstar, alpha,
                                              mlow.value, mhigh.value, sdet.to(u.Jy*u.km/u.s).value,
                                              self.footprint.to(u.deg*u.deg).value)

    
        return dmax*u.Mpc, np.float(Ndet)

    def predic_himf_res(self,mlow,mhigh,
                        alpha=-1.25,  mstar=10**9.94, phi=4.5e-3):
        """
        Predictions from HIMF for resolved sources
        Beam size and res requirements are part of object
        """
        #Effective beam has same area
        beam_eff = np.sqrt(self.beam_maj*self.beam_min)
        dmax, Nres = hit.int_himf_res(phi,mstar,alpha,mlow.value,mhigh.value,
                                            self.footprint.to(u.deg*u.deg).value,
                                            beam_eff.to(u.arcsec).value,self.nbeam)

        return dmax, Nres

    def predic_a40(self):
        """
        Predictions based on a40, including stellar masses from Shan
        """
        #read in a40 table
        table_a40 = ascii.read('/Users/adams/python/UDGs/a40_masses.csv',format='csv')
        #Want to find how many with HI masses < 10^9 also have stellar masses in that rnage
        #that is, what fraction are truly dwarfs?
        #and not gas poor more massive galaxies?
        ind_mhi9 = np.where(table_a40['logMHI'] <= 9)[0]
        ind_mstar9 = np.where(table_a40['logMstar'][ind_mhi9] <= 9.)[0]
        print(("{0:4.2f}% of galaxies with HI masses below 10^9 also have"
               " stellar masses below 10^9").format(len(ind_mstar9)/len(ind_mhi9)*100.))

        ind_mstar9 = np.where(table_a40['logMstar'] <= 9.)[0]
        ind_mhi9 = np.where(table_a40['logMHI'][ind_mstar9] > 9)[0]

        print(("{0:4.2f}% of galaxies with stellar masses below 10^9  have"
               " HI masses greater than 10^9").format(len(ind_mhi9)/len(ind_mstar9)*100.))

        

    def predic_a100(self):
        """
        Predictions based on a100
        https://ui.adsabs.harvard.edu/abs/2020AJ....160..271D/abstract
        Table 2
        """
        print("ALFALFA-based predictions not yet enabled")

    def predic_hunt(self):
        """
        Predictions based on Mstar-SFR work from L. Hunt
        This is specific to SDSS spring footprint w/in Apertif wide survey
        """
        print("Predictions from Leslie not yet enabled")
