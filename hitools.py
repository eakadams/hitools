#Useful HI tools

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Collection of HI tools that I find useful
"""

import numpy as np
import astropy.units as u

def get_nhi(sigma,sens,linewidth,beam_maj,beam_min):
    """
    Get N_HI sens at sigma level
    Inputs:
    - sigma: Number, significance level for result
    - sens: Quantity, sensitivity matched to linewidth
    - linewidth: Quantity, linewidth from N_HI calc
    - beam_fwhm: Quantity, beam fwhm
    """
    omega_beam = get_beam_area(beam_maj,beam_min)
    #hitools, hi freq
    #should maybe define globally?
    freq = 1420.405752 * u.MHz
    #convert to brightness temp
    tb = sens.to(u.K,u.brightness_temperature(freq,beam_area=omega_beam))
    #convert to N_HI
    nhi = 1.823e18 * tb * linewidth  / (u.K * u.km / u.s) / (u.cm * u.cm)
    nhi_sig = sigma*nhi
    return nhi_sig

def get_beam_area(beam_maj, beam_min):
    """
    Calculate beam area, provided FWHM of maj and min axis
    Inputs are quantities - have angular units
    """
    fwhm_to_sigma = 1. / (8 * np.log(2))**0.5
    maj_sigma = beam_maj * fwhm_to_sigma
    min_sigma = beam_min * fwhm_to_sigma
    beam_area = 2 * np.pi * min_sigma * maj_sigma
    
    return beam_area
    

def get_sens_for_res(sens,spec_res,desired_res):
    """
    Get the sensitivity for a desired spectral resolution
    Inputs:
    - sens: Quantity, sensitivey at given spec_res
    - spec_res: Quantity, spec. resolution for input sensitivity
    - desired_res: Quantity, desired spec res.
    Outputs:
    - new_sens: Quantity, sensitivity for desired_res
    """
    #first check that have same units for res
    #otherwise have to us an equivalency
    if not spec_res.unit == desired_res.unit:
        print("Units don't match")
        print("Checking if they are convertible")
        try:
            desired_res.to(spec_res.unit)
        except:
            print("Not convertible so change desired_res in freq/vel")
            desired_res = convert_chan_freq_vel(desired_res)
            #force to same units
            desired_res=desired_res.to(spec_res.unit)
            print(desired_res)

    #calculate the new sensitivity
    new_sens = sens * np.sqrt(spec_res / desired_res)

    #return the new sensitivity value
    return new_sens
        
def convert_chan_freq_vel(chan_res):
    """
    Takes an input channel resolution, either in freq or velocity
    Will return matching resolution in other dimensions
    Presumes radio definition; can consider as a param for future
    Inputs:
    chan_res (Quantity): Input resolution, with units
    Outputs:
    new_res (Quantity): Output resolution, with units, opp of input
    """
    #define conversion
    restfreq = 1420.405752 * u.MHz
    freq_to_vel = u.doppler_radio(restfreq)
    #check if in frequency; if so, get vel
    if 'Hz' in chan_res.unit.to_string():
        offsetfreq = restfreq - chan_res
        vel = restfreq.to(u.km/u.s, equivalencies = freq_to_vel)
        offsetvel = offsetfreq.to(u.km/u.s, equivalencies = freq_to_vel)
        new_res = offsetvel-vel

    else:
        print("Not frequency unit, presuming velocity")
        offsetfreq = chan_res.to(u.kHz, equivalencies = freq_to_vel)
        new_res = restfreq-offsetfreq

    return new_res
