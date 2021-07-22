#Useful HI tools

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Collection of HI tools that I find useful
"""

import numpy as np
import astropy.units as u
from astropy import constants as const
import mpmath as mp

def get_mhi(dist,sint):
    """
    Get M_HI given distance and flux integral

    Parameters
    ----------
    dist : Quantity
         Distance to the galaxy
    sint : Quantity
         Flux integral

    Returns
    -------
    mhi : Quantity
        HI mass of galaxy
    """

    factor = 2.356e5 * u.Msun * u.s  / (u.Mpc * u.Mpc * u.Jy * u.km)
    mhi = factor * dist**2 * sint

    return mhi.to(u.Msun)
    

def get_dist_res(logmhi,beam,nres):
    """
    Calculate distance at which a HI disk is resolved
    Assume M_HI-D_HI relation; don't account for scatter
    Inputs:
    - logmhi: Float, log of HI mass
    - beam: Quantity, major beam in angular units
    - nres: Float, number of beams to be considered resolved
    """
    disk_size = 10**(0.506*logmhi - 3.293) * u.kpc
    dist_res = (disk_size /
                (nres*beam )).to(u.Mpc,
                                 equivalencies=u.dimensionless_angles())
    return(dist_res)

def get_mass_resolved(dist,beam,nres):
    """
    For a given distance, find the mass which is resolved by nres beams
    Inputs:
    - Dist: Quantity, distance in Mpc
    - beam: Quantity, major beam in angular units
    - nres: Float, number of beams to be considered resolved
    """
    res_size = (dist * (nres*beam)).to(u.kpc,
                                       equivalencies=u.dimensionless_angles())
    logdisk = np.log10(res_size.value)
    logmhi = (logdisk + 3.293)/0.506
    mhi = (10**logmhi) * u.Msun

    return mhi
                                       
    

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

def get_msun_sqpc(sigma,sen,linewidth,beam_maj,beam_min):
    """
    Get column density sensivitivy in Msun/pc^2
    for specified sigma level
    Convert from N_HI, which is in units of cm^-1 (implied atoms)
    """
    nhi = get_nhi(sigma,sen,linewidth,beam_maj,beam_min)
    #put nhi into mass units, using mass of H atom
    #use atomic mass, slightly more than proton mass
    #maybe should use proton mass? const.m_p
    nhi_mass = nhi*const.u
    msun_pc2 = nhi_mass.to(u.M_sun / (u.pc * u.pc))
    return msun_pc2
    

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

def convert_vel_freq(quantity,restfreq=1420.405752 * u.MHz):
    """
    Take an input quantity that is either vel or freq and return the other
    Default for HI line
    Use radio convention (could set as param)
    """
    freq_to_vel = u.doppler_radio(restfreq)
    if 'Hz' in quantity.unit.to_string():
        new_quantity = quantity.to(u.km/u.s, equivalencies = freq_to_vel)
    else:
        #presume in vel, use try/except in case
        try:
            new_quantity = quantity.to(u.MHz,  equivalencies = freq_to_vel)
        except:
            print("Units of input not recognized")
            new_quantity = np.nan

    return new_quantity


#proper himf integration taking into account varying volume
def integrate_himf_proper(phi,mstar,alpha,mlow,mhigh,sdet,foot,dlimit=None,dmin=None):
    """
    Proper integration of the HIMF taking into account
    varying volume over which sources can be detected
    Inputs:
    - phi, mstar, alpha: HIMF parameters
    - mlow: Lower HI mass limit for integrating, solar mass units
    - mhigh: Upper HI mass limit for integrating, solar mass units
    - sdet: Integrated flux density which can be dected, Jy km/s
            This should already account for velocity width effects
    - foot: Footprint on sky in square degrees
    - dlimit: Maximum distance to consider in units of Mpc 
             (e.g., setting spatial resolution limits)
    - dmin: Minimum distance to consider, in units of Mpc
    """
    #first have to find maximum distance that can see mhigh
    maxd = np.sqrt(mhigh/(sdet*2.356e5))
    if dlimit != None:
        dmax = min(dlimit,maxd)
    else:
        dmax=maxd
    omega = foot*(np.pi/180.)**2 #survey area in steradian
    print(omega)
    if dmin != None:
        d_bins=np.arange(dmin,dmax,dmax/1000.)
    else:
        dmin = 0.
        d_bins = np.arange(0,dmax,dmax/1000.) #set up shells
    dn_full = phi*(mp.gammainc((alpha+1),mlow/mstar)-mp.gammainc((alpha+1),mhigh/mstar)) #number density for full range
    N=0.
    vol=0.
    if dmin < dmax:
        for d in d_bins:
            #print(d,omega*d**2*(dmax/1000.),vol)
            #for each shell calc number of sources
            #first have to figure out mhi_min 
            mhi_min = sdet*2.356e5*d**2
            #print(mhi_min)
            if mhi_min < mlow:
                #just use number density times volume of shell
                N = N + dn_full*omega*d**2*dmax/1000.
                #print(N)
            elif mhi_min < mhigh:
                dn = phi*(mp.gammainc((alpha+1),mhi_min/mstar)-mp.gammainc((alpha+1),mhigh/mstar)) #have to cal w/ limited lower mass
                N = N + dn*omega*d**2*dmax/1000.
            vol = vol+omega*d**2*dmax/1000.

    #print vol

    return dmax,N


def int_himf_res(phi,mstar,alpha,mlow,mhigh,foot,beam,nres):
    """
    An adapted integration of the HIMF where mass limits
    are overruled by requiring the main HI disk
    to be resolved by at least "nres" "beam"
    """
    #first, find the maximum distance to which mhigh is resolved
    #this is the upper limit
    dmax = get_dist_res(np.log10(mhigh),beam*u.arcsec,nres)
    dmax = dmax.value

    #set the survey area in steradian
    omega = foot*(np.pi/180.)**2 

    #set the bins
    d_bins = np.arange(0,dmax,dmax/1000.)

    #set numbers before iteration
    dn_full = phi*(mp.gammainc((alpha+1),mlow/mstar)-mp.gammainc((alpha+1),mhigh/mstar)) #number density for full range
    N=0.
    vol=0.

    for d in d_bins:
        #print(d,omega*d**2*(dmax/1000.),vol)
        #for each shell calc number of sources
        #first have to figure out mhi_min
        mhi_min = get_mass_resolved(d*u.Mpc,beam*u.arcsec,nres)
        mhi_min = mhi_min.value
        if mhi_min < mlow:
            #just use number density times volume of shell
            N = N + dn_full*omega*d**2*dmax/1000.
            #print(N)
        elif mhi_min < mhigh:
            dn = phi*(mp.gammainc((alpha+1),mhi_min/mstar)-mp.gammainc((alpha+1),mhigh/mstar)) #have to cal w/ limited lower mass
            N = N + dn*omega*d**2*dmax/1000.
        vol = vol+omega*d**2*dmax/1000.

    return dmax,np.float(N)


def calc_radio_lum(flux,dist):
    """
    Not strictly an HI tool but useful
    Calculate the radio luminosity given flux and dsitance
    Ignore cosmological correction because I work in low-z regime
    Maybe one day I'll worry about that sort of thing :)
    """
    lum = flux * 4 * np.pi * dist**2
    return(lum.to(u.W / u.Hz))

def calc_radio_sf(flux,dist):
    """
    Not strictly HI tool but useful
    Given a flux (1.4 GHz) and distance,
    calculate the star formation rate
    Use Bell 2003 description
    """
    lum = calc_radio_lum(flux,dist)
    lc = 6.4e21 * u.W / u.Hz
    if lum > lc:
        sf = 5.52e-22 * lum * u.Msun/u.yr * (u.Hz/u.W)
    else:
        sf = 5.52e-22 / (0.1+0.9*(lum/lc)**0.3) * lum * u.Msun/u.yr * (u.Hz/u.W)

    return sf
