#Object for individual sources

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Class for individual sources
"""

import numpy as np
import astropy.units as u
from spectral_cube import SpectralCube
from regions import EllipseSkyRegion, EllipsePixelRegion
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
try:
    from hitools.tools import get_eff_beam, get_beam_area
except:
    from tools import get_eff_beam, get_beam_area

class Source(object):
    """
    Object for working with individual sources.

    Attributes:
    -----------
    name : string
        Name of source
    ra : Float or Quantity
        Central R.A., presumed in degrees if float
    dec : Float or Quantity
        Central Dec. in degrees, presumed in degrees if float
    incl : Float or Angle/Quantity
        Inclinaiton, presumed in degrees if float
    pa : Float or Angle/Quantity
        Position angle of source, presumed in degrees if float
    radius : Float or Angle/Quantity
        Angular extent of source, presumed in arcsec if float
    mom0 : SpectralCube instance
        moment zero map of source
    
    Methods:
    --------
    get_source_ellipse(type) : Get ellipse of source based on radius, incl, pa, ra, dec. Specify observation type to include beam smearing
    """


    def __init__(self, name, ra, dec, incl, pa, major, mom0file):
        """
        Initialize the Source object
        Require name, ra, dec, inclination (incl), position angle, major and mom0 filepath
        Major is measured major axis; used for auto-masking methods
        """
        self.name = name
        if isinstance(ra, u.Quantity):
            self.ra = ra
        else:
            self.ra = ra * u.deg
        if isinstance(dec, u.Quantity):
            self.dec = dec
        else:
            self.dec = dec * u.deg
        if ( isinstance(incl, u.Quantity) or isinstance(incl, u.Angle)):
            self.incl = incl
        else:
            self.incl = incl * u.deg
        if ( isinstance(pa, u.Quantity) or isinstance(pa, u.Angle)):
            self.pa = pa
        else:
            self.pa = pa * u.deg
        if ( isinstance(major, u.Quantity) or isinstance(major, u.Angle)):
            self.major = major
        else:
            self.major = major * u.arcsec
        self.mom0 = SpectralCube.read(mom0file)

        #get source ellipse as part of initilization
        self.get_source_ellipse()

    def get_source_ellipse(self, q = 0.2, factor = 1):
        """ 
        Get an ellipse that matches source, based on major, incl and PA
        Use beam parameters (if available) to determine minor axis
        Provide disk thickness, default = 0.2
        """
        #first get naive minor, from must incl & major
        #will then decide whether disk thickness or beam is more important
        naive_minor = self.major * np.cos(self.incl.to(u.radian).value)

        #get the effective beam
        try:
            #get beam parameters. assume these are in mom0
            #want along minor axis
            eff_beam_minor = get_eff_beam(self.mom0.beam.major, self.mom0.beam.minor,
                                          self.mom0.beam.pa, (self.pa - 90*u.deg))
        except:
            print("Was not able to calculate effective beam. Setting to zero")
            eff_beam_minor = 0*u.arcsec

        #then test whether I worry about beam or disk thickness and get "true" minor
        if naive_minor < 3 * eff_beam_minor:
            print("Accounting for beam effect in minor axis")
            #get eff_beam_major as need this
            eff_beam_major = get_eff_beam(self.mom0.beam.major, self.mom0.beam.minor,
                                          self.mom0.beam.pa, self.pa)
            minor = np.sqrt( np.cos(self.incl)**2 * (self.major**2 - eff_beam_major**2) +
                             eff_beam_minor**2 )
        else:
            print(f"Accounting for a disk thickness of {q}")
            minor = self.major * np.sqrt( (1-q**2)*np.cos(self.incl)**2 + q**2)

        #then add a check that I'm not overdoing things w/ my correction
        #that is, minor should be less than major
        if minor > self.major:
            minor = self.major

        #and then define my ellipse region
        #use factor to set size of region
        center_sky = SkyCoord(self.ra, self.dec, frame='icrs')
        region = [EllipseSkyRegion(center=center_sky,height=factor * self.major,
                                   width=factor * minor, angle=self.pa)]
        #add as attribute to object
        #i should probably make this part of initialization but i'm lazt
        #and will save that for now
        self.source_region = region[0]
        self.mom0_subcube = self.mom0.subcube_from_regions(region)
        self.minor = minor

    def show_source_region(self):
        """ 
        Show the source region as a check
        """

        fig, ax = plt.subplots(1, 1, figsize=(8,6),
                           subplot_kw={'projection': self.mom0.wcs,
                                      'slices': ('x', 'y', 0)})
        ax.imshow(self.mom0[0,:,:].value)
        pixel_region = self.source_region.to_pixel(self.mom0.wcs)
        pixel_region.plot(ax = ax)

    def show_subcube(self):
        """
        Show the subcube as a check
        """
        fig, ax = plt.subplots(1, 1, figsize=(8,6),
                           subplot_kw={'projection': self.mom0_subcube.wcs,
                                      'slices': ('x', 'y', 0)})
        ax.imshow(self.mom0_subcube[0,:,:].value)

    def get_subcube_flux(self):
        """ 
        Get flux w/in subcube
        Do this by summing all pixels and multiplying by number of beams
        """
        sum_subcube = np.nansum(self.mom0_subcube[0,:,:].value)
        #npix = len(np.where(self.mom0_subcube[0,:,:] > -1)[0])
        #pixels per beam is an attribute of Spectral Cube (yay!)
        #nbeams = npix / self.mom0_subcube.pixels_per_beam
        flux = sum_subcube / self.mom0_subcube.pixels_per_beam

        return flux
                


    
