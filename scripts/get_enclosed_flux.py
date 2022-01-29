#script for getting enclosed fluxes

from __future__ import print_function

__author__ = "E.A.K. Adams"

from hitools import Source as S
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

#read the file of information need
galaxy_table = ascii.read('/Users/adams/data/pvfit/btfr_table.txt')

#initialize an array to hold enclosed HI flux values
#and also the minor axis
encl_hi_flux = np.empty(len(galaxy_table))
minor = np.empty(len(galaxy_table))

#get skycoord array
coords = SkyCoord(galaxy_table['RA'],galaxy_table['Dec'])

#now iterate through galaxies
for n, galaxy in enumerate(galaxy_table['galaxy']):
    #initialize a source object for the galaxy
    print(galaxy)
    galSource = S.Source(galaxy, coords[n].ra.deg, coords[n].dec.deg,
                         galaxy_table['inc'][n]*u.deg,
                         galaxy_table['PA'][n]*u.deg,
                         galaxy_table['DPV'][n]*u.arcsec,
                         galaxy_table['mom0path'][n])
    minor[n] = galSource.minor.to(u.arcsec).value
    #check for source name to handle flux conversion if necessary
    if galaxy[0] == 'A':
        #shield
        encl_hi_flux[n] = galSource.get_subcube_flux()
    else:
        #convert m/s to km/s
        encl_hi_flux[n] = galSource.get_subcube_flux() / 1000.

#add new columns to table
galaxy_table['Encl_HI_flux'] = encl_hi_flux
galaxy_table['minor'] = minor

#and write table out

ascii.write(galaxy_table['galaxy','DPV','minor','inc','HIflux','Encl_HI_flux'],
            '/Users/adams/data/pvfit/table_calc_encl_hi.txt', overwrite = True,
            formats = {'Encl_HI_flux': '5.2f'})

                         
    


