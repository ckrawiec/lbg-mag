import healpy as hp
import numpy as np
import sys
import os
import glob
from astropy.io import fits
from astropy.table import Table

sv_region  = '/Users/Christina/DES/data/sva1_gold_r1.0_goodregions_04_n4096.fits.gz'

red_region = '/Users/Christina/DES/data/redmagic_v6.2.15_fracdet_nside4096_map.fit.gz'
red_file   = '/Users/Christina/DES/data/redmagic_sva1_public_v6.3_faint.fits'
red_out_file  = '/Users/Christina/DES/data/random_sva1_redmagic.fits'

def main():
    redrandoms()

def getmap(region_file, ra, dec):
    #clip the randoms to masked area
    hpmap = hp.read_map(region_file)
    nside = hp.npix2nside(hpmap.size)
    theta = (90.0-dec)*np.pi/180.
    phi = ra*np.pi/180.
    pix = hp.ang2pix(nside, theta, phi)
    return hpmap, pix

def redrandoms():
    if os.path.exists(red_out_file):
        print "{} already exists".format(red_out_file)
        print "Exiting."
        exit()

    red = fits.open(red_file)[1].data

    red_ra = red['RA']
    red_dec = red['DEC']
        
    ra_spread = np.max(red_ra)-np.min(red_ra)
    dec_spread = np.max(red_dec)-np.min(red_dec)
    sindec_spread = np.max(np.sin(red_dec*np.pi/180.))-np.min(np.sin(red_dec*np.pi/180.))
    
    #random points
    #sphere 
    rra = np.min(red_ra) + np.random.rand(len(red_ra)*10)*ra_spread
    uni_sindec = np.min(np.sin(red_dec)) + np.random.rand(len(red_dec)*10)*sindec_spread
    rdec = np.arcsin(uni_sindec) * 180./np.pi
    
    mask = np.array([True]*len(rra))
    
    rmap, rpix = getmap(red_region, rra, rdec)        
    for p in set(rpix):
        frac = rmap[p]
        if frac<0:
            frac=0
    
        tot = len(mask[np.where(rpix==p)])
        mask[np.where(rpix==p)] = [True]*np.ceil(frac*tot) +[False]*np.floor((1-frac)*tot)
    
    tab = Table()
    tab['RA'] = rra[mask]
    tab['DEC'] = rdec[mask]

    sys.stderr.write('length of file = {}\nwriting...\n'.format(len(tab)))
    tab.write(red_out_file)
    sys.stderr.write('written to {}\n'.format(red_out_file))
    

    
if __name__=="__main__":
    main()
        
