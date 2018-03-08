import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import glob
from astropy.io import fits
from scipy.spatial import ckdtree
from scipy.interpolate import griddata

#params
params = {}
params['filters'] = 'GRIZ'
params['mu'] = 1.1
params['flux_cut'] = 40.
params['deep_file'] = '/Users/Christina/DES/data/y1a1_gold_dfull_cosmos.fits'
params['balrog_files'] = '/Users/Christina/DES/data/balrog_sva1_tab*_TRUTH_zp_corr_fluxes.fits'
params['sim_file_format'] = '/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab{}_SIM.fits'
params['deep_flux_column'] = 'FLUX_AUTO_{}'
params['deep_size_column'] = 'FLUX_RADIUS_I'
params['balrog_flux_column'] = 'FLUX_NOISELESS_{}'
params['balrog_size_column'] = 'HALFLIGHTRADIUS_0'


def createdeep(deep, mu, fluxcol='FLUX_AUTO_{}', sizecol='FLUX_RADIUS_I', fluxcut=40., filters='GRIZ', sizes=True, logs=True):
    #cut out low fluxes and unrealistic sizes
    flux_mask = (deep[fluxcol.format('G')]>cut) & (deep[fluxcol.format('R')]>cut) & (deep[fluxcol.format('I')]>cut) & (deep[fluxcol.format('Z')]>cut)
    size_mask = (deep[sizecol]>0.)
    
    #separate galaxies & stars using MODEST_CLASS
    gals = deep[(deep['MODEST_CLASS']==1) & flux_mask & size_mask]
    stars = deep[deep['MODEST_CLASS']==2]
    
    if logs:
        #flux vectors
        gal_vec = [np.log10(gals[fluxcol.format(f)]) for f in filters]

        #magnified flux vectors
        new_gal_vec = df_vec + np.log10(mu)
    else:
        gal_vec = [gals[fluxcol.format(f)] for f in filters]
        new_gal_vec = gal_vec * mu

    #get sizes
    gal_hlr = gethlr(gals, stars)
    if sizes:
        #add sizes to data vector
        gal_vec = np.vstack([gal_vec, gal_hlr])
        
        #add magnified sizes to magnified data vector
        new_gal_hlr = gal_hlr * np.sqrt(mu)
        new_gal_vec = np.vstack([new_gal_vec, new_gal_hlr])

    #trim data to where half light radius makes sense
    gal_data = np.array( zip(*gal_vec) )[gal_hlr>0.]
    new_gal_data = np.array( zip(*new_gal_vec) )[gal_hlr>0.]
    return gal_data, new_gal_data

def createbalrog(brog, fluxcol='FLUX_NOISELESS_{}', hlrcol='HALFLIGHTRADIUS_0',
                 fluxcut=40., filters='GRIZ', sizes=True, logs=True):
    #cut out low fluxes and unrealistic sizes
    flux_mask = (brog[fluxcol.format('G')]>0.) & (brog[fluxcol.format('R')]>0.) & (brog[fluxcol.format('I')]>0.) & (brog[fluxcol.format('Z')]>0.)
    size_mask = (brog[hlrcol]>0.)
    new_brog = brog[flux_mask & size_mask]

    if logs:
        br_vec = [np.log10(new_brog[fluxcol.format(f)]) for f in filters]
    else:
        br_vec = [new_brog[fluxcol.format(f)] for f in filters]
    if sizes:
        br_vec = np.vstack([br_vec, new_brog[hlrcol]])

    br_data = np.array( zip(*br_vec) )
    return br_data

def findmatches(br_data, deep_data, new_deep_data):
    #save match information in a dictionary
    matches = {'match radius': [],
               'magnified match radius': [],
               'index of match': [],
               'index of magnified match': []}

    #create balrog tree
    sys.stderr.write('creating tree...')
    br_tree = ckdtree.cKDTree(br_data)
    
    #query tree for fluxes
    sys.stderr.write('querying original...')
    orig_d, orig_id = br_tree.query(deep_data)
    matches['match radius'] = orig_d
        
    #query tree for magnified fluxes
    sys.stderr.write('querying magnified...')
    new_d, new_id = br_tree.query(new_deep_data)
    matched['magnified match radius'] = new_d
        
    #find balrog ids
    matches['index of match'] = orig_id
    matches['index of magnified match'] = new_id

    return matches  

def gethlr(gals, stars, sizecol='FLUX_RADIUS_I'):
    """Use median star size around galaxies & galaxy size
    from balrog to estimate true half light radius of galaxies"""
    grid_file = '/Users/Christina/DES/data/balrog/sva1/balrog_tab01_avg_star_fluxradiusi_0.1deg.fits'
    to_grid = fits.open(grid_file)[1].data
    sys.stderr.write('Grid file for average star radii:\n')
    sys.stderr.write('    {}\n'.format(grid_file))
    
    #make tree of dfull stars
    sys.stderr.write('    stars in deep table: {}\n'.format(len(stars)))
    sys.stderr.write('    galaxies in deep table: {}\n'.format(len(gals)))

    gal_pos = zip(gals['RA'], gals['DEC'])
    star_pos = zip(stars['RA'], stars['DEC'])

    #create star tree and find galaxies within 0.1 deg
    sys.stderr.write('    Creating hlr tree...')
    star_tree = ckdtree.cKDTree(star_pos)
    close = star_tree.query_ball_point(gal_pos, r=0.1) #360 arcsec
    sys.stderr.write('Done.\n')

    sys.stderr.write('    Calculating average star radii...')
    start = time.time()
    #calculate median of star radii within 0.1 deg of each galaxy
    gal_med_star_size = np.array([np.median(stars[sizecol][c]) for c in close])
    end = time.time()

    #interpolate true galaxy half light radius using balrog data
    gal_hlr = griddata(zip(to_grid['flux_radius_i'], to_grid['avg_flux_radius_i']), to_grid['hlr'],
                       zip(gals[sizecol], gal_med_star_size))
    sys.stderr.write('Done.\n')
    return gal_hlr

def getslope(table, column, mask):
    h = np.histogram(table[column][mask], bins=20)
    x_interp = np.array([np.mean(h[1][i-1:i+1]) for i in range(1,len(h[1]))])
    b, c = np.polyfit(x_interp, np.log10(h[0]), 1)
    slope = 2.5 * b
    return slope

def main(args):
    #params = parseconfigs(args[1])
    
    #open deep data file and create unlensed & lensed data vectors
    deep = fits.open(params['deep_file'])[1].data
    gal_data, new_gal_data = createdeep(deep, params['deep_flux_column'], params['deep_size_column'],
                                        params['flux_cut'], params['filters'])

    #indices of deep data
    ids = list(range(len(gal_data)))

    #cycle over all balrog tables
    balrog_matches = {}
    tabnums = []
    for tab in glob.glob(params['balrog_files']):
        tabnum = tab[41:43]
        tabnums.append(tabnum)
        sys.stderr.write('Working on balrog table {}...'.format(tabnum))

        #open tables, get object flux/size data vectors
        brog = fits.open(tab)[1].data
        br_data = createbalrog(brog, params['balrog_flux_column'], params['balrog_size_column'],
                               params['flux_cut'], params['filters'])
        these_matches = findmatches(br_data, gal_data, new_gal_data)

        #get balrog indices to match to sim tables later
        these_matches['balrog index of match'] = brog['BALROG_INDEX'][these_matches['index of match']]
        these_matches['balrog index of magnified match'] = brog['BALROG_INDEX'][these_matches['index of magnified match']]

        balrog_matches[tabnum] = these_matches

    #find closest match for each object among all tables
    key_arg = np.argmin(zip([balrog_matches[key]['match radius'] for key in balrog_matches.keys()]), axis=1)
    new_key_arg = np.argmin(zip([balrog_matches[key]['magnified match radius'] for key in balrog_matches.keys()]), axis=1)
    best_tabs = balrog_matches.keys()[key_arg]
    best_new_tabs = balrog_matches.keys()[new_key_arg]

    #save numbers of detections
    detections, new_detections = 0, 0
    for itab in range(len(tabnums)):
        #open sim catalog
        sim = fits.open(sim_file_format.format(tabnums[itab]))[1].data

        #focus on this table
        this_set = np.where(best_tabs==itab)
        this_new_set = np.where(best_new_tabs==itab)
        
        #count objects that are found (detected) in sim catalog
        sys.stderr.write('Working on sim table {}...'.format(tabnums[itab]))
        sys.stderr.write('checking for detections...')

        #detected objects whose truth fluxes match deep set
        detections += len(set(sim['BALROG_INDEX']).intersection(balrog_matches[itab]['balrog index of match'][this_set]))
        
        #detected objects whose truth fluxes match magnified deep set
        new_detections += len(set(sim['BALROG_INDEX']).intersection(balrog_matches[itab]['balrog index of magnified match'][this_new_set]))
        
    #report information
    print "Detected original matches: {}".format(detections)
    print "Detected magnified matches: {}".format(new_detections)

     
if __name__=="__main__":
    main(sys.argv)    
