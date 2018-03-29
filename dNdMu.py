import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import glob
import measuresignal as ms
from astropy.io import fits
from scipy.spatial import ckdtree
from scipy.interpolate import griddata

#params
config = {}
config['filters'] = 'GRIZ'
config['mu'] = 1.1
config['flux_cut'] = 40.
config['deep_file'] = '/Users/Christina/DES/data/y1a1_gold_dfull_cosmos.fits'
config['balrog_files'] = '/Users/Christina/DES/data/balrog_sva1_tab*_TRUTH_zp_corr_fluxes.fits'
config['sim_file_format'] = '/Users/Christina/DES/data/balrog/sva1/balrog_sva1_auto_tab{}_SIM.fits'
config['deep_flux_column'] = 'FLUX_AUTO_{}'
config['deep_size_column'] = 'FLUX_RADIUS_I'
config['balrog_flux_column'] = 'FLUX_NOISELESS_{}'
config['balrog_size_column'] = 'HALFLIGHTRADIUS_0'


def createdeep(deep, mu, fluxcol='FLUX_AUTO_{}', sizecol='FLUX_RADIUS_I',
               fluxcut=40., filters='GRIZ', sizes=True, logs=True):
    #cut out low fluxes and unrealistic sizes
    flux_mask = (deep[fluxcol.format('G')]>fluxcut) & (deep[fluxcol.format('R')]>fluxcut) & (deep[fluxcol.format('I')]>fluxcut) & (deep[fluxcol.format('Z')]>fluxcut)
    size_mask = (deep[sizecol]>0.)
    
    #separate galaxies & stars using MODEST_CLASS
    gals = deep[(deep['MODEST_CLASS']==1) & flux_mask & size_mask]
    stars = deep[deep['MODEST_CLASS']==2]
    
    if logs:
        #flux vectors
        gal_vec = [np.log10(gals[fluxcol.format(f)]) for f in filters]

        #magnified flux vectors
        new_gal_vec = gal_vec + np.log10(mu)
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
    matches['magnified match radius'] = new_d
        
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
    #dn/dmu = [N_det(mu_G) - N_det(1)] / (mu_G - 1)
    h = np.histogram(table[column][mask], bins=20)
    x_interp = np.array([np.mean(h[1][i-1:i+1]) for i in range(1,len(h[1]))])
    b, c = np.polyfit(x_interp, np.log10(h[0]), 1)
    slope = 2.5 * b
    return slope

def main(args):
    #find deep objects in true redshift range
    #find matches in balrog truth
    #find magnified matches in balrog truth
    #count detections in SIM
    #we have types for SIM detections

    #if using in code, read params
    if len(args) > 1:
        params = args
    else:
        params = config
                
    #open deep data file and create unlensed & lensed data vectors
    #DataSet(zprob_tab, data_file, id_col, output, z_col=None)
    deep = fits.open(params['deep_file'])[1].data
    gal_data, new_gal_data = createdeep(deep, params['mu'],
                                        params['deep_flux_column'],
                                        params['deep_size_column'],
                                        params['flux_cut'], params['filters'])    
    

    sys.stderr.write('Using {} objects from deep data.\n'.format(len(gal_data)))
    
    #indices of deep data
    ids = list(range(len(gal_data)))

    #cycle over all balrog tables
    balrog_matches = {}
    tabnums = []
    for tab in glob.glob(params['balrog_files'])[:2]:
        itab = tab.find('tab')
        tabnum = tab[itab+3:itab+5]
        tabnums.append(tabnum)
        sys.stderr.write('Working on balrog table {}...'.format(tabnum))

        #open tables, get object flux/size data vectors
        brog = fits.open(tab)[1].data
        br_data = createbalrog(brog, params['balrog_flux_column'],
                               params['balrog_size_column'],
                               params['flux_cut'], params['filters'])
        these_matches = findmatches(br_data, gal_data, new_gal_data)

        #get balrog indices to match to sim tables later
        these_matches['balrog index of match'] = brog['BALROG_INDEX'][these_matches['index of match']]
        these_matches['balrog index of magnified match'] = brog['BALROG_INDEX'][these_matches['index of magnified match']]

        balrog_matches[tabnum] = these_matches
        sys.stderr.write('Done.\n')

    #find closest match for each object among all tables
    match_radii = zip(*[balrog_matches[key]['match radius'] for key in balrog_matches.keys()])
    new_match_radii = zip(*[balrog_matches[key]['magnified match radius'] for key in balrog_matches.keys()])
    key_arg = np.argmin(match_radii, axis=1)
    new_key_arg = np.argmin(new_match_radii, axis=1)
    best_tabs = np.array(balrog_matches.keys())[key_arg]
    best_new_tabs = np.array(balrog_matches.keys())[new_key_arg]

    #save numbers of detections
    detections, new_detections = 0, 0
    for tabnum in tabnums:
        #open sim catalog
        sim = fits.open(params['sim_file_format'].format(tabnum))[1].data

        #focus on this table
        this_set = np.where(best_tabs==tabnum)
        this_new_set = np.where(best_new_tabs==tabnum)
        
        #count objects that are found (detected) in sim catalog
        sys.stderr.write('Working on sim table {}...'.format(tabnum))
        sys.stderr.write('checking for detections...')

        #detected objects whose truth fluxes match deep set
        truth_matches = balrog_matches[tabnum]['balrog index of match'][this_set]
        print len(truth_matches)
        detections += len(set(sim['BALROG_INDEX']).intersection(truth_matches))
        
        #detected objects whose truth fluxes match magnified deep set
        new_truth_matches = balrog_matches[tabnum]['balrog index of magnified match'][this_new_set]
        new_detections += len(set(sim['BALROG_INDEX']).intersection(new_truth_matches))

        sys.stderr.write('Done.\n')
        
    #report information
    print "Detected original matches: {}".format(detections)
    print "Detected magnified matches: {}".format(new_detections)

    k = float(new_detections - detections) / (mu - 1.)
    return k
     
if __name__=="__main__":
    main(sys.argv)    
