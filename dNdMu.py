import matplotlib.pyplot as plt
import numpy as np
import sys
import time
import glob
from astropy.io import fits
from scipy.spatial import ckdtree
from scipy.interpolate import griddata
from measuresignal import DataSet

test = True
sizes = False
logs = False
filters = 'GRIZ'
gridfile = '/Users/Christina/DES/data/balrog/sva1/balrog_tab01_avg_star_fluxradiusi_0.1deg.fits'

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
config['redshifts'] = [0.0, 9.9]

def createVecs(self, fluxcut, sizecut, modestclass=None,
               zrange=[], calchlr=False,
               magnify=False, mu=1.):
        #cut out low fluxes and unrealistic sizes
        new_data = self.data
    
        for filt in 'GRIZ':
            new_data =  new_data[new_data[self.fluxcol.format(filt)]>fluxcut]
            
        if logs:
            #flux vectors
            data_vecs = np.array(zip(*[np.log10(new_data[self.fluxcol.format(f)]) for f in filters]))
        else:
            data_vecs = np.array(zip(*[new_data[self.fluxcol.format(f)] for f in filters]))

        if sizes:
            mask = np.where(new_data[self.sizecol]>sizecut)
            new_data = new_data[mask]
            data_vecs = data_vecs[mask]
            if calchlr:
                stars = new_data[new_data['MODEST_CLASS']==2]
                #get sizes
                hlr = gethlr(self, new_data, stars)
                #trim data to where half light radius makes sense
                mask = np.where(hlr>0.)
                data_vecs = data_vecs[mask]
                new_data = new_data[mask]
                #add sizes to data vector
                unzip = zip(*data_vecs)
                if logs:
                    unzip.append(np.log10(hlr))
                else:
                    unzip.append(hlr)
                data_vecs = zip(*unzip)               
            else:
                unzip = zip(*data_vecs)
                if logs:
                    unzip.append(np.log10(new_data[self.sizecol]))
                else:
                    unzip.append(new_data[self.sizecol])
                data_vecs = zip(*unzip)

        if modestclass:
            mask = new_data['MODEST_CLASS']==modestclass
            new_data = new_data[mask]
            data_vecs = np.array(data_vecs)[mask]

        if magnify:
            #add magnified sizes to magnified data vector
            new_vecs = data_vecs
            
            if logs:
                factor = np.array([(1. + np.log10(mu)/mu)] * len(filters))
            else:
                factor = np.array([mu]*len(filters))
            
            if sizes:
                if logs:
                    factor.append(0.5 * np.log10(mu)/mu + 1.)
                else:
                    factor.append(np.sqrt(mu))
            zrange = np.array(zrange)
            zmask = ((new_data[self.zcol]>zrange.min()) & (new_data[self.zcol]<zrange.max()))
            new_vecs[zmask] = data_vecs[zmask] * factor
            
            data_vecs = new_vecs
        return np.array(data_vecs)


def findmatches(simvecs, datavecs):
    #save match information in a dictionary
    matches = {'radius': [],
               'index': []}

    #create simulation tree
    sys.stderr.write('creating tree...')
    simtree = ckdtree.cKDTree(simvecs)
    
    #query tree for nearest data
    sys.stderr.write('querying...')
    r, index = simtree.query(datavecs)
    matches['radius'] = r
        
    #indices of simulation matches
    matches['index'] = index

    return matches  

def gethlr(self, gals, stars, sizecol='FLUX_RADIUS_I'):
    """Use median star size around galaxies & galaxy size
    from balrog to estimate true half light radius of galaxies"""
  
    to_grid = fits.open(gridfile)[1].data
    
    #make tree of stars
    sys.stderr.write('Finding half-light radii...\n'.format(len(stars)))
    sys.stderr.write('    stars: {}\n'.format(len(stars)))
    sys.stderr.write('    galaxies: {}\n'.format(len(gals)))

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
    med_star_size = np.array([np.median(stars[sizecol][c]) for c in close])
    end = time.time()

    #interpolate true galaxy half light radius using balrog data
    gal_hlr = griddata(zip(to_grid['flux_radius_i'],
                           to_grid['avg_flux_radius_i']),
                       to_grid['hlr'],
                       zip(gals[sizecol], med_star_size))

    if test:
        sys.stderr.write('Making plots...')
        plt.scatter(to_grid['flux_radius_i'], to_grid['hlr'],
                    edgecolor='none', c='g', label='balrog')
        plt.scatter(gals[sizecol], gal_hlr,
                    edgecolor='none', c='b', label='galaxies')
        plt.xlabel('measured radius')
        plt.ylabel('half light radius')
        plt.xscale('log')
        plt.legend(loc='best')
        plt.savefig(self.output+'_hlrtest.png')
        plt.close()
    sys.stderr.write('Done.\n')
    return gal_hlr

def main(args):
    DataSet.createVecs = createVecs
    
    #if using in code, read params
    if len(args) > 1:
        params = args
    else:
        params = config
        sys.stderr.write('Warning: Using built-in parameters.\n')

    #open deep data file and create unlensed & lensed data vectors
    if 'deep_file' in params.keys():
        data = DataSet(params['deep_file'], output='builtin_dNdMu_deep',
                       fluxcol='FLUX_AUTO_{}',
                       sizecol='FLUX_RADIUS_I',
                       zcol='ZMINCHI2')

    elif 'Data' in params.keys():
        data = params['Data']
    else:
        sys.stderr.write('Include "Data" or "deep_file" in config. Exiting\n')
        exit()
        
    #separate galaxies using MODEST_CLASS
    gals = data.createVecs(40., 0., modestclass=1, calchlr=True)
    magnified_gals = data.createVecs(40., 0., modestclass=1, calchlr=True,
                                     zrange=params['redshifts'], mu=params['mu'],
                                     magnify=True)
    h = plt.hist(gals, histtype='step', label='original', bins=50)
    plt.hist(magnified_gals, histtype='step', linestyle='--', label='magnified', bins=h[1])
    plt.xlabel('vector elements')
    plt.legend(loc='best')
    plt.savefig(data.output+'_data_vector_hist.png')
    plt.close()
    
    sys.stderr.write('Using {} galaxies.\n'.format(len(gals.data)))
    
    #indices of galaxy data
    ids = list(range(len(gals.data)))

    #cycle over all truth tables
    truth_matches = {}
    tabnums = []
    for table_name in glob.glob(params['balrog_files']):
        itab = table_name.find('tab')
        tabnum = table_name[itab+3:itab+5]
        tabnums.append(tabnum)
        sys.stderr.write('Working on balrog table {}...'.format(tabnum))

        #open tables, get flux and/or size data vectors
        truth = DataSet(table_name, output='builtin_dNdMu_truth',
                        fluxcol='FLUX_NOISELESS_{}',
                        sizecol='HALFLIGHTRADIUS_0',
                        idcol='BALROG_INDEX')
        truth_vecs = truth.createVecs(40., 0.)

        #find matches
        original_matches  = findmatches(truth_vecs, gals)
        magnified_matches = findmatches(truth_vecs, magnified_gals)

        #get balrog indices to match to sim tables later
        truth_matches[tabnum] = {'original id': truth.data[truth.idcol][original_matches['index']],
                                 'magnified id': truth.data[truth.idcol][magnified_matches['index']],
                                 'original radius': original_matches['radius'],
                                 'magnified radius': magnified_matches['radius']}
        sys.stderr.write('Done.\n')

    #find closest match for each object among all tables
    original_radii  = zip(*[truth_matches[key]['original radius'] for key in truth_matches.keys()])
    magnified_radii = zip(*[truth_matches[key]['magnified radius'] for key in truth_matches.keys()])
    original_arg  = np.argmin(original_radii, axis=1)
    magnified_arg = np.argmin(magnified_radii, axis=1)
    original_best  = np.array(truth_matches.keys())[original_arg]
    magnified_best = np.array(truth_matches.keys())[magnified_arg]

    #save numbers of detections
    original_detections, magnified_detections = 0, 0
    for tabnum in tabnums:
        #open sim catalog
        sim = fits.open(params['sim_file_format'].format(tabnum))[1].data

        #focus on this table number
        original_set  = np.where(original_best==tabnum)
        magnified_set = np.where(magnified_best==tabnum)
        
        #count objects that are found (detected) in sim catalog
        sys.stderr.write('Working on sim table {}...'.format(tabnum))
        sys.stderr.write('checking for detections...')

        #detected objects whose truth fluxes match original set
        truth_ids = truth_matches[tabnum]['original id'][original_set]
        original_detections += len(set(sim['BALROG_INDEX']).intersection(truth_ids))
        
        #detected objects whose truth fluxes match magnified set
        truth_ids = truth_matches[tabnum]['magnified id'][magnified_set]
        magnified_detections += len(set(sim['BALROG_INDEX']).intersection(truth_ids))

        sys.stderr.write('Done.\n')
        
    #report information
    print "Detected original matches: {}".format(original_detections)
    print "Detected magnified matches: {}".format(magnified_detections)

    k = float(magnified_detections - original_detections) / (params['mu'] - 1.)
    return k


if __name__=="__main__":
    main(sys.argv)    