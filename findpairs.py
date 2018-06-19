import numpy as np
#import healpy as hp
import time
import sys
import matplotlib.pyplot as plt
import ConfigParser
import subprocess
from scipy.spatial import ckdtree
from astropy.io import fits
from astropy.table import Table

#randoms & MASKS
#sources & redmagic randoms

#radii - same units as positions
radii = np.logspace(np.log10(0.01), np.log10(0.8), 5)
#np.logspace(-2, -1, 4)

def parseconfig(config_file):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    params = {}

    #files
    params['source_file'] = config.get('I/O','source_file')
    params['lens_file']   = config.get('I/O','lens_file')
    params['source_random_file'] = config.get('I/O','source_random_file')
    params['lens_random_file']   = config.get('I/O','lens_random_file')
    params['lens_weight_file']   = config.get('I/O','lens_weight_file')
    params['output'] = config.get('I/O','output')

    #data
    params['source_ra']   = config.get('Columns', 'source_ra')
    params['source_dec']  = config.get('Columns', 'source_dec')
    params['lens_ra']     = config.get('Columns', 'lens_ra')
    params['lens_dec']    = config.get('Columns', 'lens_dec')

    #randoms
    params['source_rand_ra']     = config.get('Columns', 'source_random_ra')
    params['source_rand_dec']    = config.get('Columns', 'source_random_dec')
    params['lens_rand_ra']       = config.get('Columns', 'lens_random_ra')
    params['lens_rand_dec']      = config.get('Columns', 'lens_random_dec')

    return params

def countpairs(src, lns, rnd, n_chunk = 100):
    #radii in increasing order

    pair_dict, rpair_dict = {}, {}
    
    for ri in range(len(radii)):
        pairs, rpairs = 0, 0
        sys.stderr.write('    working on r={}\n'.format(radii[ri]))
        
        for ci in range(int(np.ceil(len(lns.positions)/float(n_chunk)))):
            start = ci*n_chunk
            end = ci*n_chunk+n_chunk
            
            pairs2  = src.tree.query_ball_point(lns.positions[start:end], r=radii[ri])
            rpairs2 = rnd.tree.query_ball_point(lns.positions[start:end], r=radii[ri])
            if ri==0:
                pairs  += np.sum(len(np.hstack(pairs2)))
                rpairs += np.sum(len(np.hstack(rpairs2)))
            else:
                pairs1  = src.tree.query_ball_point(lns.positions[start:end],
                                                    r=radii[ri-1])
                rpairs1 = rnd.tree.query_ball_point(lns.positions[start:end],
                                                    r=radii[ri-1])

                pairs  += np.sum(len(np.hstack(pairs2))) - np.sum(len(np.hstack(pairs1)))
                rpairs += np.sum(len(np.hstack(rpairs2))) - np.sum(len(np.hstack(rpairs1)))

        pair_dict[radii[ri]] = float(pairs)
        rpair_dict[radii[ri]] = float(rpairs)

    return pair_dict, rpair_dict

def mycorr(params):
        #read files
        start_d = time.time()
        sources   = DataSet(params['source_file'], params['source_ra'], params['source_dec'])
        lenses    = DataSet(params['lens_file'], params['lens_ra'], params['lens_dec'])
        rsources  = DataSet(params['source_random_file'], params['source_rand_ra'], params['source_rand_dec'])
        
        end_d = time.time()
    
        sys.stderr.write('Data initialized in {}s\n'.format(end_d-start_d))
        
        sys.stderr.write('Sources: {}\nLenses: {}\n'.format(len(sources.data), len(lenses.data)))

        #make trees
        start_tree = time.time()
        sources.initTree()
        rsources.initTree()
        end_tree = time.time()
        
        sys.stderr.write('Trees created in {}s\n'.format(end_tree-start_tree))
        
        start_q = time.time()
        
        sys.stderr.write('Starting queries...\n')
        
        #for each radius, query for all sources around lenses
        DD, DR = countpairs(sources, lenses, rsources)
        
        end_q = time.time()
        
        sys.stderr.write('Time for queries: {}s\n'.format(end_q-start_q))
        
        for k in DD.keys():
            print 'r={}: {} source-lens pairs'.format(k, DD[k])
            print '      {} random source-lens pairs'.format(DR[k])
            
        #plot pair counts
        plt.scatter(DD.keys(), DD.values(), c='b', edgecolor='none', label='sources') 
        plt.scatter(DR.keys(), DR.values(), c='r', edgecolor='none', label='random sources')
    
        plt.xscale('log')
        plt.yscale('log')
    
        #Poisson errorbars
        plt.errorbar(DD.keys(), DD.values(), yerr=np.sqrt(DD.values()), fmt='o')
        plt.errorbar(DR.keys(), DR.values(), yerr=np.sqrt(DR.values()), c='r', fmt='o')
        
        output = params['output']+'_pairs.png'
        plt.xlabel('r (deg)')
        plt.ylabel('pairs')
        plt.grid(which='both')
        plt.legend(loc='best')
        plt.savefig(output)
        sys.stderr.write('Figure saved to {}\n'.format(output))
        plt.close()
    
        #plot cross-correlations
        ratio = float(len(rsources.data)) / len(sources.data)

        w = [np.array(DD[k])/np.array(DR[k]) * ratio - 1 for k in DD.keys()]
    
        plt.scatter(DD.keys(), w,
                    c='b', edgecolor='none')
        plt.xscale('log')
        
        #Poisson errorbars
        plt.errorbar(DD.keys(), w, yerr=1./np.sqrt(DD.values()), fmt='o')
        
        plt.xlabel('r (deg)')
        plt.ylabel('w')
        
        output = params['output']+'_correlations.png'
        plt.grid(which='both')
        plt.savefig(output)
        sys.stderr.write('Figure saved to {}\n'.format(output))
        plt.close()
    
        tab = Table()
        tab['R'] = DD.keys()
        tab['w'] = w
        tab['DD'] = DD.values()
        tab['DR'] = [DR[k] for k in DD.keys()]
        tab.write(params['output']+'.fits')

def treecorr(params, units='degrees'):
    #without command line?
    #command = corr2 config
    #source/lens ra/dec column names must match
    config_file = params['output']+"_treecorr.params"
    f = open(config_file, 'w')
    f.write("""file_name = {}
file_name2 = {}
rand_file_name = {}
rand_file_name2 = {}
file_type = FITS
ra_col = {}
dec_col = {}
ra_units = {}
dec_units = {}
min_sep = 0.006
max_sep = 0.08
nbins = 5
sep_units = {}
nn_file_name = {}
verbose = 3
log_file = {}""".format(params['source_file'],
                      params['lens_file'],
                      params['random_source_file'],
                      params['random_lens_file'],
                      params['source_ra'], params['source_dec'],
                      units, units, units,
                      params['output']+'_treecorrNN.fits',
                      params['output']+'_treecorrNN.log'))
    f.close()

    #run command
    command = "/home/ckrawiec/.local/bin/corr2 "+config_file
    print command
    subprocess.call(command.split(' '))
    results = Table.read(params['output']+'_treecorrNN.fits')
    return results

class DataSet:
    def __init__(self, data_file, x_col, y_col, weight_file=None):
        self.data = fits.open(data_file)[1].data
        self.positions = np.array(zip(self.data[x_col] * np.cos(self.data[y_col]*np.pi/180.),
                                      self.data[y_col]))

        if weight_file:
            self.weights = self.getFrac(weight_file, x_col, y_col)
        else:
            self.weights = np.ones(len(self.data))

    def initTree(self):
        self.tree = ckdtree.cKDTree(self.positions)

    def getFrac(self, tab_file, ra_col, dec_col):
        hpmap = hp.read_map(tab_file, nest=True)
        nside = hp.npix2nside(hpmap.size)
    
        theta = (90.0 - self.data[dec_col])*np.pi/180.
        phi = self.data[ra_col]*np.pi/180.
        pix = hp.ang2pix(nside, theta, phi, nest=True)
        
        return hpmap[pix]
        
def main():
    params = parseconfig('findpairs.config')

    ###Use real distances
    ###might make slower, can map/change coordinates

    mycorr(params)

    treecorr(params)
    
        

if __name__=="__main__":
    main()
    
