import numpy as np
import healpy as hp
import time
import sys
import matplotlib.pyplot as plt
import ConfigParser
from scipy.spatial import ckdtree
from astropy.io import fits
from astropy.table import Table

#randoms & MASKS
#sources & redmagic randoms

#radii - same units as positions
radii = np.logspace(-2, 0, 6)

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

def countpairs(src, lns, rnd, n_chunk = 500):
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
                pairs  += np.sum(np.hstack([lns.weights[ip]*len(pairs2[ip]) for ip in range(len(pairs2))]))
                rpairs += np.sum(np.hstack([lns.weights[ip]*len(rpairs2[ip]) for ip in range(len(rpairs2))]))
            else:
                pairs1  = src.tree.query_ball_point(lns.positions[start:end],
                                                    r=radii[ri-1])
                rpairs1 = rnd.tree.query_ball_point(lns.positions[start:end],
                                                    r=radii[ri-1])

                pairs  += np.sum(np.hstack([lns.weights[ip]*len(pairs2[ip]) for ip in range(len(pairs2))]))- np.sum(np.hstack([lns.weights[ip]*len(pairs1[ip]) for ip in range(len(pairs1))]))
                rpairs += np.sum(np.hstack([lns.weights[ip]*len(rpairs2[ip]) for ip in range(len(rpairs2))]))- np.sum(np.hstack([lns.weights[ip]*len(rpairs1[ip]) for ip in range(len(rpairs1))]))

        pair_dict[radii[ri]] = pairs
        rpair_dict[radii[ri]] = rpairs

    return pair_dict, rpair_dict

class DataSet:
    def __init__(self, data_file, x_col, y_col, weight_file=None):
        self.data = fits.open(data_file)[1].data
        self.positions = np.array(zip(self.data[x_col],
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


    #read files
    start_d = time.time()
    sources   = DataSet(params['source_file'], params['source_ra'], params['source_dec'])
    lenses    = DataSet(params['lens_file'], params['lens_ra'], params['lens_dec'])
    rsources  = DataSet(params['source_random_file'], params['source_rand_ra'], params['source_rand_dec'])
    #if len(params['lens_weight_file']) > 0:
    #    rlenses   = DataSet(params['lens_file'], params['lens_ra'], params['lens_dec'],
    #                        weight_file=params['lens_weight_file'])
    end_d = time.time()

    sys.stderr.write('Data initialized in {}s\n'.format(end_d-start_d))
    
    sys.stderr.write('Sources: {}\nLenses: {}\n'.format(len(sources.data), len(lenses.data)))
    #sys.stderr.write('Source Randoms: {}\nLens Randoms: {}\n'.format(len(rsources.data), len(rlenses.data)))
    
        
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
    #rDD, rDR = countpairs(sources, rlenses, rsources)
    
    end_q = time.time()
    
    sys.stderr.write('Time for queries: {}s\n'.format(end_q-start_q))
    
    for k in DD.keys():
        print 'r={}: {} source-lens pairs'.format(k, DD[k])
        print '      {} random source-lens pairs'.format(DR[k])
        #print '      {} source-random lens pairs'.format(k, DD[k])
        #print '      {} random source-random lens pairs'.format(DR[k])


    #plot pair counts
    plt.scatter(DD.keys(), DD.values(), c='b', edgecolor='none', label='sources around lenses')
    plt.scatter(DR.keys(), DR.values(), c='r', edgecolor='none', label='source randoms around lenses')
    #plt.scatter(rDD.keys(), rDD.values(), c='c', edgecolor='none', label='sources around lens randoms')
    #plt.scatter(rDR.keys(), rDR.values(), c='m', edgecolor='none', label='source randoms around lens randoms')
    
    plt.xlabel('r (deg)')
    plt.ylabel('pairs')
    plt.xscale('log')
    output = params['output']+'_pairs.png'
    plt.grid()
    plt.legend(loc='best')
    plt.savefig(output)
    sys.stderr.write('Figure saved to {}\n'.format(output))
    plt.close()

    #plot cross-correlations
    plt.scatter(DD.keys(), np.array(DD.values())/np.array(DR.values()) - 1,
                c='b', edgecolor='none', label='lenses')
    #plt.scatter(rDD.keys(), np.array(rDD.values())/np.array(rDR.values()) - 1,
    #            c='g', edgecolor='none', label='lens randoms')
    
    plt.xlabel('r (deg)')
    plt.ylabel('w')
    plt.xscale('log')
    output = params['output']+'_correlations.png'
    plt.grid()
    plt.legend(loc='best')
    plt.savefig(output)
    sys.stderr.write('Figure saved to {}\n'.format(output))
    plt.close()

    
if __name__=="__main__":
    main()
