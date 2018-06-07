import numpy as np
import time
import sys
import matplotlib.pyplot as plt
import ConfigParser
import subprocess
import parse
from scipy.spatial import ckdtree
from astropy.io import fits
from astropy.table import Table
from dataset import DataSet

#radii - same units as positions
radii = np.logspace(np.log10(0.01), np.log10(0.5), 5)
#np.logspace(-2, -1, 4)

def countpairs(src, lns, rnd, n_chunk = 100):
    #radii in increasing order

    pair_dict, rpair_dict = {}, {}
    
    for ri in range(len(radii)):
        pairs, rpairs = 0, 0
        sys.stderr.write('    working on r={}\n'.format(radii[ri]))
        
        for ci in range(int(np.ceil(len(lns.positions)/float(n_chunk)))):
            start = ci*n_chunk
            end = ci*n_chunk+n_chunk

            #query_ball_tree
            pairs2  = src.tree.query_ball_point(lns.positions[start:end], r=radii[ri])
            rpairs2 = rnd.tree.query_ball_point(src.positions[start:end], r=radii[ri])
            if ri==0:
                pairs  += np.sum(len(np.hstack(pairs2)))
                rpairs += np.sum(len(np.hstack(rpairs2)))
            else:
                pairs1  = src.tree.query_ball_point(lns.positions[start:end],
                                                    r=radii[ri-1])
                rpairs1 = rnd.tree.query_ball_point(src.positions[start:end],
                                                    r=radii[ri-1])

                pairs  += np.sum(len(np.hstack(pairs2))) - np.sum(len(np.hstack(pairs1)))
                rpairs += np.sum(len(np.hstack(rpairs2))) - np.sum(len(np.hstack(rpairs1)))

        pair_dict[radii[ri]] = float(pairs)
        rpair_dict[radii[ri]] = float(rpairs)

    return pair_dict, rpair_dict

def mycorr(sources, lenses, rlenses, output):

        sys.stdout.write('    Sources: {}\nLenses: {}\n'.format(len(sources.data), len(lenses.data)))
        sys.stdout.flush()

        #make trees
        start_tree = time.time()
        sources.tree = ckdtree.cKDTree(sources.positions)
        rlenses.tree = ckdtree.cKDTree(rlenses.positions)
        end_tree = time.time()
        
        sys.stderr.write('    Trees created in {}s.\n'.format(end_tree-start_tree))
        
        start_q = time.time()
        sys.stderr.write('    Starting queries...\n')
        
        #for each radius, query for all sources around lenses
        DD, DR = countpairs(sources, lenses, rlenses)
        sys.stderr.write('    Done.\n')
        end_q = time.time()
        
        sys.stdout.write('    Time for pair queries: {}s.\n'.format(end_q-start_q))
        
        for k in DD.keys():
            print '    r={}: {} source-lens pairs'.format(k, DD[k])
            print '          {} random source-lens pairs'.format(DR[k])
            
        #plot pair counts
        plt.scatter(DD.keys(), DD.values(), c='b', edgecolor='none', label='sources/lenses') 
        plt.scatter(DR.keys(), DR.values(), c='r', edgecolor='none', label='sources/random lenses')
    
        plt.xscale('log')
        plt.yscale('log')
    
        #Poisson errorbars
        plt.errorbar(DD.keys(), DD.values(), yerr=np.sqrt(DD.values()), fmt='o')
        plt.errorbar(DR.keys(), DR.values(), yerr=np.sqrt(DR.values()), c='r', fmt='o')
        
        output_pairs = output+'_pairs.png'
        plt.xlabel('r (deg)')
        plt.ylabel('pairs')
        plt.grid(which='both')
        plt.legend(loc='best')
        plt.savefig(output_pairs)
        sys.stdout.write('    Pairs figure saved to {}\n'.format(output_pairs))
        plt.close()
    
        #plot cross-correlations
        ratio = float(len(rlenses.data)) / len(lenses.data)

        w = [np.array(DD[k])/np.array(DR[k]) * ratio - 1 for k in DD.keys()]
    
        plt.scatter(DD.keys(), w,
                    c='b', edgecolor='none')
        plt.xscale('log')
        
        #Poisson errorbars
        plt.errorbar(DD.keys(), w, yerr=1./np.sqrt(DD.values()), fmt='o')
        
        plt.xlabel('r (deg)')
        plt.ylabel('w')
        
        output_corr = output+'_correlations.png'
        plt.grid(which='both')
        plt.savefig(output_corr)
        sys.stdout.write('    Correlation figure saved to {}\n'.format(output_corr))
        plt.close()
    
        tab = Table()
        tab['R'] = DD.keys()
        tab['w'] = w
        tab['DD'] = DD.values()
        tab['DR'] = [DR[k] for k in DD.keys()]
        output_fits = output+'_mycorr.fits'
        sys.stdout.write('    Writing correlation results to {}.\n'.format(output_fits))
        tab.write(output_fits)

        sys.stdout.flush()

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
    command = "corr2 "+config_file
    print command
    subprocess.call(command.split(' '))
    results = Table.read(params['output']+'_treecorrNN.fits')
    return results
        
def main(config):
    params = parse.parseconfigs(config)

    #read files
    sys.stdout.write('Calculating correlations with sources from {}.\n'.format(params['source_file']))
    sys.stdout.write('                              lenses from {}.\n'.format(params['lens_file']))
    sys.stdout.write('                              random lenses from {}.\n'.format(params['random_lens_file']))
    sources   = DataSet(params['source_file'], racol='RA', deccol='DEC')
    lenses    = DataSet(params['lens_file'], racol='RA', deccol='DEC')
    random_lenses  = DataSet(params['random_lens_file'], racol='RA', deccol='DEC')

    mycorr(sources, lenses, random_lenses, params['output'])

    
    #treecorr(params)
    

if __name__=="__main__":
    main(sys.argv[1])
    
