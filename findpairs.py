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
radii = np.logspace(np.log10(0.01), np.log10(0.05), 5)
#np.logspace(-2, -1, 4)

def countpairs(src, lns, rnd, srcweights=None):
    #radii in increasing order

    annuli = {}
    
    for ri in range(len(radii)):
        annuli[radii[ri]] = {}
        
        pairs, rpairs = 0, 0
        sys.stderr.write('    working on r={}\n'.format(radii[ri]))

        #query_ball_tree
        start_tree = time.time()
        pairs2  = lns.tree.query_ball_tree(src.tree, r=radii[ri])
        rpairs2 = rnd.tree.query_ball_tree(src.tree, r=radii[ri])
        end_tree = time.time()
        sys.stderr.write('    query completed in {}s.\n'.format(end_tree-start_tree))

        p2int = [[int(pp) for pp in p] for p in pairs2]
        rp2int = [[int(rpp) for rpp in rp] for rp in rpairs2]
        if ri==0:
            pairs  += np.sum(len(np.hstack(p2int)))
            rpairs += np.sum(len(np.hstack(rp2int)))

            indices = p2int
            rindices = rp2int
        else:
            pairs1  = lns.tree.query_ball_tree(src.tree, r=radii[ri-1])
            rpairs1 = rnd.tree.query_ball_tree(src.tree, r=radii[ri-1])

            p1int = [[int(pp) for pp in p] for p in pairs1]
            rp1int = [[int(rpp) for rpp in rp] for rp in rpairs1]
            
            pairs  += np.sum(len(np.hstack(p2int))) - np.sum(len(np.hstack(p1int)))
            rpairs += np.sum(len(np.hstack(rp2int))) - np.sum(len(np.hstack(rp1int)))

            indices = [list(set(i2)-set(i1)) for i1,i2 in zip(p1int, p2int)]
            rindices = [list(set(j2)-set(j1)) for j1,j2 in zip(rp1int, rp2int)]
            
        #save pairs as sum P[indices]
        annuli[radii[ri]]['srcpairs'] = float(pairs)
        annuli[radii[ri]]['rndpairs'] = float(rpairs)
        if srcweights is not None:
            annuli[radii[ri]]['Psrcsum'] = np.sum(srcweights[[int(i) for i in np.hstack(indices)]])
            annuli[radii[ri]]['Prndsum'] = np.sum(srcweights[[int(j) for j in np.hstack(rindices)]])

    return annuli

def mycorr(sources, lenses, rlenses, output, srcweights=None):

        sys.stdout.write('    Sources: {}\nLenses: {}\n'.format(len(sources.data), len(lenses.data)))
        sys.stdout.flush()

        #make trees
        start_tree = time.time()
        sources.initTree()
        lenses.initTree()
        rlenses.initTree()
        end_tree = time.time()
        
        sys.stderr.write('    Trees created in {}s.\n'.format(end_tree-start_tree))
        
        start_q = time.time()
        sys.stderr.write('    Starting queries...\n')
        
        #for each radius, query for all sources around lenses
        annuli = countpairs(sources, lenses, rlenses, srcweights)
        sys.stderr.write('    Done.\n')
        end_q = time.time()
        
        sys.stdout.write('    Time for pair queries: {}s.\n'.format(end_q-start_q))
        
        for k in np.sort(annuli.keys()):
            print '    r={}: {} source-lens pairs'.format(k, annuli[k]['srcpairs'])
            print '          {} source-random lens pairs'.format(annuli[k]['rndpairs'])

        
            
        #plot pair counts
        r = np.sort(annuli.keys())
        DD = np.array([annuli[k]['srcpairs'] for k in r])
        DR = np.array([annuli[k]['rndpairs'] for k in r])
        plt.plot(r, DD, 'o-', c='b', markeredgecolor='none', label='sources/lenses') 
        plt.plot(r, DR, 'o-', c='r', markeredgecolor='none', label='sources/randoms')
        plt.xscale('log')
        plt.yscale('log')
    
        #Poisson errorbars
        #plt.errorbar(r, DD, yerr=np.sqrt(DD), fmt='o-')
        #plt.errorbar(r, DR, yerr=np.sqrt(DR), c='r', fmt='o-')

        #save results to table
        tab = Table()
        tab['R'] = r
        tab['DD'] = DD
        tab['DR'] = DR
        
        if srcweights is not None:
            DD_w = np.array([annuli[k]['Psrcsum'] for k in r])
            DR_w = np.array([annuli[k]['Prndsum'] for k in r])
            tab['DD_w'] = DD_w
            tab['DR_w'] = DR_w
            plt.plot(r, DD_w, 'o-',
                     c='c', markeredgecolor='none', label='sources/lenses, weighted')
            plt.plot(r, DR_w, 'o-',
                     c='m', markeredgecolor='none', label='sources/randoms, weighted') 
        
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
        w = np.array(DD)/np.array(DR) * ratio - 1
        tab['w'] = w
        plt.scatter(r, w, c='b', edgecolor='none')
        plt.xscale('log')
        
        #Poisson errorbars
        #plt.errorbar(r, w, yerr=np.sqrt(DD), fmt='o-')
        plt.xlabel('r (deg)')
        plt.ylabel('w')
        output_corr = output+'_correlations.png'
        plt.grid(which='both')
        plt.savefig(output_corr)
        sys.stdout.write('    Correlation figure saved to {}\n'.format(output_corr))
        plt.close()
    
        
        output_fits = output+'_mycorr.fits'
        sys.stdout.write('    Writing correlation results to {}.\n'.format(output_fits))
        tab.write(output_fits)

        sys.stdout.flush()

def treecorr(params, units='degrees'):
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

    mycorr(sources, lenses, random_lenses, params['output'],
           srcweights=np.array([1.]*len(sources.data)))

    
    #treecorr(params)
    

if __name__=="__main__":
    main(sys.argv[1])
    
