import numpy as np
import time
import sys
import matplotlib.pyplot as plt
import ConfigParser
import subprocess
import parse
import itertools
from scipy.spatial import ckdtree
from astropy.io import fits
from astropy.table import Table
from dataset import DataSet
from multiprocessing import Pool


#radii - same units as positions
radii = np.logspace(np.log10(0.01), np.log10(0.8), 6)

def chunkcount(ri, lns, data, weights):
    nchunks = 1000
    chunk_size = int(np.ceil(float(len(data))/nchunks))
    chunks = (zip(data['RA'], data['DEC'])[i:i+chunk_size] \
              for i in xrange(0, len(data), chunk_size))
    chunk_indices = (list(range(len(data))[i:i+chunk_size]) \
                     for i in xrange(0, len(data), chunk_size))

    pairs, rpairs = 0, 0
                     
    #query_ball_tree
    query_time = 0
    for chunk, chunk_index in itertools.izip(chunks, chunk_indices):
        start_tree = time.time()
        tree = ckdtree.cKDTree(chunk)
        pairs2 = lns.tree.query_ball_tree(tree, r=radii[ri])
            
        end_tree = time.time()
        query_time += end_tree - start_tree

        if ri!=0:
            pairs1 = lns.tree.query_ball_tree(tree, r=radii[ri-1])
            pairs += np.sum((len(item) for item in pairs2)) - \
                     np.sum((len(item) for item in pairs1))
            indices = [int(i) for i in np.hstack([list(set(i2)-set(i1)) \
                                                  for i1,i2 in zip(pairs1, pairs2)])]
        else:
            pairs  += np.sum((len(item) for item in pairs2))
            indices = [int(i) for i in np.hstack(pairs2)]

    sys.stderr.write('    source queries completed in {}s.\n'.format(query_time))
    new_weights = np.sum(weights[chunk_index][indices])
    return float(pairs), new_weights

def chunkwrapper(args):
    return chunkcount(*args)

def countpairs(src, lns, rnd, rndtype='lens',
               srcweights=None, rndweights=None,
               numthreads=1):
    #radii in increasing order
    annuli = {}

    for ri in range(len(radii)):
        annuli[radii[ri]] = {}
        sys.stderr.write('    working on r={}\n'.format(radii[ri]))

        #multiprocessing      
        annuli[radii[ri]]['srcpairs'], sum_srcweights = multi(ri, lns, src, srcweights, numthreads)
        annuli[radii[ri]]['rndpairs'], sum_rndweights = multi(ri, lns, rnd, rndweights, numthreads)
    
        if srcweights is not None:
            annuli[radii[ri]]['Psrcsum'] = sum_srcweights
        if rndweights is not None:
            annuli[radii[ri]]['Prndsum'] = sum_rndweights
            
    return annuli
    
def multi(ri, lns, dataset, weights, numthreads):
    #multiprocessing
    pool = Pool(processes=numthreads)
    n_per_process = int(np.ceil(len(dataset.data) / float(numthreads)))
    thread_chunks = (dataset.data[i:i+n_per_process] \
                     for i in xrange(0, len(dataset.data), n_per_process))
    weight_thread_chunks = (weights[i:i+n_per_process] \
                            for i in xrange(0, len(weights), n_per_process))

    results = pool.map(chunkwrapper, itertools.izip(itertools.repeat(ri), itertools.repeat(lns),
                                                    thread_chunks, weight_thread_chunks))
    chain_results = np.sum(results, axis=0)

    return chain_results

@profile
def mycorr(sources, lenses, randoms, output,
           random_type='lens', srcweights=None, rndweights=None,
           numthreads=1):
    sys.stdout.write('\n    Sources: {}\n    Lenses: {}\n'.format(len(sources.data),
                                                                  len(lenses.data)))
    sys.stdout.flush()

    if random_type=='lens':
        other_type = 'source'
        ratio = float(len(randoms.data)) / len(lenses.data)
    elif random_type=='source':
        other_type = 'lens'
        ratio = float(len(randoms.data)) / len(sources.data)
        
    #make trees
    start_tree = time.time()
    sources.initTree()
    lenses.initTree()
    randoms.initTree()
    end_tree = time.time()
        
    sys.stderr.write('    Trees created in {}s.\n'.format(end_tree-start_tree))       
    start_q = time.time()
    sys.stderr.write('    Starting queries...\n')
        
    #for each radius, query for all sources around lenses
    annuli = countpairs(sources, lenses, randoms, rndtype=random_type,
                        srcweights=srcweights, rndweights=rndweights,
                        numthreads=numthreads)
    sys.stderr.write('    Done.\n')
    end_q = time.time()
        
    sys.stdout.write('    Time for pair queries: {}s.\n'.format(end_q-start_q))
        
    for k in np.sort(annuli.keys()):
        print '    r={}: {} source-lens pairs'.format(k, annuli[k]['srcpairs'])
        print '          {} {}-random {} pairs'.format(annuli[k]['rndpairs'],
                                                       other_type,
                                                       random_type)
            
    #plot pair counts
    r = np.sort(annuli.keys())
    DD = np.array([annuli[k]['srcpairs'] for k in r])
    DR = np.array([annuli[k]['rndpairs'] for k in r])
    plt.plot(r, DD, 'o-', c='b', markeredgecolor='none', label='sources/lenses') 
    plt.plot(r, DR, 'o-', c='r', markeredgecolor='none', label='sources/randoms')
    plt.xscale('log')
    plt.yscale('log')
    
    #Poisson errorbars
    plt.errorbar(r, DD, yerr=np.sqrt(DD), fmt='o-')
    plt.errorbar(r, DR, yerr=np.sqrt(DR), c='r', fmt='o-')

    #save results to table
    tab = Table()
    tab['N_lens'] = [len(lenses.data)] * len(r)
    tab['N_source'] = [len(sources.data)] * len(r)
    tab['N_random'] = [len(randoms.data)] * len(r)
    tab['R'] = r
    tab['DD'] = DD
    tab['DR'] = DR
        
    if srcweights is not None:
        DD_w = np.array([annuli[k]['Psrcsum'] for k in r])
        tab['DD_w'] = DD_w
        plt.plot(r, DD_w, 'o-',
                 c='c', markeredgecolor='none', label='sources/lenses, weighted')
    if rndweights is not None:
        DR_w = np.array([annuli[k]['Prndsum'] for k in r])
        tab['DR_w'] = DR_w
        plt.plot(r, DR_w, 'o-',
                 c='m', markeredgecolor='none', label='sources/randoms, weighted')
    if (srcweights is not None) and (rndweights is not None):
         w_w = np.array(DD_w)/np.array(DR_w) * ratio - 1
         tab['w_w'] = w_w
        
    output_pairs = output+'_pairs.png'
    plt.xlabel('r (deg)')
    plt.ylabel('pairs')
    plt.grid(which='both')
    plt.legend(loc='best')
    plt.savefig(output_pairs)
    sys.stdout.write('    Pairs figure saved to {}\n'.format(output_pairs))
    plt.close()
    
    #plot cross-correlations
    w = np.array(DD)/np.array(DR) * ratio - 1
    tab['w'] = w
    plt.scatter(r, w, c='b', edgecolor='none')
    plt.xscale('log')
    
    #Poisson errorbars
    plt.errorbar(r, w, yerr=1./np.sqrt(DD), fmt='o-')
    plt.xlabel('r (deg)')
    plt.ylabel('w')
    output_corr = output+'_correlations.png'
    plt.grid(which='both')
    plt.savefig(output_corr)
    sys.stdout.write('    Correlation figure saved to {}\n'.format(output_corr))
    plt.close()
    
    output_fits = output+'_correlations.fits'
    sys.stdout.write('    Writing correlation results to {}.\n'.format(output_fits))
    tab.write(output_fits)

    sys.stdout.flush()
    return tab

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
    command = "/home/ckrawiec/.local/bin/corr2 "+config_file
    print command
    subprocess.call(command.split(' '))
    results = Table.read(params['output']+'_treecorrNN.fits')
    return results
        
def main(params):
    #read files
    sys.stdout.write('Calculating correlations with sources from {}.\n'.format(params['source_file']))
    sys.stdout.write('                              lenses from {}.\n'.format(params['lens_file']))
    sys.stdout.write('                              randoms ({}) from {}.\n'.format(params['random_type'],
                                                                                    params['random_file']))
    sources   = DataSet(params['source_file'], racol='RA', deccol='DEC')
    lenses    = DataSet(params['lens_file'], racol='RA', deccol='DEC')
    randoms  = DataSet(params['random_file'], racol='RA', deccol='DEC')

    if 'source_weight_column' in params:
        srcweights = sources.data[params['source_weight_column']]
        if params['random_type']=='source':
            rndweights = randoms.data[params['source_weight_column']]
        else:
            rndweights = None
    else:
        srcweights = None
        rndweights = None

    corrtab = mycorr(sources, lenses, randoms, params['output'],
                     random_type=params['random_type'],
                     srcweights=srcweights,
                     rndweights=rndweights,
                     numthreads=params['num_threads'])

    #treecorr(params)

    return corrtab
    

if __name__=="__main__":
    params = parse.parseconfigs(sys.argv[1])
    main(params)
    
