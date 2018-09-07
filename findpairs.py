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


def chunkcount(ri, lns, data, radii, weights):
    #source/lens chunks
    if radii[ri]>0.5:
        nchunks = 500
        lens_nchunks = 50
    else:
        nchunks = 2
        lens_nchunks = 2

    #size of each chunk
    chunk_size = int(np.ceil(float(len(data))/nchunks))
    lens_chunk_size = int(np.ceil(float(len(lns.data))/lens_nchunks))

    #loop over lens chunks
    query_time, pairs, weight_sum = 0, 0, 0
    for il in xrange(0, len(lns.data), lens_chunk_size):
        start_query = time.time()
        lns_chunk = ckdtree.cKDTree(zip(lns.data['RA'][il:il+lens_chunk_size],
                                        lns.data['DEC'][il:il+lens_chunk_size]))
    
        #loop over source chunks
        for i in xrange(0, len(data), chunk_size):
            #create tree from positions for this chunk
            chunk_indices = (range(i, i+chunk_size))
            tree = ckdtree.cKDTree(zip(data['RA'][i:i+chunk_size],
                                       data['DEC'][i:i+chunk_size]))

            #query the lens tree for outer radius
            pairs2 = lns_chunk.query_ball_tree(tree, r=radii[ri])
            pairs += np.sum((len(item) for item in pairs2))
 
            weight_sum += np.sum(np.sum(weights[int(d)] for d in item) for item in pairs2)
            del pairs2
            if ri!=0:
                #query the lens tree for inner radius
                pairs1 = lns_chunk.query_ball_tree(tree, r=radii[ri-1])
                del tree

                pairs -= np.sum((len(item) for item in pairs1))
#                sys.stderr.write('    pairs = {}\n'.format(pairs))
                weight_sum -= np.sum(np.sum(weights[int(d)] for d in item) for item in pairs1)
                del pairs1
#(chunk_indices[int(j)] \
 #                          for j in np.hstack(list(set(i2)-set(i1)) \
 #                                             for i1,i2 in itertools.izip(pairs1, pairs2)))
#                weight_sum += np.sum((weights[int(d)] for d in indices))

#            else:
                #if this is smallest radius, use all within "outer"
                #save sum of weights for pairs found in this chunk
#                indices = [chunk_indices[int(j)] for j in np.hstack(pairs2)]
#                weight_sum += np.sum((weights[int(d)] for d in indices))

        query_time += time.time() - start_query
        del lns_chunk

    sys.stderr.write('    queries completed in {}s.\n'.format(query_time))
    sys.stderr.write('    pairs = {}\n'.format(pairs)) 
    return float(pairs), weight_sum

def chunkwrapper(args):
    return chunkcount(*args)

def countpairs(src, lns, rnd, radii, rndtype='lens',
               srcweights=None, rndweights=None,
               numthreads=1):

    #radii in increasing order
    annuli = {}
    #sources first
    for ri in range(len(radii)):
        annuli[radii[ri]] = {}
        sys.stderr.write('    working on r={} (sources)\n'.format(radii[ri]))
        
        #multiprocessing
        if numthreads == 1:
            annuli[radii[ri]]['srcpairs'], sum_srcweights =  chunkcount(ri, lns, src.data, radii, srcweights)
            if srcweights is not None:
                annuli[radii[ri]]['Psrcsum'] = sum_srcweights            

        else:
            from multiprocessing import Pool
            annuli[radii[ri]]['srcpairs'], sum_srcweights = multi(ri, lns, src, srcweights, numthreads)            
            if srcweights is not None:
                annuli[radii[ri]]['Psrcsum'] = sum_srcweights
    del src, srcweights

    #randoms
    for ri in range(len(radii)):
        sys.stderr.write('    working on r={} (randoms)\n'.format(radii[ri]))
        if numthreads == 1:
            annuli[radii[ri]]['rndpairs'], sum_rndweights =  chunkcount(ri, lns, rnd.data, radii, rndweights)
            if rndweights is not None:
                annuli[radii[ri]]['Prndsum'] = sum_rndweights
        else:
            from multiprocessing import Pool
            annuli[radii[ri]]['rndpairs'], sum_rndweights = multi(ri, lns, rnd, rndweights, numthreads)
            if rndweights is not None:
                annuli[radii[ri]]['Prndsum'] = sum_rndweights
    del rnd, rndweights
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

def mycorr(sources, lenses, randoms, output, radii,
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
        
    start_q = time.time()
    sys.stderr.write('    Starting queries...\n')
        
    #for each radius, query for all sources around lenses
    annuli = countpairs(sources, lenses, randoms, radii, rndtype=random_type,
                        srcweights=srcweights, rndweights=rndweights,
                        numthreads=numthreads)
    sys.stderr.write('    Done.\n')
    end_q = time.time()
        
    sys.stdout.write('    Total time for pair queries: {}s.\n'.format(end_q-start_q))
        
    for k in np.sort(annuli.keys()):
        print '    r={}: {} source-lens pairs'.format(k, annuli[k]['srcpairs'])
        print '          {} {}-random {} pairs'.format(annuli[k]['rndpairs'],
                                                       other_type,
                                                       random_type)
    tab = plotnsave(sources, lenses, randoms, annuli, output,
                    random_type, srcweights, rndweights)
    return tab

def plotnsave(sources, lenses, randoms, annuli, output,
              random_type='lens', srcweights=None, rndweights=None):
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
        
    if random_type=='lens':
        other_type = 'source'
        ratio = float(len(randoms.data)) / len(lenses.data)
    elif random_type=='source':
        other_type = 'lens'
        ratio = float(len(randoms.data)) / len(sources.data)

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
    randoms   = DataSet(params['random_file'], racol='RA', deccol='DEC')

    if 'source_weight_column' in params:
        srcweights = sources.data[params['source_weight_column']]
        if params['random_type']=='source':
            rndweights = randoms.data[params['source_weight_column']]
        else:
            rndweights = None
    else:
        srcweights = None
        rndweights = None

    #radii - same units as positions
    r_min = params['r_min']
    r_max = params['r_max']
    r_bins = params['r_bins']
    radii = np.logspace(np.log10(r_min), np.log10(r_max), r_bins)

    corrtab = mycorr(sources, lenses, randoms, params['output'], radii,
                     random_type=params['random_type'],
                     srcweights=srcweights,
                     rndweights=rndweights,
                     numthreads=params['num_threads'])

    #treecorr(params)

    return corrtab
    

if __name__=="__main__":
    params = parse.parseconfigs(sys.argv[1])
    main(params)
    
