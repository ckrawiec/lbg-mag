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
from scipy.cluster.vq import kmeans, vq

jackknife = True
njackknife = 5

def chunkcount(rinner, router, lns, positions, weights):
    """
    Query tree of lens positions for sources
    in outer radius and inner radius of annulus.
    Save number of pairs and sum of weights.
    """

    #source/lens chunks
    if router>0.5:
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
            tree = ckdtree.cKDTree(positions[i:i+chunk_size])

            #query the lens tree for outer radius
            pairs2 = lns_chunk.query_ball_tree(tree, r=router)

            #add npairs(router)
            pairs += np.sum((len(item) for item in pairs2))
            weight_sum += np.sum(np.sum(weights[int(d)] for d in item) for item in pairs2)
            del pairs2
            
            if rinner>0:
                
                #query the lens tree for inner radius
                pairs1 = lns_chunk.query_ball_tree(tree, r=rinner)
                del tree

                #subtract npairs(rinner)
                pairs      -= np.sum((len(item) for item in pairs1))
                weight_sum -= np.sum(np.sum(weights[int(d)] for d in item) for item in pairs1)
                del pairs1

        query_time += time.time() - start_query
        del lns_chunk

    sys.stderr.write('    queries completed in {}s.\n'.format(query_time))
    sys.stderr.write('    pairs = {}\n'.format(pairs)) 
    return float(pairs), weight_sum

def chunkwrapper(args):
    #wrapper function for multiprocessing
    return chunkcount(*args)

def countpairs(src, lns, rnd, radii, rndtype='lens',
               srcweights=None, rndweights=None,
               numthreads=1):
    """
    Create annuli dictionary and run chunkcount
    for every inner/outer radius pair
    """
    srcpos = zip(src.data['RA'], src.data['DEC'])
    rndpos = zip(rnd.data['RA'], rnd.data['DEC'])
    
    if jackknife==True:
        jkresults = {}
        
        #separate d (src) and r (rnd) into SAME kmeans region, by position
        n_jk = njackknife
        centers, _ = kmeans(srcpos, n_jk)
        k_indices, _ = vq(srcpos, centers)

        #count DD and DR for each radius for each sample
        for k in range(n_jk):
            #radii in increasing order
            annuli = {}

            #sources first
            this_dpos = srcpos[k_indices!=k]
            this_srcweights = ???
            annuli = loopradius(annuli, 'Psrcsum', 'srcpairs'
                                radii, lns, this_dpos, this_srcweights)
            del this_dpos
    
            #again with randoms
            this_rpos = rndpos[k_indices!=k]
            annuli = loopradius(annuli, 'Prndsum', 'rndpairs',
                                radii, lns, this_rpos, rndweights)

            jkresults[k] = {}
            jkresults[k]['Psrcsum'] = [annuli[rad]['Psrcsum'] for rad in radii]
            jkresults[k]['Prndsum'] = [annuli[rad]['Prndsum'] for rad in radii]}

            #get w for this sample
            jkresults[k]['w'] = [w(annuli[rad]['Psrcsum'], annuli[rad]['Psrcsum']) \
                                 for rad in radii]
                                 
        #get jacknife estimate from average of jackknife regions
        jkresults['w'] = jk([jkresults[k]['Psrcsum'] for k in range(n_jk)],
                            [jkresults[k]['Prndsum'] for k in range(n_jk)],
                            w)
        
        #get jackknife variance from the results of jackknife regions
        jkresults['var'] = varjk([jkresults[k]['Psrcsum'] for k in range(n_jk)],
                                 [jkresults[k]['Prndsum'] for k in range(n_jk)],
                                 w)
            
    #answer with full sample
    annuli = {}

    #sources first
    srcpos = zip(src.data['RA'], src.data['DEC'])
    annuli = loopradius(annuli, 'Psrcsum', 'srcpairs'
                        radii, lns, srcpos, srcweights)
    #save space
    del src, srcpos, srcweights

    #again with randoms
    rndpos = zip(rnd.data['RA'], rnd.data['DEC'])
    annuli = loopradius(annuli, 'Prndsum', 'rndpairs',
                        radii, lns, rndpos, rndweights)
    del rnd, rndpos, rndweights

    if jackknife==True:
        return annuli, jkresults
    else:
        return annuli

def loopradius(annuli, weightkey, pairkey,
               radii, lns, positions, weights=None):

    #loop over annuli radii
    for ri in range(len(radii)):
        router = radii[ri]
        if ri==0:
            rinner = 0.
        else:
            rinner = radii[ri-1]

        #create dictionary entry and print radius
        annuli[router] = {}
        sys.stderr.write('    working on r={}\n'.format(router))
        
        #regular or multiprocessing
        if numthreads == 1:
            annuli[router][pairkey], sum_weights =  chunkcount(rinner, router,
                                                               lns, positions, weights)
            if weights is not None:
                annuli[router][weightkey] = sum_weights            

        else:
            from multiprocessing import Pool
            annuli[router][pairkey], sum_weights = multi(rinner, router,
                                                         lns, positionss, weights,
                                                         numthreads)            
            if weights is not None:
                annuli[router][weightkey] = sum_weights

    return annuli

def multi(ri, lns, dataset, weights, numthreads):
    
    #break data into chunks for each multiprocessing thread
    pool = Pool(processes=numthreads)
    n_per_process = int(np.ceil(len(dataset.data) / float(numthreads)))
    thread_chunks = (dataset.data[i:i+n_per_process] \
                     for i in xrange(0, len(dataset.data), n_per_process))
    weight_thread_chunks = (weights[i:i+n_per_process] \
                            for i in xrange(0, len(weights), n_per_process))

    #save and chain results from each thread
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

def w(DD, DR, nd, nr):
    return (float(nr)/nd) * (float(DD)/DR) - 1.

def jk(DD, DR, func):
    n = len(DD)
    nd = np.sum([DD[i] for i in range(len(DD[0]))])
    nr = np.sum([DR[i] for i in range(len(DR[0]))])
    return np.sum([func(DD[j], DR[j], nd[j], nr[j]) \
                   for j in range(n)])/float(n)

def varjk(DD, DR, func):
    n = len(DD)
    nd = np.sum([DD[i] for i in range(len(DD[0]))])
    nr = np.sum([DR[i] for i in range(len(DR[0]))])
    coeff = (n-1)/float(n)
    jk_est = jk(DD, DR, func)
    return coeff * np.sum([(func(DD[j], DR[j], nd[j], nr[j]) - jk_est)**2 for j in range(n)])

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



