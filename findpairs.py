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
radii = np.logspace(np.log10(0.01), np.log10(0.8), 6)

def countpairs(src, lns, rnd, rndtype='lens', srcweights=None, rndweights=None):
    #radii in increasing order
    annuli = {}
    
    for ri in range(len(radii)):
        annuli[radii[ri]] = {}
        
        pairs, rpairs = 0, 0
        sys.stderr.write('    working on r={}\n'.format(radii[ri]))

        #query_ball_tree
        pairs1, pairs2, rpairs1, rpairs2 = [], []

        n_chunks = 4
        chunks = [lns.tree.data[i:i+nchunks] for i in xrange(0, len(lns.tree.data), n_chunks)]
        for chunk in chunks:
            start_tree = time.time()
            lns_chunk = ckdtree.cKDTree(chunk)
            pairs2.append(lns_chunk.query_ball_tree(src.tree, r=radii[ri]))
            if rndtype == 'lens':
                #not right
                #rpairs2.append(rnd.tree.query_ball_tree(src.tree, r=radii[ri]))
            elif rndtype == 'source':
                rpairs2.append(lns_chunk.query_ball_tree(rnd.tree, r=radii[ri]))
            else:
                sys.stderr.write('Warning: random_type not specified correctly, using as lens randoms.\n')
                #rpairs2.append(rnd.tree.query_ball_tree(src.tree, r=radii[ri]))
            end_tree = time.time()
            sys.stderr.write('    query completed in {}s.\n'.format(end_tree-start_tree))

            if ri!=0:
                pairs1.append(lns_chunk.query_ball_tree(src.tree, r=radii[ri-1]))
                if rndtype == 'lens':
                    #not right
                    #rpairs1.append(rnd.tree.query_ball_tree(src.tree, r=radii[ri-1]))
                elif rndtype == 'source':
                    rpairs1.append(lns_chunk.query_ball_tree(rnd.tree, r=radii[ri-1]))
                else:
                    sys.stderr.write('Warning: random_type not specified correctly, using as lens randoms.\n')
                    #rpairs1.append(rnd.tree.query_ball_tree(src.tree, r=radii[ri-1]))

        #list of lists = len(lns.tree.data)
        p2stack = np.hstack(pairs2)
        rp2stack = np.hstack(rpairs2)
        if ri==0:
            pairs  += len(np.hstack(p2stack))
            rpairs += len(np.hstack(rp2stack))
    
            indices = [int(i) for i in np.hstack(p2stack)]
            rindices = [int(j) for j in np.hstack(rp2stack)]

        else:
            p1stack = np.hstack(pairs1)
            rp1stack = np.hstack(rpairs1)

            pairs  += len(np.hstack(p2stack)) - len(np.hstack(p1stack))
            rpairs += len(np.hstack(rp2stack)) - len(np.hstack(rp1stack))

            indices = [int(i) for i in np.hstack([list(set(i2)-set(i1)) for i1,i2 in zip(p1stack, p2stack)])]
            rindices = [int(j) for j in np.hstack([list(set(j2)-set(j1)) for j1,j2 in zip(rp1stack, rp2stack)])]
    
        #save pairs, and as sum P[indices]
        annuli[radii[ri]]['srcpairs'] = float(pairs)
        annuli[radii[ri]]['rndpairs'] = float(rpairs)
        if srcweights is not None:
            annuli[radii[ri]]['Psrcsum'] = np.sum(srcweights[indices])
        if rndweights is not None:
            annuli[radii[ri]]['Prndsum'] = np.sum(rndweights[rindices])

    return annuli

def mycorr(sources, lenses, randoms, output,
           random_type='lens', srcweights=None, rndweights=None):
    sys.stdout.write('\n    Sources: {}\n    Lenses: {}\n'.format(len(sources.data), len(lenses.data)))
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
    annuli = countpairs(sources, lenses, randoms, rndtype=random_type, srcweights=srcweights, rndweights=rndweights)
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
                     rndweights=rndweights)

    #treecorr(params)

    return corrtab
    

if __name__=="__main__":
    params = parse.parseconfigs(sys.argv[1])
    main(params)
    
