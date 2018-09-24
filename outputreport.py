import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.table import Table, join
import pickle
from astropy.io import fits

name       = 'sva1_gold_z4_6bins_5classes_r1_6bins_gr'
ntypes     = 5
zgroups = [[0.001, 1.0],
           [1.0, 2.0],
           [2.0, 3.0],
           [3.0, 4.0],
           [4.0, 9.9]]

Ntarg = 14975254
    
output_dir = '/Users/Christina/DES/magnification/output/'
outfile    = output_dir+name+'_output.pkl'
print "Using output file: ", outfile

randomSums = [3789663.46742,
              245916.197157,
              23336.0368445,
              2300.94951451,
              8494.9787729]

targetSums = [3218746.3427,
              135254.417865,
              31907.2587761,
              4727.7705871,
              28700.769912]


def checknsum(name, ntypes):
    for t in range(ntypes):
        print "Type ", t
        targtype = fits.open(output_dir+name+'_targets_type{}.fits'.format(t))[1].data
        brogtype = fits.open(output_dir+name+'_balrogsim_type{}.fits'.format(t))[1].data
        print "    Number of targets: ", len(targtype)
        print "    Number of balrogs: ", len(brogtype)
        print "    Sum of Ps:         "
        for col in targtype.columns.names:
            if col[:2] == 'P[':
                print "        Ptarg{}: {}".format(col, np.sum(targtype[col]))

        for col in brogtype.columns.names:
            if col[:2] == 'P[':
                print "        Pbrog{}: {}".format(col, np.sum(brogtype[col]))

def steradtosqarcsec(omega):
    return omega * (180./np.pi)**2. * 1.296e7
def sqdegreetosqarcsec(omega):
    return omega * 1.296e7

def annulus_area(theta, dtheta):
    #theta, dtheta in radians
    #returns solid angle in arcsec**2.
    Asterad = 2*np.pi * (np.cos(theta-dtheta/2.) - np.cos(theta+dtheta/2.))
    return steradtosqarcsec(Asterad)

def getP(output, ntypes):
    P = [] #[[P00, P01]]
    Pnames = []
    for t in range(ntypes):
        P.append([])
        Pnames.append([])
        P[-1] = [output['P_{}{}'.format(t, t2)] for t2 in range(ntypes)]
        Pnames[-1] = ['P_{}{}'.format(t, t2) for t2 in range(ntypes)]
    print "P matrix format: ", Pnames
    return np.matrix(P)

def getPprime(output, n0, ntypes):
    Pprime = []
    for t in range(ntypes):
        Pprime.append([])
        # * len(targtype)
        targtype = fits.open(output_dir+name+'_targets_type{}.fits'.format(t))[1].data
        Pprime[t] = [np.sum(targtype['P'+str(zgroups[t2])])* len(targtype)/float(n0[t2]) for t2 in range(ntypes)]
    return np.matrix(Pprime)
        
def getk(output, ntypes):
    k = []
    knames = []
    for t in range(ntypes):
        k.append([])
        knames.append([])
        k[-1] = [output['k_{}{}'.format(t, t2)] for t2 in range(ntypes)]
        knames[-1] = ['k_{}{}'.format(t, t2) for t2 in range(ntypes)]
    print "k matrix format: ", knames
    return np.matrix(k)

def getkprime(output, ntypes, inputmu=1.2):
    kprime = []
    kprimenames = []
    for t in range(ntypes):
        kprime.append([])
        kprimenames = []
        kprime[-1] = [(output['nnew_{}{}'.format(t, t2)]-output['nold_{}{}'.format(t, t2)])/(inputmu-1) for t2 in range(ntypes)]
        kprimenames[-1] = ['nnew_{}{} - nold_{}{}'.format(t, t2, t, t2) for t2 in range(ntypes)]
    print "kprime matrix format: ", kprimenames
    return np.matrix(kprime)

def getN0(output, ntypes):
    # first index = type
    N0 = np.array([output['n0_'+str(t)] for t in range(ntypes)])
    return N0

def getcorrs(fileglob):
    print fileglob
    corrs = {}
    for corr_file in glob.glob(fileglob):
        start = corr_file.find('type')
        tkey = corr_file[start+4]
        corrs[tkey] = Table.read(corr_file)
        print "Type {} correlations from {}".format(tkey, corr_file)
    return corrs

def getNw(corrs, ntypes):
    n = []
    # first index = type
    for t in range(ntypes):
        n.append(corrs[str(t)]['DD_w'])
    return np.array(n)

def getrand(corrs, ntypes):
    n = []
    for t in range(ntypes):
        n.append(corrs[str(t)]['DR_w'])
    return np.array(n)

def getareas(corrs):
    areas = []
    
    rs = range(len(corrs['0']['R']))
    for ri in rs:
        if ri==0:
            mid = corrs['0']['R'][ri]/2.
            dist = corrs['0']['R'][ri]
        else:
            mid = (corrs['0']['R'][ri] + corrs['0']['R'][ri-1])/2.
            dist = corrs['0']['R'][ri] - corrs['0']['R'][ri-1]
        areas.append(annulus_area(mid, dist)) 
    return np.sort(corrs['0']['R']), np.array(areas)

def err(nd, nr, DD, DR):
    a = DD / DR**2. #(nr/(nd*DR))**2. * DD
    b = DR * (DD/DR**2.)**2.#((nr*DD)/(nd * DR**2.))**2. * DR
    c = (DD/(nd*DR))**2. * nr
    d = ((nr*DD)/(nd**2. * DR))**2. * nd
    return np.sqrt(a+b+c+d)

def plotcorrs(corrs, n, n0, nrand, areas, ntypes, r, k):
    plt.figure(figsize=(10,6))

    rmed = getrpoints(r)

    #w with Psum and no k weighting
    #total targets and randoms
    tot_n_random = [float(corrs[str(t)]['N_random'][0]) for t in range(ntypes)]
    tot_n_tot    = [float(corrs[str(t)]['N_source'][0]) for t in range(ntypes)]
    pairs        = [np.array(corrs[str(t)]['DD_w']) for t in range(ntypes)]
    rpairs       = [np.array(corrs[str(t)]['DR_w']) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_w_nokweight')

    #w with Psum and k weighting
    tot_n_random = [float(corrs[str(t)]['N_random'][0]) for t in range(ntypes)]
    tot_n_tot    = [float(corrs[str(t)]['N_source'][0]) for t in range(ntypes)]
    pairs = [np.array(corrs[str(t)]['DD_w']) for t in range(ntypes)]
    rpairs = [np.array(corrs[str(t)]['DR_w']) for t in range(ntypes)]
    kt = [np.sum([k[t,t2] for t2 in range(ntypes)]) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_w_kweight', kweight=kt)

    #w with N and no k weighting
    tot_n_random = [float(corrs[str(t)]['N_random'][0]) for t in range(ntypes)]
    tot_n_tot    = [float(corrs[str(t)]['N_source'][0]) for t in range(ntypes)]
    pairs        = [np.array(corrs[str(t)]['DD']) for t in range(ntypes)]
    rpairs       = [np.array(corrs[str(t)]['DR']) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_nokweight')

    #w with N and k weighting
    tot_n_random = [float(corrs[str(t)]['N_random'][0]) for t in range(ntypes)]
    tot_n_tot    = [float(corrs[str(t)]['N_source'][0]) for t in range(ntypes)]
    pairs        = [np.array(corrs[str(t)]['DD']) for t in range(ntypes)]
    rpairs       = [np.array(corrs[str(t)]['DR']) for t in range(ntypes)]
    kt           = [np.sum([k[t,t2] for t2 in range(ntypes)]) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_kweight', kweight=kt)

    #w with SumPs as totals and no k weighting
    tot_n_random = [randomSums[t] for t in range(ntypes)]
    tot_n_tot    = [targetSums[t] for t in range(ntypes)]
    pairs        = [np.array(corrs[str(t)]['DD_w']) for t in range(ntypes)]
    rpairs       = [np.array(corrs[str(t)]['DR_w']) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_w_sumtots_nokweight')

    #w with SumPs as totals and no k weighting
    tot_n_random = [randomSums[t] for t in range(ntypes)]
    tot_n_tot    = [targetSums[t] for t in range(ntypes)]
    pairs        = [np.array(corrs[str(t)]['DD_w']) for t in range(ntypes)]
    rpairs       = [np.array(corrs[str(t)]['DR_w']) for t in range(ntypes)]
    kt           = [np.sum([k[t,t2] for t2 in range(ntypes)]) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_w_sumtots_kweight', kweight=kt)

    #w with SumPs as totals and no k weighting
    tot_n_random = [np.sum([randomSums[t] for t in range(ntypes)])] * ntypes
    tot_n_tot    = [np.sum([targetSums[t] for t in range(ntypes)])] * ntypes
    pairs        = [np.array(corrs[str(t)]['DD_w']) for t in range(ntypes)]
    rpairs       = [np.array(corrs[str(t)]['DR_w']) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_w_totsumtots_nokweight')

    #w with SumPs as totals and no k weighting
    tot_n_random = [np.sum([randomSums[t] for t in range(ntypes)])] * ntypes
    tot_n_tot    = [np.sum([targetSums[t] for t in range(ntypes)])] * ntypes
    pairs        = [np.array(corrs[str(t)]['DD_w']) for t in range(ntypes)]
    rpairs       = [np.array(corrs[str(t)]['DR_w']) for t in range(ntypes)]
    kt           = [np.sum([k[t,t2] for t2 in range(ntypes)]) for t in range(ntypes)]
    plotw(tot_n_tot, tot_n_random, pairs, rpairs, rmed, 'w_w_totsumtots_kweight', kweight=kt)
    
def plotw(totn, totr, pairs, rpairs, rx, pname, kweight=None):
    for t in range(ntypes):        
        ratio = float(totr[t]) / totn[t]
        print "Type ", t
        print "Nrandom = ", totr[t]
        print "Ntarget = ", totn[t]
        print "  ratio = ", ratio

        #correlation function  
        w = pairs[t]/rpairs[t] * ratio - 1.
        if kweight is not None:
            w *= kweight[t]
        plt.plot(rx[1:], [0.]*len(rx[1:]), '--', c='k')
        yerr = err(totn[t], totr[t], pairs[t], rpairs[t])
        plt.errorbar(rx[1:], w[1:], yerr=yerr[1:], fmt='o-', label=str(t))
                
    plt.legend(loc='lower right')   
    plt.ylabel('w')
    plt.xlabel('R (deg)')
    plt.xlim(0.003, 1.0)
    plt.xscale('log')
    plt.savefig('/Users/Christina/git/thesis/{}.png'.format(pname))
    plt.close()
  

def getmu(corrs, n, P, n0, k, ntypes, output):
    mu = []
    n_types = range(ntypes)
    for ri in range(len(corrs['0']['R'])):
        n_meas   = [item[ri] for item in n]
        n_unlens = [item[ri] for item in n0]
    
        #k' = dn/dmu vs. k = dlogn/dmu
        #k_prime = getkprime(output, ntypes)
       
        Y = [n_meas[t] - np.sum([P[t,t2] * n_unlens[t2] for t2 in n_types[1:]]) for t in n_types]
        #print "Y = ", Y
        
        X = [np.array([n_unlens[0]*P[t,0]] + [k[t,t2]*n_unlens[t2] for t2 in n_types[1:]]) for t in n_types]
        X = np.matrix(X)
        #print "X = ", X
        
        ans = np.linalg.lstsq(X, Y) #X*ans = Y
        print "answer({}) = {}".format(corrs['0']['R'][ri], ans[0])
        mu.append(np.array(ans[0][1:])-1)
        
    return mu

def getrpoints(r):
    rp = [0]+list(r)
    rmed = [(rp[i]+rp[i+1])/2. for i in range(len(rp)-1)]
    return rmed

def plotcoverage():
    for t in range(ntypes):
        targtype = fits.open(output_dir+name+'_targets_type{}.fits'.format(t))[1].data
        plt.plot(targtype['RA'], targtype['DEC'], 'o',
                 markeredgecolor='none', markersize=2.5, label='targets')
        brogtype = fits.open(output_dir+name+'_balrogsim_type{}.fits'.format(t))[1].data
        plt.plot(brogtype['RA'], brogtype['DEC'], 'o',
                 markeredgecolor='none', markersize=2.5, label='balrog')
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.savefig(output_dir+name+'_type{}positions.png'.format(t))
        plt.close()

def main():
    f = open(outfile, 'r')
    output = pickle.load(f)

    #plotcoverage()
    
    #total number of each type (sumP)
    n0 = getN0(output, ntypes)#/sqdegreetosqarcsec(135.)
    print "Overall number (Sum(P))"#(N/135sq.deg.)"
    for t in range(ntypes):
        print "    Type {}: {}".format(t, n0[t])

    #checknsum(name, ntypes)
    
    corrs = getcorrs(output_dir+name+'_type*_correlations.fits')
    r, areas = getareas(corrs)

    n     = getNw(corrs, ntypes)
    nrand = getrand(corrs, ntypes)
    k     = getk(output, ntypes)
    Pprime= getPprime(output, n0, ntypes)
    P     = getP(output, ntypes)

    #get mu, choose with P to use
    mu    = getmu(corrs, n, Pprime, nrand, k, ntypes, output)

    plotcorrs(corrs, n, n0, nrand, areas, ntypes, r, k)
    
    rmed = getrpoints(r)
    
    #plot magnification
    plt.figure(figsize=(8,8))
    muri = [[mu[ri][t] for ri in range(len(r))] for t in range(ntypes-1)]
    for t in range(ntypes-1):
        plt.plot(rmed[1:], muri[t][1:], 'o-', label=str(t+1))
    plt.xlabel('R ($\degree$)')
    plt.ylabel('$\mu$')
    plt.xscale('log')
    #plt.ylim(0., 10.)
    plt.yscale('log')
    plt.legend()
    plt.savefig('/Users/Christina/git/thesis/mu.png')
    plt.close()
#    print 'mu=', mu

    print "Pprime = ", Pprime
    
if __name__=="__main__":
    main()
