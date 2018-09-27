import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib as mpl
from astropy.table import Table, join
import pickle
from astropy.io import fits

mpl.rcParams['font.size']=15
name       = 'sva1_gold_z4_6bins_5classes_r1_6bins_gr_justsim_newPSTAR'
ntypes     = 5
zgroups = [[0.001, 1.0],
           [1.0, 2.0],
           [2.0, 3.0],
           [3.0, 4.0],
           [4.0, 9.9]]
    
classes = ['A','B','C','D','E']

Nlens = 94660.

output_dir = '/Users/Christina/DES/magnification/output/'
outfile    = output_dir+name+'_output.pkl'
print "Using output file: ", outfile

def err(nd, nr, DD, DR):
    D = DD / DR
    N = nr / nd
    w = D * N - 1.
    eout = np.abs(w) * np.sqrt(1./(DD/Nlens) + 1./(DR/Nlens) + 2./Nlens + 1./nd + 1./nr)
    return eout

def getcorrs(fileglob):
    print fileglob
    corrs = {}
    for corr_file in glob.glob(fileglob):
        start = corr_file.find('type')
        tkey = corr_file[start+4]
        corrs[tkey] = Table.read(corr_file)
        print "Type {} correlations from {}".format(tkey, corr_file)
    return corrs

def getk(output, ntypes):
    k = []
    for t in range(ntypes):
        k.append([])
        k[-1] = [output['k_{}{}'.format(t, t2)] for t2 in range(ntypes)]
    return np.matrix(k)

def getP(output, ntypes):
    P = [] #[[P00, P01]]
    for t in range(ntypes):
        P.append([])
        P[-1] = [output['P_{}{}'.format(t, t2)] for t2 in range(ntypes)]
    return np.matrix(P)

def getrpoints(r):
    rp = [0]+list(r)
    rmed = [(rp[i]+rp[i+1])/2. for i in range(len(rp)-1)]
    return np.array(rmed)

def wmucalc(corrs, P, k):
    ans = []
    wans = []
    werr = []
    for ri in range(len(corrs['0']['R'])):
        ratio = np.array([float(corrs[str(t)]['N_source'][0]) / corrs[str(t)]['N_random'][0]
                          for t in range(ntypes)])
        n  = np.array([corrs[str(t)]['DD'][ri] for t in range(ntypes)])
        n0 = np.array([corrs[str(t)]['DR'][ri] for t in range(ntypes)]) * ratio
        print n, n0
        wH = np.array([(n[H] - n0[H]) / n0[H] for H in range(ntypes)])
        fH1 = [P[H,0] * n[0] / (n[H] * P[0,0]) for H in range(ntypes)]
        M = np.matrix([[fH1[H]] + [k[H,G] for G in range(ntypes)[1:]] for H in range(ntypes)])

        ans.append(np.dot(M.I, wH))
        print "wH = ", wH
        print "M = ", M
        print "ans[-1] = ", ans[-1]
        print np.dot(M, ans[-1].T)
        #print "linalg ans = ", np.linalg.lstsq(M, wH)
        wans.append(wH)
        werr.append(err(np.array([float(corrs[str(t)]['N_source'][0]) for t in range(ntypes)]),
                        np.array([float(corrs[str(t)]['N_random'][0]) for t in range(ntypes)]),
                        n, n0))
                          
    return ans, wans, werr

def main():
    f = open(outfile, 'r')
    output = pickle.load(f)

    corrs = getcorrs(output_dir+name+'_type*_correlations.fits')
    r = corrs['0']['R']
    rmed = getrpoints(r)
   
    P     = getP(output, ntypes)
    k     = getk(output, ntypes)
      
    mucalc, wcalc, werr = wmucalc(corrs, P, k)

    #plot magnfication
    orig_offset = rmed*0.05
    for ri in range(len(r)):
        print "r={}, mucalc={}".format(r[ri], mucalc[ri])
    for t in range(1, ntypes):
        mut = np.array([mucalc[ri][0,t] for ri in range(len(r))])
        plt.plot(rmed[1:]+orig_offset[1:]*(t-1), mut[1:], label=t)
    plt.grid()
    plt.legend()
    plt.xscale('log')
    plt.xlim(0.0035, 1.0) 
    plt.savefig('/Users/Christina/git/thesis/{}_mucalc.png'.format(name))
    plt.close()

    #plot w
    for t in range(ntypes):
        wt = np.array([wcalc[ri][t] for ri in range(len(r))])
        wterr = np.array([werr[ri][t] for ri in range(len(r))])
        plt.errorbar(rmed[1:]+orig_offset[1:]*t, wt[1:], yerr=wterr[1:],
                     fmt='o-', label='$z_{}$={}-{}'.format(classes[t],
                                                           zgroups[t][0],
                                                           zgroups[t][1]))
    plt.grid()
    plt.legend()
    plt.xscale('log')
    plt.xlim(0.0035, 1.0) 
    plt.ylim(-0.13, 0.23)
    plt.savefig('/Users/Christina/git/thesis/{}_wcalc.png'.format(name))
    plt.close()


if __name__=="__main__":
    main()
