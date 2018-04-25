import findpairs
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
from parse import parseconfigs
from astropy.io import fits
from astropy.table import Table, vstack, join
from myutils import joincorrecttype
from matplotlib.pyplot import cm

filters = 'GRIZ'

class DataSet:
    def __init__(self, datafile,  output=None,
                 zprobtab=None, zcol=None, magcol=None, sizecol=None, fluxcol=None,
                 idcol=None, typecol=None):
        if zprobtab:
            data_tab = Table.read(datafile)
            self.data = join(zprobtab, data_tab)
            self.zranges = getzgroups(self.data.columns)
            self.probabilities = ['P'+str(zrange) for zrange in self.zranges]
        else:
            self.data = fits.open(datafile)[1].data
            #for old output version in z-prob master branch
            #self.data = joincorrecttype(data_file, zprob_file, id_col, id_col, float)
        self.output = output
        self.idcol = idcol
        if zcol:
            self.zcol = zcol
        if magcol:
            self.magcol = magcol
        if fluxcol:
            self.fluxcol = fluxcol
        if sizecol:
            self.sizecol = sizecol
        if typecol:
            self.typecol = typecol

    def assignTypes(self, zindices, belows, aboves, check=False):
        color = cm.rainbow(np.linspace(0,1,len(zindices)+1))
        
        if check:
            plt.plot(self.data[self.magcol.format('I')]-self.data[self.magcol.format('Z')],
                     self.data[self.magcol.format('R')]-self.data[self.magcol.format('I')],
                     'o', markeredgecolor='none', c=color[-1], markersize=1.5, label='all')

        self.types = {}
        #cycle through redshift groups    
        for i,c in zip(range(len(zindices)), color[:-1]):
    
            where = np.where(self.data)
            data = self.data[where]

            #cycle through this redshift group
            for j in range(len(zindices[i])):
                probs = [data[self.probabilities[zindex]] for zindex in zindices[i][j]]
    
                where = np.where((np.sum(probs, axis=0) < belows[i][j]) & (np.sum(probs, axis=0) > aboves[i][j]))
                data = data[where]

            #save data in type dictionary
            self.types[i] = data

            if check:
                #check colors for each type
                plt.plot(data[self.magcol.format('I')]-data[self.magcol.format('Z')],
                         data[self.magcol.format('R')]-data[self.magcol.format('I')],
                         '^', markeredgecolor='none', c=c, markersize=3.,
                         label='type '+str(i))

            print "Type {}: {} objects".format(i, len(data))
                
        if check:
            #save figure
            plt.legend(loc='best')
            plt.xlim(-1,4)
            plt.ylim(-1,4)
            plt.xlabel('{} - {}'.format(self.magcol.format('I'), self.magcol.format('Z')))
            plt.ylabel('{} - {}'.format(self.magcol.format('R'), self.magcol.format('I')))
            plt.savefig(self.output+'_typecolors.png')
            plt.close()

    def getSlice(zmin, zmax):
        #inclusive
        zmask = (self.data[self.zcol] >= zmin) & (self.data[self.zcol] <= zmax)
        return self.data[zmask]

    
def getzgroups(columns):
    zgroups = []
    for column in columns:
        if column[:2]=="P[":
            zgroups.append(eval(column[1:]))
    ind = np.argsort([np.min(zgrp) for zgrp in zgroups])
    out = [zgroups[i] for i in ind]

    return out

def main(args):
    #parse config file
    params = parseconfigs(args[1])

    #check existence of output files from this code
    outputs = glob.glob('{}_type*.fits'.format(params['output']))
    check_output = 0
    for output in outputs:
        if os.path.exists(output):
            sys.stderr.write('File {} already exists.\n'.format(output))
            check_output += 1
    if check_output:
        exit()
    
    #get all z-prob files
    zp_files = glob.glob(params['zp_file'])
    sys.stderr.write('Found {} z-prob files.\n'.format(len(zp_files)))
    zp_tab = Table.read(zp_files[0])
    sys.stderr.write('Reading files...')
    for zp_file in zp_files[1:]:
        zp_tab = vstack([zp_tab, Table.read(zp_file)])
    sys.stderr.write('Done.\n Found {} objects.\n'.format(len(zp_tab)))

    #create target and lens data objects
    objects = DataSet(params['data_file'],
                      zprobtab = zp_tab,
                      idcol=params['zp_id_col'],
                      output=params['output'],
                      magcol='MAG_AUTO_{}')
    lenses = fits.open(params['lens_file'])[1].data
    
    #assign types to objects from zprob_file based on their probabilities
    objects.assignTypes(params['redshift_index'],
                        params['below'], params['above'])
    
    #files from test region
    test = DataSet(params['test_data_file'],
                   zprobtab=Table.read(params['test_zp_file']),
                   idcol=params['test_data_id_col'],
                   output=params['output']+'_test',
                   zcol=params['test_z_col'],
                   magcol='MAG_AUTO_{}_d10')

    test.assignTypes(params['redshift_index'],
                     params['below'], params['above'])

    #use balrogs as randoms for correction function
    balrog = {}
    for balrog_file in glob.glob(params['sim_file_format'].format('*'))[:2]:
        itab = balrog_file.find('tab')
        tabnum = balrog_file[itab+3:itab+5]
        balrog_zp_files = glob.glob(params['balrog_zp_files'].format(tabnum))
        if len(balrog_zp_files)>0:
            balrog_zp_table = Table.read(balrog_zp_files[0])
            sys.stderr.write('Reading Balrog z-prob files for Table {}...'.format(tabnum))
            for balrog_zp_file in balrog_zp_files[1:]:
                balrog_zp_table = vstack([balrog_zp_table, Table.read(balrog_zp_file)])
            sys.stderr.write('Done.\n Found {} objects.\n'.format(len(balrog_zp_table)))

            this_balrog = DataSet(balrog_file,
                                  zprobtab=balrog_zp_table,
                                  idcol=params['balrog_id_col'],
                                  output=params['output']+'_balrog'+str(tabnum),
                                  magcol='MAG_AUTO_{}')
            this_balrog.assignTypes(params['redshift_index'],
                                    params['below'], params['above'])
            balrog[tabnum] = {}
            for objtype in objects.types:
                balrog[tabnum][objtype] = this_balrog.types[objtype]
                #######
        
    #measured quantities, vectors and matrices
    n_vec = []
    P_mat = []
    k_mat = []
    
    import dNdMu
    for objtype in objects.types:
        sys.stderr.write('Working on type {}...'.format(objtype))
        
        #write type objects to table
        this_table = '{}_type{}.fits'.format(params['output'], objtype)
        objects.types[objtype].write(this_table)

        random_table = '{}_type{}_randoms.fits'.format(params['output'], objtype)
        these_randoms = vstack([balrog[ti][objtype] for ti in balrog.keys()])
        these_randoms.write(random_table)
        
        #for each type - create parameters for treecorr 
        these_params = {}
        these_params['source_file'] = this_table
        these_params['lens_file'] = params['lens_file'] 
        these_params['random_source_file'] = random_table
        these_params['random_lens_file'] = params['random_lens_file']
        these_params['source_ra'] = 'RA'
        these_params['source_dec'] = 'DEC'
        these_params['output'] = params['output']+'_type'+str(objtype)
        
        #calculate correlation with lenses - call treecorr
        nn = findpairs.treecorr(these_params)
        n_vec.append(nn['npairs'])

        sys.stderr.write('Done.\n')
                
        P_mat.append([])
        k_mat.append([])

        #plot correlation functions
        plt.errorbar(nn['R_nom'], nn['xi'], yerr=nn['sigma_xi'],
                     fmt='o-', label='Type '+str(objtype))
        
        
        for true_objtype in objects.types:
            #"true" redshifts of true_objtype objects in test field
            z_photo = test.data[test.zcol]
            zmin = np.min(params['true_ranges'][true_objtype])
            zmax = np.max(params['true_ranges'][true_objtype])
            z_true = np.where((z_photo >= zmin) & (z_photo <= zmax))
            #ids of true_objtype objects
            ids = test.data[test.idcol][z_true]

            #"true" redshifts of objtype objects in test field
            this_z_photo = test.types[objtype][test.zcol]
            this_z_true = np.where((this_z_photo >= zmin) & (this_z_photo <= zmax))
            these_ids = test.types[objtype][test.idcol][this_z_true]

            #assigned to objtype while actually true_objtype
            both = float(len(set(ids).intersection(these_ids)))
            
            #probabilities of true/false assignment
            P_mat[-1].append(both/float(len(ids)))
                        
            #Run dNdMu with mu_G
            dNdMu_params = {'filters' : 'GRIZ',
                            'mu' : 1.5,
                            'flux_cut' : 40.,
                            'deep_file' : '/Users/Christina/DES/data/y1a1_gold_dfull_cosmos.fits',
                            'balrog_files' : '/Users/Christina/DES/data/balrog_sva1_tab*_TRUTH_zp_corr_fluxes.fits',
                            'sim_file_format' : params['sim_file_format'],
                            'deep_flux_column' : 'FLUX_AUTO_{}',
                            'deep_size_column' : 'FLUX_RADIUS_I',
                            'balrog_flux_column' : 'FLUX_NOISELESS_{}',
                            'balrog_size_column' : 'HALFLIGHTRADIUS_0',
                            'redshifts' : [zmin, zmax]}
            print dNdMu_params

            k = dNdMu.main(dNdMu_params)
            k_mat[-1].append(k)
            
    plt.legend()
    plt.ylabel('xi')
    plt.xlabel('R_nom')
    plt.ylim(-0.01, 1.0)
    plt.grid()
    plt.savefig(params['output']+'_typecorrs.png')
    plt.close()

    print "n = ", n_vec
    print "P = ", P_mat
    print "k = ", k_mat


if __name__=="__main__":
    main(sys.argv)


