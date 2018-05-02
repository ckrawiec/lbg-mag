import findpairs
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
from dataset import DataSet
from assigntypes import assignTypes
from parse import parseconfigs
from astropy.table import Table, vstack
from astropy.io import fits


filters = 'GRIZ'


def getSlice(self, zmin, zmax):
    #inclusive
    zmask = (self.data[self.zcol] >= zmin) & (self.data[self.zcol] <= zmax)
    return self.data[zmask]

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
    zp_files = glob.glob(params['zp_files'])
    sys.stderr.write('Found {} z-prob files.\n'.format(len(zp_files)))
    sys.stderr.write('Reading files...')
    zp_tabs = [Table.read(zp_file) for zp_file in zp_files]
    zp_tab = vstack(zp_tabs)
    sys.stderr.write('Done.\n Found {} objects.\n'.format(len(zp_tab)))

    #create target and lens data objects
    DataSet.assignTypes = assignTypes
    DataSet.getSlice = getSlice
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

    #use balrogs as randoms for correlation function
    balrog = {}
    for balrog_file in glob.glob(params['sim_file_format'].format('*'))[:2]:
        itab = balrog_file.find('tab')
        tabnum = balrog_file[itab+3:itab+5]
        balrog_zp_files = glob.glob(params['balrog_zp_files'].format(tabnum))
        if len(balrog_zp_files)>0:
            sys.stderr.write('Reading Balrog z-prob files for Table {}...'.format(tabnum))
            balrog_zp_tabs = [Table.read(balrog_zp_file) for balrog_zp_file in balrog_zp_files]
            balrog_zp_tab = vstack(balrog_zp_tabs)
            sys.stderr.write('Done.\n Found {} objects.\n'.format(len(balrog_zp_table)))

            this_balrog = DataSet(balrog_file,
                                  zprobtab=balrog_zp_tab,
                                  idcol=params['balrog_id_column'],
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
    DataSet.getSlice = getSlice
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
            dNdMu_params = params
            dNdMu_params['redshifts'] = [zmin, zmax]
            #k, detections = dNdMu.main(dNdMu_params)       

            #for each table, find change in detected number for this objtype
            old, new = 0, 0
            for tabnum in balrog.keys():
                #ids of objtype
                type_ids = balrog[tabnum][objtype][params['balrog_id_column']]
                #are they in detections?
           #     old += len(set(type_ids).intersection(detections[tabnum]['original matches']))
           #     new += len(set(type_ids).intersection(detections[tabnum]['magnified matches']))

            #k_mat[-1].append(float(new-old) / (params['mu'] - 1.))
            
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


