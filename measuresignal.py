import matplotlib
matplotlib.use('agg')
import findpairs
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
from dataset import DataSet, getzgroups
from assigntypes import assignTypes
from parse import parseconfigs
from astropy.table import Table, vstack
from astropy.io import fits

filters = 'GRIZ'

def calcmu(n, n0, P, k):
    diff = n - np.dot(P, n0)
    return np.dot(diff, np.linalg.inv(k))

def checkoutputs(params, outputs):
    check_output = 0
    for output in outputs:
        if 'corr' in output:
            if os.path.exists(output) and not params['overwrite_corr_files']:
                sys.stderr.write('File {} already exists and will be used.\n'.format(output))
            elif os.path.exists(output):
                sys.stderr.write('File {} already exists and will be overwritten.\n'.format(output))

        elif 'type' in output:
            if os.path.exists(output) and not params['overwrite_type_files']:
                sys.stderr.write('File {} already exists and will be used.\n'.format(output))
            elif os.path.exists(output):
                sys.stderr.write('File {} already exists and will be overwritten.\n'.format(output))
        
    if check_output:
        sys.exit('Delete files or set overwrite=True. Exiting.\n')

    return params
    

def gatherdata(params):
    #get all z-prob files
    zp_files = glob.glob(params['zp_files'])
    sys.stderr.write('Found {} z-prob files.\n'.format(len(zp_files)))
    sys.stderr.write('Reading files...')
    zp_tabs = [Table.read(zp_file) for zp_file in zp_files]
    zp_tab = vstack(zp_tabs)
    sys.stderr.write('Done.\n Found {} objects.\n'.format(len(zp_tab)))
    
    #create target and lens data objects
    objects = DataSet(params['data_file'],
                      zprobtab = zp_tab,
                      idcol=params['zp_id_col'],
                      output=params['output'],
                      magcol='MAG_AUTO_{}')
    
    #assign types to objects from zprob_file based on their probabilities
    objects.assignTypes(params['redshift_index'],
                        params['below'], params['above'])

    return objects

def main(args):
    #parse config file
    params = parseconfigs(args[1])

    #check existence of output files from this code
    outputs = glob.glob('{}_type*.fits'.format(params['output']))
    params = checkoutputs(params, outputs)

    #files from test region
    DataSet.assignTypes = assignTypes
    test = DataSet(params['test_data_file'],
                   zprobtab=Table.read(params['test_zp_file']),
                   idcol=params['test_data_id_col'],
                   output=params['output']+'_test',
                   zcol=params['test_z_col'],
                   magcol='MAG_AUTO_{}_d10')

    test.assignTypes(params['redshift_index'],
                     params['below'], params['above'])
    
    #create data or read from existing files    
    if params['overwrite_type_files']:
        objects = gatherdata(params)

    #use balrogs as randoms for correlation function
    balrog = {}
    for balrog_file in glob.glob(params['sim_file_format'].format('*')):
        itab = balrog_file.find('tab')
        tabnum = balrog_file[itab+3:itab+5]
        balrog_zp_files = glob.glob(params['balrog_zp_files'].format(tabnum))
        if len(balrog_zp_files)>0:
            sys.stderr.write('Reading Balrog z-prob files for Table {}...'.format(tabnum))
            balrog_zp_tabs = [Table.read(balrog_zp_file) for balrog_zp_file in balrog_zp_files]
            balrog_zp_tab = vstack(balrog_zp_tabs)
            sys.stderr.write('Done.\n Found {} objects.\n'.format(len(balrog_zp_tab)))

            this_balrog = DataSet(balrog_file,
                                  zprobtab=balrog_zp_tab,
                                  idcol=params['balrog_id_column'],
                                  output=params['output']+'_balrog'+str(tabnum),
                                  magcol='MAG_AUTO_{}')
            this_balrog.assignTypes(params['redshift_index'],
                                params['below'], params['above'])
            balrog[tabnum] = {}
            for objtype in test.types:
                balrog[tabnum][objtype] = this_balrog.data[this_balrog.types[objtype]]

    #measured quantities, vectors and matrices
    n_vec = []
    n0_vec = []
    P_mat = []
    k_mat = []
    
    import dNdMu
    for objtype in range(len(params['true_ranges'])):
        sys.stderr.write('Working on type {}...'.format(objtype))
        
        #write/read type objects to/from table
        this_table = '{}_type{}.fits'.format(params['output'], objtype)
        random_table = '{}_type{}_randoms.fits'.format(params['output'], objtype)
        if params['overwrite_type_files']:
            this_data = objects.data[objects.types[objtype]]
            this_data.write(this_table, overwrite=True)
            n0_vec.append(np.sum(this_data[objects.probabilities[objtype]]))
            these_randoms = vstack([balrog[ti][objtype] for ti in balrog.keys()])
            these_randoms.write(random_table, overwrite=True)
        else:
            this_data = Table.read(this_table)
            these_Ps = ['P'+str(zrange) for zrange in getzgroups(this_data.columns)]
            n0_vec.append(np.sum(this_data[these_Ps[objtype]]))
            
        this_output = params['output']+'_type'+str(objtype)+'_treecorrNN.fits'
        if (params['overwrite_corr_files']) or not os.path.exists(this_output):
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
        else:
            nn = Table.read(this_output)
            
        n_vec.append(nn['npairs'])
        sys.stderr.write('Done.\n')
                
        P_mat.append([])
        k_mat.append([])

        #plot correlation functions
        plt.errorbar(nn['R_nom'], nn['xi'], yerr=nn['sigma_xi'],
                     fmt='o-', label='Type '+str(objtype))
        
        
        for true_objtype in range(len(params['true_ranges'])):
            #"true" redshifts of true_objtype objects in test field
            z_photo = test.data[test.zcol]
            zmin = np.min(params['true_ranges'][true_objtype])
            zmax = np.max(params['true_ranges'][true_objtype])
            z_true = np.where((z_photo >= zmin) & (z_photo <= zmax))

            #ids of true_objtype objects
            ids = test.data[test.idcol][z_true]

            #"true" redshifts of objtype objects in test field
            #if params['overwrite_type_files']:
            this_z_photo = test.data[test.types[objtype]][test.zcol]
            this_z_true = np.where((this_z_photo >= zmin) & (this_z_photo <= zmax))
            these_ids = test.data[test.types[objtype]][test.idcol][this_z_true]

            #assigned to objtype while actually true_objtype
            both = float(len(set(ids).intersection(these_ids)))
            
            #probabilities of true/false assignment
            P_mat[-1].append(both/float(len(ids)))

            #Run dNdMu with mu_G
            dNdMu_params = params
            dNdMu_params['redshifts'] = [zmin, zmax]
            k, detections = dNdMu.main(dNdMu_params)       

            #for each table, find change in detected number for this objtype
            old, new = 0, 0
            for tabnum in detections.keys():
                #ids of objtype
                type_ids = balrog[tabnum][objtype][params['balrog_id_column']]
                #are they in detections?
                old += len(set(type_ids).intersection(detections[tabnum]['original matches']))
                new += len(set(type_ids).intersection(detections[tabnum]['magnified matches']))

            k_mat[-1].append(float(new-old) / (params['mu'] - 1.))
            
    plt.legend()
    plt.ylabel('xi')
    plt.xlabel('R_nom')
    #plt.ylim(-0.01, 1.0)
    plt.grid()
    plt.savefig(params['output']+'_typecorrs.png')
    plt.close()

    print "n = ", n_vec
    print "n0 = ", n0_vec
    print "P = ", P_mat
    print "k = ", k_mat

    mu = calcmu(n_vec, n0_vec, P_mat, k_mat)
    print "mu = ", mu

if __name__=="__main__":
    main(sys.argv)


