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

def assignbalrog(params):
    """
    Create dictionary of balrog tables and assigned types.
    Keys are table numbers.
    Sub-keys are object types.
    """
    balrog = {}
    #loop over balrog sim files
    for balrog_file in glob.glob(params['balrog_sim_file_format'].format('*')):
        #save table number
        itab = balrog_file.find('tab')
        tabnum = balrog_file[itab+3:itab+5]

        #gather zprob fies for this table
        balrog_zp_files = glob.glob(params['balrog_zp_files'].format(tabnum))
        this_balrog = DataSet(balrog_file,
                              zpfiles=balrog_zp_files,
                              idcol=params['balrog_id_column'],
                              output=params['output']+'_balrog'+str(tabnum),
                              magcol=params['balrog_sim_mag_column'])
        this_balrog.assignTypes(params['redshift_index'],
                                params['below'], params['above'])
        balrog[tabnum] = {}
        for objtype in this_balrog.types:
            balrog[tabnum][objtype] = this_balrog.data[this_balrog.types[objtype]]
            #write to file
            tabtype_file = '{}_balrogsimtab{}_type{}.fits'.format(params['output'],
                                                                  tabnum, objtype)
            sys.stderr.write('Writing sim types to {}\n'.format(tabtype_file))
            balrog[tabnum][objtype].write(tabtype_file, overwrite=True)

    for objtype in this_balrog.types:
        type_file = '{}_balrogsim_type{}.fits'.format(params['output'], objtype)
        type_data = vstack([balrog[ti][objtype] for ti in balrog.keys()])
        sys.stderr.write('Writing sim types to {}.\n'.format(type_file))
        type_data.write(type_file, overwrite=True)

    del balrog

def calcmu(n, n0, P, k):
    diff = n - np.dot(P, n0)
    return np.dot(diff, np.linalg.inv(k))

def checkoutputs(params):
    outputs_exist = {}
    objtypes = range(len(params['true_ranges']))

    #check type outputs
    target_files = [os.path.exists('{}_targets_type{}.fits'.format(params['output'], objtype)) for objtype in objtypes]
    outputs_exist['target types'] = np.all(target_files)

    test_files = [os.path.exists('{}_test_type{}.fits'.format(params['output'], objtype)) for objtype in objtypes]
    outputs_exist['test types'] = np.all(test_files)
    
    tabnums = []
    for sim_file glob.glob(params['balrog_sim_file_format'].format('*')):
        #save table number
        itab = sim_file.find('tab')
        tabnum.append(sim_file[itab+3:itab+5])
    base_output = params['output']+'_balrogsimtab{}_type{}.fits'
    tabtypes_exists = np.all([[os.path.exists(base_output.format(tabnum, objtype) for tabnum in tabnums] for objtype in objtypes]])
    types_exist = np.all(['{}_balrogsim_type{}.fits'.format(params['output'], objtype)] for objtype in objtypes)
    outputs_exist['balrog types'] = tabtypes_exists & types_exist

    #correlation outputs
    corr_output = params['output']+'_type{}_correlations.fits'
    corrs_exist = np.all([os.path.exists(corr_output.format(objtype)) for objtype in objtypes])
    outputs_exist['correlations'] = corrs_exist
    
    return outputs_exist


def main(config):
    #parse config file
    params = parseconfigs(config)

    #check existence of output files from this code
    outputs_exist = checkoutputs(params)

    #get test and target objects and their assigned types
    if not outputs_exist['test types']:
        test_objects = ?
    if not outputs_exist['target types']:
        target_objects = ?

    data = target_objects.data[target_objects.types[objtype]]
    data.write(data_table, overwrite=True)
    test_objects, target_objects = assigntypes.main(config)

    #get balrog objects and their assigned types and write to file
    DataSet.assignTypes = assignTypes
    if not outputs_exist['balrog types']:
        balrog = assignbalrog(params)

    #measured quantities, vectors and matrices
    n_vec = []
    n0_vec = []
    P_mat = []
    k_mat = []
    
    import dNdMu
    output_table = {}
    for objtype in range(len(params['true_ranges'])):
        sys.stderr.write('Working on type {}...\n'.format(objtype))
        
        #write/read type objects to/from table
        target_table = '{}_data_type{}.fits'.format(params['output'], objtype)
        random_table = '{}_balrogsim_type{}'.format(params['output'], objtype)
        targets = Table.read(data_table)

        Ps = ['P'+str(zrange) for zrange in getzgroups(targets.columns)]
        n0 = np.sum(targets[Ps[objtype]])
        n0_vec.append(n0)
        output_table['n0_{}'.format(objtype)] = n0

        corr_output = params['output']+'_type{}_correlations.fits'.format(objtype)
        if  not os.path.exists(corr_output):
            #for each type - create parameters for treecorr 
            these_params = {}
            these_params['source_file'] = target_table
            these_params['lens_file'] = params['lens_file'] 
            these_params['random_file'] = random_table
            these_params['random_type'] = 'source'
            these_params['weight_column'] = Ps[objtype]
            these_params['output'] = params['output']+'_type'+str(objtype)

            nn = findpairs.main(these_params)
        else:
            nn = Table.read(corr_output)

        output_table['DD{}'.format(objtype)] = nn['DD']
        output_table['DR{}'.format(objtype)] = nn['DR']
        output_table['DD_w{}'.format(objtype)] = nn['DD_w']
        output_table['DR_w{}'.format(objtype)] = nn['DR_w']
        n_vec.append(nn['DD'])
        n_weight_vec.append(nn['DD_w'])
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
            this_P = both/float(len(ids))
            output_table['P{}{}'.format(objtype, true_objtype)] = this_P
            P_mat[-1].append(this_P)

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

            this_k = float(new-old) / (params['mu'] - 1.)
            output_table['k{}{}'.format(objtype, true_objtype)] = this_k
            k_mat[-1].append(this_k)
            
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

    output_tab = Table(output_table)
    output_tab.write(params['output']+'_output.fits')

    #mu = calcmu(n_vec, n0_vec, P_mat, k_mat)
    #print "mu = ", mu

if __name__=="__main__":
    main(sys.argv[1])


