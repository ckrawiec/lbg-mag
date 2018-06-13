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
    DataSet.assignTypes = assignTypes
    if not outputs_exist['test types']:
        test_objects = DataSet(params['test_data_file'],
                               zpfiles=params['test_zp_file'],
                               idcol=params['test_data_id_col'],
                               output=params['output']+'_test',
                               zcol=params['test_z_col'],
                               magcol=params['test_mag_col'])
        test_objects.assignTypes()
        for objtype in test_objects.types:
            test_table = '{}_test_type{}.fits'.format(params['output'], objtype)
            test_data = test_objects.data[test_objects.types[objtype]]
            test_data.write(data_table, overwrite=True)
        
    if not outputs_exist['target types']:
        target_objects = DataSet(params['data_file'],
                                 zpfiles = params['zp_files'],
                                 idcol=params['zp_id_col'],
                                 output=params['output'],
                                 magcol='MAG_AUTO_{}')
        target_objects.assignTypes(params['redshift_index'],
                                   params['below'], params['above'])
        for objtype in objects.types:
            target_table = '{}_targets_type{}.fits'.format(params['output'], objtype)
            target_data = target_objects.data[target_objects.types[objtype]]
            target_data.write(target_table, overwrite=True)

    #get balrog objects and their assigned types and write to file
    if not outputs_exist['balrog types']:
        assignbalrog(params)

    #measured quantities, vectors and matrices
    n_vec = []
    n0_vec = []
    P_mat = []
    k_mat = []
    
    import dNdMu

    #test info needed for later
    z_photo = test.data[test.zcol]
    test_ids = test.data[test.idcol]
    output_table = {}
    for objtype in range(len(params['true_ranges'])):
        sys.stderr.write('Working on type {}...\n'.format(objtype))
        
        #write/read type objects to/from table
        target_table = '{}_targets_type{}.fits'.format(params['output'], objtype)
        random_table = '{}_balrogsim_type{}'.format(params['output'], objtype)
        targets = Table.read(target_table)

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

        #needed for each objtype
        type_photo_z = test.data[test.types[objtype]][test.zcol]
        type_test_ids = test.data[test.types[objtype]][test.idcol]
        type_balrog_ids = balrog[tabnum][objtype][params['balrog_id_column']]
        
        for true_objtype in range(len(params['true_ranges'])):
            zmin = np.min(params['true_ranges'][true_objtype])
            zmax = np.max(params['true_ranges'][true_objtype])

            #truth redshifts from test data
            true_z = np.where((z_photo >= zmin) & (z_photo <= zmax))

            #ids of true_objtype objects
            true_ids = test_ids[z_true]

            #objtype objects whose truth redshifts match assignment
            type_true_mask = np.where((type_photo_z >= zmin) & (type_photo_z <= zmax))
            type_true_ids = type_ids[type_true_mask]

            #assigned to objtype while actually true_objtype
            both = float(len(set(true_ids).intersection(type_true_ids)))
            
            #probabilities of true/false assignment
            P_HG = both/float(len(true_ids))
            output_table['P_{}{}'.format(objtype, true_objtype)] = P_HG
            P_mat[-1].append(P_HG)

            #Run dNdMu with mu_G
            dNdMu_params = params
            dNdMu_params['redshifts'] = [zmin, zmax]
            k, detections = dNdMu.main(dNdMu_params)       

            #for each table, find change in detected number for this objtype
            old, new = 0, 0
            for tabnum in detections.keys():
                #are objtype objects in detections?
                old += len(set(type_balrog_ids).intersection(detections[tabnum]['original matches']))
                new += len(set(type_balrog_ids).intersection(detections[tabnum]['magnified matches']))

            k_HG = float(new-old) / (params['mu'] - 1.)
            output_table['k_{}{}'.format(objtype, true_objtype)] = k_HG
            k_mat[-1].append(k_HG)

    print "n = ", n_vec
    print "n0 = ", n0_vec
    print "P = ", P_mat
    print "k = ", k_mat

    output_tab = Table(output_table)
    output_tab.write(params['output']+'_output.fits')

    mu = calcmu(n_vec, n0_vec, P_mat, k_mat)
    print "mu = ", mu

if __name__=="__main__":
    main(sys.argv[1])


