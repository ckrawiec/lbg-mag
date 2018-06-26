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

def assignbalrog(params, tabnums):
    """
    Create dictionary of balrog tables and assigned types.
    Keys are table numbers.
    Sub-keys are object types.
    """
    balrog = {}
    #loop over balrog sim files
    for tabnum in tabnums:
        #save table number
        balrog_file = params['balrog_sim_file_format'].format(tabnum)
        itab = balrog_file.find('tab')
        tabnum = balrog_file[itab+3:itab+5]

        #gather zprob fies for this table
        balrog_zp_files = params['balrog_zp_files'].format(tabnum)
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

def checkoutputs(params, tabnums):
    outputs_exist = {}
    objtypes = range(len(params['true_ranges']))

    #type outputs
    target_files = [os.path.exists('{}_targets_type{}.fits'.format(params['output'], objtype)) for objtype in objtypes]
    outputs_exist['target types'] = np.all(target_files)

    test_files = [os.path.exists('{}_test_type{}.fits'.format(params['output'], objtype)) for objtype in objtypes]
    outputs_exist['test types'] = np.all(test_files)
    
    base_output = params['output']+'_balrogsimtab{}_type{}.fits'
    tabtypes_exist = np.all([[os.path.exists(base_output.format(tabnum, objtype)) for tabnum in tabnums] for objtype in objtypes])
    types_exist = np.all([os.path.exists('{}_balrogsim_type{}.fits'.format(params['output'], objtype)) for objtype in objtypes])
    outputs_exist['balrog types'] = tabtypes_exist & types_exist

    #correlation outputs
    corr_output = params['output']+'_type{}_correlations.fits'
    corrs_exist = np.all([os.path.exists(corr_output.format(objtype)) for objtype in objtypes])
    outputs_exist['correlations'] = corrs_exist
    
    return outputs_exist


def main(config):
    #parse config file
    params = parseconfigs(config)

    #check existence of output files from this code
    tabnums = []
    for sim_file in glob.glob(params['balrog_sim_file_format'].format('*')):
        #save table number
        itab = sim_file.find('tab')
        tabnums.append(sim_file[itab+3:itab+5])
    outputs_exist = checkoutputs(params, tabnums)

    #get test and target objects and their assigned types
    DataSet.assignTypes = assignTypes

    #need test DataSet later
    test_objects = DataSet(params['test_data_file'],
                           zpfiles=params['test_zp_file'],
                           idcol=params['test_data_id_column'],
                           output=params['output']+'_test',
                           zcol=params['test_z_column'],
                           magcol=params['test_mag_column'])
    test_objects.assignTypes(params['redshift_index'],
                             params['below'], params['above'])
    if not outputs_exist['test types']:
        for objtype in test_objects.types:
            test_table = '{}_test_type{}.fits'.format(params['output'], objtype)
            test_data = test_objects.data[test_objects.types[objtype]]
            test_data.write(test_table, overwrite=True)

    #targets
    if not outputs_exist['target types']:
        target_objects = DataSet(params['source_file'],
                                 zpfiles = params['zp_files'],
                                 idcol=params['zp_id_column'],
                                 output=params['output'],
                                 magcol='MAG_AUTO_{}')
        target_objects.assignTypes(params['redshift_index'],
                                   params['below'], params['above'])
        for objtype in target_objects.types:
            target_table = '{}_targets_type{}.fits'.format(params['output'], objtype)
            target_data = target_objects.data[target_objects.types[objtype]]
            target_data.write(target_table, overwrite=True)

    #get balrog objects and their assigned types and write to file
    if not outputs_exist['balrog types']:
        assignbalrog(params, tabnums)

    #measured quantities, vectors and matrices
    n_vec = []
    n0_vec = []
    P_mat = []
    k_mat = []
    
    import dNdMu

    #test info needed for later
    test_z = test_objects.data[test_objects.zcol]
    test_ids = test_objects.data[test_objects.idcol]
    output_table = {}
    for objtype in range(len(params['true_ranges'])):
        sys.stderr.write('Working on type {}...\n'.format(objtype))
        
        #write/read type objects to/from table
        target_table = '{}_targets_type{}.fits'.format(params['output'], objtype)
        random_table = '{}_balrogsim_type{}.fits'.format(params['output'], objtype)
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
            these_params['source_weight_column'] = Ps[objtype]
            these_params['output'] = params['output']+'_type'+str(objtype)

            nn = findpairs.main(these_params)
        else:
            nn = Table.read(corr_output)

        print nn
        output_table['DD{}'.format(objtype)] = nn['DD']
        output_table['DR{}'.format(objtype)] = nn['DR']
        if 'DD_w' in nn.columns:
            output_table['DD_w{}'.format(objtype)] = nn['DD_w']
        if 'DR_w' in nn.columns:
            output_table['DR_w{}'.format(objtype)] = nn['DR_w']
        n_vec.append(nn['DD'])
        sys.stderr.write('Done.\n')
                
        P_mat.append([])
        k_mat.append([])

        #needed for each objtype
        type_test_z = test_z[test_objects.types[objtype]]
        type_test_ids = test_ids[test_objects.types[objtype]]

        #balrog ids for each table for this objtype
        type_balrog_ids = {}
        for tabnum in tabnums:
            balrog_file = params['output']+'_balrogsimtab{}_type{}.fits'.format(tabnum, objtype)
            balrog_tab = Table.read(balrog_file)
            type_balrog_ids[tabnum] = balrog_tab[params['balrog_id_column']]
        
        for true_objtype in range(len(params['true_ranges'])):
            zmin = np.min(params['true_ranges'][true_objtype])
            zmax = np.max(params['true_ranges'][true_objtype])

            #truth redshifts from test data
            true_z = np.where((test_z >= zmin) & (test_z <= zmax))

            #ids of true_objtype objects
            true_ids = test_ids[true_z]

            #objtype objects whose truth redshifts match assignment
            true_type_mask = np.where((type_test_z >= zmin) & (type_test_z <= zmax))
            true_type_ids = type_test_ids[true_type_mask]

            #assigned to objtype while actually true_objtype
            both = float(len(set(true_ids).intersection(true_type_ids)))
            
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
            output_table['k_{}{}_output'.format(objtype, true_objtype)] = k
            output_table['k_{}{}'.format(objtype, true_objtype)] = k_HG
            k_mat[-1].append(k_HG)

    print "n = ", n_vec
    print "n0 = ", n0_vec
    print "P = ", P_mat
    print "k = ", k_mat

    ######## need sum of n0 in all redshift ranges for ALL sources at each annulus
    np.save(params['output']+'_output.fits', output_table)

    #mu = calcmu(n_vec, n0_vec, P_mat, k_mat)
    #print "mu = ", mu

if __name__=="__main__":
    main(sys.argv[1])


