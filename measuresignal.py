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
                                params['below'], params['above'],
                                trueranges=params['true_ranges'])
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

    del this_balrog, type_data, balrog

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

    #targets
    if not outputs_exist['target types']:
        target_objects = DataSet(params['source_file'],
                                 zpfiles = params['zp_files'],
                                 idcol=params['zp_id_column'],
                                 output=params['output'],
                                 magcol='MAG_AUTO_{}')
        target_objects.assignTypes(params['redshift_index'],
                                   params['below'], params['above'], 
                                   trueranges=params['true_ranges'])
        for objtype in target_objects.types:
            target_table = '{}_targets_type{}.fits'.format(params['output'], objtype)
            target_data = target_objects.data[target_objects.types[objtype]]
            target_data.write(target_table, overwrite=True)

        del target_objects, target_data

    #get balrog objects and their assigned types and write to file
    if not outputs_exist['balrog types']:
        assignbalrog(params, tabnums)

    #measured quantities, vectors and matrices
    n_vec, n0_vec = [], []
    output_table = {}
    for objtype in range(len(params['true_ranges'])):
        sys.stderr.write('Working on type {}...\n'.format(objtype))
        
        #read type objects to/from table
        target_table = '{}_targets_type{}.fits'.format(params['output'], objtype)
        random_table = '{}_balrogsim_type{}.fits'.format(params['output'], objtype)

        #correlations with lenses for each type
        corr_output = params['output']+'_type{}_correlations.fits'.format(objtype)
        if  not os.path.exists(corr_output):
            #create parameters for treecorr 
            these_params = {}
            these_params['source_file'] = target_table
            these_params['random_file'] = random_table
            these_params['lens_file'] = params['lens_file'] 
            these_params['num_threads'] = params['num_threads']
            these_params['random_type'] = 'source'
            these_params['source_weight_column'] = 'P'+str(params['true_ranges'][objtype])
            these_params['output'] = params['output']+'_type'+str(objtype)

            nn = findpairs.main(these_params)
            
        else:
            nn = Table.read(corr_output)

        print nn
        n_vec.append(nn['DD'])
        sys.stderr.write('Done.\n')

        #open target table
        targets = Table.read(target_table)
        n0 = np.sum(targets[Ps[objtype]])
        n0_vec.append(n0)
        output_table['n0_{}'.format(objtype)] = n0

    k_mat, P_mat = testcompare(params, outputs_exist, tabnums, detections)
           
    import dNdMu
    #Run dNdMu with mu_G
    dNdMu_params = params
    dNdMu_params['redshifts'] = [[zmin, zmax] for zmin, zmax in params['true_ranges']]
    k, detections = dNdMu.main(dNdMu_params)
    output_table['k_output'.format(true_objtype)] = k
    print "n = ", n_vec
    print "n0 = ", n0_vec
    print "P = ", P_mat
    print "k = ", k_mat

    for key in output_table.keys():
        print "{}: {}".format(key, output_table[key])

    #save outputs to file
    import pickle
    f = open(params['output']+'_output.pkl','wb')
    pickle.dump(output_table, f)

def testcompare(params, outputs_exist, tabnums, detections):
    #test DataSet
    test_objects = DataSet(params['test_data_file'],
                           zpfiles=params['test_zp_file'],
                           idcol=params['test_data_id_column'],
                           output=params['output']+'_test',
                           zcol=params['test_z_column'],
                           magcol=params['test_mag_column'])
    test_objects.assignTypes(params['redshift_index'],
                             params['below'], params['above'],
                             trueranges=params['true_ranges'])
    
    if not outputs_exist['test types']:
        for objtype in test_objects.types:
            test_table = '{}_test_type{}.fits'.format(params['output'], objtype)
            test_data = test_objects.data[test_objects.types[objtype]]
            test_data.write(test_table, overwrite=True)
            
    #test info needed for later
    test_z = test_objects.data[test_objects.zcol]
    test_ids = test_objects.data[test_objects.idcol]
    k_mat, P_mat = []
    for objtype in range(len(params['true_ranges'])):
        #needed for each objtype
        type_test_z = test_z[test_objects.types[objtype]]
        type_test_ids = test_ids[test_objects.types[objtype]]

        #balrog ids for each table for this objtype
        type_balrog_ids = {}
        for tabnum in tabnums:
            balrog_file = params['output']+'_balrogsimtab{}_type{}.fits'.format(tabnum, objtype)
            balrog_tab = Table.read(balrog_file)
            type_balrog_ids[tabnum] = balrog_tab[params['balrog_id_column']]

        k_mat.append([])
        P_mat.append([])
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

            #for each table, find change in detected number for this objtype
            old, new = 0, 0
            for tabnum in detections.keys():
                #are objtype objects in detections?
                old += len(set(type_balrog_ids).intersection(set(detections[tabnum]['original matches'])))
                new += len(set(type_balrog_ids).intersection(set(detections[tabnum]['magnified matches'])))

            print "True type {} original/magnified matches: {}, {}".format(true_objtype, old, new)
            k_HG = (float(new)-float(old)) / (params['mu'] - 1.)
            
            output_table['k_{}{}'.format(objtype, true_objtype)] = k_HG
            k_mat[-1].append(k_HG)

    return k_mat, P_mat

if __name__=="__main__":
    main(sys.argv[1])


