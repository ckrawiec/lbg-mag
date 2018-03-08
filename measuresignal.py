import numpy as np
import matplotlib.pyplot as plt
import json
import sys
import ConfigParser
import findpairs
import glob
from astropy.io import fits
from astropy.table import Table, vstack, join
from myutils import joincorrecttype
from matplotlib.pyplot import cm

class DataSet():

    def __init__(self, zprob_tab, data_file, id_col, output):
        data_tab = Table.read(data_file)
        self.data = join(zprob_tab, data_tab)
        #for old output version in z-prob master branch
        #self.data = joincorrecttype(data_file, zprob_file, id_col, id_col, float)
        self.zranges = getzgroups(self.data.columns)
        self.probabilities = ['P'+str(zrange) for zrange in self.zranges]
        self.output = output

    def assignTypes(self, zindices, belows, aboves, check=False, magcol='MAG_AUTO_{}'):
        color = cm.rainbow(np.linspace(0,1,len(zindices)+1))
        
        if check:
            plt.plot(self.data[magcol.format('I')]-self.data[magcol.format('Z')],
                     self.data[magcol.format('R')]-self.data[magcol.format('I')],
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
                plt.plot(data[magcol.format('I')]-data[magcol.format('Z')],
                         data[magcol.format('R')]-data[magcol.format('I')],
                         '^', markeredgecolor='none', c=c, markersize=3.,
                         label='type '+str(i))

                print "Type {}: {} objects".format(i, len(data))
                
        if check:
            #save figure
            plt.legend(loc='best')
            plt.xlim(-1,4)
            plt.ylim(-1,4)
            plt.xlabel('{} - {}'.format(magcol.format('I'), magcol.format('Z')))
            plt.ylabel('{} - {}'.format(magcol.format('R'), magcol.format('I')))
            plt.savefig(self.output+'_typecolors.png')

def getzgroups(columns):
    zgroups = []
    for column in columns:
        if column[:2]=="P[":
            zgroups.append(eval(column[1:]))
    ind = np.argsort([np.min(zgrp) for zgrp in zgroups])
    out = [zgroups[i] for i in ind]

    return out


def parseconfigs(config_file):
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)

    params = {}

    #Input & Output Files
    params['zp_file'] = config.get('I/O','zp_file')
    params['data_file'] = config.get('I/O','data_file')
    params['lens_file'] = config.get('I/O','lens_file')
    params['balrog_zp_file'] = config.get('I/O','balrog_zp_file')
    params['balrog_data_file'] = config.get('I/O','balrog_data_file')
    params['output'] = config.get('I/O','output')

    params['test_zprob_file'] = config.get('I/O','test_zprob_file')
    params['test_data_file'] = config.get('I/O','test_data_file')
    params['test_k_file'] = config.get('I/O','test_k_file')
    
    #Table Column Names
    params['zp_id_col'] = config.get('Columns','zp_id_col')
    params['data_id_col'] = config.get('Columns','data_id_col')
    params['test_data_id_col'] = config.get('Columns','test_data_id_col')

    #Assignment Criteria
    params['redshift_index'] = json.loads(config.get('Assignments','redshift_index'))
    params['below'] = json.loads(config.get('Assignments','below'))
    params['above'] = json.loads(config.get('Assignments','above'))

    return params

def main():
    #parse config file
    params = parseconfigs('measuresignal.config')

    #get all z-prob files
    zp_files = glob.glob(params['zp_file'])
    sys.stderr.write('Found {} z-prob files.\n'.format(len(zp_files)))
    zp_tab = Table.read(zp_files[0])
    sys.stderr.write('Reading files...\n')
    for zp_file in zp_files[1:]:
        zp_tab = vstack([zp_tab, Table.read(zp_file)])
    sys.stderr.write('Done. Found {} objects.\n'.format(len(zp_tab)))

    #create target and lens data objects
    objects = DataSet(zp_tab, params['data_file'], params['zp_id_col'], params['output'])
    lenses = fits.open(params['lens_file'])[1].data
    
    #assign types to objects from zprob_file based on their probabilities
    objects.assignTypes(params['redshift_index'], params['below'], params['above'], check=False)

    #test files from COSMOS region
    cosmos = DataSet(Table.read(params['test_zprob_file']), params['test_data_file'],
                     params['test_data_id_col'], params['output']+'_test')
    testP = Table.read(params['test_zprob_file'])
    #testk = Table.read(params['test_k_file'])

    #measured number densities
    n_vec = []
    P_mat = []
    k_mat = []
    for objtype in objects.types:
        sys.stderr.write('Working on type {}...'.format(objtype))
        
        #for each type - calculate correlation with lenses
        these_params = {}

        #write type objects to table
        this_table = '{}_type{}.fits'.format(params['output'], objtype)
        objects.types[objtype].write(this_table)
        these_params['source_file'] = this_table
        these_params['lens_file'] = params['lens_file'] 
        these_params['random_source_file'] = params['lens_file'] ### 
        these_params['random_lens_file'] = params['lens_file'] ###
        these_params['source_ra'] = 'RA'
        these_params['source_dec'] = 'DEC'
        these_params['output'] = params['output']+'_type'+str(objtype)
        
        #call treecorr with proper config
        nn = findpairs.treecorr(these_params)
        n_vec.append(nn['npairs'])

        sys.stderr.write('Done.\n')
        #get probabilities of true/false assignment
#        P_mat.append([])
#        k_mat.append([])
#        for true_objtype in objects.types:
#            P_mat[-1].append(testP['P{}P{}'.format(objtype, true_objtype)])
#            k_mat[-1].append(testk['k{}{}'.format(objtype, true_objtype)])

    print n_vec
if __name__=="__main__":
    main()


