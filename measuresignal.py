import numpy as np
import matplotlib.pyplot as plt
import json
import ConfigParser
import findpairs
from astropy.io import fits
from myutils import joincorrecttype
from matplotlib.pyplot import cm

class DataSet():

    def __init__(self, zprob_file, data_file, id_col, output):
        self.data = joincorrecttype(data_file, zprob_file, id_col, id_col, float)
        self.zranges = getzgroups(self.data.columns)
        self.probabilities = ['P'+str(zrange) for zrange in self.zranges]
        self.output = output

    def assignTypes(self, zindices, belows, aboves, check=False):
        color = cm.rainbow(np.linspace(0,1,len(zindices)+1))
        
        if check:
            plt.plot(self.data['MAG_AUTO_I_d04']-self.data['MAG_AUTO_Z_d04'],
                     self.data['MAG_AUTO_R_d04']-self.data['MAG_AUTO_I_d04'],
                     'o', markeredgecolor='none', c=color[-1], markersize=2., label='all')

        self.types = {}      
        for i,c in zip(range(len(zindices)), color[:-1]):
    
            where = np.where(self.data)
            data = self.data[where]
            
            for j in range(len(zindices[i])):
                probs = [data[self.probabilities[zindex]] for zindex in zindices[i][j]]
    
                where = np.where((np.sum(probs, axis=0) < belows[i][j]) & (np.sum(probs, axis=0) > aboves[i][j]))
                data = data[where]

            #save data in type dictionary
            self.types[i] = data

            if check:
                #check colors for each type
                plt.plot(data['MAG_AUTO_I_d04']-data['MAG_AUTO_Z_d04'],
                         data['MAG_AUTO_R_d04']-data['MAG_AUTO_I_d04'],
                         '^', markeredgecolor='none', c=c, markersize=4.,
                         label='type '+str(i))

                print "Type {}: {} objects".format(i, len(data))
                
        if check:
            plt.legend(loc='best')
            plt.xlim(-1,4)
            plt.ylim(-1,4)
            plt.xlabel('MAG_AUTO_I_D04 - MAG_AUTO_Z_D04')
            plt.ylabel('MAG_AUTO_R_D04 - MAG_AUTO_I_D04')
            plt.savefig(self.output+'_typecolors.png')

def getzgroups(columns):
    zgroups = []
    for column in columns:
        if column[:2]=="P[":
            zgroups.append(eval(column[1:]))
    ind = np.argsort([np.min(zgrp) for zgrp in zgroups])
    return [zgroups[i] for i in ind]


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

    #Table Column Names
    params['zp_id_col'] = config.get('Columns','zp_id_col')
    params['data_id_col'] = config.get('Columns','data_id_col')

    #Assignment Criteria
    params['redshift_index'] = json.loads(config.get('Assignments','redshift_index'))
    params['below'] = json.loads(config.get('Assignments','below'))
    params['above'] = json.loads(config.get('Assignments','above'))

    return params

def main():

    params = parseconfigs('measuresignal.config')

    objects = DataSet(params['zp_file'], params['data_file'], params['zp_id_col'], params['output'])
    lenses = fits.open(params['lens_file'])[1].data
    
    #assign types to objects from zprob_file based on their probabilities
    objects.assignTypes(params['redshift_index'], params['below'], params['above'], check=True)

    #test_zprob_file must have same redshift ranges
    testP = Table.read(params['test_zprob_file'])
    testk = Table.read(params['test_k_file'])

    #measured number densities
    n_vec = []
    P_mat = []
    k_mat = []
    for objtype in objects.types:
        #for each type - calculate correlation with lenses
        these_params = {}

        #write type objects to table
        this_table = 
        these_params['source_file'] = this_table
        these_params['lens_file'] = params['lens_file']                                                                       these_params['random_source_file'] = '' 
        these_params['random_lens_file'] = params['random_lens_file']
        these_params['source_ra'] = ''
        these_params['source_dec'] = ''
        these_params['output'] = params['output']+'_type'+objtype.id
        
        #call treecorr with proper config
        nn = findpairs.treecorr(these_params)
        n_vec.append(nn[''])
        
        #get probabilities of true/false assignment
        P_mat.append([]}
        k_mat.append([])
        for true_objtype in objects.types:
            P_mat[-1].append(testP['P{}P{}'.format(objtype, true_objtype)])
            k_mat[-1].append(testk['k{}{}'.format(objtype, true_objtype)])


if __name__=="__main__":
    main()


