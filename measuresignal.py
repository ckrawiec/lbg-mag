import findpairs
import dNdMu
import numpy as np
import matplotlib.pyplot as plt
import json
import sys
import ConfigParser
import glob
from astropy.io import fits
from astropy.table import Table, vstack, join
from myutils import joincorrecttype
from matplotlib.pyplot import cm

class DataSet():

    def __init__(self, zprob_tab, data_file, id_col, output, z_col=None):
        data_tab = Table.read(data_file)
        self.data = join(zprob_tab, data_tab)
        #for old output version in z-prob master branch
        #self.data = joincorrecttype(data_file, zprob_file, id_col, id_col, float)
        self.zranges = getzgroups(self.data.columns)
        self.probabilities = ['P'+str(zrange) for zrange in self.zranges]
        self.output = output
        self.idcol = id_col
        if z_col:
            self.zcol = z_col

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
    params['balrog_data_files'] = config.get('I/O','balrog_data_files')
    params['test_zp_file'] = config.get('I/O','test_zp_file')
    params['test_data_file'] = config.get('I/O','test_data_file')
    params['output'] = config.get('I/O','output')
    
    #Table Column Names
    params['zp_id_col'] = config.get('Columns','zp_id_col')
    params['data_id_col'] = config.get('Columns','data_id_col')
    params['test_data_id_col'] = config.get('Columns','test_data_id_col')
    params['test_z_col'] = config.get('Columns','test_z_col')

    #Assignment Criteria
    params['redshift_index'] = json.loads(config.get('Assignments','redshift_index'))
    params['below'] = json.loads(config.get('Assignments','below'))
    params['above'] = json.loads(config.get('Assignments','above'))
    params['true_ranges'] = json.loads(config.get('Assignments','true_ranges'))

    return params

def main(args):
    #parse config file
    params = parseconfigs(args[1])

    #get all z-prob files
    zp_files = glob.glob(params['zp_file'])
    sys.stderr.write('Found {} z-prob files.\n'.format(len(zp_files)))
    zp_tab = Table.read(zp_files[0])
    sys.stderr.write('Reading files...')
    for zp_file in zp_files[1:]:
        zp_tab = vstack([zp_tab, Table.read(zp_file)])
    sys.stderr.write('Done.\n Found {} objects.\n'.format(len(zp_tab)))

    #create target and lens data objects
    objects = DataSet(zp_tab, params['data_file'], params['zp_id_col'], params['output'])
    lenses = fits.open(params['lens_file'])[1].data
    
    #assign types to objects from zprob_file based on their probabilities
    objects.assignTypes(params['redshift_index'], params['below'], params['above'], check=True)

    #files from test region
    test = DataSet(Table.read(params['test_zp_file']), params['test_data_file'],
                   params['test_data_id_col'], params['output']+'_test', z_col=params['test_z_col'])

    test.assignTypes(params['redshift_index'], params['below'], params['above'], check=False)
    
    #measured quantities, vectors and matrices
    n_vec = []
    P_mat = []
    k_mat = []
    for objtype in objects.types:
        sys.stderr.write('Working on type {}...'.format(objtype))
        
        #write type objects to table
        this_table = '{}_type{}.fits'.format(params['output'], objtype)
        objects.types[objtype].write(this_table)
        
        #for each type - create parameters for treecorr 
        these_params = {}
        these_params['source_file'] = this_table
        these_params['lens_file'] = params['lens_file'] 
        these_params['random_source_file'] = params['lens_file'] ### 
        these_params['random_lens_file'] = params['lens_file'] ###
        these_params['source_ra'] = 'RA'
        these_params['source_dec'] = 'DEC'
        these_params['output'] = params['output']+'_type'+str(objtype)
        
        #calculate correlation with lenses - call treecorr
        nn = findpairs.treecorr(these_params)
        n_vec.append(nn['npairs'])

        sys.stderr.write('Done.\n')
                
        P_mat.append([])
        k_mat.append([])
                 
        for true_objtype in objects.types:
            """
            get assigned types of detected objects
            dn/dmu = [N_det(mu_G) - N_det(1)] / (mu_G - 1)
            get probabilities of true/false assignment
            Using test file, true z, assign types to test targets, get P(H|G) from test templates
            Run dNdMu with mu_G
           
            k_mat[-1].append(testk['k{}{}'.format(objtype, true_objtype)])
            """
            z_photo = test.data[test.zcol]
            zmin = np.min(params['true_ranges'][true_objtype])
            zmax = np.max(params['true_ranges'][true_objtype])
            z_true = np.where((z_photo >= zmin) & (z_photo <= zmax))
            ids = test.data[test.idcol][z_true]

            this_z_photo = test.types[objtype][test.zcol]
            this_z_true = np.where((this_z_photo >= zmin) & (this_z_photo <= zmax))
            these_ids = test.types[objtype][test.idcol][this_z_true]

            #Assigned to this type while truly in this group
            both = float(len(set(ids).intersection(these_ids)))
            
            #Save number of objects called objtype when actually true_objtype
            P_mat[-1].append(both/float(len(ids)))


            #itab = tab.find('tab')
            #tabnum = tab[itab+3:itab+5]

    print n_vec
    print P_mat



if __name__=="__main__":
    main(sys.argv)


