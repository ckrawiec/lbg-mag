import numpy as np
import matplotlib.pyplot as plt
import json
import ConfigParser
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

    #need probabilities

    #assign types to objects from zprob_file based on their probabilities
    objects.assignTypes(params['redshift_index'], params['below'], params['above'], check=True)



if __name__=="__main__":
    main()

    
"""
    #get P(G|H) probabilities of type G when assigned type H from test field (d04/d10/dfull)
    test_zprob_file

    #get alphas from brog/dfull/mag - function of type G
    alpha_file

    #number of pairs as function of radii for each type H
    for each_type_H:
        findpairs(type_H, lenses, radii)

    mu = 
"""
#slope of lum func from sure ones?
#weighted lum func from assigned ones?

#choose assignment
#check original Ps

#need sources
#need measured number densities
#need number density from balrog
