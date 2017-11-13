import numpy as np
import matplotlib.pyplot as plt
from myutils import joincorrecttype

class DataSet():

    def __init__(zprob_file, data_file, id_col):
        self.data = joincorrecttype(data_file, zp_file, id_col, id_col, float)
        self.zranges = getzgroups(self.data.columns)
        self.probabilities = ['P'+str(zrange) for zrange in self.zranges]
        self.output = ''
        
    def assignTypes(self, zindices, belows, aboves, check=False):
        
        for i in range(len(zindices)):
            probs = [self.data[self.probabilities[zindex]] for zindex in zindices[i]]

            where = np.where((np.sum(probs, axis=0) < belows[i]) & (np.sum(probs, axis=0) > aboves[i]))
        
            if check:
                plt.plot(tab['MAG_AUTO_I_d04'][where]-tab['MAG_AUTO_Z_d04'][where],
                         tab['MAG_AUTO_R_d04'][where]-tab['MAG_AUTO_I_d04'][where],
                        'o', markeredgecolor='none', markersize=4., label=i)
        if check:
            plt.legend(loc='best')
            plt.xlim(0,4)
            plt.ylim(0,4)
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

    #files
    params['source_file'] = config.get('I/O','source_file')
    params['output'] = config.get('I/O','output')

    #data
    params['source_ra']   = config.get('Columns', 'source_ra')

    #randoms
    params['source_rand_ra']     = config.get('Columns', 'source_random_ra')

    return params

def main():

    params = parseconfigs('measuresignal.config')

    objects = DataSet(params['zp_file'], params['data_file'], params['id_col'])

    #need probabilities

    #assign types (H) to objects from zprob_file based on their probabilities
    objects.assignTypes(   )

    #get P(G|H) probabilities of type G when assigned type H from test field (d04/d10/dfull)
    test_zprob_file

    #get alphas from brog/dfull/mag - function of type G
    alpha_file

    #number of pairs as function of radii for each type H
    for each_type_H:
        findpairs(type_H, lenses, radii)

    mu = 

#slope of lum func from sure ones?
#weighted lum func from assigned ones?

#choose assignment
#check original Ps

#need sources
#need measured number densities
#need number density from balrog
