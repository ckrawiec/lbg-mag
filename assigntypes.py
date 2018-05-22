import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm

def assignTypes(self, zindices, belows, aboves, check=True):
    color = cm.rainbow(np.linspace(0,1,len(zindices)+1))
        
    if check:
        plt.plot(self.data[self.magcol.format('I')]-self.data[self.magcol.format('Z')],
                 self.data[self.magcol.format('R')]-self.data[self.magcol.format('I')],
                 'o', markeredgecolor='none', c=color[-1], markersize=1.5, label='all')

    self.types = {}
    
    #cycle through redshift groups    
    for i,c in zip(range(len(zindices)), color[:-1]):
        mask = np.array([True] * len(self.data))

        #cycle through this redshift group
        for j in range(len(zindices[i])):
            probs = [self.data[self.probabilities[zindex]] for zindex in zindices[i][j]]    
            where = np.array((np.sum(probs, axis=0) < belows[i][j]) & (np.sum(probs, axis=0) > aboves[i][j]))
            mask = mask & where

        #save data mask in type dictionary
        self.types[i] = mask

        if check:
            #check colors for each type
            plt.plot(self.data[mask][self.magcol.format('I')]-self.data[mask][self.magcol.format('Z')],
                     self.data[mask][self.magcol.format('R')]-self.data[mask][self.magcol.format('I')],
                     '^', markeredgecolor='none', c=c, markersize=2.5,
                     label='type '+str(i))

        print "Type {}: {} objects".format(i, len(self.data[mask]))
                
    if check:
        #save figure
        plt.legend(loc='best')
        plt.xlim(-1,4)
        plt.ylim(-1,4)
        plt.xlabel('{} - {}'.format(self.magcol.format('I'), self.magcol.format('Z')))
        plt.ylabel('{} - {}'.format(self.magcol.format('R'), self.magcol.format('I')))
        plt.savefig(self.output+'_typecolors.png')
        plt.close()

    
