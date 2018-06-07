import numpy as np
import sys
import glob
from astropy.io import fits
from astropy.table import Table, vstack, join
from myutils import joincorrecttype

def gatherfiles(zpfiles):
    #get all z-prob files
    zp_files = glob.glob(zpfiles)
    sys.stdout.write('Found {} z-prob files.\n'.format(len(zp_files)))
    sys.stdout.write('Reading files...')
    zp_tabs = [Table.read(zp_file) for zp_file in zp_files]
    zp_tab = vstack(zp_tabs)
    sys.stdout.write('Done.\n Found {} objects.\n'.format(len(zp_tab)))
    sys.stdout.flush()
    return zp_tab

def getzgroups(columns):
    zgroups = []
    for column in columns:
        if column[:2]=="P[":
            zgroups.append(eval(column[1:]))
    ind = np.argsort([np.min(zgrp) for zgrp in zgroups])
    out = [zgroups[i] for i in ind]
    return out

class DataSet:
    def __init__(self, datafile,  output=None,
                 zpfiles=None, zcol=None, magcol=None, sizecol=None, fluxcol=None,
                 idcol=None, typecol=None, racol=None, deccol=None):
        sys.stdout.write("Creating DataSet from {}.\n".format(datafile))
        if zpfiles:
            zprobtab = gatherfiles(zpfiles)
            data_tab = Table.read(datafile)
            self.data = join(zprobtab, data_tab)
            self.zranges = getzgroups(self.data.columns)
            self.probabilities = ['P'+str(zrange) for zrange in self.zranges]
        else:
            self.data = fits.open(datafile)[1].data

        sys.stdout.write("    DataSet has {} objects.\n".format(len(self.data)))
        sys.stdout.flush()
        
        self.output = output
        self.idcol = idcol
        if fluxcol:
            self.fluxcol = fluxcol
        if magcol:
            self.magcol = magcol

        if sizecol:
            self.sizecol = sizecol
        if typecol:
            self.typecol = typecol
        if zcol:
            self.zcol = zcol

        #positions
        if racol and deccol:
            self.positions = zip(self.data[racol], self.data[deccol])
