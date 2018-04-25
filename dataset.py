import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack, join
from myutils import joincorrecttype

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
                 zprobtab=None, zcol=None, magcol=None, sizecol=None, fluxcol=None,
                 idcol=None, typecol=None):
        print "Creating DataSet from {}.".format(datafile)
        if zprobtab:
            data_tab = Table.read(datafile)
            self.data = join(zprobtab, data_tab)
            self.zranges = getzgroups(self.data.columns)
            self.probabilities = ['P'+str(zrange) for zrange in self.zranges]
        else:
            self.data = fits.open(datafile)[1].data
            #for old output version in z-prob master branch
            #self.data = joincorrecttype(data_file, zprob_file, id_col, id_col, float)

        print "    DataSet has {} objects.".format(len(self.data))
        self.output = output
        self.idcol = idcol
        if zcol:
            self.zcol = zcol
        if magcol:
            self.magcol = magcol
        if fluxcol:
            self.fluxcol = fluxcol
        if sizecol:
            self.sizecol = sizecol
        if typecol:
            self.typecol = typecol
