#type S
#indices in sorted zgroups on which to cut
#Pselect is chosen P
for i in range(len(zindex)):
    Pcol = [jtab[Pcols[zi]] for zi in zindex[i]]

    where = np.where((np.sum(Pcol, axis=0) < below[i]) & (np.sum(Pcol, axis=0) > above[i]))

    plt.plot(jtab['MAG_AUTO_I_d04'][where]-jtab['MAG_AUTO_Z_d04'][where],
             jtab['MAG_AUTO_R_d04'][where]-jtab['MAG_AUTO_I_d04'][where],
             'o', markeredgecolor='none', markersize=4., label=i)
plt.legend(loc='best')
plt.xlim(0,4)
plt.ylim(0,4)
plt.show()

    
