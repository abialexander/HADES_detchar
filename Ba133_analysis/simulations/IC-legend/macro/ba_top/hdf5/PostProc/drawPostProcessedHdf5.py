# draw energy histograms from each detector

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("mplstyle.txt")


if(len(sys.argv) != 2):
    print('Usage: drawPostProcessedHdf5.py [processed.hdf5]')
    sys.exit()

filename = sys.argv[1]
plotname = filename.strip('hdf5')

df =  pd.read_hdf(filename, key="procdf")

print(df)

bins = np.arange(0, 4000, 1) #want rougly binwidth to be det resolution, e.g. 0.1keV 
plt.figure()
for det, detdf in df.groupby('detID'):
    (detdf["energy"]*1000).hist(bins=bins, histtype="step", label="detector {}".format(det))
plt.xlabel("Energy [keV]")
plt.ylabel("Counts")
plt.gca().set_xlim(0, 450)
plt.gca().set_ylim(10, plt.gca().get_ylim()[1])
plt.gca().grid(False)
plt.yscale("log")
#plt.legend(frameon=False, loc='upper right')
plt.savefig('example.png')
#plt.show()


binwidth = 0.15 #0.1 kev = rough min resolution
energies = df['energy']
energies = energies*1000
plt.figure()
plt.hist(energies, bins = np.arange(0, max(energies) + binwidth, binwidth), histtype = 'step')
plt.xlabel("Energy [keV]")
plt.ylabel("Counts")
plt.xlim(0, 450)
plt.yscale("log")
plt.savefig(plotname+'.png')

