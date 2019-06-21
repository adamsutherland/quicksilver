import imp
qs = imp.load_source('quicksilver', '/home/adam/Code/analysis/quicksilver/quicksilver.py')
import glob as glob
import pandas as pd

print
print("Processing with python")
folder = "./"
qs.G = 1.0
qs.rebound_process(folder)
qs.bary_to_central(folder)
qs.jacobi(folder)

s1 = pd.read_hdf(folder+'STAR1.hdf','jacobi',columns=["time","mass"])
s2 = pd.read_hdf(folder+'STAR2.hdf','jacobi',columns=["time","mass"])
# the number of columns written to the file is reduced
smalls = glob.glob(folder+'SM*.hdf')
smalls.sort()
for hdf in smalls:
        print hdf
        s = pd.read_hdf(hdf,'jacobi',columns=["time","mass","x","y","z","capom","inc","a"])
        s = qs.high_freq_pro(s1,s2,s,full=False)
        s.to_hdf(hdf,'jacobi',format="t")

smalls = glob.glob(folder+'PL*.hdf')
smalls.sort()
for hdf in smalls:
        print hdf
        s = pd.read_hdf(hdf,'jacobi',columns=["time","mass","x","y","z","capom","inc","a"])
        s = qs.high_freq_pro(s1,s2,s,full=False)
        s.to_hdf(hdf,'jacobi',format="t")
