# quicksilver
Python program to easily set up Mercury runs, particularly circumbinary planets.
This program uses Rachel Smullen's modified version of Mercury for circumbinary planets but will also work for single star systems and unmodified version of Mercury. 
This guide will err on the side of covering basic Mercury information in order to help those using Mercury for the first time.

Setting up runs
===========
add required files to a new directory
```
import quicksilver or import quicksilver as qs
run = quick(‘directory’)
```

Running Mercury
-----
Remember to compile the code before each run


Processing Output
==========
```
binarybary(folder)
```
You can calculate the orbit elements individually, this isn’t automatically done since it is a very computationally intensive process
```
cal_elements(folder,’coordinate system’)
```
Cicumbinary stabilty example
==========


```python
from orbitplotting import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tables
%matplotlib nbagg
```


```python
folder = "/home/adam/Projects/CBPstab/coolplot/"
a = mercury(folder)
a.primary_mass = .9
a.secondary_mass=.1
a.binary_separation=1.0
a.timestep = 10
a.output_interval = 10*365
a.stop_time = 1000*365
a.max_num_bod = 4000
for d in np.arange(1.0,4,.25):
    for theta in np.arange(0,360,10.0/d):
        #print d, theta
        a.add_small(d, theta=theta)
a.build()
```


```python
folder = "/home/adam/Projects/CBPstab/coolplot/"
hdfs = glob.glob(folder+'*SM*.hdf')
plt.figure()
plt.style.use('ggplot')
for hdf in hdfs:
    p = pd.read_hdf(hdf,"bary")
    fates = pd.read_hdf(hdf,"fate")
    fate = fates.values[0]
    if fate == 'survived':
        plt.scatter(p.x[0], p.y[0])
    if fate == 'ejected':
        plt.scatter(p.x[0], p.y[0], c='r')
    if fate == 'STAR1':
        plt.scatter(p.x[0], p.y[0], c='purple')
    if fate == 'STAR2':
        plt.scatter(p.x[0], p.y[0], c='green')
plt.axis('equal')
plt.show()
```



```python
folder = "/home/adam/Projects/CBPstab/coolplot/"
hdfs = glob.glob(folder+'*SM*.hdf')
plt.figure()
plt.style.use('ggplot')
for hdf in hdfs:
    p = pd.read_hdf(hdf,"bary")
    fates = pd.read_hdf(hdf,"fate")
    fate = fates.values[0]
    if fate == 'survived':
        plt.scatter(p.x.tail(1), p.y.tail(1))
    if fate == 'ejected':
        plt.scatter(p.x.tail(1), p.y.tail(1), c='r')
    if fate == 'STAR1':
        plt.scatter(p.x.tail(1), p.y.tail(1), c='purple')
    if fate == 'STAR2':
        plt.scatter(p.x.tail(1), p.y.tail(1), c='green')
plt.axis('equal')
plt.show()
```


Miscellaneous
-----
Why is it called quicksilver?

Please acknowledge the use of my code in any publication.

Any feedback is appreciated, especially suggestions or possible contributions.
