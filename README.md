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

Miscellaneous
-----
Why is it called quicksilver?

Please acknowledge the use of my code in any publication.

Any feedback is appreciated, especially suggestions or possible contributions.