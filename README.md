This program is designed to help with the manual classification of vesicles (or other notable features) imaged by AFM methods.
Internally the vesicles are called grains.

To use this program you must first prepare the input files.
The program requires an AFM image in the .gwy or .xyz format
It also requires that locations of the grains of interest in the .dat format as can be achieved with the use of the Gwydion software
The two files MUST have the same name (excluding the type) and be in the same folder as the program itself
See the provided Example.dat and Example.gwy
There can be more than one pair of files in the folder with EVIAN.py, in that case the pairs will be handled consecutively.
