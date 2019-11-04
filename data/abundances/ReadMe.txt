The files in this directory specify standard chemical compositions and, optionally, grains.
One of these files will be used if the command
abundances "filename.abn"
appears.

Available files
default.abn gives the default composition used if no others are specified.
Replace this file to define your own set of default abundances.

solar84.abn - abundances used by "old solar" option, used in versions 84-94
of the code.

Other .abn files - see the comments within the file for more details.

=======================================================

Format of the *.abn files

comments begin with #

The file parsing ends with a line containing a field of stars, as
****************
Everything after the field of stars is ignored, so remaining lines can be used
to document the abundances, in addition to comment lines starting with #.

The grains command can be included.  The abundance files recognize all grain command options.
This makes it possible to include grains in your mix of gas and dust.

Names of elements are followed by the abundance by number relative to hydrogen.

Hydrogen may be in specified, or may not be.  If the abundance of hydrogen is 
given and is not equal to unity the code will rescale the abundances by the entered value 
for hydrogen.  If hydrogen is not specified it is assumed to have an abundance
of 1.

All elements heavier than hydrogen are turned off before the abn files are read.
An element is turned on if it appears in the abn file.
If an element is not included in the *.abn file it will not be included in the calculation.

==========================================================

Original coding by Joshua Schlueter, NSF REU 2012 Summer
