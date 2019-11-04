This specifies optional SEDs.  One of these files will be used if the command
table SED "filename.sed" 
appears.

See Hazy1 for a complete description of how to use these files.
The following is a quick reminder.

File format:
First number is photon energy, given in Rydbergs by default.
If the keyword "units" appears on any line then it accepts 
the units option.  
All energies must be in the same units.
They must be in increasing energy order (decreasing wavelength order).

The second number is f_nu by default.  If the keyword nuFnu
appears then the flux can be given in relative nu f_nu units.

Comments begin with "#" and can occur anywhere in the file.

The list of energies and fluxes are concluded with a field of stars, 
for example, *****.
Everything after the field of stars is ignored, so can be used to document the
SED.

If the SED does not cover the full energy range of Cloudy the
missing parts of the SED are given an intensity of zero.
The "extrapolate" option tells the code to extrapolate the low-energy
part of the SED to the low-energy limit of the code using a power law.

The following SED files are now included;

akn120.sed
An AGN SED, see Hazy1 for details

cool.sed
Emission from the intracluster medium in a cool core cluster.
See Hazy1 for more details.

crab08.sed
crabdavidson.sed
Two versions of the Crab Nebula SED.
See Hazy1 for more details.

Mrk509.sed
#Mrk509 SED taken from Kaastra et al. 2011 [ A&A 534, A36 (2011)]. 

NGC5548.sed
#NGC5548 SED taken from Mehdipour et al. 2015 [A&A 575, A22 (2015)]. 

pl1.sed
A alpha = -1 power law continuum extending across the entire
energy range.

rubin.sed
A modified O star SED designed for modeling the Orion Nebula.
See Hazy1 for more details.

Trapezium.sed
An edited copy of the transmitted continuum in Trapezium/Trapezium_WMbasic.tran
To update the file, recreate Trapezium/Trapezium_WMbasic.tran, and run
update_Trapezium.pl

xdr.sed
A very simple XDR SED.
See Hazy1 for more details.

==========================================================

Original coding by Joshua Schlueter, NSF REU 2012 Summer
