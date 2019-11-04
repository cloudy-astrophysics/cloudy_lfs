Readme for Doxygen

This directory contains the file needed to create doxygen documentation for the Cloudy source.  You must have doxygen, latex, and graphviz installed for this to work.

Create the documentation with the command
doxygen Doxyfile
Open the index.html file in the html directory to view the documentation.

======================================================================

Doxygen is a source code documentation system that is widely used in open source projects.  It is available on the web at http://www.stack.nl/~dimitri/doxygen/  You will need a copy of the doxygen executable on your system to create the documentation.

graphviz
Doxygen must be able to find the graphviz.  This is used to create equations from embedded LaTex.  Download graphviz from http://www.graphviz.org/ 
If you receive the error message
> sh: dot: command not found
> Problems running dot. Check your installation!
this means that doxygen cannot find graphviz.

This directory includes the setup file "Doxyfile" that is needed to run doxygen.  The Cloudy download does not include the output documentation it generates.  To create documentation run doxygen with the command
doxygen Doxyfile
in this directory.  Doxgyen will create a new html directory below this one.  The index.html file in the html directory is the top of the documentation.

The manual for doxygen can be found at http://www.stack.nl/~dimitri/doxygen/download.html#latestman
The full description of its commands is under "Special Commands" in the "Reference Manual" section.

The document file doxygen_setup_style.txt in this directory contains some notes on how Cloudy uses doxygen.

==============================================

the file Doxyfile was created with the gui that is lauchned with the command 
doxywizard Doxyfile
this is used to set the parameters for the generated output.

Good luck,
Gary Ferland
http://www.nublado.org 
