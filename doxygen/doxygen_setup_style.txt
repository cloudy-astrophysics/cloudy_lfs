to change options in setup file (Doxyfile) do 
doxywizard Doxyfile
mode - optimize for C output

expert tab has many options 

to run doxygen and create output do
doxygen Doxyfile
in this directory

check whether INPUT parameter in Doxyfile is set to source_hot or source - this changes
source is the archived version and source_hot is the local development version

to generate the ouput doxygen also needs graphviz, another open source app.  graphviz needs
to be on the path for doxygen to find it.

Doxyfile is set up to create only html output.  other options are possible.  The html will be
in the html directory below this main directory.

Within the codebase, doxygen markup is fully contained in the headers.  All doxygen special comment blocks are using the "C style" markup.  Special comment blocks are indicated by /**.  Commands are indicated with the with the syntax \command.  All comment blocks are immediately before what they describe, with the exception of some struct members, where the comment immediately follows ("/**<" syntax).  Within special commente blocks, the \, @, &, $, #, <, >, % characters must be escaped using a preceeding \

The following commands are in use:
\file <filename> (description)  - description for a file, will appear in the output above any descriptions of items contained in the file.
\verbatim  - Doxygen outputs text enclosed in verbatim/endverbatim tags exactly (ie, preserving whitespace and newlines)
\endverbatim
\param [in|out|in,out] <parameter name> (description) - describe a parameter of a function.  Parameter name is the name of the variable and does not include type.  [in|out|in,out] is optional.
\post (description) - describe the post conditions for a function
\return (description) - describe what the function returns (in output "Returns "+description is printed)
\author Joe Blow

latex in line 
\f$ - same as in line $  - opposite is another \f$

The ideal declaration should look like the following:
=======================================================

 /**
  routine_name This is the long description of what routine_name does, appears after 
	the function name
  \brief this is routine brief description - only 1 line long
  \param iz  a description of what parameter 1 is for
  \param [in] in  a description of what parameter 2 is for
  \param [out] *out description
  \author Joe Blow
  \return explain the return value
 */ 
double routine_name(long int iz, 
  long int in ,
  double *out );

=======================================================

The manual for doxygen can be found: http://www.stack.nl/~dimitri/doxygen/download.html#latestman
The full description of these commands and all others are under "Special Commands" in the "Reference Manual" section.
