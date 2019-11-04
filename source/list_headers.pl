#!/usr/bin/perl
# program to find out the C files using the header(.h) files.

# this will be list of headers
$hfiles='headers.txt';

# List of header files  and C files using them  
$result='listfiles.list'; 

print "A list of header files will be placed in ", $hfiles , "\n";
print "All source files using each header will be placed in ", $result , "\n" ;
print "This will take a while...\n\n";

#Start of main
#system "ls *.h > $hfiles";
open(ioHEADERS,">$hfiles"); #ioHEADERS-Output file
open(LFILE,">$result");

while(defined($header=glob("*.h")))
{
  print ioHEADERS "$header\n";
}
close(ioHEADERS);

open(IHFILE,"$hfiles"); #IHFILE-Input file
while(<IHFILE>)
{
  $header=$_;
  print LFILE "$header";
  $header=~s/\n//;  
  while(defined($input=glob("*.cpp")))  #Scanning through the C files
  {
     open(CFILE,"$input");
     while(<CFILE>)
     {
       if($_=~/$header/)
       {
         print LFILE "\t$input";
       }     
     }
  }
  print LFILE "\n\n"; 
}
close(IHFILE);
close(LFILE);
#End of program

  	
