#!/bin/bash
#
# Automate the generation of a bibliography for Stout,
# and of a table showing the provenance of the species
# known to Cloudy.
#
# Chatzikos, June 14, 2018
#

echo "Gathering references..."
./db-ref-bib2json.pl -ni
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Inspecting error files..."
[ -f empty-files.txt ] && echo "> empty-files.txt" && cat empty-files.txt && echo
[ -f broken-bibtex.txt ] && echo "> broken-bibtex.txt" && cat broken-bibtex.txt && echo
[ -f unresolved-refs.txt ] && echo "> unresolved-refs.txt" && cat unresolved-refs.txt && echo
echo
[ ! -f empty-files.txt ] && [ -f broken-bibtex.txt ] && [ -f unresolved-refs.txt ] && echo "None found!"
echo "Done."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Preparing LaTeX table for Stout..."
./db-ref-json2tex.pl -e -f=s -nr=32
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Compiling Stout table PDF..."
source mktable-stout-refs-list.sh
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

if [ ! -f ../tsuite/auto/func_lines.spclab ];
then
	if [ ! -x ../source/cloudy.exe ];
	then
		if [[ $OSTYPE =~ "linux" ]];
		then
			nthreads=`nproc --all`
			nht=$( expr `grep '^flags\b' /proc/cpuinfo | tail -1 | grep ht | wc -l` + 1 )
			nthreads=`echo $nthreads / $nht | bc`
		elif [[ $OSTYPE =~ "darwin" ]];
		then
			nthreads=`sysctl -n hw.physicalcpu`
		fi

		echo "Compiling Cloudy with $nthreads threads..."
		cd ../source
		make -j $nthreads > /dev/null
		[ $? == 0 ] && echo "Done."
		cd -
	fi
	
	echo "Running tsuite/auto/func_lines.in..."
	cd ../tsuite/auto/; ../../source/cloudy.exe -r func_lines ; cd -
	[ $? != 0 ] && exit "Something went wrong.  Exiting..."
	echo; echo "Hit enter to continue..."
	read < /dev/stdin
fi

echo "Preparing table of species' provenance..."
./db-species-tex.pl -e -np=4 -nr=35 -f=s ../tsuite/auto/func_lines.spclab
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Compiling species' provenance table PDF..."
source mktable-species-db-list.sh
[ $? != 0 ] && exit "Something went wrong.  Exiting..."

echo ; echo "Stout PDF bibliography completed successfully"

echo ; echo "Cleaning up..."
source cleanup-stout-refs-list.sh
source cleanup-species-db-list.sh
echo "Done!"

echo ; echo "Copying PDFs to ../docs ..."
cp *.pdf ../docs
command="ls -trl ../docs | tail -n 2"
echo "> $command"
eval $command
echo "Done!"

echo ; echo "All done!"
