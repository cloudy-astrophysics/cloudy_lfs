#!/usr/bin/perl
#
# SUMMARY:
#
# Crawl through the atomic data base (Stout) to gather the references to the
# papers the data were obtained from.  Use flag '-ni' to non-interactively
# get the ADS links ONLY.  A JSON file that holds the data is created (or
# updated) in the database directory (e.g., data/Stout/refs.json).
#
# DESCRIPTION:
#
# The script begins by reading the default Cloudy bibliography data base in
# common/.  This is needed for updates of the data base itself, as discussed
# below.
#
# For the data base of choice (Stout), the script reads both the default and
# the 'all' masterlists (*All.ini) to learn what species are available and which
# are used by default.  Commented out entries are not processed.  A JSON file is
# searched for in the database directory (data/stout/refs.json), to obtain
# results from a previous run.  If a list of species or elements is specified on
# the command line, only these species and/or elements are processed.
#
# For each species, the script opens in succession the relevant data files,
# parses the comment section at the end of the file, and attempts to isolate the
# references.  In Stout, these are either contained between the last two fields
# of stars (i.e., lines that contain only stars) in the file, or if there is
# only one field of stars after the end of data, by the section header
# 'Reference'.  This last limitation is due to the fact that in some Stout files
# there are several fields of stars to essentially comment out unwanted sections
# of the data.
#
# Once the references have been isolated, each line is parsed in succession.
# Because there is no single format (tabs, colons, or a single space may be
# used) to distinguish the reference from any preceding text, and because
# occassionally the reference may be split across two lines, the script reports
# its best-guess list of references to the user along with the original comments
# for validation.  The user may process the list to correct a reference, or to
# provide an ADS link or bibcode, or simply to delete that reference.  Additional
# references may be entered.  This may be useful for references to private
# communications that are not picked up.
#
# Alternatively, the script may be used in non-interactive mode with the flag
# '-ni'.  In this case, it gathers only the ADS links from each file, and
# requires no user validation.  This is intended as a quick way to update the
# JSON file.  NB NB: Pre-existing references in the JSON file that are not
# ADS links are preserved by this process.  Files that contain no references
# are reported in the file 'empty-files.txt'.
#
# Note that for Stout, the script expects each reference to occupy its own line,
# or at least to be at the end of the line, separated from leading text by means
# of a tab, a colon, or an equals sign.  One reference is expected per line.
#
# Once the list is finalized, the script goes over the references one at a time,
# attempts to parse it into a list of authors, a publication year, and journal
# details, and searches for the citation in the Cloudy bibliography.  If found,
# the script proceeds to the next reference.
#
# If the reference isn't found, the script proceeds to query ADS.  ADS possesses
# 3 data bases (AST, PHYS, and PRE -- listed on the browser at the top of the
# abstract request page as 'Astronomy', 'Physics', and 'arXiv e-prints').  Each
# data base is queried successively, the results are reported to the user, and
# checked against the input reference for a match.  If a match is not found, the
# user is asked to provide a bibcode (presumably because the used filter didn't
# catch the right citation), or proceed to the next ADS data base.  If none of
# the ADS data bases provide a match, the user is asked to enter a bibcode by
# hand, or skip the reference.  In that event, the reference is reported as an
# unmatched reference in the file 'unresolved-refs.txt', in the directory from
# which the script is run.
#
# Once a bibcode is obtained, the Cloudy bibliography is queried for it, to
# make sure that the citation is truly absent from the bibliography data base.
# If so, ADS is queried anew for the BibTeX entry of the citation.  Knowledge
# of the data base where the match was found is important, as querying a
# different data base could lead to acquiring the wrong BibTeX.  (The last data
# base listed on the URL for the BibTeX request is the one that takes effect.
# This is why the original ADS query was done separately for each ADS data base:
# to preserve the origin of the match.  Note that all three data bases could be
# searched at the same time in the original query -- as on the browser.)  The
# BibTeX entry is then added to the bibliography data base, along with a unique
# human-readable cross-reference.  (The bibliography data on memory are also
# appropriately updated.)
#
# Note that in some rare occassions, the BibTeX entry may be missing some
# information.  If possible, these entries are corrected from information at
# hand, and these entries are reported in the file 'broken-bibtex.txt', in
# the directory from which the script is run.
#
# In addition, the Stout data files are appended with links to ADS abstracts,
# if they do not exist therein already.
#
# After the entire data base is parsed, the data are stored in a JSON file in
# the data base directory (e.g., data/stout/refs.json).  Subsequent runs begin
# by reading this file prior to processing the atomic data base, as discussed
# above.
#
# The error files 'empty-files.txt', 'broken-bibtex.txt', and 'unresolved-refs.txt'
# should be searched for and inspected, if they exist, after each run.  If no
# problems were found, the relevant files will not be created.
#
# This script depends on a few CPAN packages that may need to be installed:
#
#	http://search.cpan.org/dist/Astro-ADS
#	http://search.cpan.org/~makamaka/JSON-2.90/
#	http://search.cpan.org/~ambs/Text-BibTeX-0.71/
#
# You may install these by following these instructions:
#
# 1. Install cpanm to simplify the installation process, e.g.,
# 	$ cpan App::cpanminus
# 2. Install the modules one at a time, e.g.,
# 	$ cpanm JSON
# 	$ cpanm Astro::ADS
# 	$ cpanm Text::BibTeX
# 	$ cpanm LWP::Simple
#
# Usage examples:
# 1. Process the entire Stout data base, that is, all files of all species.
# 	$ ./db-ref-bib2json.pl
#
# 2. Process the energy levels in the Stout data base of Al.  All defined Al
#    ionization stages will be processed.
# 	$ ./db-ref-bib2json.pl -db=stout -ds=e al
#
# 3. Process all Stout Fe^{+} files 
# 	$ ./db-ref-bib2json.pl -db=stout -ds=e fe_2
#
# 4. Process the entire Stout data base non-interactively.
# 	$ ./db-ref-bib2json.pl -ni
#
# Chatzikos, 2016-Jan-27
# Chatzikos, 2016-Mar-28
# 	Implement non-interactive mode.
#

use warnings;
use strict;

use Cwd;
use File::Basename;
use Getopt::Long;
use JSON;

use Astro::ADS::Query;
use Text::BibTeX;
use LWP::Simple;

use BiblioToTeX;

# Immediately flush output
#
$| = 1;

#
# Global data
#
my %db_title =
(
	stout	=>	"Stout",
	chianti	=>	"Chianti",
	lamda	=>	"LAMDA",
);
my $masterlist_dir = "masterlist";
my %masterlists =
(
	stout => {
			default => "Stout.ini",
			all => "StoutAll.ini"
		},
	chianti => {
			default => "CloudyChianti.ini",
			all => "CloudyChiantiAll.ini"
		},
	lambda => {
			default => "Lamda.ini",
		},
);
my %suffix =
(
	stout => {
			energy => ".nrg",
			trans => ".tp",
			coll => ".coll",
		},
	chianti => {
			energy => ".elvlc",
			trans => ".wgfa",
			coll => ".splups",
		},
	lamda => {
			energy => ".dat"
		},
);
my %comment_sentinel =
(
	stout => '^#\s*',
	chianti => '^%\s*',
);
my %end_of_data =
(
	stout => '^\*\*+\s*$',
	chianti => '^\s*-1\s*$',
	lambda => '^! NOTE',
);

my $empty_files = "empty-files.txt";
my $broken_bibtex = "broken-bibtex.txt";
my $unresolved_ref = "unresolved-refs.txt";

my @ads_dbs = qw/ AST PHYS PRE /;

my $verbose = 1;
my $interactive;


#################################################################################
#			COMMAND-LINE INTERFACE					#
#################################################################################
sub printUsage
{
	my $program = &File::Basename::basename( $0 );
	die	"Usage:\n"
	 .	"\t$program [-h] [-f] [-ni] [-db=<db>] [-ds=<ds>] [species]\n"
	 .	"where:\n"
	 .	"\t-h      :\t show this message and exit\n" 
	 .	"\t-f      :\t force ADS query on references (not bibcodes)\n" 
	 .	"\t-ni     :\t non-interactive mode; only gather ADS links; invalidates -f option\n"
	 .	"\t-db=<db>:\t data base, one of: stout, chianti, lamda, all [OPTIONAL]\n"
	 .	"\t-ds=<ds>:\t data set, one of: e (energy), t (transitions), c (collisions) [OPTIONAL; DEFAULT: all]\n"
	 .	"\tspecies :\t species or element to parse from database [OPTIONAL];\n"
	 .	"\t         \t e.g., c_1 or c; if omitted, all database species are processed\n";
}

sub getInput
{
	my( $showUsage, $database, $dataset, $forceADSquery, $nonInteractive );

	&Getopt::Long::GetOptions
		(
			"h"	=>	\$showUsage,
			"db=s"	=>	\$database,
			"ds=s"	=>	\$dataset,
			"f"	=>	\$forceADSquery,
			"ni"	=>	\$nonInteractive,
		)
	or die "Error in command-line arguments\n";

	&printUsage()	if( defined( $showUsage ) );

	my $db = &BiblioToTeX::check_database( $database );

	if( defined( $nonInteractive ) and defined( $forceADSquery ) )
	{
		warn "Warning: both -ni and -f flags given.  Unsetting -f...\n";
		$forceADSquery = undef;
	}
	if( defined( $nonInteractive ) )
	{
		$verbose = undef;
	}

	# Set global value
	if( defined( $nonInteractive ) )
	{
		$interactive = undef;
	}
	else
	{
		$interactive = 1;
	}

	if( not defined( $dataset ) )
	{
		$dataset = 'a';
	}
	else
	{
		if( $dataset !~ m/^(e|t|c)$/ )
		{
			die "Error: Unrecognized dataset option:\t $dataset\n";
		}
	}

	my @species;
	push( @species, @ARGV );

	return	( $forceADSquery, $db, $dataset, \@species );
}


#################################################################################
#					UTILITY					#
#################################################################################
sub set_globals
{
	my $curdir = &Cwd::cwd();
	$empty_files = $curdir ."/". $empty_files;
	$broken_bibtex = $curdir ."/". $broken_bibtex;
	$unresolved_ref = $curdir ."/". $unresolved_ref;
	unlink $empty_files;
	unlink $broken_bibtex;
	unlink $unresolved_ref;
}

sub create_file
{
	my $file = shift;
	open FILE, ">$file"
	  or die "create_file: Error: Could not open file:\t $file\n";
	close FILE
	   or die "Error: Could not close file:\t $file\n";
}

sub append_file
{
	my( $filename, $contents ) = @_;

	open BIB2, ">>$filename"
	  or die "Error: Could not open: $filename\n";

	print BIB2 @$contents;

	close BIB2
	  or warn "Warning: Could not close: $filename\n";
}

sub report_verbatim
{
	my $text = shift;
	print "<<<\n";
	if( ref( $text ) eq "SCALAR" )
	{
		print $text;
	}
	elsif( ref( $text ) eq "ARRAY" )
	{
		print @$text;
	}
	print "<<<END\n";
}

sub get_response
{
	my( $shown_report ) = @_;
	my $resp;
	GET_RESP:
	{
		print "Is this correct?"
			if( defined( $shown_report ) );
		print " [y or n; enter for y]\t";
		$resp = <STDIN>;
		chomp( $resp );
		if( $resp !~ m/^(y|n|)$/i )
		{
			print "Unacceptable option: $resp\n";
			goto GET_RESP;
		}
	}
	$resp = "y"
		if( $resp =~ m/^y/i or $resp eq "" );
	$resp = "n"
		if( $resp =~ m/^n/i );
	return	$resp;
}

sub valid_data_in_hash
{
	my( $hr, @keys ) = @_;

	my $found;
	foreach my $key ( @keys )
	{
		if( exists( $$hr{$key} ) and defined( $$hr{$key} )
			and $$hr{$key} ne "" )
		{
			$found = 1;
			last;
		}
	}
	return	$found;
}

sub remove_from_array
{
	my( $array, $list_of_ind ) = @_;

	return
		if( not defined( $array ) or not @$array );
	return
		if( not defined( $list_of_ind ) or not @$list_of_ind );

	foreach my $ind ( reverse( @$list_of_ind ) )
	{
		splice( @$array, $ind, 1 );
	}
}


#################################################################################
#				AUXILIARY FILES					#
#################################################################################
sub update_aux_file
{
	my( $file, $line ) = @_;
	&create_file( $file )	if( not -e $file );
	my $contents = &BiblioToTeX::read_contents( $file );
	my $added;
	if( not grep( m/$line/, @$contents ) )
	{
		my @lines;
		push( @lines, $line );
		&append_file( $file, \@lines );
		$added = 1;
	}
	return	$added;
}


#################################################################################
#				MASTERLISTS					#
#################################################################################
sub read_this_masterlist
{
	my( $list, $this_file, $species ) = @_;
	my $contents = &BiblioToTeX::read_contents( $this_file );
	shift( @$contents );	# Get rid of magic number
	foreach my $line ( @$contents )
	{
		next	if( $line =~ m/^#/ );
		chomp( $line );
		my( $this_species, $levels ) = split( /\s+/, $line );
		next	if( $this_species =~ m/\dd$/ );
		my( $elm, $ion ) = split( /_/, $this_species );
		$$species{$this_species}{element} = $elm;
		$$species{$this_species}{ion} = $ion;
		$$species{$this_species}{list} = $list;
		#	print "$this_file\t $list\t $elm\t $ion\t <=\t $line\n"	if( $elm eq "fe" );
	}
}

sub read_masterlists
{
	my $db = shift;
	my %species;
	foreach my $list ( sort keys %{ $masterlists{$db} } )
	{
		my $this_file = $masterlist_dir ."/".  $masterlists{$db}{$list};
		&read_this_masterlist( $list, $this_file, \%species );
	}

	#	foreach my $species ( sort keys %species )
	#	{
	#		print "$species:\t". $species{$species}{list} ."\n";
	#	}
	#	die;

	return	\%species;
}


#################################################################################
#				SPECIES HASH					#
#################################################################################
sub get_species_subset
{
	my( $species_list, $all_species, $db ) = @_;

	my %species_hash;
	if( @$species_list )
	{
		foreach my $this_species ( @$species_list )
		{
			my $sp_el = $this_species;
			my $isElement;
			if( $this_species !~ '_' )
			{
				$sp_el .= '_';
				$isElement = 1;
			}
			my $found;
			foreach my $species ( sort keys %$all_species )
			{
				if( $species =~ m/^$sp_el/ )
				{
					$found = 1;
					$species_hash{$species} = $$all_species{$species};
					last
						if( not defined( $isElement ) );
					print "Species:\t$species\n"
						if( defined( $isElement ) );
				}
			}
			die	"Error: Requested species ($this_species) "
			 .	"not found in database masterlist: $db\n"
				if( not defined( $found ) );
		}
	}
	else
	{
		%species_hash = %{ $all_species };
	}

	#	foreach my $species ( sort keys %species_hash )
	#	{
	#		print "$species:\t". $species_hash{$species}{file} ."\n";
	#	}

	return	{ %species_hash };
}

sub clean_hash
{
	my $all_species = shift;

	my @keys_to_delete = qw/ inCloudyBib link ads_db /;

	foreach my $sp ( sort keys %$all_species )
	{
		my $refs = $$all_species{$sp}{ref};
		foreach my $file( sort keys( %$refs ) )
		{
			foreach my $key ( @keys_to_delete )
			{
				my @all_bibcodes;
				foreach my $ref ( @{ $$refs{$file} } )
				{
					if( exists( $$ref{$key} ) )
					{
						delete( $$ref{$key} );
					}

					delete $$ref{name}
						if( exists( $$ref{name} ) and
							( not defined( $$ref{name} )
								or $$ref{name} eq "" ) );

					push( @all_bibcodes, $$ref{bibcode} )
						if( exists( $$ref{bibcode} ) and
						    defined( $$ref{bibcode} ) );
				}

				# Remove duplicate entries
				# These may arise when both the full reference and
				# the ADS link appear in the data file
				#
				my %all_bibcodes;
				foreach my $bib ( @all_bibcodes )
				{
					next	if( not defined( $bib ) or $bib eq "" );
					if( not exists( $all_bibcodes{$bib} ) )
					{
						$all_bibcodes{$bib} = 1;
					}
					else
					{
						++$all_bibcodes{$bib};
					}
				}

				#	my $json = &JSON::to_json( $$refs{$file}, {pretty => 1} );
				#	print "BEFORE:\n$json\n";

				foreach my $bib ( sort keys %all_bibcodes )
				{
					next	if( not defined( $bib ) or $bib eq "" );
					next	if( $all_bibcodes{$bib} == 1 );

					my @indices;
					for( my $i = 0; $i < @{ $$refs{$file} }; $i++ )
					{
						push( @indices, $i )
						  if( exists( $$refs{$file}[$i]{bibcode} ) and
							( defined( $$refs{$file}[$i]{bibcode} )
							or $$refs{$file}[$i]{bibcode} eq $bib ) );
					}

					my $ind_to_drop;
					foreach my $index ( @indices )
					{
						if( $$refs{$file}[$index]{name} eq "" )
						{
							$ind_to_drop = $index;
							last;
						}
					}
					splice( @{ $$refs{$file} }, $ind_to_drop, 1 )
						if( defined( $ind_to_drop ) );
				}

				#	$json = &JSON::to_json( $$refs{$file}, {pretty => 1} );
				#	print "AFTER:\n$json\n";
				#	die	if( length( $json ) > 3 );
			}
		}
	}
}


#################################################################################
#				REFERENCES / BIBCODES				#
#################################################################################
sub replace_ands
{
	my $author = shift;
	my @nands = $$author =~ m/ and /;
	#	print $$author ."\t =>";
	for( my $i = 0; $i < @nands-1; $i++ )
	{
		$$author =~s/,* and /, /;
	}
	$$author =~ s/,* and /, & /;
	#	print $$author ."\n";
}

sub resolve_authors
{
	my( $author_list, $bt ) = @_;

	$author_list =~ s/^\s*//;
	$author_list =~ s/,\s*$//;
	&replace_ands( \$author_list );
	$$bt{author} = $author_list;

	my @authors;
	if( $author_list =~ m/et\s*al\.*/ )
	{
		$author_list =~ s/et\s*al\.*//;
		$author_list =~ s/,\s*$//;
	#	push( @authors, $author_list );
		$$bt{authors_etal} = 1;
	}
	else
	{
		if( $author_list =~ m/\& / )
		{
			if( $author_list !~ m/,\s*\&/ )
			{
				$author_list =~ s/\& /, /;
			}
			else
			{
				$author_list =~ s/,\s*\&/, /;
			}
		}
	}

	my @fields = split( /,/, $author_list );
	foreach my $f ( @fields ) { $f =~ s/(^ *| *$)//; }
	pop( @fields ) while( @fields and $fields[ $#fields ] eq '' );
	#	foreach my $f ( @fields ) { print "auth field: '$f'\n"; }
	do
	{
		my $auth = shift( @fields );
		$auth =~ s/^\s*//;
		my $givennames;
		if( $auth !~ m/\./ )
		{
			$givennames = shift( @fields );
			if( defined( $givennames ) )
			{
				if( $givennames =~ m/[A-Z]([a-zA-Z]+|&)/ )
				{
					unshift( @fields, $givennames );
				}
				elsif( $givennames =~ m/^[A-Z]\.?\s*([A-Z]\.?|)/ )
				{
					$givennames =~ s/^\s*//;
					$givennames =~ s/\s*$//;
					$auth .= ", ". $givennames;
				}
			}
		}
		else
		{
			if( $auth =~ m/\s/ )
			{
				#	print "THIS: $auth\n";
				my @fields = split( /\s+/, $auth );
				my $lastname = shift( @fields );
				if( $lastname =~ m/\./ )
				{
					my @givennames;
					push( @givennames, $lastname );
					push( @givennames, shift( @fields ) )
						while( @fields and $fields[0] =~ m/^[A-Z]\.$/ );
					$lastname = join( ' ', @fields );
					$auth = $lastname .', '. join( ' ', @givennames );
				}
				else
				{
					my $givennames = shift( @fields );
					if( $givennames !~ m/\./ )
					{
						$lastname .= " ". $givennames;
						$givennames = "";
						$givennames = join( '', @fields )
							if( @fields );
						$auth = $lastname;
						$auth .= ", ". $givennames
							if( $givennames ne "" );
					}
					else
					{
						$auth = $lastname .", ". $givennames;
					}
				}
				#	print "THAT: $auth\n";
			}
			else
			{
				# This can only happen for refs like G.J.Ferland
				my @f = split( /\./, $auth );
				$auth = pop( @f );
				$auth .= ", ";
				foreach my $n ( @f )
				{
					$auth .= "$n. ";
				}
				$auth =~ s/ $//;
			}
		}
		$auth =~ s/^\s*\&\s+//	if( $auth =~ m/\s*\&\s+/ );
		$auth =~ s/^\s+//;
		push( @authors, $auth );
	} while( @fields );

	#	foreach my $auth ( @authors ) { print "author:\t". $auth ."\n"; }
	$$bt{authors} = \@authors;
}

sub resolve_journal
{
	my( $refname, $bt ) = @_;

	my @fields = split( /,/, $refname );
	#	foreach my $f ( @fields ) { print "field: $f\n"; }

	$$bt{year} = shift( @fields );
	$$bt{journal} = shift( @fields );
	if( $$bt{journal} =~ m/\s\d+$/ )
	{
		my @f = split( /\s+/, $$bt{journal} );
		unshift( @fields, pop( @f ) );
		$$bt{journal} = join( ' ', @f );
	}
	$$bt{volume} = shift( @fields );
	$$bt{pages} = shift( @fields );
	foreach my $key ( qw/ year volume pages / )
	{
		next	if( not exists( $$bt{$key} ) or
				not defined( $$bt{$key} ) );
		$$bt{$key} =~ s/\s//g;
	}
	$$bt{journal} =~ s/^\s+//;
	$$bt{year} =~ s/[a-z]$//i;
	$$bt{pages} =~ s/\.$//;
}

sub pop_last_author
{
	my( $bt ) = @_;

	my $last_author;
	if(  not exists( $$bt{authors_etal} ) )
	{
		$last_author = pop( @{ $$bt{authors} } );
		$$bt{author} =~ s/(, |)\s*\Q$last_author\E(,|)//;
	}
	return	$last_author;
}

sub parse_reference
{
	my $ref = shift;
	print "Parsing reference:\n";
	&BiblioToTeX::print_bibentry( $ref );

	my $refname = $$ref{name};

	# Isolate author list
	my @fields = split( /\d\d\d\d/, $refname );
	shift( @fields ) while( @fields and $fields[0] eq "" );

	my %bt;
	my $author_list = shift( @fields );
	&resolve_authors( $author_list, \%bt );

	# Remove author list from reference
	$refname =~ s/$author_list//;

	&resolve_journal( $refname, \%bt );

	if( ( not defined( $bt{journal} ) or $bt{journal} eq "" )
		and $bt{volume} and $bt{pages} )
	{
		$bt{journal} = &pop_last_author( \%bt );
	}

	print "Resolved into:\n";
	&BiblioToTeX::print_bibentry( { %bt },
				\@BiblioToTeX::bibtex_basic_fields );

	return	\%bt;
}

sub resolve_bibcode
{
	my $bibcode = shift;
	my %pub;
	$pub{year} = substr( $bibcode, 0, 4, "" );
	$bibcode =~ s/\./ /g;
	my @fields = split( /\s+/, $bibcode );
	$pub{journal} = shift( @fields );
	$pub{volume} = shift( @fields );
	$pub{pages} = shift( @fields );
	$pub{pages} =~ s/[A-Z]$//
		if( defined( $pub{pages} ) );
	return	{ %pub };
}

sub form_bibcode_ending
{
	my $bt = shift;
	my $fauth_uc = uc( substr( $$bt{authors}[0], 0, 1 ) );
	my $bibcode = '';
	$bibcode .= '\.+'. $$bt{volume}
		if( exists( $$bt{volume} ) and defined( $$bt{volume} ) );
	$bibcode .= '\.+'. $$bt{pages}
		if( exists( $$bt{pages} ) and defined( $$bt{pages} ) );
	$bibcode .= $fauth_uc .'$';
	#	print "bibcode_match = $bibcode\n";
	return	$bibcode;
}


#################################################################################
#				ADS QUERIES					#
#################################################################################
sub query_ADS_one_db
{
	my( $bibcode_match, $year, $authors, $ads_db ) = @_;

	print "Querying ADS db $ads_db ...";


	my $query = Astro::ADS::Query->new(
		Authors		=>	$authors,
		AuthorLogic	=>	"AND",
		StartYear	=>	$year,
		EndYear		=>	$year );
	$query->{OPTIONS}{db_key} = $ads_db;

	my $results = $query->querydb();
	my @papers = $results->papers();

	if( not @papers )
	{
		print "\n";
		return;
	}

	print	"\nADS query of $ads_db database yielded "
	  .	@papers ." matching citations:\n";

	foreach my $paper ( @papers )
	{
		print $paper->summary( format => "ASCII" ) ."\n";
	}

	my( $bibcode, $title );
	foreach my $paper ( @papers )
	{
		if( $paper->bibcode() =~ m/$bibcode_match$/ )
		{
			$bibcode = $paper->bibcode(); 
			$title = $paper->title();
			last;
		}
	}

	if( defined( $bibcode ) and defined( $title ) )
	{
		$title = &BiblioToTeX::drop_brackets( $title, 1 );
		print "Best match:\n";
		print "\tbibcode:\t $bibcode\n";
		print "\ttitle  :\t \"$title\"\n";

		if( &get_response( 1 ) eq "n" )
		{
			$bibcode = undef;
		}
	}

	return	$bibcode;
}

sub enter_bibcode_by_hand
{
	my( $bt, $ads_db, $skip ) = @_;

	print "Please enter bibcode for:\n";
	&BiblioToTeX::print_bibentry( $bt );
	print "Please enter bibcode";
	print " (or hit enter for next ADS db)"
		if( defined( $skip ) );
	print ":\t";
	my $bibcode = <STDIN>;
	chomp( $bibcode );
	if( $bibcode eq "" and defined( $skip ) )
	{
		return;
	}
	elsif( $bibcode eq "" )
	{
		print "\t => Ignoring this reference\n";
	}
	else
	{
		if( not defined( $ads_db ) )
		{
			GET_ADS_DB:
			{
				print "Please enter ADS db:\t";
				$ads_db = <STDIN>;
				chomp( $ads_db );
				if( join( ' ', @ads_dbs ) !~ $ads_db )
				{
					print	"Illegal choice. Options are:\t"
					  .	join( ',', @ads_dbs ) ."\n";
					goto GET_ADS_DB;
				}
			}
		}
	}

	return	( $bibcode, $ads_db );
}

sub query_ADS
{
	my $bt = shift;

	my $bibcode_match = &form_bibcode_ending( $bt );

	my @authors = @{ $$bt{authors} };
	$authors[0] = "^". $authors[0];
	#	print "'@authors'\n";

	my( $bibcode, $ads_db );
	foreach my $this_ads_db ( @ads_dbs )
	{
		$ads_db = $this_ads_db;
		$bibcode = &query_ADS_one_db( $bibcode_match, $$bt{year},
						\@authors, $this_ads_db );
		last
			if( defined( $bibcode ) );

		( $bibcode, $ads_db ) =
			&enter_bibcode_by_hand( $bt, $ads_db, 1 );

		last
			if( defined( $bibcode ) );

		$bibcode = $ads_db = undef;
	}

	( $bibcode, $ads_db ) = &enter_bibcode_by_hand( $bt, $ads_db )
		if( not defined( $bibcode ) );

	#	print "Got bibcode:\t $bibcode\n";
	#	die;

	return	( $bibcode, $ads_db );
}

sub get_ads_bibtex
{
	my( $bibcode, $ads_db ) = @_;
	my $req = "http://adsabs.harvard.edu/cgi-bin/nph-bib_query"
		. "?bibcode=$bibcode"
		. "&data_type=BIBTEX"
		. "&db_key=$ads_db"
		. "&nocookieset=1";
	#	print "$req\n";
	return	&LWP::Simple::get( $req );
}


#################################################################################
#				ACQUIRE REFERENCES				#
#################################################################################
sub preprocess_refs
{
	my $refs = shift;

	#
	# Drop references that appear more than once
	#
	my @indices;
	for( my $i = 0; $i < @$refs-1; $i++ )
	{
		my $found;
		foreach my $ind ( @indices )
		{
			$found = 1
				if( $ind == $i );
		}
		next	if( defined( $found ) );
		my $key = "name";
		   $key = "bibcode"
			if( not defined( $$refs[$i]{name} )
				or $$refs[$i]{name} eq "" );
		for( my $j = $i+1; $j < @$refs; $j++ )
		{
			next
				if( not defined( $$refs[$j]{$key} )
					or $$refs[$j]{$key} eq "" );
			push( @indices, $j )
				if( $$refs[$i]{$key} eq $$refs[$j]{$key} );
		}
	}
	@indices = sort @indices;

	&remove_from_array( $refs, \@indices );
}

sub report_refs
{
	my $refs = shift;

	print "\nNumber of refs acquired:\t". @$refs ."\n";
	for(my $iref = 0; $iref < @$refs; $iref++ )
	{
		printf "%d:\t %s",
			$iref+1, $$refs[ $iref ]{name};
		printf "\t%s", $$refs[ $iref ]{link}
			if( exists( $$refs[ $iref ]{link} ) and
			   defined( $$refs[ $iref ]{link} ) );
		print "\n";
	}
}
sub enter_ref_by_hand
{
	my $ref = shift;

	print	"Please enter one of:\n"
	  .	"    d to delete entry\n"
	  .	" or r to enter a reference (as in, 'r: Ferland et al 2013, RMxAA, 49, 137'):\n"
	  .	" or l to enter a link (as in: 'l: http://adsabs.harvard.edu/abs/2007ApJ...654.1171A'):\n"
	  .	" or b to enter a bibcode (as in: 'b: 2007ApJ...654.1171A'):\t";
	my $read_ref = <STDIN>;
	chomp( $read_ref );
	$read_ref =~ s/^\s*//;
	if( $read_ref =~ m/^r:/ )
	{
		$read_ref =~ s/^r:\s*//;
		$$ref{name} = $read_ref;
		delete( $$ref{link} )
			if( exists( $$ref{link} ) );
	}
	elsif( $read_ref =~ m/^l:/ )
	{
		$read_ref =~ s/^l:\s*//;
		$$ref{link} = $read_ref;
		delete( $$ref{name} )
			if( exists( $$ref{name} ) );
	}
	elsif( $read_ref =~ m/^b:/ )
	{
		$read_ref =~ s/^b:\s*//;
		$$ref{bibcode} = $read_ref;
		delete( $$ref{name} )
			if( exists( $$ref{name} ) );
	}
	else
	{
		delete $$ref{name}	if( exists( $$ref{name} ) );
		delete $$ref{link}	if( exists( $$ref{link} ) );
		delete $$ref{bibcode}	if( exists( $$ref{bibcode} ) );
	}
}

sub confirm_ref_postparse_all
{
	my( $iref, $ref ) = @_;

	printf "%d:", $iref;
	if( $$ref{name} ne "" )
	{
		print "\t$$ref{name}\n";
	}
	else
	{
		print "\t$$ref{link}\n";
	}

	return
		if( &get_response( 1 ) =~ m/^y$/i );

	return	&enter_ref_by_hand( $ref );
}

sub confirm_refs
{
	my( $filename, $contents_orig, $refs ) = @_;

	print "\nProcessing File:\t $filename\n";
	print "=" x42 ."\n";
	print "Original Comments: ";
	&report_verbatim( $contents_orig );

	&report_refs( $refs );

	if( &get_response( 1 ) eq 'n' )
	{
		my @indices;
		for(my $iref = 0; $iref < @$refs; $iref++ )
		{
			&confirm_ref_postparse_all( $iref+1, $$refs[ $iref ] );
			push( @indices, $iref )
			  if( not exists( $$refs[ $iref ]{name} ) and
			      not exists( $$refs[ $iref ]{link} ) );
		}
		&remove_from_array( $refs, \@indices );

		print "\nAdd more references?";
		my @added_refs;
		while( &get_response() eq 'y' )
		{
			my %ref;
			&enter_ref_by_hand( \%ref );
			push( @added_refs, { %ref } )
				if( defined( &valid_data_in_hash( \%ref,
						qw/name link bibcode/ ) ) );
			print "\nAdd more references?";
		}
		push( @$refs, @added_refs );
		#	my $json = &JSON::to_json( $refs, {pretty => 1} );
		#	print "BEFORE:\n$json\n";
	}
}

sub parse_stout_comments
{
	my( $filename, $contents_orig ) = @_;

	my @contents = @$contents_orig;
	if( grep( /Reference/, @contents ) )
	{
		shift( @contents ) while ( @contents and
						$contents[0] !~ m/Reference/ );
	}

	my @refs;
	foreach my $line ( @contents )
	{
		chomp( $line );

		$line =~ s/^#Reference://;
		$line =~ s/$comment_sentinel{stout}//;
		my $ref = $line;

		if( $line =~ m/\t\w+$/ )
		{
			( undef, $ref ) = split( /\t+/, $line );
		}
		elsif( $line =~ m/:\s*/ )
		{
			my @ref = split( ':', $line );
			shift( @ref )	# Get rid of text, as in "# text: ref"
				if( $ref[0] !~ m/http/ );
			$ref = join( ':', @ref );
		}
		elsif( $line =~ m/=/ )
		{
			( undef, $ref ) = split( '=', $line );
		}

		#	print "ISIT:\t". defined( $ref )."\n";
		if( not defined( $ref ) or $ref eq "" )
		{
			$ref = $line;
		}
		else
		{
			$ref =~ s/^\s*\"//;
			$ref =~ s/\"\s*$//;
			#	print "'$ref'\n";
		}
		#	elsif( defined( $ref ) )
		#	{
		#		print "ref= '$ref'\n";
		#	}
		#	print	"line:\t $line\n"
		#	  .	"ref :\t $ref\n";

		$ref =~ s/^\s*//;
		$ref =~ s/\s*$//;

		my %ref;
		if( $ref =~ m#http://adsabs# )
		{
			#	print "ref:\t '$ref'\n";
			my( $name, $http ) = split( 'http', $ref );
			#	print "name= '$name'\t http= '$http'\n";
			$ref{name} = $name;
			$ref{link} = 'http'. $http;
			my @f = split( '/', $ref{link} );
			$ref{bibcode} = pop( @f );
		}
		elsif( defined( $interactive ) )
		{
			if( $ref =~ m/NIST/ )
			{
				$ref{name} = $ref;
			}
			else
			{
				# Citations employ at least 3 commas
				#
				my $ncommas = () = $ref =~ m/,/g;
				$ref{name} = $ref
					if( $ncommas >= 3 );
			}
		}
		push( @refs, \%ref )
			if( %ref );
	}

	return	\@refs;
}

sub parse_chianti_comments
{
	my( $filename, $contents_orig ) = @_;

	my @contents = @$contents_orig;

	my @refs;
	foreach my $line ( @contents )
	{
		chomp( $line );
		$line =~ s/$comment_sentinel{chianti}//;
		my $ref = $line;
		if( $line =~ m/:\s/ )
		{
			( undef, $ref ) = split( ':', $line );
		}

		#	print "ISIT:\t". defined( $ref )."\n";
		if( not defined( $ref ) or $ref eq "" )
		{
			$ref = $line;
		}
		else
		{
			$ref =~ s/^\s*\"//;
			$ref =~ s/\"\s*$//;
			#	print "'$ref'\n";
		}
		#	elsif( defined( $ref ) )
		#	{
		#		print "ref= '$ref'\n";
		#	}
		#	print	"line:\t $line\n"
		#	  .	"ref :\t $ref\n";

		$ref =~ s/^\s*//;

		my %ref;
		if( $ref =~ m#http://adsabs# )
		{
			#	print "ref:\t '$ref'\n";
			my( $name, $http ) = split( 'http', $ref );
			#	print "name= '$name'\t http= '$http'\n";
			$ref{name} = $name;
			$ref{link} = "http" . $http;
			my @f = split( '/', $ref{link} );
			$ref{bibcode} = pop( @f );
		}
		elsif( $ref =~ m/NIST/ )
		{
			$ref{name} = $ref;
		}
		else
		{
			# Citations employ at least 3 commas
			#
			my $ncommas = () = $ref =~ m/,/g;
			$ref{name} = $ref
				if( $ncommas >= 3 );
		}
		push( @refs, \%ref )
			if( %ref );
	}

	#
	# Replace long NIST references with version numbers
	#
	my @indices;
	for( my $i = 0; $i < @refs; $i++ )
	{
		#	print "NAME:\t". $refs[$i]{name} ."\n";
		if( $i < @refs-1
			and $refs[$i]{name} =~
			'Martin, W.C., Sugar, J., Musgrove, A., Dalton, G.R., 1995,'
			and $refs[$i+1]{name} =~ 'NIST Database for Atomic Spectroscopy'
			and $refs[$i+1]{name} =~ 'Version 1.0'
			and $refs[$i+1]{name} =~ 'NIST Standard Reference Database 61'
			)
		{
			$refs[$i]{name} = 'NIST v1.0';
			$i++;
			push( @indices, $i );
		}
		elsif( $refs[$i]{name} =~
			'Fuhr, J.R., et al., "NIST Atomic Spectra Database" Ver. 2.0, March 1999, NIST Physical Reference Data' )
		{
			$refs[$i]{name} = 'NIST v2.0';
		}
		elsif( $refs[$i]{name} =~ 'Version 3.0.3 of the NIST online database' )
		{
			$refs[$i]{name} = 'NIST v3.0.3';
		}
	}
	#	die "@indices\n";

	&remove_from_array( \@refs, \@indices );

	return	\@refs;
}

sub get_bibcodes_update_biblio
{
	my( $forceADSquery, $species, $refs ) = @_;

	return	if( not defined( $refs) or not @$refs );

	print "\nAcquiring bibcodes:\n"
		if( defined( $verbose ) );
	my @bibcodes;
	for( my $i = 0; $i < @$refs; $i++ )
	{
		my $this_ref = $$refs[ $i ];
		print "\nRef ". ($i+1) .":\n"
			if( defined( $verbose ) );
		if( exists( $$this_ref{bibcode} ) and
				defined( $$this_ref{bibcode} ) )
		{
			print "Have bibcode => $$this_ref{bibcode}\n"
				if( defined( $verbose ) );
			&add_ref_to_bibliography( $this_ref );
			push( @bibcodes, $$this_ref{bibcode} );
		}
		elsif( $$this_ref{name} !~ m/NIST/i and
			$$this_ref{name} !~ m/Chianti/i and
			$$this_ref{name} !~ m/(unpublished|private comm)/ )
		{
			my $bt = &parse_reference( $this_ref );
			if( not defined( $forceADSquery ) )
			{
				$$this_ref{bibcode} =
					&BiblioToTeX::get_bibcode_from_Cloudy_bib( $bt );
				if( defined( $$this_ref{bibcode} ) )
				{
					$$this_ref{inCloudyBib} = 1;
				}
				else
				{
					( $$this_ref{bibcode}, $$this_ref{ads_db} ) =
						&query_ADS( $bt );
				}
			}
			else
			{
				( $$this_ref{bibcode}, $$this_ref{ads_db} ) =
					&query_ADS( $bt );
				$$this_ref{inCloudyBib} = 1
					if( defined(
						&BiblioToTeX::get_bibcode_from_Cloudy_bib( $bt ) ) );
			}

			if( defined( $$this_ref{bibcode} ) and
				$$this_ref{bibcode} ne "" )
			{
				push( @bibcodes, $$this_ref{bibcode} );
				&add_ref_to_bibliography( $this_ref )
					if( not exists( $$this_ref{inCloudyBib} ) );
			}
			elsif( $$this_ref{bibcode} eq "" )
			{
				my $line = sprintf( "%s:\t%s\n",
						$species,
						$$this_ref{name} );
				my $added = &update_aux_file( $unresolved_ref, $line );
				print "\t=> Added unresolved ref to: $unresolved_ref\n"
					if( defined( $added ) );
			}
		}
		elsif( defined( $verbose ) )
		{
			print	"$$this_ref{name}\n";
		}
	}

	return	\@bibcodes;
}

sub fix_broken_record
{
	my( $bibcode, $bibdef ) = @_;

	my $bt = &resolve_bibcode( $bibcode );

	my $entry = &BiblioToTeX::convert_bibtex_to_hash( $bibdef );

	my @missing_records;
	foreach my $field ( qw/journal year volume pages/ )
	{
		if( ( not exists( $$entry{$field} ) or not defined( $$entry{$field} ) )
			and exists( $$bt{$field} ) and defined( $$bt{$field} ) )
		{
			&BiblioToTeX::insert_bibtex_field( $bibdef, $field, $$bt{$field} );
			push( @missing_records, $field );
		}
	}

	if( @missing_records )
	{
		my $line = sprintf( "%15s\t=>\t%s\n",
				$bibcode,
				join( ",", @missing_records ) );
		my $added = &update_aux_file( $broken_bibtex, $line );
		print "ADS BibTeX def missing records:\t $line\n"
		 if( defined( $added ) );
	}
}

sub create_crossref
{
	my $bibdef = shift;

	my $entry = &BiblioToTeX::convert_bibtex_to_hash( $bibdef );

	my( $lastname, $givennames ) =
		&BiblioToTeX::resolve_author_name( $$entry{authors}[0] );
	$lastname =~ s/\s//g;
	$lastname =~ s/\W//g;
	my $crossref = $lastname . $$entry{year};
	if( defined( &BiblioToTeX::get_crossref_abbrv( $crossref ) ) )
	{
		my $this_crossref;
		my $ichar = 0;
		do
		{
			$this_crossref = $crossref . chr( $ichar + ord('b') );
			$ichar++;
		}
		while( $ichar < 26 and
			defined( &BiblioToTeX::get_crossref_abbrv( $this_crossref ) ) );
		$crossref = $this_crossref;
	}

	my $crossref_bibdef =
	sprintf(
		"@%s{%s,\n"
	   .	"  crossref = \"%s\"\n"
	   .	"}\n",
		$entry->type,
	   	$crossref,
		$entry->key
		);
	print "\t=> Created BibTeX crossref key:\t $crossref\n";

	return	( $crossref_bibdef, $entry->key, $crossref );
}

sub add_ref_to_bibliography
{
	my( $ref ) = @_;

	return	if( exists( $$ref{inCloudyBib} ) );

	my $bt = &resolve_bibcode( $$ref{bibcode} );
	#	print "$$bt{year}\n";

	if( defined( &BiblioToTeX::is_bibcode_in_Cloudy_bib(
					$$ref{bibcode}, $$bt{year} ) ) )
	{
		print	"BibTeX for $$ref{bibcode}"
		  .	" already in Cloudy biblio\n"
		   if( defined( $interactive ) );
		return;
	}

	my $resp;
	if( exists( $$ref{ads_db} ) and defined( $$ref{ads_db} ) )
	{
		$resp = &get_ads_bibtex( $$ref{bibcode}, $$ref{ads_db} );
	}
	else
	{
		foreach my $ads_db ( @ads_dbs )
		{
			$resp = &get_ads_bibtex( $$ref{bibcode}, $ads_db );
			if( defined( $resp ) and $resp ne "" )
			{
				print 	"ADS database queried: $ads_db,"
				  .	" for: $$ref{bibcode}, got:\n";
				&report_verbatim( $resp );
				if( $resp =~ m/$$ref{bibcode}/ )
				{
					$$ref{ads_db} = $ads_db;
					last;
				}
				else
				{
					if( &get_response( 1 ) eq "y" )
					{
						$$ref{ads_db} = $ads_db;
						last;
					}
				}
			}
		}
	}

	if( not defined( $resp ) )
	{
		print "Did not find $$ref{bibcode}\n";
	}
	else
	{
		my @bibdef;
		foreach my $line ( split( "\n", $resp ) )
		{
			next if $line =~ /^\s*$/;
			next if $line =~ /^(Retrieved|Query)/;
			push( @bibdef, "$line\n" );
		}
		print "\nAcquired BibTeX: ";
		&report_verbatim( \@bibdef );

		&fix_broken_record( $$ref{bibcode}, \@bibdef );

		my( $crossref_bibdef, $bibcode, $crossref ) =
			&create_crossref( \@bibdef );

		&BiblioToTeX::add_crossref_abbrv( $crossref, $bibcode );

		my @bibdefs_new;
		push( @bibdefs_new, $crossref_bibdef );
		push( @bibdefs_new, @bibdef );
		push( @bibdefs_new, "\n" );
		&append_file( $BiblioToTeX::bibliography, \@bibdefs_new );
		print	"\t=> Added bibtex entry for $$ref{bibcode} to "
		  .	"$BiblioToTeX::bibliography\n";
		&BiblioToTeX::update_biblio_hash( \@bibdef );
	}
}

sub update_datafile
{
	my( $bibcodes, $contents, $filename ) = @_;
	return	if( not defined( $bibcodes ) or not @$bibcodes );

	my @lines;
	foreach my $bibcode ( @$bibcodes )
	{
		my $bibcode_esc = $bibcode;
		$bibcode_esc =~ s/\./\\./g;
		#	print "bibcode = $bibcode\t bibcode_esc = $bibcode_esc\n";
		if( not grep( /$bibcode_esc/, @$contents ) )
		{
			my $line = "# http://adsabs.harvard.edu/abs/" . $bibcode ."\n";
			push( @lines, $line );
		}
	}
	return	if( not @lines );

	&append_file( $filename, \@lines );

	print "\t=> Updated datafile: $filename\n";
}

sub get_number_of_fields_of_stars
{
	my $contents = shift;
	my $nlines_stars = 0;
	foreach my $line ( @$contents )
	{
		$nlines_stars++
			if( $line =~ m/^\*\*+ *$/ );
	}
	#	print "$nlines_stars\n";
	return	$nlines_stars;
}

sub get_data
{
	my( $db, $contents ) = @_;
	my $nlines_stars = &get_number_of_fields_of_stars( $contents );
	#	print "$nlines_stars\n";
	my( $ilines_stars, @data ) = ( 0 );
	do
	{
		push( @data, shift( @$contents ) )
			while( defined( $$contents[0] ) and
					$$contents[0] !~ $end_of_data{$db} );
		shift( @$contents );	# get rid of field of stars
		$ilines_stars++;
	}
	while( $ilines_stars < $nlines_stars );
	return	\@data;
}

sub get_stout_refs
{
	my( $db, $contents ) = @_;

	my( $ilines_stars, @refs ) = ( 0 );
	do
	{
		@refs = ();
		unshift( @refs, pop( @$contents ) )
			while( defined( $$contents[ -1 ] ) and
					$$contents[ -1 ] !~ $end_of_data{$db} );
		pop( @$contents );
		$ilines_stars++;
		#$ilines_stars = 2
		#	if( @refs and $refs[0] =~ 'Reference' );
	}
	while( $ilines_stars < 1 );

	pop( @$contents )
		if( defined( $$contents[-1] ) and
			$$contents[-1] =~ $end_of_data{$db} );

	return	\@refs;
}

sub report_empty_files
{
	my( $filename, $data, $refs ) = @_;

	#	print "DATA:\t >>>@$data<<<\n";
	#	print "REFS:\t >>>@$refs<<<\n";

	my $added;
	if( defined( $data ) and not @$data and defined( $refs ) and not @$refs and $filename !~ m/.coll$/ )
	{
		$added = &update_aux_file( $empty_files, "$filename: NO data, NO refs\n" );
	}
	elsif( defined( $data ) and @$data and defined( $refs ) and not @$refs )
	{
		$added = &update_aux_file( $empty_files, "$filename: has data, NO refs\n" );
	}
	elsif( defined( $data ) and not @$data and defined( $refs ) and @$refs and $filename !~ m/.coll$/ )
	{
		$added = &update_aux_file( $empty_files, "$filename: NO data, has refs\n" );
	}
	print "\t=> Reported missing data / refs to:\t $empty_files\n"
		if( defined( $added ) );
}

sub get_file_references
{
	my( $forceADSquery, $species, $db, $filename ) = @_;

	if( 0 )
	{
		print "pwd:\t". &Cwd::cwd() ."\n";
		print "filename= '$filename'\n";
		print "db = '$db'\n";
	}
	my $contents = &BiblioToTeX::read_contents( $filename );
	my $data;
	if( $db eq 'stout' )
	{
		shift( @$contents );	# Get rid of magic number
		my $refs = &get_stout_refs( $db, $contents );
		$data = $contents;
		$contents = $refs;
		&report_empty_files( "../data/stout/".$filename, $data, $refs );
	}
	else
	{
		$data = &get_data( $db, $contents );
	}
	return
		if( not @$data );
	#	print "data:\t". @$data ."\n";
	#	print "rest:\t". @$contents ."\n";

	my $refs;
	if( $db eq "stout" )
	{
		$refs = &parse_stout_comments( $filename, $contents );
	}
	elsif( $db eq "chianti" )
	{
		$refs = &parse_chianti_comments( $filename, $contents );
	}
	elsif( $db eq "lamda" )
	{
		$refs = &parse_lamda_comments( $filename, $contents );
	}


	if( defined( $interactive ) )
	{
		#	&report_refs( $refs );
		&preprocess_refs( $refs );
		#	&report_refs( $refs );
		&confirm_refs( $filename, $contents, $refs );
		&report_refs( $refs );
	}

	my $bibcodes = &get_bibcodes_update_biblio( $forceADSquery,
				$species, $refs );
	&update_datafile( $bibcodes, $contents, $filename )
		if( defined( $interactive ) and $db eq "stout" );

	return	$refs;
}

sub pick_datatypes
{
	my( $db, $ds ) = @_;

	my @datatypes;
	if( $ds eq 'a' )
	{
		@datatypes = sort keys %{ $suffix{$db} };
	}
	else
	{
		if( $ds eq 'e' )
		{
			push( @datatypes, "energy" );
		}
		elsif( $ds eq 't' )
		{
			push( @datatypes, "trans" );
		}
		elsif( $ds eq 'c' )
		{
			push( @datatypes, "coll" );
		}
	}

	return	\@datatypes;
}

sub get_stored_ref
{
	my( $file_ref, $ref_type, $refs_arr ) = @_;

	my $this_hr;
	foreach my $hr ( @{ $refs_arr } )
	{
		if( exists( $$hr{$ref_type} ) and
			$$hr{$ref_type} eq $$file_ref{$ref_type} )
		{
			$this_hr = $hr;
			last;
		}
	}

	return	$this_hr;
}

sub update_refs_data
{
	my( $refs_data, $sp, $datatype, $file_refs ) = @_;

	return
		if( not defined( $file_refs ) );

	for( my $iref = 0; $iref < @$file_refs; $iref++ )
	{
		my $ref = $$file_refs[ $iref ];

		my $this_hr;
		if( exists( $$ref{bibcode} ) and $$ref{bibcode} ne "" )
		{
			#	print "ref bibcode =\t $$ref{bibcode}\n";
			$this_hr = &get_stored_ref( $ref, 'bibcode',
					$$refs_data{$sp}{ref}{$datatype} );
		}
		elsif( exists( $$ref{name} ) and $$ref{name} ne "" )
		{
			#	print "ref name=\t $$ref{name}\n";
			$this_hr = &get_stored_ref( $ref, 'name',
					$$refs_data{$sp}{ref}{$datatype} );
		}
		#	print defined( $this_hr ) ? "yes\n" : "no\n";

		if( not defined( $this_hr ) )
		{
			push( @{ $$refs_data{$sp}{ref}{$datatype} }, $ref );
		}
	}

	if( defined( $interactive ) )
	{
		# Now do the reverse: prune all stored data that are not in the
		# current version of the file
		my @rm_index;
		#	print @{ $$refs_data{$sp}{ref}{$datatype} }."\n";
		for( my $iref = 0; $iref < @{ $$refs_data{$sp}{ref}{$datatype} }; $iref++ )
		{
			my $ref = $$refs_data{$sp}{ref}{$datatype}[ $iref ];

			my $this_hr;
			if( exists( $$ref{bibcode} ) and $$ref{bibcode} ne "" )
			{
				#	print 'bibcode=  '. $$ref{bibcode} ."\n";
				$this_hr = &get_stored_ref( $ref, 'bibcode', $file_refs );
			}
			elsif( exists( $$ref{name} ) and $$ref{name} ne "" )
			{
				#	print 'name=  '. $$ref{name} ."\n";
				$this_hr = &get_stored_ref( $ref, 'name', $file_refs );
			}
			#	print defined( $this_hr ) ? "yes\n" : "no\n";

			if( not defined( $this_hr ) )
			{
				push( @rm_index, $iref );
			}
		}

		&remove_from_array( $$refs_data{$sp}{ref}{$datatype}, \@rm_index );
	}
}

sub get_references
{
	my( $forceADSquery, $db, $ds, $species_hash ) = @_;

	my $datatypes = &pick_datatypes( $db, $ds );

	foreach my $sp ( sort keys %$species_hash )
	{
		if( defined( $interactive ) )
		{
			print "\nProcess species: $sp?";
			my $resp = &get_response();
			next	if( $resp eq 'n' );
		}
		else
		{
			print "Processing species:\t $sp\n";
		}

		my $basename = $$species_hash{$sp}{element} ."/". $sp ."/". $sp;
		foreach my $datatype ( @$datatypes )
		{
			my $suff = $suffix{$db}{$datatype};
			my $file = $basename . $suff;
			if( defined( $interactive ) )
			{
				print "\nProcess file: $file?";
				my $resp = &get_response();
				next	if( $resp eq 'n' );
			}

			my $file_refs = &get_file_references(
					$forceADSquery, $sp, $db, $file );
			&update_refs_data( $species_hash, $sp, $datatype, $file_refs );
		}
	}
}


#################################################################################
#				MAIN PROGRAM					#
#################################################################################

my( $forceADSquery, $db, $ds, $species_list ) = &getInput();

&Astro::ADS::Query::ads_mirror( 'adsabs.harvard.edu' );

&BiblioToTeX::set_globals();
&set_globals();

&BiblioToTeX::load_cloudy_bibliography();

foreach my $db ( @$db )
{
	my $curdir = &Cwd::cwd();

	my $dbdir = $BiblioToTeX::data_dir ."/$db";
	chdir $dbdir
	   or die "Error: Could not chdir to db: $dbdir\n";

	my $all_species = &read_masterlists( $db );
	#	my @all = keys %$all_species;
	#	die "$dbdir:\t @all\n";

	&BiblioToTeX::load_json( $all_species );
	#	my @all = keys %$all_species;
	#	die "$dbdir:\t @all\n";

	my $species_hash = &get_species_subset( $species_list, $all_species, $db );
	#	my @all = keys %$species_hash;
	#	die "$dbdir:\t @all\n";

	&get_references( $forceADSquery, $db, $ds, $species_hash );
	&clean_hash( $all_species );
	&BiblioToTeX::store_json( $all_species );

	chdir $curdir
	   or die "Error: Could not chdir back to: $curdir\n";
}

print "For missing data / refs, see:\t $empty_files\n"
 if( -s $empty_files );
print "For broken bibtex records, see:\t $broken_bibtex\n"
 if( -s $broken_bibtex );
print "For unresolved refs, see:\t $broken_bibtex\n"
 if( -s $unresolved_ref );
