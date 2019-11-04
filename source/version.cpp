/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "date.h"
#include "version.h"
#include "service.h"

static const int CLD_MAJOR = 13;
static const int CLD_MINOR = 0;
static const int CLD_PATCH = 0;

#ifdef SVN_REVISION
static const char* svn_revision = SVN_REVISION;
#else
static const char* svn_revision = "rev_not_set";
#endif

static const char* cUrl = "$HeadURL: svn://svn.nublado.org/cloudy/trunk/source/version.cpp $";

t_version::t_version()
{
	static const char chMonth[12][4] =
		{ "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

	// first analyze the URL to determine where we live, the chVersion string is derived from that
	// the code below is based on the following naming scheme:
	//
	// /branches/c08_branch -- release branch, all bug fixes are submitted here 
	//
	// /tags/develop/c08.00_rc1 -- release candidates go here 
	//
	// /tags/release/c08.00 -- first official release
	// /tags/release/c08.01 -- first bug-fix rollup, etc... 
	//
	// /tags/patch_versions/c08.00_pl00 -- identical to /tags/release/c08.00
	// /tags/patch_versions/c08.00_pl01 -- first patch update, etc... 
	//
	// /trunk -- this will be labeled as "experimental"
	// /tags/stable -- ditto
	// /branches/* -- ditto, note that "*" can be anything except c??_branch

	vector<string> Part;
	string Url = cUrl;
	Split( Url, "/", Part, SPM_RELAX );
	if( Part.size() >= 3 )
	{
		// the last two parts are "source" and "version.cpp $", we don't need them...
		// the one before is the relevant identifier (e.g. "trunk", "newmole", "c08.01")
		// for conciseness we will refer to it below as the "branch"
		string Branch = Part[Part.size()-3];

		bool lgReleaseTag = ( Url.find("/tags/release/") != string::npos );
		bool lgPatchTag = ( Url.find("/tags/patch_versions/") != string::npos );
		bool lgDevelopTag = ( Url.find("/tags/develop/") != string::npos );
		// this expects a branch name like "c08_branch"
		lgReleaseBranch = ( Url.find("/branches/") != string::npos &&
				    Branch.size() == 10 && Branch[0] == 'c' &&
				    Branch.find("_branch") != string::npos );

		lgRelease = ( lgReleaseTag || lgPatchTag );

		// determine if this is a beta version
		string::size_type ptr;
		if( lgDevelopTag && ( ptr = Branch.find( "_rc" ) ) != string::npos )
			// this expects a branch name like "c08.00_rc1"
			sscanf( Branch.substr( ptr+3 ).c_str(), "%ld", &nBetaVer );
		else
			nBetaVer = 0;

		// beta versions are always generated from the release branch...
		if( nBetaVer > 0 )
			lgReleaseBranch = true;

		int nMajorLevel=0, nMinorLevel=0, nPatchLevel=0;

		if( lgReleaseBranch || lgRelease || nBetaVer > 0 )
		{
			// this expects a branch name starting with "c08"
			sscanf( Branch.substr(1,2).c_str(), "%d", &nMajorLevel );
			if( nMajorLevel != CLD_MAJOR )
				fprintf( ioQQQ, "PROBLEM - CLD_MAJOR mismatch, please check version.cpp\n" );
		}

		if( lgRelease || nBetaVer > 0 )
		{
			// this expects a branch name starting with "c08.01"
			sscanf( Branch.substr(4,2).c_str(), "%d", &nMinorLevel );
			if( nMinorLevel != CLD_MINOR )
				fprintf( ioQQQ, "PROBLEM - CLD_MINOR mismatch, please check version.cpp\n" );
		}
		
		if( lgPatchTag )
		{
			// this expects a branch name like "c08.01_pl02"
			sscanf( Branch.substr(9,2).c_str(), "%d", &nPatchLevel );
			if( nPatchLevel != CLD_PATCH )
				fprintf( ioQQQ, "PROBLEM - CLD_PATCH mismatch, please check version.cpp\n" );
			// c08.00_pl00 is identical to release c08.00, so pass it off as the latter...
			if( nPatchLevel == 0 )
				lgReleaseTag = true;
		}

		string pps = ( isdigit(svn_revision[0]) ) ? "r" : "";

		ostringstream oss;
		if( lgReleaseTag )
			// this expects a branch name like "c08.01"
			oss << Branch.substr(1,5);
		else if( lgPatchTag )
			// this expects a branch name like "c08.01_pl02"
			oss << Branch.substr(1,5) << " (patch level " << nPatchLevel << ")";
		else if( nBetaVer > 0 )
			// this expects a branch name like "c08.00_rc1"
			oss << Branch.substr(1,5) << " beta " << nBetaVer << " (prerelease)";
		else if( lgReleaseBranch )
			// this expects a branch name like "c08_branch"
			oss << "(" << Branch << ", " << pps << svn_revision << ", prerelease)";
		else
			// the branch name can be anything except "c??_branch"
			oss << "(" << Branch << ", " << pps << svn_revision << ", experimental)";
		chVersion = oss.str();
	}
	else
	{
		// create a default version string in case HeadURL was not expanded

		/* is this a release branch? */
		lgReleaseBranch = false;
		/* is this a release version? */
		lgRelease = false;

		/* is this a beta version?  0 for no
		 * if this is non-zero then lgRelease above should be false */
		nBetaVer = 0;

		ostringstream oss;
		if( lgRelease )
		{
			oss << setfill('0') << setw(2) << CLD_MAJOR << "." << setw(2) << CLD_MINOR;
			if( CLD_PATCH > 0 )
				oss << " (patch level " << CLD_PATCH << ")";
		}
		else if( nBetaVer > 0 )
		{
			oss << setfill('0') << setw(2) << CLD_MAJOR << "." << setw(2) << CLD_MINOR;
			oss << " beta " << nBetaVer << " (prerelease)";
		}
		else
		{
			oss << setfill('0') << setw(2) << YEAR%100 << "." << setw(2) << MONTH+1;
			oss << "." << setw(2) << DAY;
		}
		chVersion = oss.str();
	}

	ostringstream oss2;
	oss2 << setfill('0') << setw(2) << YEAR%100 << chMonth[MONTH] << setw(2) << DAY;
	chDate = oss2.str();

	char mode[8];
	if( sizeof(int) == 4 && sizeof(long) == 4 && sizeof(long*) == 4 )
		strncpy( mode, "ILP32", sizeof(mode) );
	else if( sizeof(int) == 4 && sizeof(long) == 4 && sizeof(long*) == 8 )
		strncpy( mode, "IL32P64", sizeof(mode) );
	else if( sizeof(int) == 4 && sizeof(long) == 8 && sizeof(long*) == 8 )
		strncpy( mode, "I32LP64", sizeof(mode) );
	else if( sizeof(int) == 8 && sizeof(long) == 8 && sizeof(long*) == 8 )
		strncpy( mode, "ILP64", sizeof(mode) );
	else
		strncpy( mode, "UNKN", sizeof(mode) );

	bool flag[2];
	flag[0] = ( cpu.i().min_float()/2.f > 0.f );
	flag[1] = ( cpu.i().min_double()/2. > 0. );

	/* now generate info on how we were compiled, including compiler version */
	ostringstream oss3;
	oss3 << "Cloudy compiled on " << __DATE__ << " in OS " << __OS << " using the ";
	oss3 << __COMP << " " << __COMP_VER << " compiler. Mode " << mode << ", ";
	oss3 << "denormalized float: " << TorF(flag[0]) << " double: " << TorF(flag[1]) << ".";
	chInfo = oss3.str();
}
