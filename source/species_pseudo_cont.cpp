/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveSpecies generate output for the save species command */
#include "cddefines.h"
#include "continuum.h"
#include "trace.h"
#include "save.h"
#include "radius.h"
#include "generic_state.h"
#include "mole.h"
#include "species.h"
#include "lines.h"
#include "lines_service.h"



t_pseudo_cont pseudoContDef;


enum intenType
{
	UNSET = 0,
	INWARD = 1,
	OUTWARD = 2,
	TOTAL = 3
};

STATIC string getIntenTypeStr( const int ipContType )
{
	DEBUG_ENTRY( "getIntenTypeStr()" );

	switch( ipContType )
	{
		case UNSET:
			return "Not Set";
		case INWARD:
			return "Inward";
		case OUTWARD:
			return "Outward";
		case TOTAL:
			return "Total";
		default:
			TotalInsanity();
	}
}

void getSpecies( const string &speciesLabel, genericState &species )
{
	DEBUG_ENTRY( "getSpecies()" );

	vector<genericState> v = matchGeneric( speciesLabel, false );
	if( v.size() != 1 )
	{
		fprintf( ioQQQ, "Error: Incorrect number of matches (%d) for species '%s'\n",
				int(v.size()), speciesLabel.c_str() );
		cdEXIT( EXIT_FAILURE );
	}
	//	printf( "sp= '%s'\n", v[0].label().c_str() );
	species = v[0];
}



/*==============================================================================*/
/*				   BASE CLASS					*/
/*==============================================================================*/
class band_cont
{
protected:
	string speciesLabel;
	long int nBins;
	genericState species;
	vector<realnum> inten_inward,
			inten_outward;
public:
	band_cont() : nBins(0) {}
	band_cont(const band_cont&) = default;
	band_cont& operator= (const band_cont&) = default;
	band_cont(band_cont&&) = default;
	band_cont& operator= (band_cont&&) = default;
	virtual ~band_cont() = default;
	string label() const { return speciesLabel; }
	long bins() const { return nBins; }
protected:
	bool check_index( const long ibin ) const
	{
		return	(ibin < 0 || ibin >= nBins) ? false : true;
	}
	virtual void check_index_fatal( const long ibin ) const = 0;
	virtual void sumBand( double *sumOutward, double *sumInward ) const = 0;
public:
	void accumulate( bool lgReset, double dVeffAper );
	virtual realnum getWl( const long ibin ) const = 0;
	double getInten( const long ibin, const int ipContType ) const;
};

void band_cont::accumulate( bool lgReset = true, double dVeffAper = 1.0 )
{
	DEBUG_ENTRY( "band_cont::accumulate()" );

	// initialize
	if( lgReset )
	{
		for( long ibin = 0; ibin < nBins; ibin++ )
		{
			inten_inward[ ibin ] = 0.;
			inten_outward[ ibin ] = 0.;
		}
	}

	// integrate
	vector<double> sumOutward( nBins ), sumInward( nBins );

	sumBand( get_ptr(sumOutward), get_ptr(sumInward) );

	for( long ibin = 0; ibin < nBins; ibin++ )
	{
		inten_inward[ ibin ] += realnum(sumInward[ ibin ] * dVeffAper);
		inten_outward[ ibin ] += realnum(sumOutward[ ibin ] * dVeffAper);
	}
}

double band_cont::getInten( const long ibin, const int ipContType ) const
{
	DEBUG_ENTRY( "band_cont::getInten()" );

	check_index_fatal( ibin );

	double inten = 0.;
	switch( ipContType )
	{
		case INWARD:
			inten = inten_inward[ ibin ];
			break;
		case OUTWARD:
			inten = inten_outward[ ibin ];
			break;
		case TOTAL:
			inten = inten_inward[ ibin ] + inten_outward[ ibin ];
			break;
		default:
			fprintf( ioQQQ, "Error: Illegal continuum type: %d\n",
					ipContType );
			cdEXIT( EXIT_FAILURE );
			break;
	}

	return inten;
}


/*==============================================================================*/
/*				PSEUDO CONTINUA					*/
/*==============================================================================*/
class pseudo_cont : public band_cont
{
private:
	double wlLo,
		wlHi;
	double log_wlLo,
		log_step;
	vector<realnum> wl;
public:
	void setup( string &label, double wlo, double whi, long nb );
private:
	void sumBand( double *sumOutward, double *sumInward ) const;
	void check_index_fatal( const long ibin ) const
	{
		if( ! check_index( ibin ) )
		{
			fprintf( ioQQQ, "Error: Pseudo-continuum bin (%ld) "
					"for species '%s' "
					"out of range (0, %ld)\n",
					ibin,
					speciesLabel.c_str(),
					nBins-1
				);
			cdEXIT( EXIT_FAILURE );
		}
	}
public:
	realnum getWl( const long ibin ) const
	{
		check_index_fatal( ibin );
		double mid_wl = 0.5 * ( wl[ ibin ] + wl[ ibin+1 ] );
		return realnum( AnuUnit( realnum( RYDLAM/ mid_wl ) ) );
	}
};

void pseudo_cont::setup( string &label, double wlo, double whi, long nb )
{
	DEBUG_ENTRY( "pseudo_cont::setup()" );

	speciesLabel = label;
	wlLo = wlo;
	wlHi = whi;
	nBins = nb;
	wl.resize( nBins + 1 );
	inten_inward.resize( nBins );
	inten_outward.resize( nBins );

	log_step = log10( wlHi / wlLo ) / nBins;
	log_wlLo = log10( wlLo );
	for( long ibin = 0; ibin < nBins+1; ++ibin )
	{
		// bounds of cells in Angstroms
		wl[ ibin ] = exp10( log_wlLo + ibin*log_step );
		//	fprintf( stdout, "ibin= %ld\t wl= %g\n",
		//			ibin, wl[ ibin ] );
	}

	getSpecies( speciesLabel, species );
	//	printf( "sp= '%s'\n", species.label().c_str() );
}

void pseudo_cont::sumBand( double *sumOutward, double *sumInward ) const
{
	DEBUG_ENTRY( "pseudo_cont::sumBand()" );

	for( long i=0; i<nBins; ++i )
	{
		sumOutward[i] = 0.;
		sumInward[i] = 0.;
	}

	if( density( species ) <= SMALLFLOAT )
		return;

	if( species.sp->lines == NULL )
	{
		fprintf( ioQQQ,
			"WARNING: Species '%s' does not have any data for 'save species continuum'.\n",
			species.label().c_str() );
		return;
	}

	double log_step_inv = 1. / log_step;

	for( TransitionProxy::iterator tr = species.sp->lines->begin();
		tr != species.sp->lines->end(); ++tr )
	{
		if( (*tr).WLAng() < wlLo || (*tr).WLAng() > wlHi )
			continue;

		double bin = (log10( (*tr).WLAng() ) - log_wlLo) * log_step_inv;
		long ibin = long( bin );
		if( ! check_index( ibin ) )
			continue;

		{
			enum { DEBUG_BAND = false };
			if( DEBUG_BAND && (*tr).Emis().xIntensity() > 0. )
			{
				fprintf( ioQQQ,
					"WLAng= %g\t wl(i)= %g\t wl(i+1)= %g\t "
					"bin= %g\t ibin= %ld\t inten= %g\t "
					"frac_inw= %g\n",
					(*tr).WLAng(),
					wl[ibin], wl[ibin+1],
					bin, ibin,
					(*tr).Emis().xIntensity(),
					(*tr).Emis().FracInwd() );
			}
		}

		sumOutward[ ibin ] += (*tr).Emis().xIntensity() *
			MAX2( 0., 1-(*tr).Emis().FracInwd() );
		sumInward[ ibin ] += (*tr).Emis().xIntensity() *
			(*tr).Emis().FracInwd();
	}
}
/*============================================================================*/

static vector< pseudo_cont > PseudoCont;

STATIC void getPseudoIndex( const string &speciesLabel,
				vector<pseudo_cont>::iterator &this_it )
{
	DEBUG_ENTRY( "getPseudoIndex()" );

	this_it = PseudoCont.end();

	for( vector<pseudo_cont>::iterator it = PseudoCont.begin();
		it != PseudoCont.end(); ++it )
	{
		//	fprintf( ioQQQ, "'%s'\t '%s'\n",
		//		speciesLabel.c_str(),
		//		(*it).label().c_str() );
		if( speciesLabel == (*it).label() )
		{
			this_it = it;
			break;
		}
	}
}

STATIC long getAdjPseudoIndex( const string &speciesLabel )
{
	DEBUG_ENTRY( "getAdjPseudoIndex()" );

	long index = -1;
	for( size_t jps = 0; jps < save.setPseudoCont.size(); jps++ )
	{
		//	fprintf( ioQQQ,
		//		"'%s'\t '%s'\n",
		//		speciesLabel.c_str(),
		//		save.setPseudoCont[ jps ].c_str() );
		if( speciesLabel == save.setPseudoCont[ jps ].speciesLabel )
		{
			index = long( jps );
		}
	}

	return index;
}

STATIC void getPseudoWlRange( const string &speciesLabel,
			double &wlLo, double &wlHi, long &nBins )
{
	DEBUG_ENTRY( "getPseudoWlRange()" );

	wlLo = pseudoContDef.wlLo;
	wlHi = pseudoContDef.wlHi;
	nBins = pseudoContDef.nBins;

	long index = getAdjPseudoIndex( speciesLabel );

	if( index >= 0 && index < long(save.setPseudoCont.size()) )
	{
		wlLo = save.setPseudoCont[ index ].wlLo;
		wlHi = save.setPseudoCont[ index ].wlHi;
		nBins = save.setPseudoCont[ index ].nBins;
	}

	//	fprintf( stdout,
	//		"'%s'\t xLamLow = %g\t xLamHi = %g\t nBins = %ld\n",
	//		speciesLabel.c_str(), wlLo, wlHi, nBins );
}

STATIC void PseudoContCreate( long ips )
{
	DEBUG_ENTRY( "PseudoContCreate()" );

	double wlLo,
		wlHi;
	long int nBins;

	getPseudoWlRange( save.contSaveSpeciesLabel[ ips ],
		wlLo, wlHi, nBins );

	PseudoCont[ ips ].setup( save.contSaveSpeciesLabel[ ips ],
		wlLo, wlHi, nBins );
	return;
}

void SpeciesPseudoContCreate()
{
	DEBUG_ENTRY( "PseudoContCreate()" );

	if( PseudoCont.size() != 0 )
		return;

	PseudoCont.resize( save.contSaveSpeciesLabel.size() );

	for( size_t ips = 0; ips < save.contSaveSpeciesLabel.size(); ips++ )
	{
		PseudoContCreate( ips );
	}
}

/*============================================================================*/

void SpeciesPseudoContAccum()
{
	DEBUG_ENTRY( "PseudoContAccum()" );

	if( LineSave.ipass != 1 )
		return;

	for( vector<pseudo_cont>::iterator it = PseudoCont.begin();
		it != PseudoCont.end(); ++it )
	{
		(*it).accumulate( nzone == 1, radius.dVeffAper );
	}
}

/*============================================================================*/

STATIC long resolveSpecType( const char *saveSpec )
{
	DEBUG_ENTRY( "resolveSpecType()" );

	long ipContType = UNSET;
	if( strcmp( saveSpec, "CONi" ) == 0 )
	{
		ipContType = INWARD;
	}
	else if( strcmp( saveSpec, "CONo" ) == 0 )
	{
		ipContType = OUTWARD;
	}
	else if( strcmp( saveSpec, "CONt" ) == 0 )
	{
		ipContType = TOTAL;
	}
	return ipContType;
}

void SaveSpeciesPseudoCont( const long ipPun, const string &speciesLabel )
{
	DEBUG_ENTRY( "SaveSpeciesPseudoCont()" );

	vector<pseudo_cont>::iterator it;
	getPseudoIndex( speciesLabel, it );
	if( it == PseudoCont.end() )
	{
		fprintf( ioQQQ,
			"Error: Species continuum data unmatched for species '%s'\n",
			speciesLabel.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	if( save.punarg[ipPun][0] == 1 )
	{
		long ipContType = resolveSpecType( save.chSaveArgs[ ipPun ] );

		// row format used for grids - one time print of wavelength grid first
		if( save.lgSaveHeader(ipPun) )
		{
			// one time print of continuum header
			fprintf( save.params[ipPun].ipPnunit,
				"#%s ", speciesLabel.c_str() );
			switch( ipContType )
			{
				case INWARD:
					fprintf( save.params[ipPun].ipPnunit, "inward" );
					break;
				case OUTWARD:
					fprintf( save.params[ipPun].ipPnunit, "outward" );
					break;
				case TOTAL:
					fprintf( save.params[ipPun].ipPnunit, "total" );
					break;
				default:
					TotalInsanity();
			}
			fprintf( save.params[ipPun].ipPnunit,
				": Energy\t%s Int[erg cm-2 s-1]\n",
				getIntenTypeStr( ipContType ).c_str() );

			for( long ibin = 0; ibin < (*it).bins(); ibin++ )
			{
				if( ibin > 0 )
					fprintf( save.params[ipPun].ipPnunit, "\t" );
				fprintf( save.params[ipPun].ipPnunit,
					"%.5e",
					(*it).getWl( ibin ) );
			}
			fprintf( save.params[ipPun].ipPnunit, "\n");
			save.SaveHeaderDone(ipPun);
		}

		/* give intensities */
		//	fprintf( save.params[ipPun].ipPnunit, "%s\t",
		//		getIntenTypeStr( ipContType ).c_str() );
		for( long ibin = 0; ibin < (*it).bins(); ibin++ )
		{
			fprintf( save.params[ipPun].ipPnunit, "%e",
				(*it).getInten( ibin, ipContType ));
			if( ibin < (*it).bins()-1 )
			{
				fprintf( save.params[ipPun].ipPnunit, "\t" );
			}
		}
		fprintf( save.params[ipPun].ipPnunit, "\n");
	 }
	 else if( save.punarg[ipPun][0] == 2 )
	 {
		// the default, four columns, wl, total, inward, outward
		// one time print of header
		if( save.lgSaveHeader(ipPun) )
		{
			fprintf( save.params[ipPun].ipPnunit,
				"#%s ", speciesLabel.c_str() );
			fprintf( save.params[ipPun].ipPnunit,
				"wl\ttotal\tinward\toutward\n" );
			save.SaveHeaderDone(ipPun);
		}
		//	printf("ips= %ld\tbins= %ld\n", ips, PseudoCont[ ips ].bins());
		for( long ibin = 0; ibin < (*it).bins(); ibin++ )
		{
			fprintf( save.params[ipPun].ipPnunit, "%.5f",
				(*it).getWl( ibin ) );
			fprintf( save.params[ipPun].ipPnunit, "\t%e",
				(*it).getInten( ibin, TOTAL ) );
			fprintf( save.params[ipPun].ipPnunit, "\t%e",
				(*it).getInten( ibin, INWARD ) );
			fprintf( save.params[ipPun].ipPnunit, "\t%e",
				(*it).getInten( ibin, OUTWARD ) );
			fprintf( save.params[ipPun].ipPnunit, "\n");
		}
	}
}



/*==============================================================================*/
/*					BANDS					*/
/*==============================================================================*/
class bands_file
{
private:
	string fileBands;
	vector<realnum> prt_wl,
			wlLo,
			wlHi;
	long nBands;
public:
	void setup( const string &fname )
	{
		fileBands = fname;
	}
	string bandFilename() const
	{
		return fileBands;
	}
	bool load();
	long get_nBands() const
	{
		return nBands;
	}
	realnum getWl( const long iband ) const
	{
		return double( prt_wl[ iband ] );
	}
	realnum getWlLo( const long iband ) const
	{
		return wlLo[ iband ];
	}
	realnum getWlHi( const long iband ) const
	{
		return wlHi[ iband ];
	}
};

bool bands_file::load()
{
	DEBUG_ENTRY( "bands_file::load()" );

	/* get FeII band data */
	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " BandsCreate opening %s:", fileBands.c_str() );
	}

	FILE *ioDATA = open_data( fileBands, "r" );

	/* now count how many bands are in the file */
	nBands = 0;

	string chLine;

	if( fileBands == "FeII_bands.ini" )
	{ 
		/* first line is a version number and does not count */
		if( !read_whole_line( chLine, ioDATA ) )
		{
			fprintf( ioQQQ, " BandsCreate could not read first line of %s.\n",
				fileBands.c_str() );
			return false;
		}
	}

	while( read_whole_line( chLine, ioDATA ) )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
			++nBands;
	}

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " BandsCreate could not rewind %s.\n",
			fileBands.c_str() );
		return false;
	}

	prt_wl.resize( nBands );
	wlLo.resize( nBands );
	wlHi.resize( nBands );

	/* first line is a version number - now confirm that it is valid */
	if( fileBands == "FeII_bands.ini" )
	{
		if( !read_whole_line( chLine, ioDATA ) )
		{
			fprintf( ioQQQ, " BandsCreate could not read first line of %s.\n",
				fileBands.c_str() );
			return false;
		}
		{
			long i = 1;
			const long int iyr = 9, imo=6 , idy = 11;
			long iyrread, imoread , idyread;
			bool lgEOL;
			iyrread = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			imoread = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			idyread = (long)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);

			if(( iyrread != iyr ) ||
			  (  imoread != imo ) ||
			  (  idyread != idy ) )
			{
				fprintf( ioQQQ, 
					" PROBLEM BandsCreate: the version of %s is not the "
					"current version.\n",
					fileBands.c_str() );
				fprintf( ioQQQ, 
					"         BandsCreate: I expected the magic numbers %li %li %li "
					"but found %li %li %li.\n", 
					iyr, imo , idy ,iyrread, imoread , idyread  );
				return false;
			}
		}
	}

	/* now read in data again, but save it this time */
	long k = 0;
	while( read_whole_line( chLine, ioDATA ) )
	{
		bool lgEOL;

		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
		{
			long i = 1;
			prt_wl[k] = (realnum)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this band column 1.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine.c_str() );
				return false;
			}
			wlLo[k] = (realnum)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this band column 2.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine.c_str() );
				return false;
			}
			wlHi[k] = (realnum)FFmtRead(chLine.c_str(),&i,chLine.length(),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this band column 3.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine.c_str() );
				return false;
			}
			/* fprintf(ioQQQ,
				" band data %f %f %f \n",
				prt_wl[k], wlLo[k],wlHi[k] ); */
			++k;
		}
	}
	/* now validate this incoming data */
	for( long i=0; i < nBands; ++i )
	{
		/* make sure all are positive */
		if( prt_wl[i] <=0. || wlLo[i] <=0. || wlHi[i] <=0. )
		{
			fprintf( ioQQQ, " band %li has none positive entry.\n",i );
			return false;
		}
		/* make sure bands bounds are in correct order, shorter - longer wavelength*/
		if( wlLo[i] >= wlHi[i] )
		{
			fprintf( ioQQQ, " band %li has improper bounds.\n" ,i );
			return false;
		}
	}

	fclose(ioDATA);

	/* return success */
	return true;
}

/*==============================================================================*/
static vector<bands_file> Bands;

STATIC void findBandsFile( const string &filename,
				vector<bands_file>::iterator &this_it )
{
	DEBUG_ENTRY( "findBandsFile()" );

	this_it = Bands.end();
	for( vector<bands_file>::iterator it = Bands.begin();
		it != Bands.end(); ++it )
	{
		if( (*it).bandFilename() == filename )
		{
			this_it = it;
			break;
		}
	}
}

STATIC void addBandsFile( const string &filename,
				vector<bands_file>::iterator &it )
{
	DEBUG_ENTRY( "addBandsFile()" );

	findBandsFile( filename, it );

	if( it == Bands.end() )
	{
		bands_file b_tmp;
		b_tmp.setup( filename );
		b_tmp.load();
		Bands.push_back( b_tmp );
		it = Bands.end();
		--it;	// iterator pointer to last vector element
	}
}



/*==============================================================================*/
/*				SPECIES BAND EMISSION				*/
/*==============================================================================*/
class species_bands : public band_cont
{
private:
	string bandLabel,
		comment;
	vector<bands_file>::iterator bands_it;
public:
	void setup( const string &splab, vector<bands_file>::iterator it )
	{
		speciesLabel = splab;
		getSpecies( splab, species );
		//	printf("species: '%s'\n", species.label().c_str());

		bands_it = it;
		nBins = (*bands_it).get_nBands();
		inten_inward.resize( nBins );
		inten_outward.resize( nBins );

		for( long iband = 0; iband < nBins; ++iband )
		{
			inten_inward[ iband ] = 0.;
			inten_outward[ iband ] = 0.;
		}
	}
private:
	void check_index_fatal( const long iband ) const
	{
		if( ! check_index( iband ) )
		{
			fprintf( ioQQQ, "Error: Band (%ld) "
					"for species '%s' "
					"from file '%s' "
					"out of range (0, %ld)\n",
					iband,
					speciesLabel.c_str(),
					(*bands_it).bandFilename().c_str(),
					nBins-1
				);
			cdEXIT( EXIT_FAILURE );
		}
	}
	void sumBand( double *sumOutward, double *sumInward ) const;
public:
	void insert();
public:
	string bandFilename() const { return (*bands_it).bandFilename(); }
	realnum getWl( const long iband ) const
	{
		check_index_fatal( iband );
		return (*bands_it).getWl( iband );
	}
	realnum getWlLo( const long iband ) const
	{
		check_index_fatal( iband );
		return (*bands_it).getWlLo( iband );
	}
	realnum getWlHi( const long iband ) const
	{
		check_index_fatal( iband );
		return (*bands_it).getWlHi( iband );
	}
};

void species_bands::sumBand( double *sumOutward, double *sumInward ) const
{
	DEBUG_ENTRY( "species_bands::sumBand()" );

	for( long i=0; i<nBins; ++i )
	{
		sumOutward[i] = 0.;
		sumInward[i] = 0.;
	}

	if( density( species ) <= SMALLFLOAT )
		return;

	if( species.sp->lines == NULL )
	{
		fprintf( ioQQQ,
			"WARNING: Species '%s' does not have any data for 'save species bands'.\n",
			species.label().c_str() );
		return;
	}

	for( TransitionProxy::iterator tr = species.sp->lines->begin();
		tr != species.sp->lines->end(); ++tr )
	{
		for( long iband = 0; iband < nBins; ++iband )
		{
			if( (*tr).WLAng() >= getWlLo( iband ) &&
				(*tr).WLAng() < getWlHi( iband ) )
			{
				sumOutward[ iband ] += (*tr).Emis().xIntensity() *
					MAX2( 0., 1-(*tr).Emis().FracInwd() );
				sumInward[ iband ] += (*tr).Emis().xIntensity() *
					(*tr).Emis().FracInwd();
			}
		}
	}
}

void species_bands::insert()
{
	DEBUG_ENTRY( "species_bands::insert()" );

	if( bandLabel.length() == 0 || comment.length() == 0 )
	{
		string spectralLabel;
		chemical_to_spectral( speciesLabel, spectralLabel );
		//	printf("spectralLabel = '%s'\n", spectralLabel.c_str());
		bandLabel = spectralLabel + "b";
		comment = spectralLabel + " emission in bands defined in " +
				(*bands_it).bandFilename();
	}

	for( long iband = 0; iband < nBins; iband++ )
	{
		long ipnt;
		PntForLine( double( getWl( iband ) ), bandLabel.c_str(), &ipnt );
		lindst( inten_inward[ iband ] + inten_outward[ iband ],
			getWl( iband ),
			bandLabel.c_str(),
			ipnt, 't', false,
			(" total " + comment ).c_str() );
		lindst( inten_inward[ iband ],
			getWl( iband ),
			"InwdBnd",
			ipnt, 't', false,
			(" inward " + comment ).c_str() );
	}
}
/*============================================================================*/

static vector<species_bands> SpecBands;

STATIC void getSpecBandsIndex( const string &speciesLabel, const string &fileBands,
				 vector<species_bands>::iterator &this_it )
{
	DEBUG_ENTRY( "getSpecBandsIndex()" );

	this_it = SpecBands.end();

	for( vector<species_bands>::iterator it = SpecBands.begin();
		it != SpecBands.end(); ++it )
	{
		if( speciesLabel == (*it).label() &&
			fileBands == (*it).bandFilename() )
		{
			this_it = it;
			break;
		}
	}
}

/*============================================================================*/

void SpeciesBandsCreate()
{
	DEBUG_ENTRY( "SpeciesBandsCreate()" );

	// Already initialized
	if( SpecBands.size() != 0 )
		return;

	for( vector<save_species_bands>::iterator it = save.specBands.begin();
		it != save.specBands.end(); ++it )
	{
		vector<bands_file>::iterator b_it;
		addBandsFile( (*it).filename, b_it );

		species_bands sb_tmp;
		sb_tmp.setup( (*it).speciesLabel, b_it );
		SpecBands.push_back( sb_tmp );
	}
}

/*============================================================================*/

void SpeciesBandsAccum()
{
	DEBUG_ENTRY( "SpeciesBandsAccum()" );

	for( vector<species_bands>::iterator it = SpecBands.begin();
		it != SpecBands.end(); ++it )
	{
		(*it).accumulate();
		(*it).insert();
	}
}

/*============================================================================*/


void SaveSpeciesBands( const long ipPun, const string &speciesLabel,
			const string &fileBands )
{
	DEBUG_ENTRY( "SaveSpeciesBands()" );

	vector<species_bands>::iterator it;
	getSpecBandsIndex( speciesLabel, fileBands, it );
	if( it == SpecBands.end() )
	{
		fprintf( ioQQQ,
			"Error: Species band data unmatched for species "
			"'%s' and bands from file '%s'\n",
			speciesLabel.c_str(), fileBands.c_str() );
		cdEXIT( EXIT_FAILURE );
	}

	if( save.lgSaveHeader(ipPun) )
	{
		// one time print of header
		fprintf( save.params[ipPun].ipPnunit, "#Wl(A)\t Intensity: Total\t Inward\t Outward\n" );
		save.SaveHeaderDone(ipPun);
	}

	for( long iband = 0; iband < (*it).bins(); iband++ )
	{
		fprintf( save.params[ipPun].ipPnunit, "%g",
			(*it).getWl( iband ) );
		fprintf( save.params[ipPun].ipPnunit, "\t%e",
			(*it).getInten( iband, TOTAL ) );
		fprintf( save.params[ipPun].ipPnunit, "\t%e",
			(*it).getInten( iband, INWARD ) );
		fprintf( save.params[ipPun].ipPnunit, "\t%e",
			(*it).getInten( iband, OUTWARD ) );
		fprintf( save.params[ipPun].ipPnunit, "\n" );
	}
}
