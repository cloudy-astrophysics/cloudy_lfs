/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseMetal parse parameters on metal command */
#include "cddefines.h"
#include "input.h"
#include "optimize.h"
#include "grainvar.h"
#include "called.h"
#include "abund.h"
#include "parser.h"

void ParseMetal(Parser &p)
{
	DEBUG_ENTRY( "ParseMetal()" );

	/* parse the metals command */
	abund.lgAbnSolar = false;	
	if( p.nMatch("DEPL") )
	{
		/* keyword depletion is present
		 * deplete by set of scale factors */
		abund.lgDepln = true;
		for( long int i=0; i < LIMELM; i++ )
		{
			abund.depset[i] = abund.Depletion[i];
		}
	}
	else
	{
		/* metal depletion factor, if negative then it is the log */
		abund.ScaleMetals = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("the depletion factor");

		/* sort out whether log */
		bool lgLogOn = false;
		if( p.nMatch(" LOG") )
		{
			lgLogOn = true;
		}
		else if( p.nMatch("LINE") )
		{
			lgLogOn = false;
		}

		double dmlog;
		if( abund.ScaleMetals <= 0. || lgLogOn )
		{
			dmlog = abund.ScaleMetals;
			abund.ScaleMetals = exp10(abund.ScaleMetals);
		}
		else
		{
			dmlog = log10(abund.ScaleMetals);
		}

		/* option to vary grain abundance as well */
		bool lgGrains;
		if( p.nMatch("GRAI") )
		{
			lgGrains = true;
			gv.GrainMetal = abund.ScaleMetals;
		}
		else
		{
			lgGrains = false;
			gv.GrainMetal = 1.;
		}

		/* vary option */
		if( optimize.lgVarOn )
		{
			strcpy( optimize.chVarFmt[optimize.nparm], "METALS= %f LOG" );
			if( lgGrains )
				strcat( optimize.chVarFmt[optimize.nparm], " GRAINS" );
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			optimize.vparm[0][optimize.nparm] = (realnum)dmlog;
			optimize.vincr[optimize.nparm] = 0.5;
			optimize.nvarxt[optimize.nparm] = 1;
			++optimize.nparm;
		}
	}
	return;
}
