/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "abund.h"
t_abund abund;

void t_abund::zero()
{
	DEBUG_ENTRY( "t_abund::zero()" );
	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		/* these are depletion scale factors */
		depset[nelem] = 1.;
		/*begin sanity check */
		if( depset[nelem] == 0. )
		{
			fprintf( ioQQQ, " ZERO finds insane abundance or depletion.\n" );
			fprintf( ioQQQ, " atomic number=%6ld abundance=%10.2e depletion=%10.2e\n", 
			  nelem, solar[nelem], depset[nelem] );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
		/*end sanity check */
	}

	/* typical ISM depletion factors, subjective mean of Cowie and Songaila
	 * and Jenkins 
	 * */
	Depletion[0] = 1.;
	Depletion[1] = 1.;
	Depletion[2] = .16f;
	Depletion[3] = .6f;
	Depletion[4] = .13f;
	Depletion[5] = 0.4f;
	Depletion[6] = 1.0f;
	Depletion[7] = 0.6f;
	Depletion[8] = .3f;
	Depletion[9] = 1.f;
	Depletion[10] = 0.2f;
	Depletion[11] = 0.2f;
	Depletion[12] = 0.01f;
	Depletion[13] = 0.03f;
	Depletion[14] = .25f;
	Depletion[15] = 1.0f;
	Depletion[16] = 0.4f;
	Depletion[17] = 1.0f;
	Depletion[18] = .3f;
	Depletion[19] = 1e-4f;
	Depletion[20] = 5e-3f;
	Depletion[21] = 8e-3f;
	Depletion[22] = 6e-3f;
	Depletion[23] = 6e-3f;
	Depletion[24] = 5e-2f;
	Depletion[25] = 0.01f;
	Depletion[26] = 0.01f;
	Depletion[27] = 0.01f;
	Depletion[28] = .1f;
	Depletion[29] = .25f;

	lgDepln = false;
	ScaleMetals = 1.;
	
}
