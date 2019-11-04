/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atom_levelN compute an arbitrary N level atom */
#include "cddefines.h"
#include "atoms.h"
#include "physconst.h"
#include "thermal.h"
#include "conv.h"
#include "phycon.h"
#include "dense.h"
#include "trace.h"
#include "newton_step.h"
#include "dynamics.h"
#include "vectorize.h"
#include "prt.h"

// #define EIGEN
#ifdef EIGEN
extern "C" {
	void dgeev_(const char*, const char*,int*,double*,int*,double*,
					double*,
					double*,int*,double*,int*,double*,int*,int*,int,int);
}
#endif

// Set collision rate from collision strength
void setCollRate::operator()(
	long nlev,
	double TeInverse,
	double **col_str,
	const double ex[],	
	const double g[], 
	double **CollRate
)
{
	resize(nlev);
	for( long ihi=1; ihi < nlev; ++ihi )
	{
		double fac = dsexp((ex[ihi]-ex[ihi-1])*TeInverse);
		for( long ilo=0; ilo < ihi-1; ++ilo )
		{
			excit[ihi][ilo] = fac*excit[ihi-1][ilo];
		}
		excit[ihi][ihi-1] = fac;
	}
	
	if( trace.lgTrace && trace.lgTrLevN )
	{
		fprintf( ioQQQ, " coll str\n" );
		fprintf( ioQQQ, "  hi  lo" );
		for( long ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( long ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( long ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (col_str)[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " excit, te=%10.2e\n", phycon.te );
		fprintf( ioQQQ, "  hi  lo" );
		
		for( long ilo=0; ilo < (nlev-1); ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );
		
		for( long ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( long ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", excit[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}
	}
	
	for( long ihi=1; ihi < nlev; ihi++ )
	{
		for( long ilo=0; ilo < ihi; ilo++ )
		{
			/* this should be a collision strength */
			ASSERT( (col_str)[ihi][ilo]>= 0. );
			/* this is deexcitation rate */
			(CollRate)[ihi][ilo] = (col_str)[ihi][ilo]/g[ihi]*dense.cdsqte;
			/* this is excitation rate */
			(CollRate)[ilo][ihi] = (CollRate)[ihi][ilo]*g[ihi]/g[ilo]*
				excit[ihi][ilo];
		}
	}	
}

void Atom_LevelN::operator()(
	long int nLevelCalled,
	/* ABUND is total abundance of species, used for nth equation
	 * if balance equations are homogeneous */
	realnum abund, 
	/* G(nlev) is stat weight of levels */
	const vector<double>& g, 
	/* EX(nlev) is excitation potential of levels, deg K or wavenumbers
	 * 0 for lowest level, all are energy rel to ground NOT d(ENER) */
	const vector<double>& ex, 
	/* this is 'K' for ex[] as Kelvin deg, is 'w' for wavenumbers */
	char chExUnits,
	/* populations [cm-3] of each level as deduced here */
	vector<double>& pops, 
	/* departure coefficient, derived below */
	vector<double>& depart,
	/* net transition rate, A * esc prob, s-1 */
	multi_arr<double,2>& AulEscp, 
	/* AulDest[ihi][ilo] is destruction rate, trans from ihi to ilo, A * dest prob,
	 * asserts confirm that [ihi][ilo] is zero */
	multi_arr<double,2>& AulDest, 
	/* AulPump[ilo][ihi] is pumping rate from lower to upper level (s^-1), (hi,lo) must be zero  */
	multi_arr<double,2>& AulPump, 
	/* collision rates (s^-1), evaluated here and returned for cooling by calling function,
	 * unless following flag is true.  If true then calling function has already filled
	 * in these rates.  CollRate[i][j] is rate from i to j */
	const multi_arr<double,2>& CollRate,
	/* this is an additional creation rate from continuum, normally zero, units cm-3 s-1 */
	const vector<double>& source,
	/* this is an additional destruction rate to continuum, normally zero, units s-1 */
	const vector<double>& sink,
	/* total cooling and its derivative, set here but nothing done with it*/
	double *cooltl, 
	double *coolder, 
	/* string used to identify species in case of error */
	const char *chLabel, 
	/* flag to print matrices input to solvers */
	const bool lgPrtMatrix,
	/* nNegPop flag indicating what we have done
	 * positive if negative populations occurred
	 * zero if normal calculation done
	 * negative if too cold, matrix not solved, since highest level had little excitation
	 * (for some atoms other routine will be called in this case) */
	int *nNegPop,
	/* true if populations are zero, either due to zero abundance of very low temperature */
	bool *lgZeroPop,
	/* option to print debug information */
	bool lgDeBug,
	/* option to do atom in LTE, can be omitted, default value false */
	bool lgLTE,
	/* array that will be set to the cooling per line, can be omitted, default value NULL */
	multi_arr<double,2> *Cool,
	/* array that will be set to the cooling derivative per line, can be omitted, default value NULL */
	multi_arr<double,2> *dCooldT,
	/* total excitation rate out of ground state */
	double *grnd_excit)
{
	long int level, 
	  ihi, 
	  ilo, 
	  j; 

	int32 ner;

	double cool1,
	  TeInverse,
	  TeConvFac;

	DEBUG_ENTRY( "atom_levelN()" );

	*nNegPop = -1;
	if (grnd_excit != NULL) *grnd_excit = 0.0;

	/* >>chng 05 dec 14, units of ex[] can be Kelvin (old default) or wavenumbers */
	if( chExUnits=='K' )
	{
		/* ex[] is in temperature units - this will multiply ex[] to
		 * obtain Boltzmann factor */
		TeInverse = 1./phycon.te;
		/* this multiplies ex[] to obtain energy difference between levels */
		TeConvFac = 1.;
	}
	else if( chExUnits=='w' )
	{
		/* ex[] is in wavenumber units */
		TeInverse = 1./phycon.te_wn;
		TeConvFac = T1CM;
	}
	else
		TotalInsanity();

	avx_ptr<double> arg(1,nLevelCalled), excit_gnd(1,nLevelCalled);
	for( ihi=1; ihi < nLevelCalled; ++ihi )
		arg[ihi] = -(ex[ihi]-ex[0])*TeInverse;
	vexp( arg.ptr0(), excit_gnd.ptr0(), 1, nLevelCalled );

	long int nlev = nLevelCalled;
	// decrement number of levels until we have positive excitation rate,
	while( nlev > 1 && excit_gnd[nlev-1]+AulPump[0][nlev-1] < SMALLFLOAT &&
	       source[nlev-1] == 0. )
	{
		pops[nlev-1] = 0.;
		depart[nlev-1] = 0.;
		--nlev;
	}

	/* exit if zero abundance or all population in ground */
	ASSERT( abund>= 0. );
	if( abund == 0. || nlev==1 )
	{
		*cooltl = 0.;
		*coolder = 0.;
		/* says calc was ok */
		*nNegPop = 0;
		*lgZeroPop = true;

		pops[0] = abund;
		depart[0] = 1.;
		for( level=1; level < nlev; level++ )
		{
			pops[level] = 0.;
			depart[level] = 0.;
		}
		return;
	}

#	ifndef NDEBUG
	/* excitation temperature of lowest level must be zero */
	ASSERT( ex[0] == 0. );

	for( ihi=0; ihi < nlev; ihi++ )
	{
		for( ilo=0; ilo < nlev; ilo++ )
		{
			/* AulDest[ilo][ihi] - so that spontaneous transitions only proceed from high energy to low
			 * AulEscp[ilo][ihi] - so that spontaneous transitions only proceed from high energy to low */
			ASSERT( ihi == ilo || (AulDest)[ihi][ilo] >= 0. );
			ASSERT( ihi == ilo || (AulEscp)[ihi][ilo] >= 0 );
		}
	}

	for( ihi=1; ihi < nlev; ihi++ )
	{
		for( ilo=0; ilo < ihi; ilo++ )
		{
			ASSERT( (AulPump)[ihi][ilo] >= 0. );
		}
	}
#	endif

	if( lgDeBug || (trace.lgTrace && trace.lgTrLevN) )
	{
		fprintf( ioQQQ, " atom_levelN trace printout for atom=%s with tot abund %e \n", chLabel, abund);
		fprintf( ioQQQ, " AulDest\n" );

		fprintf( ioQQQ, "  hi  lo" );

		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (AulDest)[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " A*esc\n" );
		fprintf( ioQQQ, "  hi  lo" );
		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (AulEscp)[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " AulPump\n" );

		fprintf( ioQQQ, "  hi  lo" );
		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (AulPump)[ilo][ihi] );
			}
			fprintf( ioQQQ, "\n" );
		}

		fprintf( ioQQQ, " coll rate\n" );
		fprintf( ioQQQ, "  hi  lo" );
		for( ilo=0; ilo < nlev-1; ilo++ )
		{
			fprintf( ioQQQ, "%4ld      ", ilo );
		}
		fprintf( ioQQQ, "      \n" );

		for( ihi=1; ihi < nlev; ihi++ )
		{
			fprintf( ioQQQ, "%3ld", ihi );
			for( ilo=0; ilo < ihi; ilo++ )
			{
				fprintf( ioQQQ, "%10.2e", (CollRate)[ihi][ilo] );
			}
			fprintf( ioQQQ, "\n" );
		}
	}

	/* we will predict populations */
	*lgZeroPop = false;

	if( !lgLTE )
	{
		bvec.resize(nlev);
		amat.alloc(nlev,nlev);

		/* rate equations equal zero */
		amat.zero();

		/* eqns for destruction of level
		 * AulEscp[iho][ilo] are A*esc p, CollRate is coll excit in direction
		 * AulPump[low][high] is excitation rate from low to high */
		for( level=0; level < nlev; level++ )
		{
			/* leaving level to below */
			for( ilo=0; ilo < level; ilo++ )
			{
				double rate = (CollRate)[level][ilo] + (AulEscp)[level][ilo] +
					(AulDest)[level][ilo] + (AulPump)[level][ilo];
				amat[level][level] += rate;
				amat[level][ilo] = -rate;
			}

			/* leaving level to above */
			for( ihi=level + 1; ihi < nlev; ihi++ )
			{
				double rate = (CollRate)[level][ihi] + (AulPump)[level][ihi];
				amat[level][level] += rate;
				amat[level][ihi] = -rate;
			}
		}
		if (grnd_excit != NULL) 	*grnd_excit = amat[0][0];

		if( (dynamics.lgAdvection || dynamics.lgTimeDependentStatic) &&
			iteration > dynamics.n_initial_relax )
		{
			double slowrate = -1.0;
			long slowlevel = -1;
			// level == 0 is intentionally excluded... 
			for (level=1; level < nlev; ++level)
			{
				double rate = dynamics.timestep*amat[level][level];
				if (rate < slowrate || slowrate < 0.0)
				{
					slowrate = rate;
					slowlevel = level;
				}
			}
			if ( slowrate < 3.0 )
			{
				static map<string,double> failures;
				ostringstream oss;
				oss << chLabel << "[level=" << slowlevel << ']';
				string str=oss.str();
				if (failures.find(str) == failures.end())
				{
					failures[str] = slowrate;
					// >>chng 16 feb 18: demoted the following message from a PROBLEM
					//   to a CAUTION to prevent a flood of messages in minor.txt...
					fprintf(ioQQQ," CAUTION in atom_levelN: " 
							  "non-equilibrium appears important for %s, "
							  "rate * timestep is %g\n",
							  str.c_str(), slowrate);
				}
				else 
				{
					if ( slowrate < failures[str])
						failures[str] = slowrate;
				}
			}
		}

		double maxdiag = 0.;
		for( level=0; level < nlev; level++ )
		{
			// source terms from elsewhere
			bvec[level] = source[level];
			if( amat[level][level] > maxdiag )
				maxdiag = amat[level][level]; 
			amat[level][level] += sink[level];
		}

		// If there are no sources or sinks, then one of the rows of the
		// linear system is linearly dependent on the others so there is
		// no matrix inverse.  If the sources are sufficiently small,
		// the answer will be numerically ill-conditioned.
		// 
		// For strictly zero sources, this may be remedied by applying a
		// conservation constraint instead of one of the redundant rows.
		// To handle the broader case, we here add this conservation
		// constraint scaled to switch in smoothly as the condition
		// number of the matrix becomes poorer (and hence the accuracy
		// of the linear system becomes poor).
		//
		// The condition number estimate here is quick but very crude.
		//
		// Need to specify smallest condition number before we decide
		// that conservation will get a better answer than the raw
		// linear solver.  Numerical Recipes (3rd edition, S2.6.2)
		// suggests that ~1e-14 might be appropriate for the solution of
		// matrices in double precision.

		const double CONDITION_SCALE = 1e-10;
		double conserve_scale = maxdiag*CONDITION_SCALE;
		ASSERT( conserve_scale > 0.0 );

		for( level=0; level<nlev; ++level )
		{
			amat[level][0] += conserve_scale;
		}
		bvec[0] += abund*conserve_scale;

		if( lgDeBug )
		{
			fprintf(ioQQQ," amat matrix follows:\n");
			for( level=0; level < nlev; level++ )
			{
				for( j=0; j < nlev; j++ )
				{
					fprintf(ioQQQ," %.4e" , amat[level][j]);
				}
				fprintf(ioQQQ,"\n");
			}
			fprintf(ioQQQ," Vector follows:\n");
			for( j=0; j < nlev; j++ )
			{
				fprintf(ioQQQ," %.4e" , bvec[j] );
			}
			fprintf(ioQQQ,"\n");
		}

#ifdef EIGEN
		if (lgDeBug)
		{
			multi_arr<double,2,C_TYPE> mcol(nlev,nlev),mrad(nlev,nlev),
				vl(nlev,nlev),vr(nlev,nlev);
			mcol.zero();
			mrad.zero();

			for( level=0; level < nlev; level++ )
			{
				/* leaving level to below */
				for( ilo=0; ilo < level; ilo++ )
				{
					double rate = (*CollRate)[level][ilo];
					mcol[level][level] += rate;
					mcol[level][ilo] = -rate;
					rate = (*AulEscp)[level][ilo] +
						(*AulDest)[level][ilo] + (*AulPump)[level][ilo];
					mrad[level][level] += rate;
					mrad[level][ilo] = -rate;
				}
				/* leaving level to above */
				for( ihi=level + 1; ihi < nlev; ihi++ )
				{
					double rate = (*CollRate)[level][ihi];
					mcol[level][level] += rate;
					mcol[level][ihi] = -rate;
					rate = (*AulPump)[level][ihi];
					mrad[level][level] += rate;
					mrad[level][ihi] = -rate;
				}
			}
			multi_arr<double,2,C_TYPE> mtst(nlev,nlev);
			for( ilo=0; ilo < nlev; ilo++ )
			{
				for( ihi=0; ihi < nlev; ++ihi)					
				{					
					mtst[ilo][ihi] = mcol[ilo][ihi];
				}
			}
			const char job[2] = "V";
			int lwork = 4*nlev,info,inlev=nlev;
			vector<double> wr(nlev),wi(nlev),work(lwork);
			dgeev_(job,job,&inlev,get_ptr(mtst.vals()),&inlev,get_ptr(wr),
					 get_ptr(wi),get_ptr(vl.vals()),&inlev,
					 get_ptr(vr.vals()),&inlev,get_ptr(work),&lwork,&info,1,1);
			fprintf(ioQQQ,"info %d\nW=",info);
			for( level=0; level < nlev; level++ )
			{
				fprintf(ioQQQ,"%ld %15.8g+%15.8gi\n",level,wr[level],wi[level]);
			}

			// Reconstruct collision matrix
			// mcol = u w v^t
			// mcol[ihi][ilo] = sum(level,vl[level][ilo]*w[level]*vr[level][ihi])
			// ( u^t mrad v + w ) v^t b = 0

			// Normalize left eigen-vectors so L & R are inverses
			for( level=0; level < nlev; level++ )
			{
				double tot=0.0;
				for( ilo=0; ilo < nlev; ilo++ )
				{
					tot += vl[level][ilo]*vr[level][ilo];
				}
				tot = 1.0/tot;
				for( ilo=0; ilo < nlev; ilo++ )
				{
					vl[level][ilo] *= tot;
				}
			}

			mtst.zero();
			for( ilo=0; ilo < nlev; ilo++ )
			{
				for( ihi=0; ihi < nlev; ihi++ )
				{
					for( level=0; level < nlev; level++ )
					{
						mtst[ilo][ihi] += mrad[level][ihi] * vr[ilo][level]; // or mcol to test...
					}
				}
			}
			multi_arr<double,2,C_TYPE> mtst1(nlev,nlev);
			mtst1.zero();
			for( ilo=0; ilo < nlev; ilo++ )
			{
				for( ihi=0; ihi < nlev; ihi++ )
				{
					for( level=0; level < nlev; level++ )
					{
						mtst1[ilo][ihi] += vl[ihi][level]*mtst[ilo][level];
					}
				}
			}
			for( level=0; level < nlev; level++ )
			{
				fprintf(ioQQQ,"%ld %15.8g %15.8g\n",level,wr[level],mtst1[level][level]);
			}
		}
#endif

		if( lgPrtMatrix )
		{
			prt.matrix.prtRates( nlev, amat, bvec );
		}

		ner = solve_system(amat.vals(), bvec, nlev, NULL);

		if( ner != 0 )
		{
			fprintf( ioQQQ, " PROBLEM atom_levelN: dgetrs finds singular or ill-conditioned matrix for %s\n",
				chLabel);
			cdEXIT(EXIT_FAILURE);
		}

		/* set populations */
		for( level=0; level < nlev; level++ )
		{
			/* save bvec into populations */
			pops[level] = bvec[level];
		}

		/* now find total energy exchange rate, normally net cooling and its 
		 * temperature derivative */
		*cooltl = 0.;
		*coolder = 0.;
		for( ihi=1; ihi < nlev; ihi++ )
		{
			for( ilo=0; ilo < ihi; ilo++ )
			{
				/* this is now net cooling rate [erg cm-3 s-1] */
				cool1 = (pops[ilo]*(CollRate)[ilo][ihi] - pops[ihi]*(CollRate)[ihi][ilo])*
					(ex[ihi] - ex[ilo])*BOLTZMANN*TeConvFac;
				*cooltl += cool1;
				if( Cool != NULL )
					(*Cool)[ihi][ilo] = cool1;
				/* derivative wrt temperature - use Boltzmann factor relative to ground */
				/* >>chng 03 aug 28, use real cool1 */
				double dcooldT1 = cool1*( (ex[ihi] - ex[0])*thermal.tsq1 - thermal.halfte );
				*coolder += dcooldT1;
				if( dCooldT != NULL )
					(*dCooldT)[ihi][ilo] = dcooldT1;
			}
		}

		/* fill in departure coefficients */
		if( pops[0] > SMALLFLOAT && excit_gnd[nlev-1] > SMALLFLOAT )
		{
			/* >>chng 00 aug 10, loop had been from 1 and 0 was set to total abundance */
			depart[0] = 1.;
			for( ihi=1; ihi < nlev; ihi++ )
			{
				/* these are off by one - lowest index is zero */
				depart[ihi] = (pops[ihi]/pops[0])*(g[0]/g[ihi])/excit_gnd[ihi];
			}
		}

		else
		{
			/* >>chng 00 aug 10, loop had been from 1 and 0 was set to total abundance */
			for( ihi=0; ihi < nlev; ihi++ )
			{
				/* Boltzmann factor or abundance too small, set departure coefficient to zero */
				depart[ihi] = 0.;
			}
			depart[0] = 1.;
		}
	}
	else
	{
		/* this branch calculates LTE level populations */

		/* get the value of the partition function and the derivative */
		valarray<double> help1(nlev), help2(nlev);
		double partfun = help1[0] = g[0];
		double dpfdT = help2[0] = 0.; /* help2 is d(help1)/dT */
		for( level=1; level < nlev; level++ )
		{
			help1[level] = g[level]*excit_gnd[level];
			partfun += help1[level];
			help2[level] = help1[level]*(ex[level]-ex[0])*TeInverse/phycon.te;
			dpfdT += help2[level];
		}

		/* calculate the level populations and departure coefficients,
		 * as well as the derivative of the level pops */
		valarray<double> dndT(nlev);
		for( level=0; level < nlev; level++ )
		{
			pops[level] = abund*help1[level]/partfun;
			dndT[level] = abund*(help2[level]*partfun - dpfdT*help1[level])/pow2(partfun);
			depart[level] = 1.;
		}

		/* now find the net cooling from the atom */
		*cooltl = 0.;
		*coolder = 0.;
		for( ihi=1; ihi < nlev; ihi++ )
		{
			for( ilo=0; ilo < ihi; ilo++ )
			{
				double deltaE = (ex[ihi] - ex[ilo])*BOLTZMANN*TeConvFac;
				double cool1 = (pops[ihi]*((AulEscp)[ihi][ilo] + (AulPump)[ihi][ilo]) -
						pops[ilo]*(AulPump)[ilo][ihi])*deltaE;
				*cooltl += cool1;
				if( Cool != NULL )
					(*Cool)[ihi][ilo] = cool1;
				double dcooldT1 = (dndT[ihi]*((AulEscp)[ihi][ilo] + (AulPump)[ihi][ilo]) -
						   dndT[ilo]*(AulPump)[ilo][ihi])*deltaE;
				*coolder += dcooldT1;
				if( dCooldT != NULL )
					(*dCooldT)[ihi][ilo] = dcooldT1;
			}
		}
	}

	/* sanity check for valid solution - non negative populations */
	*nNegPop = 0;
	/* the limit we allow the fractional population to go below zero before announcing failure. */
	const double poplimit = 1.0e-10;
	for( level=0; level < nlev; level++ )
	{
		if( pops[level] < 0. )
		{
			if( fabs(pops[level]/abund) > poplimit )
			{
				//nNegPop >= 1 leads to a failure
				*nNegPop = *nNegPop + 1;
			}
			else
			{
				pops[level] = SMALLFLOAT;
				//fprintf( ioQQQ, "\n PROBLEM Small negative populations were found in atom = %s . "
				//		"The problem was ignored and the negative populations were set to SMALLFLOAT",chLabel );
			}
		}
	}

	if( *nNegPop!=0 )
	{
		ASSERT( *nNegPop>=1 );
		// *nNegPop 0 is valid solution, nNegPop> 1 negative populations found
		fprintf( ioQQQ, "\n PROBLEM atom_levelN found negative population, nNegPop=%i, atom=%s lgSearch=%c\n",
			*nNegPop , chLabel , TorF( conv.lgSearch ) );

		for( level=0; level < nlev; level++ )
		{
			fprintf( ioQQQ, "%10.2e", pops[level] );
		}

		fprintf( ioQQQ, "\n" );
		for( level=0; level < nlev; level++ )
		{
			pops[level] = (double)MAX2(0.,pops[level]);
		}
	}

	if(  lgDeBug || (trace.lgTrace && trace.lgTrLevN) )
	{
		fprintf( ioQQQ, "\n atom_leveln absolute population   " );
		for( level=0; level < nlev; level++ )
		{
			fprintf( ioQQQ, " %10.2e", pops[level] );
		}
		fprintf( ioQQQ, "\n" );

		fprintf( ioQQQ, " departure coefficients" );
		for( level=0; level < nlev; level++ )
		{
			fprintf( ioQQQ, " %10.2e", depart[level] );
		}
		fprintf( ioQQQ, "\n\n" );
	}

#	ifndef NDEBUG
	/* these were reset to non zero values by the solver, but we will
	 * assert that they are zero (for safety) when routine reenters so must
	 * set to zero here, but only if asserts are in place */
	for( ilo=0; ilo < nlev; ilo++ )
	{
		for( ihi=ilo+1; ihi < nlev; ihi++ )
		{
			/* zero ots destruction rate */
			AulDest[ilo][ihi] = 0.;
		}
	}
	for( ihi=1; ihi < nlev; ihi++ )
	{
		for( ilo=0; ilo < ihi; ilo++ )
		{
			/* both AulDest and AulPump (low, hi) are not used, should be zero */
			AulPump[ihi][ilo] = 0.;
		}
	}
	for( ilo=0; ilo < nlev; ilo++ )
	{
		for( ihi=ilo+1; ihi < nlev; ihi++ )
		{
			AulEscp[ilo][ihi] = 0;
		}
	}
#	endif

	// -1 return had meant too cool and no populations determined, no longer used
	// since we now decrement until populations can be determined
	ASSERT( *nNegPop>=0 );

	return;
}
