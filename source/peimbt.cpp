/* This file is part of Cloudy and is copyright (C)1978-2019 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "peimbt.h"
t_peimbt peimbt;

void t_peimbt::zero()
{
	DEBUG_ENTRY( "t_peimbt::zero()" );
	tsqden = 1e7;
}
