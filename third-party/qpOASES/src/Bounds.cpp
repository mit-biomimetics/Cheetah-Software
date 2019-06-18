/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file src/Bounds.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the Bounds class designed to manage working sets of
 *	bounds within a QProblem.
 */


#include <qpOASES/Bounds.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	B o u n d s
 */
Bounds::Bounds( ) : SubjectTo( )
{
}


/*
 *	B o u n d s
 */
Bounds::Bounds( int_t _n ) : SubjectTo( _n )
{
	init( _n );
}


/*
 *	B o u n d s
 */
Bounds::Bounds( const Bounds& rhs ) : SubjectTo( rhs )
{
	copy( rhs );
}


/*
 *	~ B o u n d s
 */
Bounds::~Bounds( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
Bounds& Bounds::operator=( const Bounds& rhs )
{
	if ( this != &rhs )
	{
		clear( );
		SubjectTo::operator=( rhs );
		copy( rhs );
	}

	return *this;
}



/*
 *	i n i t
 */
returnValue Bounds::init(	int_t _n
							)
{
	if ( _n < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	clear( );

	if ( _n >= 0 )
	{
		freee.init( _n );
		fixed.init( _n );
	}

	return SubjectTo::init( _n );
}



/*
 *	s e t u p B o u n d
 */
returnValue Bounds::setupBound(	int_t number, SubjectToStatus _status
								)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Add bound index to respective index list. */
	switch ( _status )
	{
		case ST_INACTIVE:
			if ( this->addIndex( this->getFree( ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		case ST_LOWER:
			if ( this->addIndex( this->getFixed( ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		case ST_UPPER:
			if ( this->addIndex( this->getFixed( ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
			break;

		default:
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A l l F r e e
 */
returnValue Bounds::setupAllFree( )
{
	return setupAll( ST_INACTIVE );
}


/*
 *	s e t u p A l l L o w e r
 */
returnValue Bounds::setupAllLower( )
{
	return setupAll( ST_LOWER );
}


/*
 *	s e t u p A l l U p p e r
 */
returnValue Bounds::setupAllUpper( )
{
	return setupAll( ST_UPPER );
}


/*
 *	m o v e F i x e d T o F r e e
 */
returnValue Bounds::moveFixedToFree( int_t number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of fixed variables to that of free ones. */
	if ( this->removeIndex( this->getFixed( ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( this->addIndex( this->getFree( ),number,ST_INACTIVE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	m o v e F r e e T o F i x e d
 */
returnValue Bounds::moveFreeToFixed(	int_t number, SubjectToStatus _status
										)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of free variables to that of fixed ones. */
	if ( this->removeIndex( this->getFree( ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( this->addIndex( this->getFixed( ),number,_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	f l i p F i x e d
 */
returnValue Bounds::flipFixed( int_t number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( status != 0 )
		switch (status[number])
		{
			case ST_LOWER: status[number] = ST_UPPER; break;
			case ST_UPPER: status[number] = ST_LOWER; break;
			default: return THROWERROR( RET_MOVING_BOUND_FAILED );
		}

	return SUCCESSFUL_RETURN;
}


/*
 *	s w a p F r e e
 */
returnValue Bounds::swapFree(	int_t number1, int_t number2
								)
{
	/* consistency check */
	if ( ( number1 < 0 ) || ( number1 >= n ) || ( number2 < 0 ) || ( number2 >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Swap index within indexlist of free variables. */
	return this->swapIndex( this->getFree( ),number1,number2 );
}


/*
 *	s h i f t
 */
returnValue Bounds::shift(	int_t offset )
{
	int_t i;

	/* consistency check */
	if ( ( offset == 0 ) || ( n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > n/2 ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( ( n % offset ) != 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );


	/* 1) Shift types and status. */
	for( i=0; i<n-offset; ++i )
	{
		setType( i,getType( i+offset ) );
		setStatus( i,getStatus( i+offset ) );
	}

	/* 2) Construct shifted index lists of free and fixed variables. */
	Indexlist shiftedFreee( n );
	Indexlist shiftedFixed( n );

	for( i=0; i<n; ++i )
	{
		switch ( getStatus( i ) )
		{
			case ST_INACTIVE:
				if ( shiftedFreee.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_LOWER:
				if ( shiftedFixed.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_UPPER:
				if ( shiftedFixed.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			default:
				return THROWERROR( RET_SHIFTING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	freee = shiftedFreee;
	fixed = shiftedFixed;

	return SUCCESSFUL_RETURN;
}


/*
 *	r o t a t e
 */
returnValue Bounds::rotate( int_t offset )
{
	int_t i;

	/* consistency check */
	if ( ( offset == 0 ) || ( offset == n ) || ( n <= 1 ) )
		return SUCCESSFUL_RETURN;

	if ( ( offset < 0 ) || ( offset > n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );


	/* 1) Rotate types and status. */
	SubjectToType*   typeTmp   = new SubjectToType[offset];
	SubjectToStatus* statusTmp = new SubjectToStatus[offset];

	for( i=0; i<offset; ++i )
	{
		typeTmp[i] = getType( i );
		statusTmp[i] = getStatus( i );
	}

	for( i=0; i<n-offset; ++i )
	{
		setType( i,getType( i+offset ) );
		setStatus( i,getStatus( i+offset ) );
	}

	for( i=n-offset; i<n; ++i )
	{
		setType( i,typeTmp[i-n+offset] );
		setStatus( i,statusTmp[i-n+offset] );
	}

	delete[] statusTmp; delete[] typeTmp;

	/* 2) Construct shifted index lists of free and fixed variables. */
	Indexlist rotatedFreee( n );
	Indexlist rotatedFixed( n );

	for( i=0; i<n; ++i )
	{
		switch ( getStatus( i ) )
		{
			case ST_INACTIVE:
				if ( rotatedFreee.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_LOWER:
				if ( rotatedFixed.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_UPPER:
				if ( rotatedFixed.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			default:
				return THROWERROR( RET_ROTATING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	freee = rotatedFreee;
	fixed = rotatedFixed;

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue Bounds::print( )
{
	if ( n == 0 )
		return SUCCESSFUL_RETURN;

	#ifndef __SUPPRESSANYOUTPUT__

	char myPrintfString[MAX_STRING_LENGTH];

	int_t nFR = getNFR( );
	int_t nFX = getNFX( );

	int_t* FR_idx;
	getFree( )->getNumberArray( &FR_idx );

	int_t* FX_idx;
	getFixed( )->getNumberArray( &FX_idx );

	snprintf( myPrintfString,MAX_STRING_LENGTH,"Bounds object comprising %d variables (%d free, %d fixed):\n",(int)n,(int)nFR,(int)nFX );
	myPrintf( myPrintfString );

	REFER_NAMESPACE_QPOASES print( FR_idx,nFR,"free " );
	REFER_NAMESPACE_QPOASES print( FX_idx,nFX,"fixed" );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue Bounds::clear( )
{
	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue Bounds::copy(	const Bounds& rhs
							)
{
	freee = rhs.freee;
	fixed = rhs.fixed;

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p A l l
 */
returnValue Bounds::setupAll( SubjectToStatus _status )
{
	int_t i;

	/* 1) Place unbounded variables at the beginning of the index list of free variables. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_UNBOUNDED )
		{
			if ( setupBound( i,_status ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 2) Add remaining (i.e. bounded but possibly free) variables to the index list of free variables. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_BOUNDED )
		{
			if ( setupBound( i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 3) Place implicitly fixed variables at the end of the index list of free variables. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_EQUALITY )
		{
			if ( setupBound( i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	/* 4) Moreover, add all bounds of unknown type. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_UNKNOWN || getType( i ) == ST_DISABLED )
		{
			if ( setupBound( i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_BOUND_FAILED );
		}
	}

	return SUCCESSFUL_RETURN;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
