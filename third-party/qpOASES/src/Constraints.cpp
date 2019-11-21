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
 *	\file src/Constraints.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the Constraints class designed to manage working sets of
 *	constraints within a QProblem.
 */


#include <qpOASES/Constraints.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	C o n s t r a i n t s
 */
Constraints::Constraints( ) : SubjectTo( )
{
}


/*
 *	C o n s t r a i n t s
 */
Constraints::Constraints( int_t _n ) : SubjectTo( _n )
{
	init( _n );
}


/*
 *	C o n s t r a i n t s
 */
Constraints::Constraints( const Constraints& rhs ) : SubjectTo( rhs )
{
	copy( rhs );
}


/*
 *	~ C o n s t r a i n t s
 */
Constraints::~Constraints( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
Constraints& Constraints::operator=( const Constraints& rhs )
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
returnValue Constraints::init(	int_t _n
								)
{
	if ( _n < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	clear( );

	if ( _n >= 0 )
	{
		active.init(   _n );
		inactive.init( _n );
	}

	return SubjectTo::init( _n );
}



/*
 *	s e t u p C o n s t r a i n t
 */
returnValue Constraints::setupConstraint(	int_t number, SubjectToStatus _status
											)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Add constraint index to respective index list. */
	switch ( _status )
	{
		case ST_INACTIVE:
			if ( this->addIndex( this->getInactive( ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
			break;

		case ST_LOWER:
			if ( this->addIndex( this->getActive( ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
			break;

		case ST_UPPER:
			if ( this->addIndex( this->getActive( ),number,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
			break;

		default:
			return THROWERROR( RET_INVALID_ARGUMENTS );
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t u p A l l I n a c t i v e
 */
returnValue Constraints::setupAllInactive( )
{
	return setupAll( ST_INACTIVE );
}


/*
 *	s e t u p A l l L o w e r
 */
returnValue Constraints::setupAllLower( )
{
	return setupAll( ST_LOWER );
}


/*
 *	s e t u p A l l U p p e r
 */
returnValue Constraints::setupAllUpper( )
{
	return setupAll( ST_UPPER );
}


/*
 *	m o v e A c t i v e T o I n a c t i v e
 */
returnValue Constraints::moveActiveToInactive( int_t number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of active constraints to that of inactive ones. */
	if ( this->removeIndex( this->getActive( ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( this->addIndex( this->getInactive( ),number,ST_INACTIVE ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	m o v e I n a c t i v e T o A c t i v e
 */
returnValue Constraints::moveInactiveToActive(	int_t number, SubjectToStatus _status
												)
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	/* Move index from indexlist of inactive constraints to that of active ones. */
	if ( this->removeIndex( this->getInactive( ),number ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	if ( this->addIndex( this->getActive( ),number,_status ) != SUCCESSFUL_RETURN )
		return THROWERROR( RET_MOVING_BOUND_FAILED );

	return SUCCESSFUL_RETURN;
}


/*
 *	f l i p F i x e d
 */
returnValue Constraints::flipFixed( int_t number )
{
	/* consistency check */
	if ( ( number < 0 ) || ( number >= n ) )
		return THROWERROR( RET_INDEX_OUT_OF_BOUNDS );

	if ( status != 0 )
		switch (status[number])
		{
			case ST_LOWER: status[number] = ST_UPPER; break;
			case ST_UPPER: status[number] = ST_LOWER; break;
			default: return THROWERROR( RET_MOVING_CONSTRAINT_FAILED );
		}

	return SUCCESSFUL_RETURN;
}


/*
 *	s h i f t
 */
returnValue Constraints::shift( int_t offset )
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
	Indexlist shiftedActive( n );
	Indexlist shiftedInactive( n );

	for( i=0; i<n; ++i )
	{
		switch ( getStatus( i ) )
		{
			case ST_INACTIVE:
				if ( shiftedInactive.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_LOWER:
				if ( shiftedActive.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			case ST_UPPER:
				if ( shiftedActive.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_SHIFTING_FAILED );
				break;

			default:
				return THROWERROR( RET_SHIFTING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	active = shiftedActive;
	inactive = shiftedInactive;

	return SUCCESSFUL_RETURN;
}


/*
 *	r o t a t e
 */
returnValue Constraints::rotate( int_t offset )
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
	Indexlist rotatedActive( n );
	Indexlist rotatedInactive( n );

	for( i=0; i<n; ++i )
	{
		switch ( getStatus( i ) )
		{
			case ST_INACTIVE:
				if ( rotatedInactive.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_LOWER:
				if ( rotatedActive.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			case ST_UPPER:
				if ( rotatedActive.addNumber( i ) != SUCCESSFUL_RETURN )
					return THROWERROR( RET_ROTATING_FAILED );
				break;

			default:
				return THROWERROR( RET_ROTATING_FAILED );
		}
	}

	/* 3) Assign shifted index list. */
	active = rotatedActive;
	inactive = rotatedInactive;

	return SUCCESSFUL_RETURN;
}


/*
 *	p r i n t
 */
returnValue Constraints::print( )
{
	if ( n == 0 )
		return SUCCESSFUL_RETURN;

	#ifndef __SUPPRESSANYOUTPUT__

	char myPrintfString[MAX_STRING_LENGTH];

	int_t nIAC = getNIAC( );
	int_t nAC  = getNAC( );

	int_t* IAC_idx;
	getInactive( )->getNumberArray( &IAC_idx );

	int_t* AC_idx;
	getActive( )->getNumberArray( &AC_idx );

	snprintf( myPrintfString,MAX_STRING_LENGTH,"Constraints object comprising %d constraints (%d inactive, %d active):\n",(int)n,(int)nIAC,(int)nAC );
	myPrintf( myPrintfString );

	REFER_NAMESPACE_QPOASES print( IAC_idx,nIAC,"inactive" );
	REFER_NAMESPACE_QPOASES print( AC_idx, nAC, "active  " );

	#endif /* __SUPPRESSANYOUTPUT__ */

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue Constraints::clear( )
{
	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue Constraints::copy(	const Constraints& rhs
								)
{
	active   = rhs.active;
	inactive = rhs.inactive;

	return SUCCESSFUL_RETURN;
}



/*
 *	s e t u p A l l
 */
returnValue Constraints::setupAll( SubjectToStatus _status )
{
	int_t i;

	/* 1) Place unbounded constraints at the beginning of the index list of inactive constraints. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_UNBOUNDED )
		{
			if ( setupConstraint( i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}

	/* 2) Add remaining (i.e. "real" inequality) constraints to the index list of inactive constraints. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_BOUNDED )
		{
			if ( setupConstraint( i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}

	/* 3) Place implicit equality constraints at the end of the index list of inactive constraints. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_EQUALITY )
		{
			if ( setupConstraint( i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}

	/* 4) Moreover, add all constraints of unknown type. */
	for( i=0; i<n; ++i )
	{
		if ( getType( i ) == ST_UNKNOWN || getType( i ) == ST_DISABLED )
		{
			if ( setupConstraint( i,_status ) != SUCCESSFUL_RETURN )
				return THROWERROR( RET_SETUP_CONSTRAINT_FAILED );
		}
	}


	return SUCCESSFUL_RETURN;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
