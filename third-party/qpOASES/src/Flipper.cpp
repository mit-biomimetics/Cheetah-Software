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
 *	\file src/Flipper.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the Flipper class designed to manage working sets of
 *	constraints and bounds within a QProblem.
 */


#include <qpOASES/Flipper.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	F l i p p e r
 */
Flipper::Flipper( )
{
	R = 0;
	Q = 0;
	T = 0;
	
	init( );
}


/*
 *	F l i p p e r
 */
Flipper::Flipper(	uint_t _nV,
					uint_t _nC
					)
{
	R = 0;
	Q = 0;
	T = 0;
	
	init( _nV,_nC );
}


/*
 *	F l i p p e r
 */
Flipper::Flipper( const Flipper& rhs )
{
	R = 0;
	Q = 0;
	T = 0;

	copy( rhs );
}


/*
 *	~ F l i p p e r
 */
Flipper::~Flipper( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
Flipper& Flipper::operator=( const Flipper& rhs )
{
	if ( this != &rhs )
	{
		clear( );
		copy( rhs );
	}

	return *this;
}



/*
 *	i n i t
 */
returnValue Flipper::init(	uint_t _nV,
							uint_t _nC
							)
{
	clear( );

	nV = _nV;
	nC = _nC;

	return SUCCESSFUL_RETURN;
}



/*
 *	g e t
 */
returnValue Flipper::get(	Bounds* const _bounds,
							real_t* const _R,
							Constraints* const _constraints,
							real_t* const _Q,
							real_t* const _T 
							) const
{
	if ( _bounds != 0 )
		*_bounds = bounds;

	if ( _constraints != 0 )
		*_constraints = constraints;

	if ( ( _R != 0 ) && ( R != 0 ) )
		memcpy( _R,R, nV*nV*sizeof(real_t) );

	if ( ( _Q != 0 ) && ( Q != 0 ) )
		memcpy( _Q,Q, nV*nV*sizeof(real_t) );

	if ( ( _T != 0 ) && ( T != 0 ) )
		memcpy( _T,T, getDimT()*sizeof(real_t) );

	return SUCCESSFUL_RETURN;
}


/*
 *	s e t
 */
returnValue Flipper::set(	const Bounds* const _bounds,
							const real_t* const _R,
							const Constraints* const _constraints,
							const real_t* const _Q,
							const real_t* const _T
							)
{
	if ( _bounds != 0 )
		bounds = *_bounds;

	if ( _constraints != 0 )
		constraints = *_constraints;

	if ( _R != 0 )
	{
		if ( R == 0 )
			R = new real_t[nV*nV];

		memcpy( R,_R, nV*nV*sizeof(real_t) );
	}

	if ( _Q != 0 )
	{
		if ( Q == 0 )
			Q = new real_t[nV*nV];

		memcpy( Q,_Q, nV*nV*sizeof(real_t) );
	}

	if ( _T != 0 )
	{
		if ( T == 0 )
			T = new real_t[getDimT()];

		memcpy( T,_T, getDimT()*sizeof(real_t) );
	}

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue Flipper::clear( )
{
	if ( R != 0 )
	{
		delete[] R;
		R = 0;
	}
	
	if ( Q != 0 )
	{
		delete[] Q;
		Q = 0;
	}
	
	if ( T != 0 )
	{
		delete[] T;
		T = 0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue Flipper::copy(	const Flipper& rhs
							)
{
	return set( &(rhs.bounds),rhs.R, &(rhs.constraints),rhs.Q,rhs.T );
}


uint_t Flipper::getDimT( ) const
{
	if ( nV > nC )
		return nC*nC;
	else
		return nV*nV;
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
