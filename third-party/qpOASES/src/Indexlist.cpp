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
 *	\file src/Indexlist.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
 *	\version 3.2
 *	\date 2007-2017
 *
 *	Implementation of the Indexlist class designed to manage index lists of
 *	constraints and bounds within a QProblem_SubjectTo.
 */


#include <qpOASES/Indexlist.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	I n d e x l i s t
 */
Indexlist::Indexlist( )
{
	number = 0;
	iSort  = 0;

	init( );
}


/*
 *	I n d e x l i s t
 */
Indexlist::Indexlist( int_t n )
{
	number = 0;
	iSort  = 0;

	init( n );
}


/*
 *	I n d e x l i s t
 */
Indexlist::Indexlist( const Indexlist& rhs )
{
	copy( rhs );
}


/*
 *	~ I n d e x l i s t
 */
Indexlist::~Indexlist( )
{
	clear( );
}


/*
 *	o p e r a t o r =
 */
Indexlist& Indexlist::operator=( const Indexlist& rhs )
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
returnValue Indexlist::init(	int_t n
								)
{
	if ( n < 0 )
		return THROWERROR( RET_INVALID_ARGUMENTS );

	clear( );

	length = 0;
	physicallength = n;

	if ( n > 0 )
	{
		number = new int_t[n];
		iSort  = new int_t[n];
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t N u m b e r A r r a y
 */
returnValue Indexlist::getNumberArray( int_t** const numberarray ) const
{
	if (numberarray == 0)
		return THROWERROR( RET_INVALID_ARGUMENTS );

	*numberarray = number;

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t I S o r t A r r a y
 */
returnValue Indexlist::getISortArray( int_t** const iSortArray ) const
{
	*iSortArray = iSort;

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t I n d e x
 */
int_t Indexlist::getIndex( int_t givennumber ) const
{
	int_t index = findInsert(givennumber);
	return number[iSort[index]] == givennumber ? iSort[index] : -1;
}


/*
 *	a d d N u m b e r
 */
returnValue Indexlist::addNumber( int_t addnumber )
{
	if ( length >= physicallength )
		return THROWERROR( RET_INDEXLIST_EXCEEDS_MAX_LENGTH );

	int_t i, j;
	number[length] = addnumber;
	j = findInsert(addnumber);
	for (i = length; i > j+1; i--)
		iSort[i] = iSort[i-1];
	iSort[j+1] = length;
	++length;

	return SUCCESSFUL_RETURN;
}


/*
 *	r e m o v e N u m b e r
 */
returnValue Indexlist::removeNumber( int_t removenumber )
{
	int_t i;
	int_t idx = findInsert( removenumber );
	int_t iSidx = iSort[idx];

	/* nothing to be done if number is not contained in index set */
	if ( number[iSidx] != removenumber )
		return SUCCESSFUL_RETURN;

	/* update sorted indices iSort first */
	for (i = 0; i < length; i++)
		if (iSort[i] > iSidx) iSort[i]--;
	for (i = idx+1; i < length; i++)
		iSort[i-1] = iSort[i];

	/* remove from numbers list */
	for( i=iSidx; i<length-1; ++i )
		number[i] = number[i+1];
	number[length-1] = -1;

	--length;

	return SUCCESSFUL_RETURN;
}


/*
 *	s w a p N u m b e r s
 */
returnValue Indexlist::swapNumbers( int_t number1, int_t number2 )
{
	int_t index1 = findInsert( number1 );
	int_t index2 = findInsert( number2 );

	/* consistency check */
	if ( ( number[iSort[index1]] != number1 ) || ( number[iSort[index2]] != number2 ) )
		return THROWERROR( RET_INDEXLIST_CORRUPTED );

	int_t tmp;
	/* swap numbers */
	tmp = number[iSort[index1]];
	number[iSort[index1]] = number[iSort[index2]];
	number[iSort[index2]] = tmp;
	/* swap sorting indices */
	tmp = iSort[index1];
	iSort[index1] = iSort[index2];
	iSort[index2] = tmp;

	return SUCCESSFUL_RETURN;
}



/*****************************************************************************
 *  P R O T E C T E D                                                        *
 *****************************************************************************/

/*
 *	c l e a r
 */
returnValue Indexlist::clear( )
{
	if ( iSort != 0 )
	{
		delete[] iSort;
		iSort = 0;
	}

	if ( number != 0 )
	{
		delete[] number;
		number = 0;
	}

	return SUCCESSFUL_RETURN;
}


/*
 *	c o p y
 */
returnValue Indexlist::copy(	const Indexlist& rhs
								)
{
	int_t i;

	length = rhs.length;
	physicallength = rhs.physicallength;

	if ( rhs.number != 0 )
	{
		number = new int_t[physicallength];
		for( i=0; i<physicallength; ++i )
			number[i] = rhs.number[i];
		iSort = new int_t[physicallength];
		for( i=0; i<physicallength; ++i )
			iSort[i] = rhs.iSort[i];
	}
	else
	{
		number = 0;
		iSort = 0;
	}

	return SUCCESSFUL_RETURN;
}

int_t Indexlist::findInsert(int_t i) const
{
	/* quick check if index can be appended */
	if (length == 0 || i < number[iSort[0]]) return -1;
	if (i >= number[iSort[length-1]]) return length-1;

	/* otherwise, perform bisection search */
	int_t fst = 0, lst = length-1, mid;

	while (fst < lst - 1)
	{
		mid = (fst + lst) / 2;
		if (i >= number[iSort[mid]]) fst = mid;
		else lst = mid;
	}

	return fst;
}

END_NAMESPACE_QPOASES


/*
 *	end of file
 */
