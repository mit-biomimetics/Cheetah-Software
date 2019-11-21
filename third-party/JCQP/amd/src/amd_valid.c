/* ========================================================================= */
/* === AMD_valid =========================================================== */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* Check if a column-form matrix is valid or not.  The matrix A is
 * n_row-by-n_col.  The row indices of entries in column j are in
 * Ai [Ap [j] ... Ap [j+1]-1].  Required conditions are:
 *
 *	n_row >= 0
 *	n_col >= 0
 *	nz = Ap [n_col] >= 0	    number of entries in the matrix
 *	Ap [0] == 0
 *	Ap [j] <= Ap [j+1] for all j in the range 0 to n_col.
 *      Ai [0 ... nz-1] must be in the range 0 to n_row-1.
 *
 * If any of the above conditions hold, AMD_INVALID is returned.  If the
 * following condition holds, AMD_OK_BUT_JUMBLED is returned (a warning,
 * not an error):
 *
 *	row indices in Ai [Ap [j] ... Ap [j+1]-1] are not sorted in ascending
 *	    order, and/or duplicate entries exist.
 *
 * Otherwise, AMD_OK is returned.
 *
 * In v1.2 and earlier, this function returned TRUE if the matrix was valid
 * (now returns AMD_OK), or FALSE otherwise (now returns AMD_INVALID or
 * AMD_OK_BUT_JUMBLED).
 */

#include "amd_internal.h"

GLOBAL Int AMD_valid
(
    /* inputs, not modified on output: */
    Int n_row,		/* A is n_row-by-n_col */
    Int n_col,
    const Int Ap [ ],	/* column pointers of A, of size n_col+1 */
    const Int Ai [ ]	/* row indices of A, of size nz = Ap [n_col] */
)
{
    Int nz, j, p1, p2, ilast, i, p, result = AMD_OK ;

    if (n_row < 0 || n_col < 0 || Ap == NULL || Ai == NULL)
    {
	return (AMD_INVALID) ;
    }
    nz = Ap [n_col] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	AMD_DEBUG0 (("column 0 pointer bad or nz < 0\n")) ;
	return (AMD_INVALID) ;
    }
    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	AMD_DEBUG2 (("\nColumn: "ID" p1: "ID" p2: "ID"\n", j, p1, p2)) ;
	if (p1 > p2)
	{
	    /* column pointers must be ascending */
	    AMD_DEBUG0 (("column "ID" pointer bad\n", j)) ;
	    return (AMD_INVALID) ;
	}
	ilast = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    AMD_DEBUG3 (("row: "ID"\n", i)) ;
	    if (i < 0 || i >= n_row)
	    {
		/* row index out of range */
		AMD_DEBUG0 (("index out of range, col "ID" row "ID"\n", j, i));
		return (AMD_INVALID) ;
	    }
	    if (i <= ilast)
	    {
		/* row index unsorted, or duplicate entry present */
		AMD_DEBUG1 (("index unsorted/dupl col "ID" row "ID"\n", j, i));
		result = AMD_OK_BUT_JUMBLED ;
	    }
	    ilast = i ;
	}
    }
    return (result) ;
}
