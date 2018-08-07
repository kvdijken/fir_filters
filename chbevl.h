#ifndef _chbevl_h_
#define _chbevl_h_

/*                                                     chbevl.c

   Copyright (c) 2001, 2002 Enthought, Inc.
   All rights reserved.
  
   Copyright (c) 2003-2017 SciPy Developers.
   All rights reserved.
   
       Evaluate Chebyshev series



   SYNOPSIS:

   int N;
   double x, y, coef[N], chebevl();

   y = chbevl( x, coef, N );



   DESCRIPTION:

   Evaluates the series

          N-1
           - '
    y  =   >   coef[i] T (x/2)
           -            i
          i=0

   of Chebyshev polynomials Ti at argument x/2.

   Coefficients are stored in reverse order, i.e. the zero
   order term is last in the array.  Note N is the number of
   coefficients, not the order.

   If coefficients are for the interval a to b, x must
   have been transformed to x -> 2(2x - b - a)/(b-a) before
   entering the routine.  This maps x from (a, b) to (-1, 1),
   over which the Chebyshev polynomials are defined.

   If the coefficients are for the inverted interval, in
   which (a, b) is mapped to (1/b, 1/a), the transformation
   required is x -> 2(2ab/x - b - a)/(b-a).  If b is infinity,
   this becomes x -> 4a/x - 1.



   SPEED:

   Taking advantage of the recurrence properties of the
   Chebyshev polynomials, the routine requires one more
   addition per loop than evaluating a nested polynomial of
   the same degree.

*/

double chbevl( double x, double array[], int n );

#endif
