#ifndef _i0_h_
#define _i0_h_

/*
   From scipy code at https://github.com/scipy/scipy/blob/master/scipy/special/cephes/i0.h

   Copyright (c) 2001, 2002 Enthought, Inc.
   All rights reserved.
  
   Copyright (c) 2003-2017 SciPy Developers.
   All rights reserved.
   

*/

/*                                                     i0.c

       Modified Bessel function of order zero



   SYNOPSIS:

   double x, y, i0();

   y = i0( x );



   DESCRIPTION:

   Returns modified Bessel function of order zero of the
   argument.

   The function is defined as i0(x) = j0( ix ).

   The range is partitioned into the two intervals [0,8] and
   (8, infinity).  Chebyshev polynomial expansions are employed
   in each interval.



   ACCURACY:

                        Relative error:
   arithmetic   domain     # trials      peak         rms
      IEEE      0,30        30000       5.8e-16     1.4e-16

*/
/*              i0e.c

    Modified Bessel function of order zero,
    exponentially scaled



   SYNOPSIS:

   double x, y, i0e();

   y = i0e( x );



   DESCRIPTION:

   Returns exponentially scaled modified Bessel function
   of order zero of the argument.

   The function is defined as i0e(x) = exp(-|x|) j0( ix ).



   ACCURACY:

                        Relative error:
   arithmetic   domain     # trials      peak         rms
      IEEE      0,30        30000       5.4e-16     1.2e-16
   See i0().

*/
double i0(double x);

#endif
