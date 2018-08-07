/*
 * Translated from Python code at https://github.com/scipy/scipy/blob/master/scipy/signal/fir_filter_design.py

   Copyright (c) 2001, 2002 Enthought, Inc.
   All rights reserved.
  
   Copyright (c) 2003-2017 SciPy Developers.
   All rights reserved.
   
 
 */

#include <math.h>
#include "fir_filters.h"
#include "windows.h"

/*
 * 
 */
float sinc( float x ) {
  if( x == 0.0 )
    return 1.0;
  else
    return sin(M_PI * x) / (M_PI * x);
}


/*
 * 
 */
int firwin_Kaiser_bandpass(
    float *h, 
    int numtaps, 
    float cutoff_L, 
    float cutoff_H, 
    float beta,
    bool scale, 
    float fs )
{
    bool result = true;

#ifdef DEBUG
    Serial.println( "firwin_Kaiser_bandpass" );
    Serial.print( "numtaps = " );
    Serial.println( numtaps );
    Serial.print( "cutoff_L = " );
    Serial.println( cutoff_L );
    Serial.print( "cutoff_H = " );
    Serial.println( cutoff_H );
    Serial.print( "beta = " );
    Serial.println( beta );
    Serial.print( "scale = " );
    Serial.println( scale );
    Serial.print( "sampling frequency = " );
    Serial.println( fs );
#endif

    // nyq = 0.5 * _get_fs(fs, nyq)
    float nyq = 0.5 * fs;
#ifdef DEBUG
    Serial.print( "Nyquist frequency = " );
    Serial.println( nyq );
#endif

    // cutoff = np.atleast_1d(cutoff) / float(nyq)
    float cutoff[4] = { 0, cutoff_L/nyq, cutoff_H/nyq, 0 };
#ifdef DEBUG
    Serial.print( "cutoff_L/nyq = " );
    Serial.println( cutoff_L/nyq, 6 );
    Serial.print( "cutoff_H/nyq = " );
    Serial.println( cutoff_H/nyq, 6 );
#endif

    // Check for invalid input.
    // if cutoff.ndim > 1:
    //     raise ValueError("The cutoff argument must be at most "
    //                      "one-dimensional.")
    // if cutoff.size == 0:
    //     raise ValueError("At least one cutoff frequency must be given.")
    // if cutoff.min() <= 0 or cutoff.max() >= 1:
    //     raise ValueError("Invalid cutoff frequency: frequencies must be "
    //                      "greater than 0 and less than fs/2.")
    // if np.any(np.diff(cutoff) <= 0):
    //     raise ValueError("Invalid cutoff frequencies: the frequencies "
    //                      "must be strictly increasing.")

    // if width is not None:
        // A width was given.  Find the beta parameter of the Kaiser window
        // and set `window`.  This overrides the value of `window` passed in.
    //     atten = kaiser_atten(numtaps, float(width) / nyq)
    //     beta = kaiser_beta(atten)
    //     window = ('kaiser', beta)
    // KvD: for now, assume width was given (in Python it could be None, or optional)
/*    
 *     
 float kaiser_beta;
    {
      float atten = kaiser_atten( numtaps, width / nyq );
      float beta = kaiser_beta( atten );
      window = WINDOW_KAISER;
      kaiser_beta = beta;
    }
*/

    // pass_nyquist = bool(cutoff.size & 1) ^ pass_zero
    // KvD ^=xor, ^=and
//    bool pass_nyquist = false;  // (remember this is a bandpass)
    
    // if pass_nyquist and numtaps % 2 == 0:
    //     raise ValueError("A filter with an even number of coefficients must "
    //                      "have zero response at the Nyquist frequency.")

    // Insert 0 and/or 1 at the ends of cutoff so that the length of cutoff
    // is even, and each pair in cutoff corresponds to passband.
    // cutoff = np.hstack(([0.0] * pass_zero, cutoff, [1.0] * pass_nyquist))
    // KvD: This has been done on declaring cutoff

    // `bands` is a 2D array; each row gives the left and right edges of
    // a passband.
    // bands = cutoff.reshape(-1, 2)
    float bands[1][2] = { {cutoff[1], cutoff[2]} }; 
                       //{cutoff[0], cutoff[1]} };
#ifdef DEBUG
      Serial.print( "bands[0][0] = " );
      Serial.println( bands[0][0], 6 );
      Serial.print( "bands[0][1] = " );
      Serial.println( bands[0][1], 6 );
//      Serial.print( "bands[1][0] = " );
//      Serial.println( bands[1][0], 6 );
//      Serial.print( "bands[1][1] = " );
//      Serial.println( bands[1][1], 6 );
#endif
    

    // Build up the coefficients.
    // alpha = 0.5 * (numtaps - 1)
    float alpha = 0.5 * (numtaps - 1);
    // m = np.arange(0, numtaps) - alpha
    float m[numtaps];
    for( int i=0; i<numtaps; i++ )
      m[i] = i - alpha;
    // h = 0
#ifdef DEBUG
    for( int j=0; j<numtaps; j++ ) {
      Serial.print( "m[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( m[j], 6 );
    }
#endif
//    float h[numtaps]; // KvD: delivered externally
    for( int i=0; i<numtaps; i++ )
      h[i] = 0;
#ifdef DEBUG
    for( int j=0; j<numtaps; j++ ) {
      Serial.print( "h[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( h[j], 6 );
    }
#endif
    // for left, right in bands:
    int numBands = 1;
    for( int i=0; i<numBands; i++ )
    {
      float left = bands[i][0];
      float right = bands[i][1];
#ifdef DEBUG
      Serial.print( "left = " );
      Serial.println( left, 6 );
      Serial.print( "right = " );
      Serial.println( right, 6 );
#endif
    //     h += right * sinc(right * m)
    //     h -= left * sinc(left * m)
      for( int j=0; j<numtaps; j++ ) {
#ifdef DEBUG
      Serial.print( "right * m[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( right * m[j], 6 );
      
      Serial.print( "left * m[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( left * m[j], 6 );
      
      Serial.print( "sinc( right * m[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( sinc( right * m[j] ), 6 );
      
      Serial.print( "sinc( left * m[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( sinc( left * m[j] ), 6 );

      Serial.println();
#endif
        h[j] += right * sinc( right * m[j] );
        h[j] -= left * sinc( left * m[j] );
      }
    };
#ifdef DEBUG
    for( int j=0; j<numtaps; j++ ) {
      Serial.print( "h[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( h[j], 6 );
    }
#endif

    // Get and apply the window function.
    // from .signaltools import get_window
    // win = get_window(window, numtaps, fftbins=False)
    // h *= win
    float w[numtaps];
    // KvD: for now, only use Kaiser window
    result = result && kaiser( beta, numtaps, w );
    for( int i=0; i<numtaps; i++ )
      h[i] *= w[i];
#ifdef DEBUG
    Serial.println( "Kaiser window:" );
    for( int j=0; j<numtaps; j++ ) {
      Serial.print( "w[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( w[j], 6 );
    }
    Serial.println( "Unscaled coefficients:" );
    for( int j=0; j<numtaps; j++ ) {
      Serial.print( "h[" );
      Serial.print( j );
      Serial.print( "] = " );
      Serial.println( h[j] * 1000000, 6 );
    }
#endif
    // Now handle scaling if desired.
    // if scale:
    //     # Get the first passband.
    //     left, right = bands[0]
    //     if left == 0:
    //         scale_frequency = 0.0
    //     elif right == 1:
    //         scale_frequency = 1.0
    //     else:
    //         scale_frequency = 0.5 * (left + right)
    //     c = np.cos(np.pi * m * scale_frequency)
    //     s = np.sum(h * c)
    //     h /= s
    if( scale )
    {
      // Get the first passband
      float left = bands[0][0];
      float right = bands[0][1];
      float scale_frequency;
      if( left == 0 )
        scale_frequency = 0.;
      else if ( right == 1 )
        scale_frequency == 1.;
      else
        scale_frequency = 0.5 * (left + right);
#ifdef DEBUG
      Serial.print("Scaling frequency = ");
      Serial.println(scale_frequency,6);
#endif
      
      float s = 0.;
      for( int i=0; i<numtaps; i++ ) {
        float c = cos( M_PI * m[i] * scale_frequency );
        s += h[i] * c;
      }
      for( int i=0; i< numtaps; i++ )
        h[i] /= s;
#ifdef DEBUG
      Serial.println( "Scaled coefficients:" );
      for( int j=0; j<numtaps; j++ ) {
        Serial.print( "h[" );
        Serial.print( j );
        Serial.print( "] = " );
        Serial.println( h[j] * 1000000, 6 );
      }
#endif
    }

    return result;
}

