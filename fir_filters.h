#ifndef _fir_filters_h_
#define _fir_filters_h_


/**
 *     firwin_Kaiser_bandpass
 *     
 *     Computes the coefficients for a bandpass filter using a Kaiser window.
 *     
 *     
    Parameters
    ----------
    h : array[numtaps] with the computed taps
    numtaps : int
        Length of the filter (number of coefficients, i.e. the filter
        order + 1).
    cutoff_L : Low cutoff frequency of filter (expressed in the same units as `nyq`)
    cutoff_H : High cutoff frequency of filter (expressed in the same units as `nyq`)
    beta : Beta value for the Kaiser window
    scale : bool, optional
        Set to True to scale the coefficients so that the frequency
        response is exactly unity at a certain frequency.
        That frequency is either:

        - 0 (DC) if the first passband starts at 0 (i.e. pass_zero
          is True)
        - `nyq` (the Nyquist frequency) if the first passband ends at
          `nyq` (i.e the filter is a single band highpass filter);
          center of first passband otherwise

    fs : The sampling frequency of the signal.  Each frequency in `cutoff`
         must be between 0 and ``fs/2``.  Default is 2.

    Returns
    -------
    h : (numtaps,) ndarray
        Coefficients of length `numtaps` FIR filter.


 */
 
int firwin_Kaiser_bandpass(
    float *h, 
    int numtaps, 
    float cutoff_L, 
    float cutoff_H, 
    float beta,
    bool scale, 
    float fs );

#endif
