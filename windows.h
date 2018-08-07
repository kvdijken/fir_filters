#ifndef _windows_h_
#define _windows_h_


/*
 * Calculates a Kaiser window for use in filter construction.
 * 
 * Parameters:
 * beta: Shape parameter, determines trade-off between main-lobe width and
 *       side lobe level. As beta gets large, the window narrows.
 * M: number of points in the output window
 * window: pointer to an array of n float's
 * 
 * Returns:
 * The window, with the maximum value normalized to 1 (though the value 1
 * does not appear if `n` is even).
 */
 
int kaiser( float beta, int M, float *window );

#endif
