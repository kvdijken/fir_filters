# fir_filters

This library allows the computation of the coefficients for a FIR (Finite Impulse Response) filter. It is a direct translation
of the firwin code in the Python scipy library. The code for firwin was translated to C, as well as some other library functions 
in the scipy library. The scipy repository is at https://github.com/scipy/scipy/tree/master/scipy. The Python code is still
present in the C sources. This has been done to allow checking the validity of the code.

Not all functionality in firwin is translated. The firwin code allows all kind of FIR filters to be constructed: (multi) bandpass,
notch, lowpass and highpass filters, with all kind of windows. The code in this library only allows construction of (single) bandpass
filters using a Kaiser window. Including other filters and windows should not be so difficult.

Parameters which this library understand are
- numtaps: number of taps to be used
- cutoff_L: lower cutoff frequency
- cutoff_H: higher cutoff frequency
- beta: the beta value for the Kaiser window
- scale: 'true' if the coefficients should be scaled, 'false' if not
- fs: sampling frequency

What still needs to be done:
- testing: This code seems to be working, but it is still very new. There may still be bugs in it.
- error checking. At this moment the code is still _very_ rough, and does not perform any error checking
- include other filters
- include other windows

This code has been tested on a Teensy 3.6, compiled from within an Arduino environment.

Disclaimer: I am neither an expert on digital signal processing, nor a Python programmer.
