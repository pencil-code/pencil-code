sine-Gordon
===========

## Maintainer:

Axel Brandenburg <brandenb/nordita[dot]org>

## Added:

10-Oct-2024

## Status:

succeeds

## Recommended resolution:

128x1x1

## Comments:

reproduces the analytic solution
sep-6-2025/alberto: modified to use klein_gordon module

## Links:
* http://pencil-code.nordita.org/samples/1d-tests/sine-Gordon

## References:

*  [wikipedia](https://en.wikipedia.org/wiki/Sine-Gordon_equation)

## Suggested experiments:

* Try with lower precision: "DERIV=deriv_2nd" or higher precision: "DERIV=deriv_10th".
* Try with longer or shorter time steps and measure the speed of the wave.
* Try with larger resolution; determine the maximum permissible length of the time step dt < C<sub>CFL</sub> dx/v
