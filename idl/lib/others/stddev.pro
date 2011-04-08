;$Id: //depot/idl/IDL_71/idldir/lib/stddev.pro#1 $
;
; Copyright (c) 1997-2009, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       STDDEV
;
; PURPOSE:
;       This function computes the stddev of an N-element vector. 
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = stddev(X)
;
; INPUTS:
;       X:      An N-element vector of type integer, float or double.
;
; KEYWORD PARAMETERS:
;       DOUBLE: IF set to a non-zero value, computations are done in
;               double precision arithmetic.
;       NAN:    If set, treat NaN data as missing.
;
; EXAMPLE:
;       Define the N-element vector of sample data.
;         x = [65, 63, 67, 64, 68, 62, 70, 66, 68, 67, 69, 71, 66, 65, 70]
;       Compute the standard deviation.
;         result = stddev(x)
;       The result should be:
;       8.16292
;
; PROCEDURE:
;       STDDEV calls the IDL function MOMENT.
;
; REFERENCE:
;       APPLIED STATISTICS (third edition)
;       J. Neter, W. Wasserman, G.A. Whitmore
;       ISBN 0-205-10328-6
;
; MODIFICATION HISTORY:
;       Written by:  GSL, RSI, August 1997
;-
FUNCTION STDDEV, X, Double = Double, NaN = NaN

  ON_ERROR, 2

  RETURN, sqrt((moment( X, Double=Double, Maxmoment=2, NaN = NaN ))[1])
END
