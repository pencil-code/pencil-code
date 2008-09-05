;$Id$
;
; Copyright (c) 1997-1999, Research Systems, Inc.  All rights reserved.
;       Unauthorized reproduction prohibited.
;+
; NAME:
;       MEAN
;
; PURPOSE:
;       This function computes the mean of an N-element vector. 
;
; CATEGORY:
;       Statistics.
;
; CALLING SEQUENCE:
;       Result = MEAN(X)
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
;         result = MEAN(x)
;       The result should be:
;       66.7333
;
; PROCEDURE:
;       MEAN calls the IDL function MOMENT.
;
; REFERENCE:
;       APPLIED STATISTICS (third edition)
;       J. Neter, W. Wasserman, G.A. Whitmore
;       ISBN 0-205-10328-6
;
; MODIFICATION HISTORY:
;       Written by:  GSL, RSI, August 1997
;-
FUNCTION MEAN, X, Double = Double, NaN = nan
  COMPILE_OPT IDL2,HIDDEN
  ON_ERROR, 2

  RETURN, (moment( X, Double=Double, Maxmoment=1, NaN=nan ))[0]
END
