function gaussian, xi, parms, pderiv, DOUBLE=double
;+
; NAME:
;       GAUSSIAN
; PURPOSE:
;       Compute the 1-d Gaussian function and optionally the derivative
; EXPLANATION:
;       Compute the 1-D Gaussian function and optionally the derivative 
;       at an array of points.
;
; CALLING SEQUENCE:
;       y = gaussian( xi, parms,[ pderiv ])
;
; INPUTS:
;       xi = array, independent variable of Gaussian function.
;
;       parms = parameters of Gaussian, 2, 3 or 4 element array:
;               parms[0] = maximum value (factor) of Gaussian,
;               parms[1] = mean value (center) of Gaussian,
;               parms[2] = standard deviation (sigma) of Gaussian.
;               (if parms has only 2 elements then sigma taken from previous
;               call to gaussian(), which is stored in a  common block).
;               parms[3] = optional, constant offset added to Gaussian.
; OUTPUT:
;       y -  Function returns array of Gaussian evaluated at xi.    Values will
;            be floating pt. (even if xi is double) unless the /DOUBLE keyword
;            is set.
;
; OPTIONAL INPUT:
;       /DOUBLE - set this keyword to return double precision for both
;             the function values and (optionally) the partial derivatives.
; OPTIONAL OUTPUT:
;       pderiv = [N,3] or [N,4] output array of partial derivatives,
;               computed only if parameter is present in call.
;
;               pderiv[*,i] = partial derivative at all xi absisca values
;               with respect to parms[i], i=0,1,2,[3].
;
;
; EXAMPLE:
;       Evaulate a Gaussian centered at x=0, with sigma=1, and a peak value
;       of 10 at the points 0.5 and 1.5.   Also compute the derivative
;
;       IDL> f = gaussian( [0.5,1.5], [10,0,1], DERIV )
;       ==> f= [8.825,3.25].   DERIV will be a 2 x 3 array containing the
;       numerical derivative at the two points with respect to the 3 parameters.
; 
; COMMON BLOCKS:
;       None
; HISTORY:
;       Written, Frank Varosi NASA/GSFC 1992.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Use machar() for machine precision, added /DOUBLE keyword,
;       add optional constant 4th parameter    W. Landsman   November 2001
;-
  On_error,2
  common gaussian, sigma

  if N_params() LT 2 then begin
        print,'Syntax - y = GAUSSIAN( xi, parms,[ pderiv, /DOUBLE ])'
        print,'         parms[0] = maximum value (factor) of Gaussian'
        print,'         parms[1] = mean value (center) of Gaussian'
        print,'         parms[2] = standard deviation (sigma) of Gaussian'
        print,'         parms[3] = optional constant to be added to Gaussian'
        return, -1
  endif

  common gaussian, sigma

        Nparmg = N_elements( parms )
        npts = N_elements(xi) 
        ptype = size(parms,/type)
        if (ptype LE 3) or (ptype GE 12) then parms = float(parms)
        if (Nparmg GE 3) then sigma = parms[2]

        double = keyword_set(DOUBLE)
        if double then $       ;Double precision?
            gauss = dblarr( npts ) else $
            gauss = fltarr( npts )
 
        z = ( xi - parms[1] )/sigma
        zz = z*z

; Get smallest value expressible on computer.   Set lower values to 0 to avoid
; floating underflow
        minexp = alog((machar(DOUBLE=double)).xmin)     
 
        w = where( zz LT -2*minexp, nw )
        if (nw GT 0) then gauss[w] = exp( -zz[w] / 2 )

        if N_params() GE 3 then begin

                if double then $ 
                pderiv = dblarr( npts, Nparmg ) else $
                pderiv = fltarr( npts, Nparmg )
                fsig = parms[0] / sigma

                pderiv[0,0] = gauss
                pderiv[0,1] = gauss * z * fsig

                if (Nparmg GE 3) then  pderiv[0,2] = gauss * zz * fsig
                if (Nparmg GE 4) then  pderiv[0,3] = replicate(1, npts)
           endif

 if Nparmg LT 4 then return, parms[0] * gauss else $
                     return, parms[0] * gauss + parms[3]
 end
