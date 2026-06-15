PRO FLOW_REB_CUBE, DATA, DIM, STEP, STD, REB, VX, VY, BOXCAR=boxcar, ADIF=adif, CORR=corr, QFIT2=qfit2, CROSSD=crossd, SILENT=silent
;+
; NAME:
;       FLOW_REB_CUBE
;
; PURPOSE:
;       Compute flow maps.
;
; INPUTS:
;       DATA    3D data cube e.g. [x,y,time] with DIM=3
;       DIM     number of the time dimension (1-3)
;       STEP    distance between reference and life image
;       STD     width for smoothing window
;       REB     factor to "congrid" the flow map
;
; KEYWORDS:
;       BOXCAR  if set, a boxcar window of width STD is used. Hence, STD must be an odd number
;       ADIF    uses an absolute differences algorithm
;       CORR    uses a multiplicative algorithm. Default is the sum of square of the local differences
;       QFIT2   uses 9 points fitting procedure
;       CROSSD  uses cross derivative interpolation formulae
;       SILENT  if set, suppresses displaying the size of the array at the terminal
;
; OUTPUTS:
;       VX,VY   proper motion map: X and Y components
;
; SIDE EFFECTS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; PROCEDURE:
;       It uses November's method of shifting BOTH images. The defaults for
;       the different methods are 1) square differences for the matching,
;       2) gaussian convolution for the smoothing, and 3) FIVEPOINT for the
;       subpixel extrema finding procedure.
;
; EXAMPLE:
;       Let DATA be an array of image frames with time in the DIM dimension.
;       For an image resolution of 0.2" per pixel, we want to compute a flow
;       map using a FWHM of 4" for the smoothing gaussian window:
;
;       IDL> FLOW_REB_CUBE, data, dim, step, 8.5, reb, vx, vy, /silent
;
;       FWHM of gaussian window is 8.5 * .2 * 2.355 = 4 arcsec.
;
; REFERENCES:
;       November, L.J. and Simon, G.W.: 1988, Ap.J., 333, 427
;       Darvann, T.: 1991, Master's Thesis, University of Oslo.
;
; MODIFICATION HISTORY:
;       Written by Roberto Luis Molowny Horas, May 1992.
;       Keywords added in July 1992.
;       Adapted by Philippe-A. Bourdin, June 2010.
;
;-
ON_ERROR, 2

	s1 = size(data)
	CASE dim OF
		1 : BEGIN
			nim = s1[1]
			xdim = s1[2]
			ydim = s1[3]
			a = reform(data[0,*,*], xdim, ydim)
		END
		2 : BEGIN
			nim = s1[2]
			xdim = s1[1]
			ydim = s1[3]
			a = reform(data[*,0,*], xdim, ydim)
		END
		3 : BEGIN
			nim = s1[3]
			xdim = s1[1]
			ydim = s1[2]
			a = reform(data[*,*,0], xdim, ydim)
		END
	ENDCASE
	n = nim - step

	aa = congrid(a, xdim/reb, ydim/reb, /interp)
	s = size(aa)

	IF KEYWORD_SET(silent) THEN silent = 0

	; cummulative correlation function
	cc = FLTARR(s[1], s[2], 3, 3)

	FOR k = 0,n-1 DO BEGIN

		CASE dim OF
			1 : BEGIN
				a = reform(data[k,*,*], xdim, ydim)
				b = reform(data[k+step,*,*], xdim, ydim)
			END
			2 : BEGIN
				a = reform(data[*,k,*], xdim, ydim)
				b = reform(data[*,k+step,*], xdim, ydim)
			END
			3 : BEGIN
				a = reform(data[*,*,k], xdim, ydim)
				b = reform(data[*,*,k+step], xdim, ydim)
			END
		ENDCASE

		aa = congrid(a, xdim/reb, ydim/reb, /interp)
		bb = congrid(b, xdim/reb, ydim/reb, /interp)

		; remove mean
		aa = aa - TOTAL(aa) / s[4]
		bb = bb - TOTAL(bb) / s[4]

		FOR i = -1,1 DO BEGIN
			FOR j = -1,1 DO BEGIN
				; method selection:
				CASE 1 OF
					KEYWORD_SET(adif): BEGIN
						; absolute differences
						cc[0,0,i+1,j+1] = cc[*,*,i+1,j+1] + ABS(SHIFT(aa,i,j) - SHIFT(bb,-i,-j))
					END
					KEYWORD_SET(corr): BEGIN
						; cross products
						cc[0,0,i+1,j+1] = cc[*,*,i+1,j+1] + SHIFT(aa,i,j) * SHIFT(bb,-i,-j)
					END
					ELSE: BEGIN
						; square differences
						dumb = SHIFT(aa,i,j) - SHIFT(bb,-i,-j)
						; is this really faster than (...)^2 in modern versions of IDL?
						cc[0,0,i+1,j+1] = cc[*,*,i+1,j+1] + dumb*dumb
						dumb = 0
					END
				ENDCASE
			END
		END
		aa[*,*] = 0.
		bb[*,*] = 0.
	ENDFOR

	; treat the edges
	cc[0,0,0,0] = cc[1,*,*,*]
	cc[0,0,0,0] = cc[*,1,*,*]
	cc[s[1]-1,0,0,0] = cc[s[1]-2,*,*,*]
	cc[0,s[2]-1,0,0] = cc[*,s[2]-2,*,*]

	std_reb = std/reb
	IF KEYWORD_SET(boxcar) THEN BEGIN
		; boxcar smoothing
		FOR i = 0,2 DO BEGIN
			FOR j = 0,2 DO BEGIN
				cc[0,0,i,j] = SMOOTHE(cc[*,*,i,j], std_reb)
			END
		END
	END ELSE BEGIN
		; Gausian convolution
		FOR i = 0,2 DO BEGIN
			FOR j = 0,2 DO BEGIN
				cc[0,0,i,j] = GCONVOL(cc[*,*,i,j], std_reb)
			END
		END
	END

	CASE 1 OF
		KEYWORD_SET(qfit2): BEGIN
			; 9-ppoint fitting
			QFIT2, cc, vx, vy
		END
		KEYWORD_SET(crossd): BEGIN
			; cross derivative
			CROSSD, cc, vx, vy
		END
		ELSE: BEGIN
			; default
			FIVEPOINT, cc, vx, vy
		END
	ENDCASE
	cc = 0

	; compute the original size and values of the flow map
	vx = congrid(vx, xdim, ydim, /interp) * reb
	vy = congrid(vy, xdim, ydim, /interp) * reb

	; scale the result
	vx = 2 * vx
	vy = 2 * vy
END

