function convolve, image, psf, FT_PSF=psf_FT, FT_IMAGE=imFT, NO_FT=noft, $
			CORRELATE=correlate, AUTO_CORRELATION=auto
;+
; NAME:
;	CONVOLVE
; PURPOSE:
;	Convolution of an image with a Point Spread Function (PSF)
; EXPLANATION:
;	The default is to compute the convolution using a product of 
;	Fourier transforms (for speed).
;
; CALLING SEQUENCE:
;
;	imconv = convolve( image1, psf, FT_PSF = psf_FT )
;  or:
;	correl = convolve( image1, image2, /CORREL )
;  or:
;	correl = convolve( image, /AUTO )
;
; INPUTS:
;	image = 2-D array (matrix) to be convolved with psf
;	psf = the Point Spread Function, (size < or = to size of image).
;
; OPTIONAL INPUT KEYWORDS:
;
;	FT_PSF = passes out/in the Fourier transform of the PSF,
;		(so that it can be re-used the next time function is called).
;	FT_IMAGE = passes out/in the Fourier transform of image.
;
;	/CORRELATE uses the conjugate of the Fourier transform of PSF,
;		to compute the cross-correlation of image and PSF,
;		(equivalent to IDL function convol() with NO rotation of PSF)
;
;	/AUTO_CORR computes the auto-correlation function of image using FFT.
;
;	/NO_FT overrides the use of FFT, using IDL function convol() instead.
;		(then PSF is rotated by 180 degrees to give same result)
; METHOD:
;	When using FFT, PSF is centered & expanded to size of image.
; HISTORY:
;	written, Frank Varosi, NASA/GSFC 1992.
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
	sp = size( psf_FT )  &  sif = size( imFT )
	sim = size( image )  &  sc = sim/2  &  npix = N_elements( image )

	if (sim[0] NE 2) OR keyword_set( noft ) then begin
		if keyword_set( auto ) then begin
			message,"auto-correlation only for images with FFT",/INF
			return, image
		  endif else if keyword_set( correlate ) then $
				return, convol( image, psf ) $
			else	return, convol( image, rotate( psf, 2 ) )
	   endif

	if (sif[0] NE 2) OR (sif[sif[0]+1] NE 6) OR $
	   (sif[1] NE sim[1]) OR (sif[2] NE sim[2]) then imFT = FFT( image,-1 )

	if keyword_set( auto ) then $
	 return, shift( npix*float( FFT( imFT*conj( imFT ),1 ) ), sc[1],sc[2] )

	if (sp[0] NE 2) OR (sp[sp[0]+1] NE 6) OR $
	   (sp[1] NE sim[1]) OR (sp[2] NE sim[2]) then begin
		sp = size( psf )
		if (sp[0] NE 2) then begin
			message,"must supply PSF matrix (2nd arg.)",/INFO
			return, image
		   endif
		Loc = ( sc - sp/2 ) > 0		;center PSF in new array,
		s = (sp/2 - sc) > 0	   ;handle all cases: smaller or bigger
		L = (s + sim-1) < (sp-1)
		psf_FT = complexarr( sim[1], sim[2] )
		psf_FT[ Loc[1], Loc[2] ] = psf[ s[1]:L[1], s[2]:L[2] ]
		psf_FT = FFT( psf_FT, -1, /OVERWRITE )
	   endif

	if keyword_set( correlate ) then $
		conv = npix * float( FFT( imFT * conj( psf_FT ), 1 ) ) $
	  else	conv = npix * float( FFT( imFT * psf_FT, 1 ) )

	sc = sc + (sim MOD 2)	;shift correction for odd size images.

return, shift( conv, sc[1], sc[2] )
end
