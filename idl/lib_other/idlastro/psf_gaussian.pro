function psf_gaussian, parameters, NPIXEL=npix, NDIMENSION=ndim, FWHM=fwhm,  $
                        DOUBLE = double, CENTROID=cntrd, ST_DEV=st_dev,  $
                        XY_CORREL=xy_corr, NORMALIZE=normalize
;+
; NAME:
;       PSF_GAUSSIAN
;
; PURPOSE:
;       Create a 1-d, 2-d, or 3-d Gaussian with specified FWHM, center 
; EXPLANATION:
;       Return a point spread function having Gaussian profiles,
;       as either a 1D vector, a 2D image, or 3D volumetric-data.
;
; CALLING SEQUENCE:
;       psf = psf_Gaussian( NPIXEL=, FWHM= , CENTROID = 
;                     [ /DOUBLE, /NORMALIZE, ST_DEV=,  NDIMEN= ] ) 
; or:
;       psf = psf_Gaussian( parameters, NPIXEL = ,NDIMEN = )
;
; REQUIRED INPUT KEYWORD:
;       NPIXEL = number pixels for each dimension, specify as an array,
;               or just one number to make all sizes equal.
;
; OPTIONAL KEYWORDS:
;       CENTROID = floating scalar or vector giving position of  PSF center.    
;               default is exact center of requested vector/image/volume.
;               The number of elements in CENTROID should equal the number of
;               dimensions.    **The definition of Centroid was changed in
;               March 2002, and now an integer defines the center of a pixel.**
;
;       /DOUBLE  = If set, then the output array is computed in double precision
;               the default is to return a floating point array.
;
;       FWHM = the desired Full-Width Half-Max (pixels) in each dimension,
;               specify as an array, or single number to make all the same.
;
;       NDIMEN = integer dimension of result: either 1 (vector), 2 (image), or 
;                3 (volume), default = 2 (an image result).
;
;       /NORMALIZE causes resulting PSF to be normalized so Total( psf ) = 1.
;
;       ST_DEV = optional way to specify width by standard deviation param.
;                Ignored if FWHM is specified.
;
;       XY_CORREL = scalar between 0 and 1 specifying correlation coefficient
;               Use this keyword, for example, to specify an elliptical 
;               Gaussian oriented at an angle to the X,Y axis.   Only valid
;               for 2-dimensional case.
;
;
; INPUTS (optional):
;
;       parameters = an NDIMEN by 3 array giving for each dimension:
;                       [ maxval, center, st_dev ],  overrides other keywords.
;
; EXAMPLE:
;       (1) Create a 31 x 31 array containing a normalized centered Gaussian 
;       with an X FWHM = 4.3 and a Y FWHM = 3.6
;
;       IDL> array = PSF_GAUSSIAN( Npixel=31, FWHM=[4.3,3.6], /NORMAL )
;
;       (2) Create a 50 pixel 1-d Gaussian vector with a maximum of 12, 
;          centered at  pixel 23 with a sigma of 19.2
;
;       IDL> psf = psf_gaussian([12,23,19.2],npixel=50)
; EXTERNAL CALLS:
;       function Gaussian()
; NOTES:
;       To improve speed, floating underflow exceptions are suppressed (using 
;       the MASK=32  keyword of CHECK_MATH() rather than being flagged.
;
; HISTORY:
;       Written, Frank Varosi NASA/GSFC 1991.
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Suppress underflow messages, add DOUBLE keyword. **Modified centroid
;       definition so integer position is pixel center** W. Landsman March 2002
;       Allow use of the ST_DEV (not STDEV) keyword W. Landsman Nov. 2002
;-
        On_error,2

        if (N_params() LT 1 ) and $
            not (keyword_set( FWHM) or keyword_set(ST_DEV)) then begin
                print,'Syntax - psf = PSF_GAUSSIAN( parameters, NPIXEL = )'
                print, $
       'or       psf = PSF_GAUSSIAN( FWHM = ,ST_DEV = ,NPIXEL = ,[CENTROID = ])'
                return, -1
        endif

        sp = size( parameters )
        if sp[0] EQ 1 then begin               ;Vector supplied?
                ndim = 1
                factor = parameters[0]
                cntrd = parameters[1]
                st_dev = parameters[2] 
         endif  else  if (sp[0] GE 1) then begin    ;Ndimen x 3 array supplied?
                 ndim = sp[1]
                 factor = total( parameters[*,0] )/float( ndim )
                cntrd = parameters[*,1]
                st_dev = parameters[*,2]
           endif

        double = keyword_set(double)
        if double then idltype = 5 else idltype = 4
        if N_elements( ndim ) NE 1 then ndim=2
        ndim = ndim>1

        if N_elements( npix ) LE 0 then begin
                message,"must specify size of result with NPIX=",/INFO
                return,(-1)
          endif else if N_elements( npix ) LT ndim then $
                        npix = replicate( npix[0], ndim )

        if (N_elements( cntrd ) LT ndim) AND (N_elements( cntrd ) GT 0) then $
                        cntrd = replicate( cntrd[0], ndim )

        if N_elements( cntrd ) LE 0 then cntrd=(npix-1)/2. 
        if N_elements( fwhm ) GT 0 then begin 
               st_dev = fwhm/( 2.0d* sqrt( 2.0d* aLog(2.0d) ) )
               if not double then st_dev  = float(st_dev)
        endif 

        if N_elements( st_dev ) LE 0 then begin
                message,"must specify ST_DEV= or FWHM=",/INFO
                return,(-1)
          endif

        if N_elements( st_dev ) LT ndim then $
                        st_dev = replicate( st_dev[0], ndim )

        CASE ndim OF

        1: BEGIN
                x = findgen( npix[0] ) - cntrd[0]
                psf = gaussian( x, [1,0,st_dev] )
             END

        2: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]

                if N_elements( xy_corr ) EQ 1 then begin
                        sigfac = 1 / (2. * st_dev^2 )
                        y2 = sigfac[1] * y^2
                        x1 = sigfac[0] * x
                        yc = y * ( xy_corr/(st_dev[0]*st_dev[1]) )
                        for j=0,npix[1]-1 do begin
                                zz = x * (yc[j] + x1) + y2[j]
                                w = where( zz LT 86, nw )
                                if (nw GT 0) then psf[w,j] = exp( -zz[w] )
                          endfor
                  endif else begin
                        psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE=double )
                        psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE=double )
                        error = check_math(/print, MASK=32)
                        save_except = !EXCEPT & !EXCEPT = 0
                        for j=0,npix[1]-1 do psf[0,j] = psfx * psfy[j]
                        error = check_math(MASK=32)    ;Clear floating underflow
                        !EXCEPT = save_except  
                   endelse
             END

        3: BEGIN
                psf = make_array( DIM=npix[0:ndim-1], TYPE = idltype )
                x = make_array( npix[0], /INDEX, TYPE=idltype ) - cntrd[0]
                y = make_array( npix[1], /INDEX, TYPE=idltype ) - cntrd[1]
                z = make_array( npix[2], /INDEX, TYPE=idltype ) - cntrd[2]
                psfx = gaussian( x, [ 1, 0, st_dev[0] ], DOUBLE = double )
                psfy = gaussian( y, [ 1, 0, st_dev[1] ], DOUBLE = double)
                psfz = gaussian( z, [ 1, 0, st_dev[2] ], DOUBLE = double )
                error = check_math(MASK=32,/PRINT)
                save_except = !EXCEPT & !EXCEPT = 0
                for k=0,npix[2]-1 do begin
                    for j=0,npix[1]-1 do psf[0,j,k] = psfx * psfy[j] * psfz[k]
                 endfor
                 error = check_math(MASK=32)
                 !EXCEPT = save_except  
             END

        ENDCASE

        if keyword_set( normalize ) then return, psf/total( psf )

        if N_elements( factor ) EQ 1 then begin
                if (factor NE 1) then return,factor*psf else return,psf
           endif else return, psf
end
