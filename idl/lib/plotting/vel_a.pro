;;;;;;;;;;;;;;;;;;;;;;
;;;   vel_a.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;
;;;  Author: axel, wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   8-Jul-1999
;;;  Version: vel_a 1.2
;;;  CVS $Revision: 1.6 $
;;;  Based on: vel.pro,v 1.4 1997/01/15 03:11:50 idl Exp,
;;;  Description: A clone of IDL's vel allowing for
;;;    a) X and Y arguments and the corresponding {X,Y}RANGE,
;;;       {X,Y}STYLE (like velovect, contour, ..)
;;;    b) a SEED argument specifying the random seed
;;;    c) the keywords NOERASE and OVERPLOT
;;;    If called without X and Y arguments, the (slightly bizarre)
;;;    behaviour of IDL's `vel' is reproduced.
;;;    All unknown keywords are passed over to PLOTS (which might not
;;;    be very clever, but makes vel_a,..,color=120 work)

function vel_mybi_a,a,x,y
on_error,2                      ;Return to caller if an error occurs
sizea=size(a)
nx=sizea[1]
i=long(x)+nx*long(y)
q=y-long(y)
p=x-long(x)
q1 = 1.-q
p1 = 1.-p

; Weighting factors were wrong for a(i+1) & a(i+nx), switched them.

aint=p1*q1*a[i] + p*q1*a[i+1] + q*p1*a[i+nx] + p*q*a[i+nx+1]
return,aint
end

PRO ARRHEAD_A,X

ON_ERROR,2                      ;Return to caller if an error occurs
theta = 30 * !radeg
TANT = TAN(THETA)
NP=3.0
SCAL=8.

SX=SIZE(X)
N=SX[2]

boundp = (x[*,N-4,0] eq 0) or (x[*,N-4,0] eq 1) or (x[*,N-4,1] eq 0) or (x[*,N-4,1] eq 1) 
bigl=sqrt((X[*,N-4,0]-X[*,N-5,0])^2+(X[*,N-4,1]-X[*,N-5,1])^2)
wbigl=where((bigl ne 0.0) and (not boundp))
wnbigl=where((bigl eq 0.0) or boundp, count)

LL  = SCAL*TANT*bigl[wbigl]/NP

DX = LL*(X[wbigl,N-4,1]-X[wbigl,N-5,1])/BIGL[wbigl]
DY = LL*(X[wbigl,N-4,0]-X[wbigl,N-5,0])/BIGL[wbigl]

XM = X[wbigl,N-4,0]-(SCAL-1)*(X[wbigl,N-4,0]-X[wbigl,N-5,0])/NP
YM = X[wbigl,N-4,1]-(SCAL-1)*(X[wbigl,N-4,1]-X[wbigl,N-5,1])/NP

X[wbigl,N-3,0] = XM-DX
X[wbigl,N-2,0] = X[wbigl,N-4,0]
X[wbigl,N-1,0] = XM+DX

X[wbigl,N-3,1] = YM+DY
X[wbigl,N-2,1] = X[wbigl,N-4,1]
X[wbigl,N-1,1] = YM-DY

if count ge 1 then begin  ;No head for 0 length
	X[wnbigl,N-3,0] = x[wnbigl,n-4,0]
	X[wnbigl,N-2,0] = X[wnbigl,n-4,0]
	X[wnbigl,N-1,0] = X[wnbigl,n-4,0]
	
	X[wnbigl,N-3,1] = X[wnbigl,N-4,1]
	X[wnbigl,N-2,1] = X[wnbigl,N-4,1]
	X[wnbigl,N-1,1] = X[wnbigl,N-4,1]
	endif

return
END

function arrows_a,u,v,n,length, $
                nsteps=nsteps,seed=seed, $
                dx=dx,dy=dy
on_error,2                      ;Return to caller if an error occurs
su=size(u)
nx=su[1]
ny=su[2]

if (n_elements(dx) le 0) then dx = 1.
if (not keyword_set(dy)) then dy = 1.
lmax=sqrt(max(u^2+v^2))		;Max vector length
lthx=dy*length/lmax/nsteps
lthy=dx*length/lmax/nsteps
xt=randomu(seed,n)		;Starting position
yt=randomu(seed,n)
x=fltarr(n,nsteps+3,2)
x[*,0,0]=xt
x[*,0,1]=yt
for i=1L,nsteps-1 do begin
  xt[*]=(nx-1)*x[*,i-1,0]
  yt[*]=(ny-1)*x[*,i-1,1]
  ut=vel_mybi_a(u,xt,yt)
  vt=vel_mybi_a(v,xt,yt)
  x[*,i,0]=x[*,i-1,0]+ut*lthx
  x[*,i,1]=x[*,i-1,1]+vt*lthy
  ;; Reset the points located outside the drawing window:
  bad=where(x[*,i,0] lt 0)
  if (bad[0] ge 0) then begin
    ;; Interpolate to get the correct point on the boundary
    x[bad,i,1] = x[bad,i-1,1] - $
        x[bad,i-1,0]*(x[bad,i,1]-x[bad,i-1,1])/(x[bad,i,0]-x[bad,i-1,0])
    x[bad,i,0] = 0
  endif
  bad=where(x[*,i,0] gt 1)
  if (bad[0] ge 0) then begin
    ;; Interpolate to get the correct point on the boundary
    x[bad,i,1] = x[bad,i-1,1] + $
        (1-x[bad,i-1,0])*(x[bad,i,1]-x[bad,i-1,1])/(x[bad,i,0]-x[bad,i-1,0])
    x[bad,i,0] = 1
  endif
  bad=where(x[*,i,1] lt 0)
  if (bad[0] ge 0) then begin
    ;; Interpolate to get the correct point on the boundary
    x[bad,i,0] = x[bad,i-1,0] - $
        x[bad,i-1,1]*(x[bad,i,0]-x[bad,i-1,0])/(x[bad,i,1]-x[bad,i-1,1])
    x[bad,i,1] = 0
  endif
  bad=where(x[*,i,1] gt 1)
  if (bad[0] ge 0) then begin
    ;; Interpolate to get the correct point on the boundary
    x[bad,i,0] = x[bad,i-1,0] + $
        (1-x[bad,i-1,1])*(x[bad,i,0]-x[bad,i-1,0])/(x[bad,i,1]-x[bad,i-1,1])
    x[bad,i,1] = 1
  endif
end
ARRHEAD_A,X
return,x
end


PRO VEL_A,U,W,xx,yy, $
          LENGTH=length, XMAX=xmax, $
          NVECS=nvecs, NSTEPS=nsteps, $
          TITLE=title, OVERPLOT=overplot, $
          SEED=seed, NOERASE=noerase, $
          COLOR=color, $
          _EXTRA=extra
;+
; NAME:
;	VEL
;
; PURPOSE:
;	Draw a velocity (flow) field with arrows following the field 
;	proportional in length to the field strength.  Arrows are composed 
;	of a number of small segments that follow the streamlines.
;
; CATEGORY:
;	Graphics, two-dimensional.
;
; CALLING SEQUENCE:
;	VEL, U, V
;
; INPUTS:
;	U:	The X component at each point of the vector field.  This 
;		parameter must be a 2D array.
;
;	V:	The Y component at each point of the vector field.  This 
;		parameter must have the same dimensions as U.
;
; KEYWORD PARAMETERS:
;	NVECS:	The number of vectors (arrows) to draw.  If this keyword is
;		omitted, 200 vectors are drawn.
;
;	XMAX:	X axis size as a fraction of Y axis size.  The default is 1.0.
;		This argument is ignored when !p.multi is set.
;
;	LENGTH:	The length of each arrow line segment expressed as a fraction 
;		of the longest vector divided by the number of steps.  The 
;		default is 0.1.
;
;	NSTEPS:	The number of shoots or line segments for each arrow.  The
;		default is 10.
;
;	TITLE:	A string containing the title for the plot.
;	
; OUTPUTS:
;	No explicit outputs.  A velocity field graph is drawn on the current
;	graphics device.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	A plot is drawn on the current graphics device.
;
; RESTRICTIONS:
;	none
;
; PROCEDURE:
;	NVECS random points within the (u,v) arrays are selected.
;	For each "shot" the field (as bilinearly interpolated) at each
;	point is followed using a vector of LENGTH length, tracing
;	a line with NSTEPS segments.  An arrow head is drawn at the end.
;
; MODIFICATION HISTORY:
;	Neal Hurlburt, April, 1988.
;	12/2/92	- modified to handle !p.multi (jiy-RSI)
;       7/12/94 HJM - Fixed error in weighting factors in function
;                     vel_mybi() which produced incorrect velocity vectors.
;       sometimes - ab added xx and yy arguments and introduces
;                   keywords OVERPLOT, NOERASE and SEED
;       8/7/99  - wd improved the xx and yy argument handling and fixed
;                 a bug for different x- and y-axis ranges.
;       31/8/99 - wd switched loop indices to type long integer
;-

on_error,2                      ;Return to caller if an error occurs
if n_elements(Nvecs) le 0 then nvecs=200
if n_elements(nsteps) le 0 then nsteps = 10
if n_elements(length) le 0 then length=.1
if n_elements(title) le 0 then title='Velocity Field'
if n_elements(seed) le 0 then seed=3.11
if n_elements(color) le 0 then begin
  if !d.name eq 'PS' then color=0 else color=255
endif
;
if n_elements(xx) le 0 then begin ; Mimic the old behaviour
  X=ARROWS_A(U,W,Nvecs,LENGTH, nsteps=nsteps,seed=seed)
  if (!p.multi[1] eq 0 and !p.multi[1] eq 0) then begin
    if (n_elements(xmax) eq 0) then xmax = 1.0
    IF XMAX GT 1. THEN position=[0.20,(0.5-0.30/XMAX),0.90,(0.5+0.40/XMAX)]$
    else position=[(0.5-0.30*XMAX),0.20,(0.5+0.40*XMAX),0.90]
    plot,[0,1,1,0,0],[0,0,1,1,0],title=title,pos=position
  endif else begin
    plot,[0,1,1,0,0],[0,0,1,1,0],title=title
  endelse
  FOR I=0L,Nvecs-1 DO PLOTS,X[I,*,0],X[I,*,1],_EXTRA=extra
endif else begin
  x1=min(xx) & dxx=max(xx)-x1
  y1=min(yy) & dyy=max(yy)-y1
  if keyword_set(noerase) then begin
    if not keyword_set(overplot) then plot,x1+dxx*[0,1,1,0,0],y1+dyy*[0,0,1,1,0],$
        title=title,/nodata,/noerase,color=color
  end else begin
    if not keyword_set(overplot) then plot,x1+dxx*[0,1,1,0,0],y1+dyy*[0,0,1,1,0],$
        title=title,/nodata,color=color
  end
  X=ARROWS_A(U,W,Nvecs,LENGTH, nsteps=nsteps,seed=seed,dx=dxx,dy=dyy)
  FOR I=0L,Nvecs-1 DO PLOTS,x1+dxx*X[I,*,0],y1+dyy*X[I,*,1],noclip=0,color=color
endelse
RETURN
end


PRO VEL,ux,uy,xx,xy, $
          _EXTRA=extra
;; Redefine VEL as a wrapper to the above defined function `vel_a'. It
;; is better to call `vel_a' in an IDL program, since then you are
;; sure what you get.
  vel_a, ux,uy,xx,xy, _EXTRA=extra
  RETURN
end
