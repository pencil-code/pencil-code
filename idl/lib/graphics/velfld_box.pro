PRO velfld_box,verts,lmax=lmax,face=face,thick=thick,zof=zof
;+
; NAME:        VELFLD
; PURPOSE:     Draw a velocity (flow) field with the length of arrows
;              following the field and proportional in length to the
;              field strength.  Arrows are composed of a number of small
;              segments, following the streamlines.
; CATEGORY:    2d/3d graphics
; CALLING SEQUENCE:
;              VELFLD, ax, ay, xmax, ymax
; INPUTS:
;        ax = Vector field X component at each point.  Must be a 2d array.
;        ay = Y component at each point, must be same dimensions as u.
;        verts = (3,8) array holding positions of 8 corners of box
; KEYWORD PARAMETERS:
;        Nvecs = # of vectors (arrows) to draw, if omitted 200 vectors
;                  are drawn.
;        Length = length of each step expressed as a fraction of the longest
;                  vector divided by the number of steps.  Default = 0.1.
;        Nsteps = # of shoots or line segments for each arrow, default = 10.
;        lmax = scaling factor for arrows, default=0.1
;        face = plane (in 3d) in which to draw field.
;               0 for z=1 plane, 4 for z=0 plane, 1 for x plane, 2 for y plane.
; OUTPUTS:
;        A velocity field graph.
; COMMON BLOCKS:
;        ps_com.
; SIDE EFFECTS:
;        A plot is drawn.
; RESTRICTIONS:
;        Position is fixed within plot window.
; PROCEDURE:
;        Nvecs random points within the (ax,ay) arrays are selected.
;        For each "shot" the field (as bilinearly interpolated) at each
;        point is followed using a vector of "length" length, tracing
;        a line with Nsteps segments.  An arrow is drawn at the end.
; MODIFICATION HISTORY:
;        Neal Hurlburt, April, 1988. Amended for 3d, April 1992, dpb.
;-

        common ps_com, dname,xthick,ythick,pthick,pfont,ps_file,xwsize,ywsize

if n_elements(face) le 0 then face=0
if n_elements(thick) le 0 then thick=1
if n_elements(zof) le 0 then zof=0

xw0=verts(0,0)
xw1=verts(0,1)
xwsize=xw1-xw0
yw0=verts(1,0)
yw1=verts(1,2)
ywsize=yw1-yw0
zw0=verts(2,0)
zw1=verts(2,4)
zwsize=zw1-zw0

if face eq 0 then begin
  xbox=[xw0,xw1,xw1,xw0,xw0]
  ybox=[yw0,yw0,yw1,yw1,yw0]
  zbox=replicate(zw1,5)
endif else if face eq 1 then begin
    xbox=[xw0,xw1,xw1,xw0,xw0]
    zbox=[zw0,zw0,zw1,zw1,zw0]
    ybox=replicate(yw0,5)
  endif else if face eq 2 then begin
      zbox=[zw0,zw1,zw1,zw0,zw0]
      ybox=[yw0,yw0,yw1,yw1,yw0]
      xbox=replicate(xw0,5)
  endif else if face eq 4 then begin
       xbox=[xw0,xw1,xw1,xw0,xw0]
       ybox=[yw0,yw0,yw1,yw1,yw0]
       zbox=replicate(zw0,5)-zof
    endif

plots,xbox,ybox,zbox,/t3d,/noclip,color=!p.color
RETURN
end
