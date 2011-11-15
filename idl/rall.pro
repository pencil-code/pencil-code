;;;;;;;;;;;;;;;;;;;;
;;;   rall.pro   ;;;
;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   09-Sep-2001
;;;  $Id$
;;;
;;;  Description:
;;;   Read data from all processors and combine them into one array
;;;   for each variable.
;;;   Overwrites nx, ny, etc, from start.pro, thus you need to run
;;;   start.pro after this if you want to continue working on data
;;;   from individual processors (e.g. if you want to run r.pro).

function param2
COMPILE_OPT HIDDEN 
; Dummy to keep IDL from complaining. The real param2() routine will be
; compiled below
end

;
;  need to run start first: check whether this has been done
;
if (n_elements(started) le 0) then begin
  message, "You need to run start.pro first: use `.rnew start'"
endif

;
;  Read startup parameters
;
@readstartpars

;
;  read global sizes
;
nprocx=0L & nprocy=0L & nprocz=0L
close,1
openr,1,datatopdir+'/'+dimfile
readf,1,mx,my,mz,nvar,naux
readf,1,prec
readf,1,nghostx,nghosty,nghostz
readf,1,nprocx,nprocy,nprocz
close,1
;
ncpus = nprocx*nprocy*nprocz
;
;  read local sizes
;
datadir=datatopdir+'/proc0'
mxloc=0L & myloc=0L & mzloc=0L
;
close,1
openr,1,datadir+'/'+dimfile
readf,1,mxloc,myloc,mzloc
close,1
;
nxloc=mxloc-2*nghostx
nyloc=myloc-2*nghosty
nzloc=mzloc-2*nghostz
;
;  read data
;
default, dataopdir, 'data'
default, varfile, 'var.dat'
default, ldustvelocity, 0
default, ldustdensity, 0
;
x=fltarr(mx)*ONE & y=fltarr(my)*ONE & z=fltarr(mz)*ONE
xloc=fltarr(mxloc)*ONE & yloc=fltarr(myloc)*ONE & zloc=fltarr(mzloc)*ONE
;
dx_1=fltarr(mx)*zero & dy_1=fltarr(my)*zero & dz_1=fltarr(mz)*zero
dx_1_loc=fltarr(mxloc)*zero & dy_1_loc=fltarr(myloc)*zero & dz_1_loc=fltarr(mzloc)*zero
dx_tilde=fltarr(mx)*zero & dy_tilde=fltarr(my)*zero & dz_tilde=fltarr(mz)*zero
dx_tilde_loc=fltarr(mxloc)*zero & dy_tilde_loc=fltarr(myloc)*zero & dz_tilde_loc=fltarr(mzloc)*zero

;
;  Read data
;
varcontent=pc_varcontent(QUIET=quiet)
totalvars=(size(varcontent))[1]

; Prepare for read
readstring=''
content=''
for i=0L,totalvars-1L do begin
  readstring = readstring + ',' + varcontent[i].idlvarloc
  content    = content + ', ' + varcontent[i].variable
  ; Initialise variable
  if (varcontent[i].variable eq 'UNKNOWN') then $
           message, 'Unknown variable at position ' + str(i)  $
                                    + ' needs declaring in pc_varcontent.pro', /INFO   
  if (execute(varcontent[i].idlvar+'='+varcontent[i].idlinit,0) ne 1) then $
           message, 'Error initialising ' + varcontent[i].variable $
                                    +' - '+ varcontent[i].idlvar, /INFO
  if (execute(varcontent[i].idlvarloc+'='+varcontent[i].idlinitloc,0) ne 1) then $
           message, 'Error initialising ' + varcontent[i].variable $
                                    +' - '+ varcontent[i].idlvarloc, /INFO
;If it's a vector quantity skip the required number of elements
  i=i+varcontent[i].skip
end

content = strmid(content,2)
;
for i=0,ncpus-1 do begin        ; read data from individual files
  tag='proc'+str(i)
  datadir=datatopdir+'/proc'+strtrim(i,2)
  ; read processor position
  dummy=''
  ipx=0L &ipy=0L &ipz=0L
  close,1
  openr,1,datadir+'/'+dimfile
  readf,1, dummy
  readf,1, dummy
  readf,1, dummy
  readf,1, ipx,ipy,ipz
  ; read data
  if ((i eq 0) and (quiet le 2)) then begin
    print,'File '+varfile+' contains: ', content
    print, FORMAT='(A,$)', "Reading: "
  endif

  close,1
  openr,1, datadir+'/'+varfile, /F77
  if (quiet le 2) then print, FORMAT='(A," ",$)', tag
  if (execute('readu,1'+readstring) ne 1) then $
      message, 'Error reading: ' + 'readu,1'+readstring
    ;
    ;  read deltay in case of shear
    ;
  if (lshear) then begin
    readu,1, t, xloc, yloc, zloc, dx, dy, dz, deltay
  end else begin
    readu,1, t, xloc, yloc, zloc
  end

  close,1
  gridfile=datadir+'/'+'grid.dat'
  if (any(lequidist eq 0)) then begin
    openr,1,gridfile,/F77
    point_lun,1,pos
    readu,1, dx_1_loc,     dy_1_loc,     dz_1_loc
    readu,1, dx_tilde_loc, dy_tilde_loc, dz_tilde_loc
    close,1
  endif else begin
    ;
    ;  Ensure we don't use these values
    ;
    dx_1_loc = dx_1_loc*!values.f_nan
    dy_1_loc = dy_1_loc*!values.f_nan
    dz_1_loc = dz_1_loc*!values.f_nan
    dx_tilde_loc = dx_tilde_loc*!values.f_nan
    dy_tilde_loc = dy_tilde_loc*!values.f_nan
    dz_tilde_loc = dz_tilde_loc*!values.f_nan
  endelse

  ;
  ;  Don't overwrite ghost zones of processor to the left (and
  ;  accordingly in y and z direction makes a difference on the
  ;  diagonals)
  ;
  if (ipx eq 0) then begin
    i0x=ipx*nxloc & i1x=i0x+mxloc-1
    i0xloc=0 & i1xloc=mxloc-1
  endif else begin
    i0x=ipx*nxloc+nghostx & i1x=i0x+mxloc-1-nghostx
    i0xloc=nghostx & i1xloc=mxloc-1
  endelse
  ;
  if (ipy eq 0) then begin
    i0y=ipy*nyloc & i1y=i0y+myloc-1
    i0yloc=0 & i1yloc=myloc-1
  endif else begin
    i0y=ipy*nyloc+nghosty & i1y=i0y+myloc-1-nghosty
    i0yloc=nghosty & i1yloc=myloc-1
  endelse
  ;
  if (ipz eq 0) then begin
    i0z=ipz*nzloc & i1z=i0z+mzloc-1
    i0zloc=0 & i1zloc=mzloc-1
  endif else begin
    i0z=ipz*nzloc+nghostz & i1z=i0z+mzloc-1-nghostz
    i0zloc=nghostz & i1zloc=mzloc-1
  endelse
  ;
  x[i0x:i1x] = xloc[i0xloc:i1xloc]
  y[i0y:i1y] = yloc[i0yloc:i1yloc]
  z[i0z:i1z] = zloc[i0zloc:i1zloc]
  ;
  dx_1[i0x:i1x] = dx_1_loc[i0xloc:i1xloc]
  dy_1[i0y:i1y] = dy_1_loc[i0yloc:i1yloc]
  dz_1[i0z:i1z] = dz_1_loc[i0zloc:i1zloc]
  ;
  dx_tilde[i0x:i1x] = dx_tilde_loc[i0xloc:i1xloc]
  dy_tilde[i0y:i1y] = dy_tilde_loc[i0yloc:i1yloc]
  dz_tilde[i0z:i1z] = dz_tilde_loc[i0zloc:i1zloc]

  for iv=0L,totalvars-1L do begin
    if (varcontent[iv].variable eq 'UNKNOWN') then continue
    cmd =   varcontent[iv].idlvar $
          + "[i0x:i1x,i0y:i1y,i0z:i1z,*]=" $
          + varcontent[iv].idlvarloc $
          +"[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]"         
    if (execute(cmd) ne 1) then $
        message, 'Error combining data for ' + varcontent[iv].variable         
  ; For vector quantities skip the required number of elements
    iv=iv+varcontent[iv].skip
  endfor

endfor
if (quiet le 2) then print

;
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
rr = sqrt(xx^2+yy^2+zz^2)

; reset datadir to more reasonable default
datadir=datatopdir+'/proc0'
;
;  reset boundary values and nx,ny,nz for (full) physical domain (not
;  sub-domain)
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
l1=3 & l2=mx-4 & l12=l1+indgen(nx)
m1=3 & m2=my-4 & m12=m1+indgen(ny)
n1=3 & n2=mz-4 & n12=n1+indgen(nz)
;
;  fix z3=ztop which was local top in start.pro
;
if (lgravz) then begin
  ztop=z[n2] & z3=ztop
endif
;
;  Summarize data
;
@varcontent_stats
;
;  free memory
;
undefine, uu_loc
undefine, lnrho_loc
undefine, ss_loc
undefine, aa_loc
undefine, lncc_loc
;
read_all = 1                    ; marker for r.pro

end

; End of file
