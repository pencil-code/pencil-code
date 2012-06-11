;
;  Note: this routine carries a more general name than what it deserves
;  at the moment. At the moment this routine does a specific task while
;  looping through all horizontal processor layers (calculate bm and uxbtestm).
;  In future, the special features of this routine could be put into
;  include files (but this would obscure the debugging of this initial phase).
;
;  15-jun-08/axel: adapted from rallxy.pro
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   rallxy_zloop.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: nils/ adapted from rall.pro 20-nov-2003    
;;;  $Id$
;;;
;;;  Description:
;;;   Read data from all processors in a xy-plane (select plane by 
;;;   setting 'target_ipz' (default=0)) and combine them into one array
;;;   for each variable.
;;;   Overwrites nx, ny, etc, from start.pro, thus you need to run
;;;   start.pro after this if you want to continue working on data
;;;   from individual processors (e.g. if you want to run r.pro).
;
function param2
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
;  define output arrays based on full mz here:
;
nzgrid=mz-2*nghostz
njtest=ntestfield/3
zzz=fltarr(nzgrid/nprocz,nprocz)
bm=fltarr(nzgrid/nprocz,nprocz,3)
uxbtestm=fltarr(nzgrid/nprocz,nprocz,3,njtest)

;
;  must set mz=mzloc such that varcontent defines the correct arrays
;  allow this value to be read in interactively if equal to default.
;
mz=mzloc
;
;  begin z-loop
;
for target_ipz=0,nprocz-1 do begin

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
x = fltarr(mx) & y = fltarr(my) & z = fltarr(mz)
xloc = fltarr(mxloc) & yloc = fltarr(myloc) & zloc = fltarr(mzloc)

;
;  Read data
;
varcontent=pc_varcontent()
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
                                    + ' needs declaring in varcontent.pro', /INFO   
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
    print,'File contains: '+content
    print, FORMAT='(A,$)', "Reading: "
  endif

  if (ipz eq target_ipz) then begin
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
      i0z=0 & i1z=i0z+mzloc-1
      i0zloc=0 & i1zloc=mzloc-1
      ;
      x[i0x:i1x] = xloc[i0xloc:i1xloc]
      y[i0y:i1y] = yloc[i0yloc:i1yloc]
      z[i0z:i1z] = zloc[i0zloc:i1zloc]

      for iv=0L,totalvars-1L do begin
          cmd =   varcontent[iv].idlvar $
            + "[i0x:i1x,i0y:i1y,i0z:i1z,*]=" $
            + varcontent[iv].idlvarloc $
            +"[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*]"         
          if (execute(cmd) ne 1) then $
            message, 'Error combining data for ' + varcontent[iv].variable         
          ; For vector quantities skip the required number of elements
          iv=iv+varcontent[iv].skip
      endfor
  endif
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
l1=3 & l2=mx-4
m1=3 & m2=my-4
n1=3 & n2=mz-4
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
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

;
;  bottom of z-loop
;
print,'bottom of z-loop'
;
;  mean field
;
zzz(*,target_ipz)=z(n1:n2)
bm(*,target_ipz,*)=haverv((curl(aa))(l1:l2,m1:m2,n1:n2,*))
aatest=reform(aatest,mx,my,mz,3,njtest)
for jtest=0,njtest-1 do begin
  uxbtestm(*,target_ipz,*,jtest)=haverv((reform(cross(uu,curl(aatest(*,*,*,*,jtest)))))(l1:l2,m1:m2,n1:n2,*))
endfor
help,target_ipz,bm,uxbtestm
;
endfor
zzz=reform(zzz,nzgrid)
bm=reform(bm,nzgrid,3)
uxbtestm=reform(uxbtestm,nzgrid,3,njtest)
;
;  save result
;
save,file=varfile+'.uxbtestm',bm,uxbtestm,t,zzz
;
end

; End of file
