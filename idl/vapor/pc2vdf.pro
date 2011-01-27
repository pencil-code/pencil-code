; pc2vdf, program to create vdf files
; vdf is the format required by vapor for 3d visualization
; last modified 14/07/2010: Guerrero

pro pc2vdf,tinit=tinit,tend=tend,rhof=rhof,sph=sph,variables=variables,beq=beq

; KEYWORD TINIT: INITIAL VARFILE
; KEYWORD TEND: FINAL VARFILE
; KEYWORD RHOF: COMPUTES DENSITY PERTURBATIONS INSTEAD OF DENSITY 
; KEYWORD SPH: READ PC VARFILES IN SPHERICAL COORDS.
; KEYWORD BEQ: DEFINES A NORMALIZATION VALUE FOR THE MAGNETIC FIELD
;e.g.:pc2vdf,tinit=2,tend=2,variebles=['uu','bb','lnrho'],/rhof,beq=0.26

default, tinit, '0'
default, tend, '1'
default,variables,['uu','lnrho']
default,beq,1.


k = strmatch(variables,'bb')
mag = where(k EQ 1)

if (mag GE 0) then begin
   if (keyword_set(sph)) then begin
      varnames=['rho_sph','vr','vt','vp','br','bt','bp','Eb_sph']
   endif else begin
      varnames=['rho','vx','vy','vz','bx','by','bz','Eb']
   endelse
endif else begin
   if (keyword_set(sph)) then begin
      varnames=['rho_sph','vr','vt','vp']
   endif else begin
      varnames=['rho','vx','vy','vz']
   endelse
endelse
;
cd, current=c
readpath=c
spawn,'mkdir VAPOR'
writepath=c+'/VAPOR'
cd, writepath
if (tinit EQ 0) then begin
  varfile1 = 'var.dat'
endif else begin
  varfile1 = 'VAR'+strtrim(tinit,2)
endelse
;
cd, readpath 
      ;
      OPENR,1, readpath+'data/proc0/'+varfile1, ERROR = err
      close,1
      pc_read_grid,o=g,/trimxyz,/quiet
;
cd, writepath
;	$Id: WriteVDF.pro,v 1.6 2008/09/03 20:51:55 clynejp Exp $
;
;	This example shows how to create a data set containing 
;	a single variable and at a single time step
;
; 
;	Dimensions of the data volumes - all volumes in a data set must be 
;	of the same dimension
;
cd, writepath
nx = (size(g.x))[1]   ;   (size(test.uu))[1]-6
ny = (size(g.y))[1]   ;   (size(test.uu))[2]-6
nz = (size(g.z))[1]   ;  (size(test.uu))[3]-6

if (keyword_set(sph)) then dim = [nz,ny,nx] else dim = [nx,ny,nz]
if (keyword_set(rhof)) then rho_f = fltarr(nx,ny,nz)
;
;
; The number of coarsened approximations to create. A value of 0
; indicates that the data should not be transformed, i.e. the
; data will only exist at it's native resolution. A value of 1
; indicates that a single coarsening should be applied, and so on.
;
num_levels = 1
;	Create a new VDF metadata object of the indicated dimension and 
;	transform
;	level. vdf_create() returns a handle for future operations on 
;	the metadata object.
;
mfd = vdf_create(dim,num_levels)

;	Set the maximum number of timesteps in the data set. Note, valid data 
;	set may contain less than the maximum number of time steps, 
;	but not more
;
timesteps = tend-tinit+1
vdf_setnumtimesteps, mfd,timesteps
;
;	Set the names of the variables the data set will contain. In this case,
;
numvar = size(varnames)
numvar = numvar(1)
print,'number of variables to write =', numvar
;
vdf_setvarnames, mfd, varnames
;
; Set the extents of the volume in user-defined
; physical coordinates. Note, as the aspect ratio of
; the user-defined coordinates do not match that of
; the volume resolution (512x512x512), the volume
; will be stretched when rendered. I.e. the spacing
; between the Z coordinate samples is defined to
; be twice that of X or Y.
;
if (keyword_set(sph)) then begin 
   vdf_setcoordtype,mfd,'spherical'
; order for spherical coordinates phi(z),theta(y),r(x) 
   extents = [-180+(180.*g.z[0]/!pi),  -90+(180.*g.y[0]/!pi), g.x[0], $ 
              -180+(180.*g.z[nz-1]/!pi),-90+(180.*g.y[ny-1]/!pi), g.x[nx-1]]
   print,extents
;new variables
   rho_sph = fltarr(nz,ny,nx)
   vr = fltarr(nz,ny,nx)
   vt = fltarr(nz,ny,nx)
   vp = fltarr(nz,ny,nx)
   print,'size of v = ', size(vr)
   br = fltarr(nz,ny,nx)
   bt = fltarr(nz,ny,nx)
   bp = fltarr(nz,ny,nx)
   Eb_sph = fltarr(nz,ny,nx)
   last_dim = nx-1
endif else begin
   extents = [g.x[0], g.y[0], g.z[0], g.x[nx-1], g.y[ny-1], g.z[nz-1]]
   last_dim = nz-1
endelse
;
vdf_setextents, mfd, extents
;
;	Store the metadata object in a file for subsequent use
;
vdffile = 'work.vdf'
vdf_write, mfd, vdffile
;
;	Destroy the metadata object. We're done with it.
;
vdf_destroy, mfd
;
;	At this point we've defined all of the metadata associated with 
;	the test data set. Now we're ready to begin populating the data
;	set with actual data. 
;
;
;	Create a "buffered write" data transformation object. The data
;	transformation object will permit us to write (transform) raw
;	data into the data set. The metadata for the data volumes is
;	obtained from the metadata file we created previously. I.e.
;	'vdffile' must contain the path to a previously created .vdf
;	file. The 'vdc_bufwritecreate' returns a handle, 'dfd', for 
;	subsequent ;	operations.

; Prepare the data set for writing. We need to identify the time step
; and the name of the variable that we wish to store
;
;
itt = 0
for it = tinit,tend do begin
vars = variables
fac=1
 varfile = 'VAR'+strtrim(string(format="(i2)",it),1)  
 OPENR,1, readpath+'data/proc0/'+varfile, ERROR = err
      close,1
  cd, readpath
     print,'Reading: ',varfile
     print,'vars = ',vars
     pc_read_var,varfile=varfile,obj=test,/trimall,variables=vars,/magic,/quiet

     rho=exp(test.lnrho)
;
    if (keyword_set(rhof)) then begin
     rhoy = total(rho[*,*,*],2)/ny
     rhoxy = total(rhoy[*,*],1)/nx
     for k = 0,nz-1 do rho_f[*,*,k] = rho[*,*,k] - rhoxy[k]
     rho = rho_f
    endif
;
    ss=[0,ny/2,0]
    rho = shift(rho,ss)
    v1 =shift(test.uu[*,*,*,0],ss)
    v2 =shift(test.uu[*,*,*,1],ss)  
    v3 =shift(test.uu[*,*,*,2],ss)
  if (mag GT 0) then begin
      b1 = shift(test.bb[*,*,*,0],ss)/beq
      b2 = shift(test.bb[*,*,*,1],ss)/beq
      b3 = shift(test.bb[*,*,*,2],ss)/beq
      Eb = alog(b1^2 + b2^2 + b2^2)
  endif
; reordering
    if (keyword_set(sph)) then begin
       for l = 0, nx-1 do begin
          for m =0, ny-1 do begin
             for n =0, nz-1 do begin
                rho_sph(n,m,l) = rho(l,m,n)
                vr(n,m,l) = v1(l,m,n)
                vt(n,m,l) = v2(l,m,n)
                vp(n,m,l) = v3(l,m,n) 
                if (mag GT 0) then begin
                   br(n,m,l) = b1(l,m,n)
                   bt(n,m,l) = b2(l,m,n)
                   bp(n,m,l) = b3(l,m,n)
                   Eb_sph(n,m,l) = Eb(l,m,n)
                endif
             endfor
          endfor
       endfor
    endif
;    
    print,'changing to previous directory'
    cd, writepath
    cd, current=c
    print,c

 for nvar = 0, numvar-1 do begin
      dfd = vdc_bufwritecreate(vdffile)
      vdc_openvarwrite, dfd, itt, varnames(nvar), -1
; Write (transform) the volume to the data set one slice at a time

      if (keyword_set(sph)) then begin
         case nvar of
            0: tmp_var = rho_sph
            1: tmp_var = vr
            2: tmp_var = vt
            3: tmp_var = vp
            4: tmp_var = br
            5: tmp_var = bt
            6: tmp_var = bp
            7: tmp_var = Eb_sph
         endcase
      endif else begin
         case nvar of
            0: tmp_var = rho
            1: tmp_var = v1
            2: tmp_var = v2
            3: tmp_var = v3
            4: tmp_var = b1
            5: tmp_var = b2
            6: tmp_var = b3
            7: tmp_var = Eb
         endcase
      endelse
      for z = 0, last_dim do begin
         vdc_bufwriteslice, dfd, float(tmp_var[*,*,z])
      endfor
      
; Close the currently opened variable/time-step. We're done writing
; to it
      vdc_closevar, dfd
   endfor
  itt++

endfor
;	Destroy the "buffered write" data transformation object. 
;	We're done with it.
cd, readpath
vdc_bufwritedestroy,dfd
end
