; pc2vdf, program to create vdf files
; vdf is the format required by vapor for 3d visualization
; last modified 14/07/2010: Guerrero

pro pc2vdf,tinit=tinit,tend=tend,rhof=rhof,sph=sph,variables=variables

; KEYWORD TINIT: INITIAL VARFILE
; KEYWORD TEND: FINAL VARFILE
; KEYWORD RHOF: COMPUTES DENSITY PERTURBATIONS INSTEAD OF DENSITY 
; KEYWORD SPH: READ PC VARFILES IN SPHERICAL COORDS.

default, tinit, '0'
default, tend, '1'
default,variables,['uu','lnrho']

k = strmatch(variables,'aa')
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
  
     print,'Reading:   ',varfile1
     pc_read_var,variables=['uu'],varfile=varfile1,obj=test,/quiet
     pc_read_grid,o=g,/trimxyz
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
nx = (size(test.uu))[1]-6
ny = (size(test.uu))[2]-6
nz = (size(test.uu))[3]-6
l1=3
m1=3
n1=3
l2=l1+nx-1
m2=m1+ny-1
n2=n1+nz-1

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
;   extents = [-180,-60.,0.6,180,60.,1.03]
   print,extents
;new variables
   rho_sph = fltarr(nz,ny,nx)
   vr = fltarr(nz,ny,nx)
   vt = fltarr(nz,ny,nx)
   vp = fltarr(nz,ny,nx)
   print,size(vr)
   br = fltarr(nz,ny,nx)
   bt = fltarr(nz,ny,nx)
   bp = fltarr(nz,ny,nx)
   Eb_sph = fltarr(nz,ny,nx)
   last_dim = nx-1
endif else begin
   extents = [g.x[0], g.y[0], g.z[0], g.x[nx-1], g.y[ny-1], g.z[nz-1]]
   lastdim = nz-1
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
fac=1
 varfile = 'VAR'+strtrim(string(format="(i3)",it),1)  
 OPENR,1, readpath+'data/proc0/'+varfile, ERROR = err
      close,1
  cd, readpath
     print,'Reading:   ',varfile
     pc_read_var,variables=variables,varfile=varfile,obj=test,/quiet   
     rho=exp(test.lnrho[l1:l2,m1:m2,n1:n2])
;
    if (keyword_set(rhof)) then begin
     rhoy = total(rho[*,*,*],2)/ny
     rhoxy = total(rhoy[*,*],1)/nx
     for k = 0,nz-1 do rho_f[*,*,k] = rho[*,*,k] - rhoxy[k]
     rho = rho_f
    endif
;
    v1 =test.uu[l1:l2,m1:m2,n1:n2,0]
    v2 =test.uu[l1:l2,m1:m2,n1:n2,1]  
    v3 =test.uu[l1:l2,m1:m2,n1:n2,2]
    if (mag GT 0) then begin
       bb = curl(test.aa)
       b1 =bb[l1:l2,m1:m2,n1:n2,0]
       b2 =bb[l1:l2,m1:m2,n1:n2,1]
       b3 =bb[l1:l2,m1:m2,n1:n2,2]
       Eb = b1^2 + b2^2 + b2^2 
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
;            4: tmp_var = br
;            5: tmp_var = bt
;            6: tmp_var = bp
;            7: tmp_var = Eb_sph
         endcase
      endif else begin
         case nvar of
            0: tmp_var = rho
            1: tmp_var = v1
            2: tmp_var = v2
            3: tmp_var = v3
;            4: tmp_var = b1
;            5: tmp_var = b2
;            6: tmp_var = b3
;            7: tmp_var = Eb
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
