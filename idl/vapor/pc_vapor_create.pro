pro pc_vapor_create,tinit,tend

readpath='./'
writepath='./Vapor'

cd, writepath, CURRENT=readpath
varfile1 = 'VAR'+strtrim(tinit,2)

;filename=run+varfile

cd, readpath, CURRENT=writepath
      ;
      OPENR,1, readpath+'data/proc0/'+varfile1, ERROR = err
      close,1
  
     print,'Reading:   ',varfile1
     pc_read_var,variables=['uu'],varfile=varfile1,obj=test,/quiet

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
dim = [nx,ny,nz]

;
;
; The number of coarsened approximations to create. A value of 0
; indicates that the data should not be transformed, i.e. the
; data will only exist at it's native resolution. A value of 1
; indicates that a single coarsening should be applied, and so on.
;
num_levels = 2

;
;	Create a new VDF metadata object of the indicated dimension and 
;	transform
;	level. vdf_create() returns a handle for future operations on 
;	the metadata object.
;
mfd = vdf_create(dim,num_levels)

;
;	Set the maximum number of timesteps in the data set. Note, valid data 
;	set may contain less than the maximum number of time steps, 
;	but not more
;
timesteps = tend-tinit
vdf_setnumtimesteps, mfd,timesteps

;
;	Set the names of the variables the data set will contain. In this case,
;	only a single variable will be present, "ml"
;
varnames = ['rho','vx','vy','vz','ox','oy','oz']
numvar = size(varnames)
numvar = numvar(1)
print,numvar

vdf_setvarnames, mfd, varnames
;
; Set the extents of the volume in user-defined
; physical coordinates. Note, as the aspect ratio of
; the user-defined coordinates do not match that of
; the volume resolution (512x512x512), the volume
; will be stretched when rendered. I.e. the spacing
; between the Z coordinate samples is defined to
; be twice that of X or Y.

extents = [0, 0, 0, 1, 1, 1]
vdf_setextents, mfd, extents

;
;	Store the metadata object in a file for subsequent use
;
vdffile = 'shear_2.vdf'
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
;


; Prepare the data set for writing. We need to identify the time step
; and the name of the variable that we wish to store
;

for it = 0, tend-tinit-1 do begin
fac=1
itt=it+tinit
 varfile = 'VAR'+strtrim(itt*fac,2)
  OPENR,1, readpath+'data/proc0/'+varfile, ERROR = err
      close,1
  cd, readpath
     print,'Reading:   ',varfile
     pc_read_var,variables=['uu','lnrho'],varfile=varfile,obj=test,/quiet   

    rho=exp(test.lnrho[l1:l2,m1:m2,n1:n2])
    v1 =test.uu[l1:l2,m1:m2,n1:n2,0]
    v2 =test.uu[l1:l2,m1:m2,n1:n2,1]  
    v3 =test.uu[l1:l2,m1:m2,n1:n2,2]
    oo=curl(test.uu)
    o1 =oo[l1:l2,m1:m2,n1:n2,0]
    o2 =oo[l1:l2,m1:m2,n1:n2,1]
    o3 =oo[l1:l2,m1:m2,n1:n2,2]
cd, writepath
 for nvar = 0, numvar-1 do begin
      
      dfd = vdc_bufwritecreate(vdffile)
      vdc_openvarwrite, dfd, it, varnames(nvar), -1
; Write (transform) the volume to the data set one slice at a time
      case nvar of
         
         0: tmp_var = rho
         1: tmp_var = v1
         2: tmp_var = v2
         3: tmp_var = v3
         4: tmp_var = o1
         5: tmp_var = o2
         6: tmp_var = o3
      endcase
      
      for z = 0, dim[2]-1 do begin
         vdc_bufwriteslice, dfd, float(tmp_var[*,*,z])
      endfor
      
; Close the currently opened variable/time-step. We're done writing
; to it
      vdc_closevar, dfd
   endfor
endfor
plot,  oo[l1:l2,m2/2,n1:n2,1]
;	Destroy the "buffered write" data transformation object. 
;	We're done with it.
vdc_bufwritedestroy, dfd


cd, readpath

end
