;
; $Id$
;
;   Read pvar.dat, or other PVAR file
;
pro pc_read_pvar, object=object, varfile=varfile_, datadir=datadir, ivar=ivar, $
    npar_max=npar_max, stats=stats, quiet=quiet, swap_endian=swap_endian, $
    rmv=rmv, irmv=irmv, trmv=trmv, oldrmv=oldrmv, solid_object=solid_object, $
    theta_arr=theta_arr, savefile=savefile
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
;  Defaults.
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
default, stats, 1
default, quiet, 0
default, rmv, 0
default, oldrmv, 0
default, savefile, 1
;
if (n_elements(ivar) eq 1) then begin
  default,varfile_,'PVAR'
  varfile=varfile_+strcompress(string(ivar),/remove_all)
endif else begin
  default,varfile_,'pvar.dat'
  varfile=varfile_
endelse
;
;  Set rmv if solid_objects
;
if (keyword_set(solid_object)) then rmv=1
;
;  Get necessary dimensions.
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
pc_read_pdim, obj=pdim, datadir=datadir, /quiet
pc_set_precision, dim=dim, /quiet
;
;  Check if file exists.
;
dummy=file_search('./data/param2.nml', COUNT=countfile)
;
; Check if we are inserting particles continuously
;
if (countfile gt 0) then begin
  pc_read_param, object=param2, /param2, datadir=datadir, quiet=quiet
  linsert_particles_continuously=param2.linsert_particles_continuously
endif else begin
  pc_read_param, object=param, datadir=datadir, quiet=quiet
  linsert_particles_continuously=param.linsert_particles_continuously
endelse
;
;  Derived dimensions.
;
mpvar=pdim.mpvar
npar=pdim.npar
default, npar_max, npar & if (npar_max gt npar) then npar_max=npar
ncpus=dim.nprocx*dim.nprocy*dim.nprocz
;
;  Time and grid parameters.
;
t=zero
x=fltarr(dim.mx)*one & y=fltarr(dim.my)*one & z=fltarr(dim.mz)*one
dx=zero &  dy=zero &  dz=zero
;
;  Read processor dimensions if reading just one processor.
;
if (ncpus gt 1) then begin
  pc_read_dim, obj=procdim, datadir=datadir, proc=0, /quiet
endif else begin
  procdim=dim
endelse
xloc=fltarr(procdim.mx)*one & yloc=fltarr(procdim.my)*one & zloc=fltarr(procdim.mz)*one
;
;  Read variable indices from index.pro
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
openr, 1, datadir+'/index.pro'
line=''
while (not eof(1)) do begin
  readf, 1, line, format='(a)'
  if (execute(line) ne 1) then $
    message, 'There was a problem with index.pro', /INF
endwhile
close, 1
;
;  Define structure for data
;
varcontent=replicate( $
    {varcontent_all_par, $
    variable   : 'UNKNOWN', $
    idlvar     : 'dummy', $
    idlinit    : 'fltarr(npar)*one', $
    skip       : 0}, $
    mpvar+1)
;
INIT_SCALAR  = 'fltarr(npar)*one'
INIT_3VECTOR = 'fltarr(npar,3)*one'
;
;  Go through all possible particle variables
;
default, ixp, 0
varcontent[ixp].variable = 'Particle position (xx)'
varcontent[ixp].idlvar   = 'xx'
varcontent[ixp].idlinit  = INIT_3VECTOR
varcontent[ixp].skip     = 2
;
default, ivpx, 0
varcontent[ivpx].variable = 'Particle velocity (vv)'
varcontent[ivpx].idlvar   = 'vv'
varcontent[ivpx].idlinit  = INIT_3VECTOR
varcontent[ivpx].skip     = 2
;
default, iap, 0
varcontent[iap].variable = 'Particle radius (a)'
varcontent[iap].idlvar   = 'a'
varcontent[iap].idlinit  = INIT_SCALAR
;
default, inpswarm, 0
varcontent[inpswarm].variable = 'Particle internal number (npswarm)'
varcontent[inpswarm].idlvar   = 'npswarm'
varcontent[inpswarm].idlinit  = INIT_SCALAR
;
default, irhopswarm, 0
varcontent[irhopswarm].variable = 'Particle mass density (rhopswarm)'
varcontent[irhopswarm].idlvar   = 'rhopswarm'
varcontent[irhopswarm].idlinit  = INIT_SCALAR
;
varcontent[0].variable    = 'UNKNOWN'
varcontent[0].idlvar      = 'UNKNOWN'
varcontent[0].idlinit     = '0.'
varcontent[0].skip        = 0
;
;  Put variable names in array
;
variables = (varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)[1:*]
;
;  Define arrays from contents of varcontent
;
totalvars = mpvar
for iv=1L,totalvars do begin
  if (varcontent[iv].variable eq 'UNKNOWN') then $
      message, 'Unknown variable at position ' + str(iv) $
      + ' needs declaring in pc_read_pvar.pro', /INF
  if (execute(varcontent[iv].idlvar+'='+varcontent[iv].idlinit,0) ne 1) then $
      message, 'Error initialising ' + varcontent[iv].variable $
      +' - '+ varcontent[iv].idlvar, /INFO
  iv=iv+varcontent[iv].skip
endfor
;
;  Define arrays for temporary storage of data.
;
array=fltarr(npar_max,totalvars)*one
ipar=lonarr(npar)
tarr=fltarr(ncpus)*one
t=zero
npar_loc=0L
if (rmv) then begin
  ipar_rmv=lonarr(npar)
  npar_rmv=0L
  trmv    =fltarr(npar)*one
endif
;
;  Loop over processors.
;
for i=0,ncpus-1 do begin
;
  if (not keyword_set(quiet)) then $
      print,'Loading chunk ', strtrim(str(i+1)), ' of ', $
      strtrim(str(ncpus)), ' (', $
      strtrim(datadir+'/proc'+str(i)+'/'+varfile), ')...'
;
  filename=datadir+'/proc'+strtrim(i,2)+'/'+varfile 
;
;  Check if file exists.
;
  dummy=file_search(filename, COUNT=countfile)
  if (not countfile gt 0) then begin
    print, 'ERROR: cannot find file '+ filename
    stop
  endif
;
;  Get a unit number and open file.
;
  get_lun, file
  close, file
  openr, file, filename, /f77, swap_endian=swap_endian
;
;  Read the number of particles at the local processor together with their
;  global index numbers.
;
  readu, file, npar_loc
;
;  Read particle data (if any).
;
  if (npar_loc ne 0) then begin
;    
    ipar_loc=lonarr(npar_loc)
    readu, file, ipar_loc
;
;  Register particle indices for later check if all particles have been read.
;  
    for k=0L,npar_loc-1 do begin
      ipar[ipar_loc[k]-1]=ipar[ipar_loc[k]-1]+1
    endfor
;
;  Read local processor data.
;
    array_loc=fltarr(npar_loc,mpvar)*one
    readu, file, array_loc
;
;  Put local processor data into proper place in global data array
;
    for k=0L,npar_loc-1 do begin
      if (ipar_loc[k] le npar_max) then array[ipar_loc[k]-1,*]=array_loc[k,*]
    endfor
;
  endif
;
;  Read time and grid data.
;
  readu, file, t, xloc, yloc, zloc, dx, dy, dz
;
;  Close time.
;
  close, file
  free_lun, file
;
;  Read information about removed particles.
;
  if (rmv) then begin
    filename1=datadir+'/proc'+strtrim(i,2)+'/rmv_ipar.dat'
    filename2=datadir+'/proc'+strtrim(i,2)+'/rmv_par.dat'
    file_exists=file_test(filename1)
    if (file_exists) then begin
;
;  Read indices and removal times of removed particles.
;
      get_lun, file1 & close, file1
      openr, file1, filename1
      ipar_rmv_loc=0L & t_rmv_loc=0.0*one
;
;  Find out how many particles were removed at this processor. This is done
;  by counting lines until eof, without actually using the read data.
;
      nrmv=0L
      while not (eof(file1)) do begin
        readf, file1, ipar_rmv_loc, t_rmv_loc
        nrmv=nrmv+1
      endwhile
      close, file1
;
;  Read indices and times into array. The index is read as a real number,
;  but converted to integer afterwards.
;
      array_loc=fltarr(2,nrmv)*one
      openr, file1, filename1
      readf, file1, array_loc
      close, file1 & free_lun, file1
      ipar_rmv_loc=reform(long(array_loc[0,*]))
      t_rmv_loc   =reform(     array_loc[1,*])
;
;  Read positions, velocities, etc.
;
      get_lun, file2 & close, file2
      openr, file2, filename2, /f77
      array_loc=fltarr(mpvar)*one
      k=0L
      ipar_rmv_loc_dummy=0L
      for k=0,nrmv-1 do begin
        if (oldrmv) then begin
          readu, file2, ipar_rmv_loc_dummy, array_loc
        endif else begin
          readu, file2, array_loc
        endelse
        if (t_rmv_loc[k] lt t) then begin
          array[ipar_rmv_loc[k]-1,*]=array_loc
          ipar_rmv[ipar_rmv_loc[k]-1]=ipar_rmv[ipar_rmv_loc[k]-1]+1
          trmv[ipar_rmv_loc[k]-1]=t_rmv_loc[k]
        endif
      endfor
      close, file2 & free_lun, file2
    endif
  endif
;
;  Create global x, y and z arrays from local ones.
;
  if (ncpus gt 1) then begin
    pc_read_dim, object=procdim, datadir=datadir, proc=i, /quiet
;
    if (procdim.ipx eq 0L) then begin
      i0x=0L
      i1x=i0x+procdim.mx-1L
      i0xloc=0L
      i1xloc=procdim.mx-1L
    endif else begin
      i0x=procdim.ipx*procdim.nx+procdim.nghostx
      i1x=i0x+procdim.mx-1L-procdim.nghostx
      i0xloc=procdim.nghostx & i1xloc=procdim.mx-1L
    endelse
;
    if (procdim.ipy eq 0L) then begin
      i0y=0L
      i1y=i0y+procdim.my-1L
      i0yloc=0L
      i1yloc=procdim.my-1L
    endif else begin
      i0y=procdim.ipy*procdim.ny+procdim.nghosty
      i1y=i0y+procdim.my-1L-procdim.nghosty
      i0yloc=procdim.nghosty
      i1yloc=procdim.my-1L
    endelse
;
    if (procdim.ipz eq 0L) then begin
      i0z=0L
      i1z=i0z+procdim.mz-1L
      i0zloc=0L
      i1zloc=procdim.mz-1L
    endif else begin
      i0z=procdim.ipz*procdim.nz+procdim.nghostz
      i1z=i0z+procdim.mz-1L-procdim.nghostz
      i0zloc=procdim.nghostz
      i1zloc=procdim.mz-1L
    endelse
;
    x[i0x:i1x] = xloc[i0xloc:i1xloc]
    y[i0y:i1y] = yloc[i0yloc:i1yloc]
    z[i0z:i1z] = zloc[i0zloc:i1zloc]
;
  endif else begin
;
    x=xloc
    y=yloc
    z=zloc
;
  endelse
  tarr[i]=t
;
endfor
;
;  Give indices of removed particles and the time of removal.
;
if (rmv) then begin
  irmv=where(ipar_rmv eq 1)
  if (irmv[0] ne -1) then begin
    trmv=trmv[irmv]
    npar_rmv=n_elements(irmv)
  endif
endif
;
;  Trim x, y and z arrays.
;
x=x[dim.l1:dim.l2]
y=y[dim.m1:dim.m2]
z=z[dim.n1:dim.n2]
;
;  Put data into sensibly named arrays.
;
for iv=1L,mpvar do begin
  res=varcontent[iv].idlvar+'=array[*,iv-1:iv-1+varcontent[iv].skip]'
  if (execute(res,0) ne 1) then $
    message, 'Error putting data into '+varcontent[iv].idlvar+' array'
  iv=iv+varcontent[iv].skip
endfor
;
;  Check if all particles found exactly once.
;  Allow for particles not beeing found if particles are beeing 
;  inserted continuously.
;
if (not keyword_set(quiet)) then begin
  if (linsert_particles_continuously) then begin
    if (rmv) then begin
    if ( (max(ipar+ipar_rmv) gt 1) ) then begin
      print, 'Warning: Some particles found more'
      print, 'than once in snapshot files.'
      print, 'Particle number---No. of occurences'
      for i=0,npar-1 do begin
        if ( (ipar[i]+ipar_rmv[i] gt 1) ) then print, i, ipar[i], ipar_rmv[i]
      endfor
    endif
    endif else begin
    if ( (max(ipar) gt 1)) then begin
      print, 'Warning: Some particles found more'
      print, 'than once in snapshot files.'
      print, 'Particle number---No. of occurences'
      for i=0,npar-1 do begin
        if ( ipar[i] gt 1 ) then print, i, ipar[i]
      endfor
    endif
  endelse  
  endif else begin
  if (rmv) then begin
    if ( (max(ipar+ipar_rmv) ne 1) or (min(ipar+ipar_rmv) ne 1) ) then begin
      print, 'Warning: Some particles not found at all or found more'
      print, 'than once in snapshot files.'
      print, 'Particle number---No. of occurences'
      for i=0,npar-1 do begin
        if ( (ipar[i]+ipar_rmv[i] ne 1) ) then print, i, ipar[i], ipar_rmv[i]
      endfor
    endif
  endif else begin
    if ( (max(ipar) ne 1) or (min(ipar) ne 1) ) then begin
      print, 'Warning: Some particles not found at all or found more'
      print, 'than once in snapshot files.'
      print, 'Particle number---No. of occurences'
      for i=0,npar-1 do begin
        if ( ipar[i] ne 1 ) then print, i, ipar[i]
      endfor
    endif
  endelse
  endelse
endif
;
;  Check if times are consistent between processors.
;
if (min(tarr) ne max(tarr)) then begin
  print, 'The time of the snapshot is inconsistent among the processors!'
  print, 'min(t), max(t)=', min(tarr), max(tarr)
endif
;
;  Print out total number of particles.
;
if (not quiet) then begin
  if (rmv) then begin
    print, ''
    print, 'Found '+strtrim(n_elements(where(ipar eq 1)),2)+' particles ' + $
        'and '+strtrim(npar_rmv,2)+' removed particles'
  endif else begin
    print, ''
    print, 'Found '+strtrim(n_elements(where(ipar eq 1)),2)+' particles'
  endelse
endif
;
;  Print out time.
;
if (not quiet) then begin
  print, ''
  print, 't = ', mean(tarr)
endif
;
;  If requested print a summary.
;
if ( keyword_set(stats) and (not quiet) ) then begin
  print, ''
  print, 'VARIABLE SUMMARY:'
  print, '  name              minval          maxval          mean            stddev'
  for iv=0,n_elements(variables)-1 do begin
    command='size=size('+variables[iv]+') & isvector=size[0] eq 2'
    result=execute(command)
;  Vector.
    if (isvector) then begin
      for ivec=0,2 do begin
        command='minval=min('+variables[iv]+'[*,ivec])'
        result=execute(command)
        command='maxval=max('+variables[iv]+'[*,ivec])'
        result=execute(command)
        command='meanval=mean('+variables[iv]+'[*,ivec])'
        result=execute(command)
        if (npar ge 2) then begin
          command='stdval=stddev('+variables[iv]+'[*,ivec])'
          result=execute(command)
        endif else begin
          stdval=0.0d
        endelse
        if (ivec eq 0) then ind='x'
        if (ivec eq 1) then ind='y'
        if (ivec eq 2) then ind='z'
        if (!prompt eq 'GDL> ') then begin  ; GDL does not understand A-10
          print, ' '+variables[iv]+'_'+ind, '-->', $
              minval, maxval, meanval, stdval, format='(A10,A5,4e15.6)'
        endif else begin
          print, ' '+variables[iv]+'_'+ind, '-->', $
              minval, maxval, meanval, stdval, format='(A-10,A5,4e15.6)'
        endelse
      endfor
    endif else begin
;  Scalar.
      command='minval=min('+variables[iv]+')'
      result=execute(command)
      command='maxval=max('+variables[iv]+')'
      result=execute(command)
      command='meanval=mean('+variables[iv]+')'
      result=execute(command)
      if (npar ge 2) then begin
        command='stdval=stddev('+variables[iv]+')'
        result=execute(command)
      endif else begin
        stdval=0.0d
      endelse
      if (!prompt eq 'GDL> ') then begin  ; GDL does not understand A-10
        print, ' '+variables[iv], '-->', $
            minval, maxval, meanval, stdval, format='(A10,A5,4e15.6)'
      endif else begin
        print, ' '+variables[iv], '-->', $
            minval, maxval, meanval, stdval, format='(A-10,A5,4e15.6)'
      endelse
    endelse
  endfor
endif
;
;  Put data and parameters in object.
;
makeobject="object = create_struct(name=objectname," + $
    "['t','x','y','z','dx','dy','dz'," + $
    arraytostring(variables,quote="'",/noleader) + "]," + $
    "t,x,y,z,dx,dy,dz," + $
    arraytostring(variables,/noleader) + ")"
if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /info
  undefine, object
endif
;
end
