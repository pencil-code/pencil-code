;
; $Id$
;+
;   Read pvar.dat, or other PVAR file
;-
;  1-mar-17/MR: added new keyword parameters "single" and "array".
;               /single would force to store the data in single precision.
;               array=<name> would force output of the data in an array <name>
;               instead in "object". This saves half the RAM space as no duplication
;               of the (internally anyway used) array into individual variables is performed.
;               "object" contains then the names and starting positions of the variables in array.
;
pro pc_read_pvar, object=object, varfile=varfile_, datadir=datadir, ivar=ivar, $
    npar_max=npar_max, stats=stats, quiet=quiet, swap_endian=swap_endian, $
    rmv=rmv, irmv=irmv, trmv=trmv, oldrmv=oldrmv, $
    solid_object=solid_object, theta_arr=theta_arr, savefile=savefile, help=help, $
    proc=proc, ipar=ipar, trimxyz=trimxyz, id_proc=id_proc, single=single, array=array, $
    lnostoponnopart=lnostoponnopart              
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl

if (keyword_set(help)) then begin
  doc_library, 'pc_read_pvar'
  return
endif
;
;  Defaults.
;
datadir = pc_get_datadir(datadir)
default, stats, 1
default, quiet, 0
default, rmv, 0
default, oldrmv, 0
default, savefile, 1
default, proc, -1
default, trimxyz, 1
default, id_proc, 0
default, single, 0
default, lnostoponnopart, 0
if (keyword_set(lnostoponnopart)) then lnostoponnopart=1

objout = not arg_present(array)
;
if (n_elements(ivar) eq 1) then begin
  default, varfile_, 'PVAR'
  varfile = varfile_ + strcompress (string (ivar), /remove_all)
  if (file_test (datadir+'/allprocs/'+varfile[0]+'.h5')) then varfile += '.h5'
endif else begin
  default_varfile = 'pvar.dat'
  if (file_test (datadir+'/allprocs/pvar.h5')) then default_varfile = 'pvar.h5'
  default, varfile_, default_varfile
  varfile = varfile_
endelse
;
; Load HDF5 varfile if requested or available.
;
  if (strmid (varfile, strlen(varfile)-3) eq '.h5') then begin
    message, "pc_read_pvar: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    pc_read_grid, object=grid, dim=dim, param=param, datadir=datadir, /quiet
    t = pc_read ('time', file=varfile, datadir=datadir,single=single)
    quantities = h5_content('part')
    num_quantities = n_elements (quantities)
    distribution = pc_read ('proc/distribution')
    object = { t:t, x:grid.x, y:grid.y, z:grid.z, dx:grid.dx, dy:grid.dy, dz:grid.dz, distribution:distribution, npar_found:total(distribution, /preserve_type) }
    if (proc ge 0) then begin
      start = 0
      if (proc ge 1) then start = total (object.distribution[0:proc-1])
      count = object.distribution[proc]
    end
    for pos = 0, num_quantities-1 do begin
      quantity = quantities[pos]
      if (strlowcase (quantity) eq 'id') then quantity = 'ipar'
      object = create_struct (object, quantity, pc_read ('part/'+quantities[pos], start=start, count=count), single=single)
    end
    h5_close_file
    return
  end
;
;  Set rmv if solid_objects
;
if (keyword_set(solid_object)) then rmv=1
;
;  Get necessary dimensions.
;
if (proc eq -1) then  $
  pc_read_dim, obj=dim, datadir=datadir, /quiet $
else $
  pc_read_dim, obj=dim, datadir=datadir, proc=proc, /quiet

pc_read_pdim, obj=pdim, datadir=datadir, /quiet
;
; Check if we are inserting particles continuously
;
if (file_test('./data/param2.nml')) then begin
  pc_read_param, object=param2, /param2, datadir=datadir, quiet=quiet, single=single
  linsert_particles_continuously=param2.linsert_particles_continuously
endif else begin
  pc_read_param, object=param, datadir=datadir, quiet=quiet, single=single
  linsert_particles_continuously=param.linsert_particles_continuously
endelse
;
;  Derived dimensions.
;
mpvar=pdim.mpvar
mpaux=pdim.mpaux
if (proc ne -1) then begin
  filename=datadir+'/proc'+strtrim(proc,2)+'/'+varfile
;
;  Check if file exists.
;
  if (not file_test(filename)) then begin
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
;  Read the number of particles at the chosen processor.
;
  npar=0L
  readu, file, npar
;
  close, file
  free_lun, file
endif else begin
  npar=pdim.npar
endelse
default, npar_max, npar & if (npar_max gt npar) then npar_max=npar
ncpus=dim.nprocx*dim.nprocy*dim.nprocz
;
;  Time and grid parameters.
;
t=zero
x=make_array(dim.mx, type=type_idl) & y=make_array(dim.my, type=type_idl) & z=make_array(dim.mz, type=type_idl)
dx=zero &  dy=zero &  dz=zero
;
;  Read processor dimensions if reading just one processor.
;
if (ncpus gt 1) then begin
  pc_read_dim, obj=procdim, datadir=datadir, proc=0, /quiet
endif else begin
  procdim=dim
endelse
xloc=make_array(procdim.mx, type=type_idl) & yloc=make_array(procdim.my, type=type_idl) & zloc=make_array(procdim.mz, type=type_idl)
;
;  Read particle variable indices from particle_index.pro
;
datadir = pc_get_datadir(datadir)
openr, lun, datadir+'/particle_index.pro', /get_lun
line=''
while (not eof(lun)) do begin
  readf, lun, line, format='(a)'
  if (execute(line) ne 1) then $
      message, 'There was a problem with particle_index.pro', /INF
endwhile
close, lun
free_lun, lun
;
;  Define structure for data
;
INIT_SCALAR  = single ? 'fltarr(npar)' : 'make_array(npar, type=type_idl)'
INIT_3VECTOR = single ? 'fltarr(npar,3)' : 'make_array(npar,3, type=type_idl)'
;
varcontent=replicate( $
    {varcontent_all_par, $
    variable   : 'UNKNOWN', $
    idlvar     : 'dummy', $
    idlinit    : INIT_SCALAR, $
    skip       : 0}, $
    mpvar+mpaux+1)
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
default, ivpx_cart, 0
varcontent[ivpx_cart].variable = 'Particle velocity cartesian (vvp_cart)'
varcontent[ivpx_cart].idlvar   = 'vv_cart'
varcontent[ivpx_cart].idlinit  = INIT_3VECTOR
varcontent[ivpx_cart].skip     = 2
;
default, idXp1, 0
varcontent[idXp1].variable = 'Particle Separation (dXp)'
varcontent[idXp1].idlvar   = 'dXp'
varcontent[idXp1].idlinit  = INIT_3VECTOR
varcontent[idXp1].skip     = 2
;
default, idVp1, 0
varcontent[idVp1].variable = 'Particle Velocity Difference(dVp)'
varcontent[idVp1].idlvar   = 'dVp'
varcontent[idVp1].idlinit  = INIT_3VECTOR
varcontent[idVp1].skip     = 2
;
;default, idXpo1, 0
;varcontent[idXpo1].variable = 'Old Particle Separation (dXpo)'
;varcontent[idXpo1].idlvar   = 'dXpo'
;varcontent[idXpo1].idlinit  = INIT_3VECTOR
;varcontent[idXpo1].skip     = 2
;
default, ibpx, 0
varcontent[ibpx].variable = 'Particle magnetic field (bp)'
varcontent[ibpx].idlvar   = 'bp'
varcontent[ibpx].idlinit  = INIT_3VECTOR
varcontent[ibpx].skip     = 2
;
default, iap, 0
varcontent[iap].variable = 'Particle radius (a)'
varcontent[iap].idlvar   = 'a'
varcontent[iap].idlinit  = INIT_SCALAR
;
default, icaustics, 0
varcontent[icaustics].variable = 'Number of caustics'
varcontent[icaustics].idlvar   = 'caustics'
varcontent[icaustics].idlinit  = INIT_SCALAR
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
default, ivpxt, 0
varcontent[ivpxt].variable = 'Particle turbulent speed (vpxt)'
varcontent[ivpxt].idlvar   = 'vpxt'
varcontent[ivpxt].idlinit  = INIT_SCALAR
;
default, ivpzt, 0
varcontent[ivpzt].variable = 'Particle turbulent speed (vpzt)'
varcontent[ivpzt].idlvar   = 'vpzt'
varcontent[ivpzt].idlinit  = INIT_SCALAR
;
default, iaps, 0
varcontent[iaps].variable = 'Particle sink radius (aps)'
varcontent[iaps].idlvar   = 'aps'
varcontent[iaps].idlinit  = INIT_SCALAR
;
default, imp, 0
varcontent[imp].variable = 'Particle mass (mp)'
varcontent[imp].idlvar   = 'mp'
varcontent[imp].idlinit  = INIT_SCALAR
;
default, born, 0
varcontent[born].variable = 'Time of birth (tb)'
varcontent[born].idlvar   = 'born'
varcontent[born].idlinit  = INIT_SCALAR
;
default, iTp, 0
varcontent[iTp].variable = 'Particle temperature (Tp)'
varcontent[iTp].idlvar   = 'Tp'
varcontent[iTp].idlinit  = INIT_SCALAR
;
default, incol, 0
varcontent[incol].variable = 'Number of collisions'
varcontent[incol].idlvar   = 'ncol'
varcontent[incol].idlinit  = INIT_SCALAR
;
default, inucl_Se, 0
varcontent[inucl_Se].variable = 'Supersaturation at nucleation'
varcontent[inucl_Se].idlvar   = 'nucl_Se'
varcontent[inucl_Se].idlinit  = INIT_SCALAR
;
default, iint_Se, 0
varcontent[iint_Se].variable = 'Mass-integrated supersaturation'
varcontent[iint_Se].idlvar   = 'int_Se'
varcontent[iint_Se].idlinit  = INIT_SCALAR
;
default, inucl_T, 0
varcontent[inucl_T].variable = 'Temperature at nucleation'
varcontent[inucl_T].idlvar   = 'nucl_Temp'
varcontent[inucl_T].idlinit  = INIT_SCALAR
;
default, iint_T, 0
varcontent[iint_T].variable = 'Mass-integrated temperature'
varcontent[iint_T].idlvar   = 'int_Temp'
varcontent[iint_T].idlinit  = INIT_SCALAR
;
default, inucl_mix_frac, 0
varcontent[inucl_mix_frac].variable = 'Mixture fraction at nucleation'
varcontent[inucl_mix_frac].idlvar   = 'nucl_mixfrac'
varcontent[inucl_mix_frac].idlinit  = INIT_SCALAR
;
default, iint_mix_frac, 0
varcontent[iint_mix_frac].variable = 'Mass-integrated mixture fraction'
varcontent[iint_mix_frac].idlvar   = 'int_mixfrac'
varcontent[iint_mix_frac].idlinit  = INIT_SCALAR
;
default, iup11, 0
varcontent[iup11].variable = 'grad uu at particle(up11)'
varcontent[iup11].idlvar   = 'up11'
varcontent[iup11].idlinit  = INIT_SCALAR
;
default, iup12, 0
varcontent[iup12].variable = 'grad uu at particle(up12)'
varcontent[iup12].idlvar   = 'up12'
varcontent[iup12].idlinit  = INIT_SCALAR
;
;
default, iup13, 0
varcontent[iup13].variable = 'grad uu at particle(up13)'
varcontent[iup13].idlvar   = 'up13'
varcontent[iup13].idlinit  = INIT_SCALAR
;
default, iup21, 0
varcontent[iup21].variable = 'grad uu at particle(up21)'
varcontent[iup21].idlvar   = 'up21'
varcontent[iup21].idlinit  = INIT_SCALAR
;
default, iup22, 0
varcontent[iup22].variable = 'grad uu at particle(up22)'
varcontent[iup22].idlvar   = 'up22'
varcontent[iup22].idlinit  = INIT_SCALAR
;
default, iup23, 0
varcontent[iup23].variable = 'grad uu at particle(up23)'
varcontent[iup23].idlvar   = 'up23'
varcontent[iup23].idlinit  = INIT_SCALAR
;
default, iup31, 0
varcontent[iup31].variable = 'grad uu at particle(up31)'
varcontent[iup31].idlvar   = 'up31'
varcontent[iup31].idlinit  = INIT_SCALAR
;
default, iup32, 0
varcontent[iup32].variable = 'grad uu at particle(up32)'
varcontent[iup32].idlvar   = 'up32'
varcontent[iup32].idlinit  = INIT_SCALAR
;
default, iup33, 0
varcontent[iup33].variable = 'grad uu at particle(up33)'
varcontent[iup33].idlvar   = 'up33'
varcontent[iup33].idlinit  = INIT_SCALAR
;  Check if there is other pvar data written by the special module.
;
file_special=datadir+'/index_special_particles.pro'
exist_specialvar=file_test(file_special)
if (exist_specialvar eq 1) then begin
  openr, lun, file_special, /get_lun
  line=''
  while (not eof(lun)) do begin
    readf, lun, line
    str_tmp=strsplit(line," ",/extract)
    str=str_tmp[0] & istr=fix(str_tmp[1])
    if (istr gt 0) then begin
      varcontent[istr].variable   = strtrim(str,2)
      varcontent[istr].idlvar     = strtrim(str,2)
      varcontent[istr].idlinit    = INIT_SCALAR
    endif
  endwhile
  close, lun
  free_lun, lun
endif
;
varcontent = varcontent[1:*]
;
;  Put variable names in array
;
variables = (varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)
;
;  Define arrays from contents of varcontent
;
totalvars = mpvar+mpaux
for iv=0L,totalvars-1L do begin
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
ipar=lonarr(npar)
;
; Define array for processor id storage
;
if (id_proc) then begin
  iproc=replicate(-1,npar)
  variables=[variables,'iproc']
endif
;
tarr=make_array(ncpus, type=type_idl)
t=zero
npar_loc=0L
npar=0L
if (rmv) then begin
  npar=pdim.npar
  ipar_rmv=lonarr(npar)
  npar_rmv=0L
  trmv    =make_array(npar, type=type_idl)
endif
;
;  Read from single processor.
;
if (proc ne -1) then begin
  filename=datadir+'/proc'+strtrim(proc,2)+'/'+varfile
  if (not keyword_set(quiet)) then $
      print, 'Loading data from processor ', strtrim(str(proc))
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
  readu, file, npar
;
;  Read particle data (if any).
;
  if (npar ne 0) then begin
;
    ipar=lonarr(npar)
    readu, file, ipar
;
;  Read local processor data.
;
    array=make_array(npar,mpvar+mpaux, type=type_idl)

    readu, file, array
;
  endif
;
;  Read time and grid data.
;
  readu, file, t, xloc, yloc, zloc, dx, dy, dz
  x=xloc
  y=yloc
  z=zloc
;
;  Close time.
;
  close, file
  free_lun, file
;
endif else begin
;
  array=make_array(npar_max,totalvars, type=single ? 4 : type_idl)
;
;  Loop over processors.
;
  for i=0,ncpus-1 do begin
;
    if (not keyword_set(quiet)) then $
        print, 'Loading chunk ', strtrim(str(i+1)), ' of ', $
        strtrim(str(ncpus)), ' (', $
        strtrim(datadir+'/proc'+str(i)+'/'+varfile), ')...'
;
    filename=datadir+'/proc'+strtrim(i,2)+'/'+varfile
;
;  Check if file exists.
;
    if (not file_test(filename)) then begin
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
      if (n_elements(iipar) eq 0) then iipar=ipar_loc else iipar = [iipar, ipar_loc]
;
;  Register particle indices for later check if all particles have been read.
;
      for k=0L,npar_loc-1 do begin
        ipar[ipar_loc[k]-1]=ipar[ipar_loc[k]-1]+1
      endfor
;
;  Read local processor data.
;
      array_loc=make_array(npar_loc,mpvar+mpaux, type=type_idl)
      readu, file, array_loc
;
;  Put local processor data into proper place in global data array
;
      for k=0L,npar_loc-1 do begin
         if (ipar_loc[k] le npar_max) then array[ipar_loc[k]-1,*]=array_loc[k,*]
;
;  Assign processor number so that we can trace where the particles
;  ended up after migrating
;
         if (id_proc) then begin
           if (ipar_loc[k] le npar_max) then iproc[ipar_loc[k]-1]=i
         endif
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
      if (file_test(filename1)) then begin
;
;  Read indices and removal times of removed particles.
;
        openr, lun, filename1, /get_lun
        ipar_rmv_loc=0L
        t_rmv_loc=zero
;
;  Find out how many particles were removed at this processor. This is done
;  by counting lines until eof, without actually using the read data.
;
        nrmv=0L
        dummy=''
        while not (eof(lun)) do begin
          readf, lun, dummy
          nrmv=nrmv+1
        endwhile
        close, lun
        free_lun, lun
        nfields=n_elements(strsplit(dummy,' '))
;
;  Read indices and times into array. The index is read as a real number,
;  but converted to integer afterwards.
;
        array_loc=make_array(nfields,nrmv, type=type_idl)
        get_lun, file1 & close, file1
        openr, file1, filename1
        readf, file1, array_loc
        close, file1 & free_lun, file1
        ipar_rmv_loc=reform(long(array_loc[0,*]))
        t_rmv_loc   =reform(     array_loc[1,*])
        if (nfields eq 3) then begin
          ipar_sink_rmv_loc=reform(long(array_loc[2,*]))
          if (n_elements(array_sink) eq 0) then $
              array_sink=make_array(npar_max,totalvars, type=single ? 4 : type_idl)
          array_sink_loc=make_array(mpvar+mpaux, type=type_idl)
        endif
;
;  Read positions, velocities, etc.
;
        get_lun, file2 & close, file2
        openr, file2, filename2, /f77
        array_loc=make_array(mpvar+mpaux, type=type_idl)
        k=0L
        ipar_rmv_loc_dummy=0L
        for k=0,nrmv-1 do begin
          if (oldrmv) then begin
            readu, file2, ipar_rmv_loc_dummy, array_loc
          endif else begin
            if (nfields eq 3) then begin
              readu, file2, array_loc, array_sink_loc
            endif else begin
              readu, file2, array_loc
            endelse
          endelse
;
;  Unfortunately the time of removal is written in ASCII in ipar_rmv. This can
;  lead to round off errors when comparing t to r_rmv. We therefore determine
;  whether t is single or double precision and allow variations on the last
;  significant digit of t_rmv.
;
          if (n_elements(epsi) eq 0) then begin
            itype=size(t,/type)
            if (itype eq 4) then epsi=1.0d-6
            if (itype eq 5) then epsi=1.0d-12
          endif
;
          if (t_rmv_loc[k] le t*(1.0d + epsi)) then begin
            array[ipar_rmv_loc[k]-1,*]=array_loc
            if (nfields eq 3) then begin
               array_sink[ipar_rmv_loc[k]-1,*]=array_sink_loc
            endif
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
endelse
;
;  Trim x, y and z arrays.
;
if (trimxyz) then begin
  x=x[dim.l1:dim.l2]
  y=y[dim.m1:dim.m2]
  z=z[dim.n1:dim.n2]
endif
;
;  Put data into sensibly named arrays.
;
if objout then begin
  for iv=0L,mpvar+mpaux-1L do begin
    res=varcontent[iv].idlvar+'=array[*,iv:iv+varcontent[iv].skip]'
    if (execute(res,0) ne 1) then $
      message, 'Error putting data into '+varcontent[iv].idlvar+' array'
    iv=iv+varcontent[iv].skip
  endfor
  undefine, array
endif
;
;  Check if all particles found exactly once.
;  Allow for particles not being found if particles are being
;  inserted continuously.
;
if (proc eq -1 and not keyword_set(quiet)) then begin
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
if (proc eq -1 and (min(tarr) ne max(tarr))) then begin
  print, 'The time of the snapshot is inconsistent among the processors!'
  print, 'min(t), max(t)=', min(tarr), max(tarr)
endif
;
;  Print out total number of particles.
;
if (not quiet) then begin
  if (proc eq -1) then begin
    if (rmv) then begin
      print, ''
      print, 'Found '+strtrim(n_elements(where(ipar eq 1)),2)+' particles ' + $
          'and '+strtrim(npar_rmv,2)+' removed particles'
    endif else begin
      print, ''
      print, 'Found '+strtrim(n_elements(where(ipar eq 1)),2)+' particles'
    endelse
  endif else begin
    print, ''
    print, 'Processor snapshot contains '+strtrim(npar,2)+' particles'
  endelse
endif
;
;  Print out time.
;
if (not quiet) then begin
  print, ''
  if (proc eq -1) then print, 't = ', mean(tarr) else $
      print, 't =', t
endif
;
;  If requested print a summary.
;
if objout then $
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
npar_found=n_elements(where(ipar eq 1))+0L
if objout then $
  makeobject="object = create_struct(name=objectname," + $
             "['t','x','y','z','dx','dy','dz','npar_found','ipar'," + $
             arraytostring(variables,quote="'",/noleader) + "]," + $
             "t,x,y,z,dx,dy,dz,npar_found,iipar," + $
             arraytostring(variables,/noleader) + ")" $
else begin
;
;  Turning skips into variable lengths
;
  for iv=0,mpvar+mpaux-1L do begin
    if varcontent[iv].idlvar ne 'dummy' then begin
      iskip = varcontent[iv].skip
      varcontent[iv].skip = iskip+1
      iv += iskip
    endif
  endfor

  makeobject="object = create_struct(name=objectname," + $
             "['t','x','y','z','dx','dy','dz','npar_found','ipar'," + $
             "'varnames', 'varlens']," + $
             "t,x,y,z,dx,dy,dz,npar_found,iipar," + $
             "varcontent.idlvar, varcontent.skip )"
endelse

if (execute(makeobject) ne 1) then begin
  message, 'Error: building of object failed, but data locally available as t,x,y,z,dx,dy,dz,npar_found,iipar'+arraytostring(varcontent.idlvar)+'.', /info
  undefine, object
  if (not lnostoponnopart) then begin
     stop
  endif
endif
;
end
