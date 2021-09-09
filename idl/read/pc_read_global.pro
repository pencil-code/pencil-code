;
;  $Id$
;+
;  Same as pc_read_var, but for global variables.
;-
pro pc_read_global,                                                  $
    object=object, varfile=varfile, variables=variables, tags=tags,  $
    validate_variables=validate_variables, trimall=trimall,          $
    nameobject=nameobject, allprocs=allprocs, old_format=old_format, $
    dim=dim, grid=grid, param=param, datadir=datadir, proc=proc,     $
    stats=stats, nostats=nostats, quiet=quiet, help=help,            $
    swap_endian=swap_endian, varcontent=varcontent,                  $
    scalar=scalar, run2D=run2D, single=single

COMPILE_OPT IDL2,HIDDEN
;
; Use common block belonging to derivative routines etc. so we can
; set them up properly.
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_coords, coord_system
;
; If asked for, show some help.
;
  if (keyword_set(help)) then begin
    doc_library, 'pc_read_global'
    return
  endif
;
; Default settings
;
  default, validate_variables, 1
  default, trimall, 0
  default, allprocs, 0
  default, single, 0
;
; Default data directory
;
  datadir = pc_get_datadir(datadir)
;
; Name and path of varfile to read
;
  if (file_test (datadir+'/global.h5')) then default, varfile, 'global.h5'
  default, varfile, 'global.dat'
;
;  Read meta data.
;
  default, varcontent, pc_varcontent_global(datadir=datadir,dim=dim, $
                       param=param,quiet=quiet,scalar=scalar,run2D=run2D,single=single)
;
; Load HDF5 varfile if requested or available.
;
  if (not keyword_set (old_format) and strmid (varfile, strlen(varfile)-3) eq '.h5') then begin

    message, "pc_read_global: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    
    quantities = varcontent[*].idlvar
    num_quantities = n_elements (quantities)
    for pos = 0, num_quantities-1 do begin
      quantity = quantities[pos]
      if (quantity eq 'dummy') then continue
      if (quantity eq 'gg') then begin
        quantity = [ 'global_ggx', 'global_ggy', 'global_ggz' ]
      end else if (quantity eq 'glnTT') then begin
        quantity = [ 'glnTx', 'glnTy', 'glnTz' ]
      end
      if (pos eq 0) then begin
        object = create_struct (quantities[pos], pc_read (quantity, file=varfile, datadir=datadir, trimall=trimall, processor=proc, dim=dim, single=single))
      end else begin
        object = create_struct (object, quantities[pos], pc_read (quantity, trimall=trimall, processor=proc, dim=dim, single=single))
      end
    end
    h5_close_file
    return
  end
;
; Get necessary dimensions quietly
;
  pc_read_dim, object=dim, datadir=datadir, proc=proc, /quiet
  pc_read_param, object=param, dim=dim, datadir=datadir, /quiet, single=single
;
; We know from start.in whether we have to read 2-D or 3-D data.
;
  default, run2D, 0
  if (param.lwrite_2d) then run2D=1
;
; Set coord_system
;
  coord_system=param.coord_system
;
; Call pc_read_grid to make sure any derivative stuff is correctly set in the
; common block. Don't need the data for anything though.
;
  if (is_defined(grid) or (allprocs eq 1)) then $
    procgrid=grid $  ; no validity check of dim
  else $
    pc_read_grid, obj=procgrid, dim=dim, datadir=datadir, param=param, proc=proc, swap_endian=swap_endian, /quiet, single=single
;
; Read dimensions (global)...
;
  if ((n_elements(proc) eq 1) or (allprocs eq 1)) then $
    procdim=dim $    ; no validity check of dim
  else $
    pc_read_dim, object=procdim, datadir=datadir, proc=0, /quiet
;
; Should ghost zones be returned?
;
  if (trimall) then trimxyz=1
;
; Local shorthand for some parameters
;
  nx=dim.nx
  ny=dim.ny
  nz=dim.nz
  nw=nx*ny*nz
  mx=dim.mx
  my=dim.my
  mz=dim.mz
  mw=mx*my*mz
  l1=dim.l1
  l2=dim.l2
  m1=dim.m1
  m2=dim.m2
  n1=dim.n1
  n2=dim.n2
  nghostx=dim.nghostx
  nghosty=dim.nghosty
  nghostz=dim.nghostz
  mvar=dim.mvar
  mxloc=procdim.mx
  myloc=procdim.my
  mzloc=procdim.mz
;
; Number of processors over which to loop.
;
  if ((n_elements(proc) eq 1) or (allprocs eq 1)) then begin
    nprocs=1
  endif else begin
    nprocs=dim.nprocx*dim.nprocy*dim.nprocz
  endelse
;
;  Set up variable/tag lists
;
  totalvars=(size(varcontent))[1]
;
  if (n_elements(variables) ne 0) then begin
    if (keyword_set(additional)) then begin
      filevars=(varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)
      variables=[filevars,variables]
      if (n_elements(tags) ne 0) then begin
        tags=[filevars,tags]
      endif
    endif
  endif else begin
    default,variables,(varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)
  endelse
;
; Default tags are set equal to the variables.
;
  default, tags, variables
;
; Sanity check for variables and tags
;
  if (n_elements(variables) ne n_elements(tags)) then begin
    message, 'ERROR: variables and tags arrays differ in size'
  endif
;
; Prepare for read (build read command)
;
  res=''
  content=''
  for iv=0L,totalvars-1L do begin
    if (n_elements(proc) eq 1) then begin
      res=res+','+varcontent[iv].idlvar
    endif else begin
      res=res+','+varcontent[iv].idlvarloc
    endelse
    content=content+', '+varcontent[iv].variable
;
; Initialise read buffers
;
    if (varcontent[iv].variable eq 'UNKNOWN') then $
        message, 'Unknown variable at position ' + str(iv)  $
        + ' needs declaring in varcontent.pro', /info
    if (execute(varcontent[iv].idlvar+'='+varcontent[iv].idlinit,0) ne 1) then $
        message, 'Error initialising ' + varcontent[iv].variable $
        +' - '+ varcontent[iv].idlvar, /info
    if (n_elements(proc) ne 1) then begin
      if (execute(varcontent[iv].idlvarloc+'='+varcontent[iv].idlinitloc,0) ne 1) then $
          message, 'Error initialising ' + varcontent[iv].variable $
          +' - '+ varcontent[iv].idlvarloc, /info
    endif
;
; For vector quantities skip the required number of elements of the f array
;
    iv=iv+varcontent[iv].skip
  endfor
;
; Display information about the files contents
;
  content = strmid(content,2)
  if ( not keyword_set(quiet) ) then $
      print,'File '+varfile+' contains: ', content
;
; Loop over processors
;
  for i=0,nprocs-1 do begin
    ; Build the full path and filename
    if (allprocs eq 2) then begin
      filename=datadir+'/proc'+str(i)+'/'+varfile
    endif else if (allprocs eq 1) then begin
      filename=datadir+'/allprocs/'+varfile
    endif else if (n_elements(proc) eq 1) then begin
      filename=datadir+'/proc'+str(proc)+'/'+varfile
    endif else begin
      filename=datadir+'/proc'+str(i)+'/'+varfile
      if (not keyword_set(quiet)) then $
          print, 'Loading chunk ', strtrim(str(i+1)), ' of ', $
          strtrim(str(nprocs)), ' (', $
          strtrim(datadir+'/proc'+str(i)+'/'+varfile), ')...'
      ; Read processor box dimensions
      pc_read_dim, object=procdim, datadir=datadir, proc=i, /quiet
    endelse
;
; Check for existence and read the data.
;
    if (not file_test(filename)) then begin
      message, 'ERROR: cannot find file '+ filename
    endif
;
;  Don't overwrite ghost zones of processor to the left (and
;  accordingly in y and z direction makes a difference on the
;  diagonals)
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
; Open a varfile and read some data!
;
    openr, lun, filename, /f77, swap_endian=swap_endian, /get_lun
    if (execute('readu, lun '+res) ne 1) then $
        message, 'Error reading: ' + 'readu, lun ' + res
;
    if (n_elements(proc) ne 1) then begin
;
; Loop over variables.
;
      for iv=0L,totalvars-1L do begin
;
; For 2-D run with lwrite_2d=T we only need to read 2-D data.
;
        if (keyword_set(run2D)) then begin
          if (nx eq 1) then begin
; 2-D run in (y,z) plane.
            cmd =   varcontent[iv].idlvar $
                + "[dim.l1,i0y:i1y,i0z:i1z,*,*]=" $
                + varcontent[iv].idlvarloc $
                +"[i0yloc:i1yloc,i0zloc:i1zloc,*,*]"
          endif else if (ny eq 1) then begin
; 2-D run in (x,z) plane.
            cmd =   varcontent[iv].idlvar $
                + "[i0x:i1x,dim.m1,i0z:i1z,*,*]=" $
                + varcontent[iv].idlvarloc $
                +"[i0xloc:i1xloc,i0zloc:i1zloc,*,*]"
          endif else begin
; 2-D run in (x,y) plane.
            cmd =   varcontent[iv].idlvar $
                + "[i0x:i1x,i0y:i1y,dim.n1,*,*]=" $
                + varcontent[iv].idlvarloc $
                +"[i0xloc:i1xloc,i0yloc:i1yloc,*,*]"
          endelse 
        endif else begin
;
; Regular 3-D run.
;        
          cmd =   varcontent[iv].idlvar $
              + "[i0x:i1x,i0y:i1y,i0z:i1z,*,*]=" $
              + varcontent[iv].idlvarloc $
              +"[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*,*]"
        endelse
        if (execute(cmd) ne 1) then $
            message, 'Error combining data for ' + varcontent[iv].variable
;
; For vector quantities skip the required number of elements
;
        iv=iv+varcontent[iv].skip
      endfor
;
    endif
;
    close, lun
    free_lun, lun
  endfor
;
; Tidy memory a little
;
  if (n_elements(proc) ne 1) then begin
    for iv=0L,totalvars-1L do begin
      undefine, varcontent[iv].idlvarloc
    endfor
  endif
;
; Check variables one at a time and skip the ones that give errors.
; This way the program can still return the other variables, instead
; of dying with an error. One can turn off this option in order to
; reduce execution time.
;
  if (validate_variables) then begin
    skipvariable=make_array(n_elements(variables),/INT,value=0)
    for iv=0,n_elements(variables)-1 do begin
      res=execute(tags[iv]+'='+variables[iv])
      if (not res) then begin
        if (not keyword_set(quiet)) then $
            print,"% Skipping: "+tags[iv]+" -> "+variables[iv]
        skipvariable[iv]=1
      endif
    endfor
    if (min(skipvariable) ne 0) then return
    if (max(skipvariable) eq 1) then begin
      variables=variables[where(skipvariable eq 0)]
      tags=tags[where(skipvariable eq 0)]
    endif
  endif
;
; Save changs to the variables array (but don't include the effect of /TRIMALL)
;
  variables_in=variables
;
; Remove ghost zones if requested.
;
  if (keyword_set(trimall)) then variables = 'pc_noghost('+variables+',dim=dim)'
;
; Make structure out of the variables.
;
  makeobject = "object = " + "CREATE_STRUCT(name=nameobject,[" + $
      arraytostring(tags,QUOTE="'",/noleader) + "]" + $
      arraytostring(variables) + ")"
;
; Execute command to make the structure.
;
  if (execute(makeobject) ne 1) then begin
    message, 'ERROR evaluating variables: '+makeobject
    undefine, object
  endif
;
; If requested print a summary (actually the default - unless being quiet.)
;
  if (keyword_set(stats) or $
     (not (keyword_set(nostats) or keyword_set(quiet)))) then begin
    if (not keyword_set(quiet)) then print, ''
    if (not keyword_set(quiet)) then print, 'VARIABLE SUMMARY:'
    pc_object_stats, object, dim=dim, trim=trimall, quiet=quiet
  endif
;
end
