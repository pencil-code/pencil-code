;
;
;
pro pc_read_global, t=t,                                          $
            object=object, varfile=varfile_, ASSOCIATE=ASSOCIATE, $
            variables=variables,tags=tags, MAGIC=MAGIC,           $
            TRIMXYZ=TRIMXYZ, TRIMALL=TRIMALL,                     $
            nameobject=nameobject,                                $
            dim=dim,param=param,                                  $
            ivar=ivar,                                            $
            datadir=datadir,proc=proc,ADDITIONAL=ADDITIONAL,      $
            nxrange=nxrange,nyrange=nyrange,nzrange=nzrange,      $
            STATS=STATS,NOSTATS=NOSTATS,QUIET=QUIET,HELP=HELP,    $
            SWAP_ENDIAN=SWAP_ENDIAN,varcontent=varcontent,        $
            scalar=scalar,run2D=run2D

COMPILE_OPT IDL2,HIDDEN
;
; Use common block belonging to derivative routines etc. so we can
; set them up properly.
;
  common cdat,x,y,z,mx,my,mz,nw,ntmax,date0,time0
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
  common pc_precision, zero, one
;
; Default data directory
;
; default, datadir, 'data'
  IF (not keyword_set(datadir)) THEN datadir='data'
;
; Name and path of varfile to read
;
  default,varfile_,'global.dat'
  varfile=varfile_
;
; Get necessary dimensions QUIETly
;
  if (n_elements(dim) eq 0) then pc_read_dim,object=dim,datadir=datadir,proc=proc,/quiet
  if (n_elements(param) eq 0) then pc_read_param,object=param,dim=dim,datadir=datadir,/QUIET
;
; Call pc_read_grid to make sure any derivative stuff is correctly set in the common block
; Don't need the data fro anything though
;
  pc_read_grid,dim=dim,datadir=datadir,param=param,/QUIET,SWAP_ENDIAN=SWAP_ENDIAN
;
; Read problem dimensions (global)
;
  if (n_elements(proc) eq 1) then begin
    procdim=dim
  endif else begin
    pc_read_dim,object=procdim,datadir=datadir,proc=0,/QUIET
  endelse
;
; and check pc_precision is set for all Pencil Code tools
;
  pc_set_precision,dim=dim,quiet=quiet
;
; Should ghost zones be returned?
;
  if (keyword_set(TRIMALL)) then begin
    TRIMXYZ=1L
  endif else begin
    TRIMALL=0
  endelse
;
; Local shorthand for some parameters
;
  nx=dim.nx
  ny=dim.ny
  nz=dim.nz
  nw=dim.nx*dim.ny*dim.nz
  mx=dim.mx
  my=dim.my
  mz=dim.mz
  mvar=dim.mvar
  precision=dim.precision
  mxloc=procdim.mx
  myloc=procdim.my
  mzloc=procdim.mz
;
; Number of procs overwhich to loop
;
  if (n_elements(proc) eq 1) then nprocs=1 else nprocs = dim.nprocx*dim.nprocy*dim.nprocz
;
; Initialize / set default returns for ALL variables
;
  t=zero
  x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
  dx=zero &  dy=zero &  dz=zero & deltay=zero

  if (n_elements(proc) ne 1) then begin
    xloc=fltarr(procdim.mx)*one & yloc=fltarr(procdim.my)*one & zloc=fltarr(procdim.mz)*one
  endif
;
;  Read meta data and set up variable/tag lists
;
  default,varcontent,pc_varcontent_global(datadir=datadir,dim=dim, $
                         param=param,quiet=quiet,scalar=scalar,run2D=run2D)
  totalvars=(size(varcontent))[1]-1L
;
  if n_elements(variables) ne 0 then begin
    VALIDATE_VARIABLES=1
    if keyword_set(ADDITIONAL) then begin
      filevars=(varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)[1:*]
      variables=[filevars,variables]
      if n_elements(tags) ne 0 then begin
        tags=[filevars,tags]
      endif
    endif
  endif else begin
    default,variables,(varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)[1:*]
  endelse
  default,tags,variables
;
; Sanity check variables and tags
;
  if (n_elements(variables) ne n_elements(tags)) then begin
    message, 'ERROR: variables and tags arrays differ in size'
  endif
;
; Apply "magic" variable transformations for derived quantities
;
  if keyword_set(MAGIC) then pc_magic_var,variables,tags,param=param, $
      datadir=datadir
;
; Get a free unit number
;
  GET_LUN, file
;
; Prepare for read (build read command)
;
  res=''
  content=''
  for iv=1L,totalvars do begin
    if (n_elements(proc) eq 1) then begin
      res     = res + ',' + varcontent[iv].idlvar
    endif else begin
      res     = res + ',' + varcontent[iv].idlvarloc
    endelse
    content = content + ', ' + varcontent[iv].variable
;
; Initialise read buffers
;
    if (varcontent[iv].variable eq 'UNKNOWN') then $
             message, 'Unknown variable at position ' + str(iv)  $
                      + ' needs declaring in varcontent.pro', /INFO
    if (execute(varcontent[iv].idlvar+'='+varcontent[iv].idlinit,0) ne 1) then $
             message, 'Error initialising ' + varcontent[iv].variable $
                                      +' - '+ varcontent[iv].idlvar, /INFO
    if (n_elements(proc) ne 1) then begin
      if (execute(varcontent[iv].idlvarloc+'='+varcontent[iv].idlinitloc,0) ne 1) then $
               message, 'Error initialising ' + varcontent[iv].variable $
                                      +' - '+ varcontent[iv].idlvarloc, /INFO
    endif
;
; For vector quantities skip the required number of elements of the f array
;
    iv=iv+varcontent[iv].skip
  end
;
; Display information about the files contents
;
  content = strmid(content,2)
  IF ( not keyword_set(QUIET) ) THEN print,'File '+varfile+' contains: ', content
;
; Loop over processors
;
  for i=0,nprocs-1 do begin
    if (n_elements(proc) eq 1) then begin
      ; Build the full path and filename
      filename=datadir+'/proc'+str(proc)+'/'+varfile
    endif else begin
      filename=datadir+'/proc'+str(i)+'/'+varfile
      if (not keyword_set(QUIET)) then $
          print, 'Loading chunk ', strtrim(str(i+1)), ' of ', $
          strtrim(str(nprocs)), ' (', $
          strtrim(datadir+'/proc'+str(i)+'/'+varfile), ')...'
      pc_read_dim,object=procdim,datadir=datadir,proc=i,/QUIET
    endelse
    ; Check for existance and read the data
    dummy=findfile(filename, COUNT=countfile)
    if (not countfile gt 0) then begin
      message, 'ERROR: cannot find file '+ filename
    endif
;
; Setup the coordinates mappings from the processor
; to the full domain.
;
    if (n_elements(proc) eq 1) then begin
    endif else begin
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
    endelse
;
; Open a varfile and read some data!
;
    close,file
    openr,file, filename, /F77,SWAP_ENDIAN=SWAP_ENDIAN
    if (execute('readu,file'+res) ne 1) then $
        message, 'Error reading: ' + 'readu,'+str(file)+res
;
    if (n_elements(proc) ne 1) then begin

      for iv=1L,totalvars do begin
        if (varcontent[iv].variable eq 'UNKNOWN') then continue
  ;DEBUG: tmp=execute("print,'Minmax of "+varcontent[iv].variable+" = ',minmax("+varcontent[iv].idlvarloc+")")
        if (not keyword_set(run2D)) then begin
          ; classical 3D-run (x,y,z)
          cmd =   varcontent[iv].idlvar $
              + "[i0x:i1x,i0y:i1y,i0z:i1z,*,*]=" $
              + varcontent[iv].idlvarloc $
              +"[i0xloc:i1xloc,i0yloc:i1yloc,i0zloc:i1zloc,*,*]"
        endif else begin
          if (ny eq 1) then begin
            ; 2D-run in plane (x,z)
            cmd =   varcontent[iv].idlvar $
                + "[i0x:i1x,i0z:i1z,*,*]=" $
                + varcontent[iv].idlvarloc $
                +"[i0xloc:i1xloc,i0zloc:i1zloc,*,*]"
           endif else begin
             ; 2D-run in plane (x,y)
             cmd =   varcontent[iv].idlvar $
                 + "[i0x:i1x,i0y:i1y,*,*]=" $
                 + varcontent[iv].idlvarloc $
                 +"[i0xloc:i1xloc,i0yloc:i1yloc,*,*]"
           endelse
        endelse
        if (execute(cmd) ne 1) then $
            message, 'Error combining data for ' + varcontent[iv].variable

        ; For vector quantities skip the required number of elements
        iv=iv+varcontent[iv].skip
      endfor

    endif

    close,file
    FREE_LUN,file

  endfor
;
; Tidy memory a little
;
  if (n_elements(proc) ne 1) then begin
    for iv=1L,totalvars do begin
      undefine,varcontent[iv].idlvarloc
    endfor
  endif
;
  if (keyword_set(VALIDATE_VARIABLES)) then begin
    skipvariable=make_array(n_elements(variables),/INT,value=0)
    for iv=0,n_elements(variables)-1 do begin
    ;  res1=execute("testvariable=n_elements("+variables[iv]+")")
      res=execute(tags[iv]+'='+variables[iv])
      if (not res) then begin
        if (not keyword_set(QUIET)) then print,"% Skipping: "+tags[iv]+" -> "+variables[iv]
        skipvariable[iv]=1
      endif
    endfor
    testvariable=0
    if (min(skipvariable) ne 0) then begin
      return
    endif
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
; Trim x, y and z if requested
;
  if (keyword_set(TRIMXYZ)) then begin
    xyzstring="x[dim.l1:dim.l2],y[dim.m1:dim.m2],z[dim.n1:dim.n2]"
  endif else begin
    xyzstring="x,y,z"
  endelse

  if keyword_set(TRIMALL) then begin
  ;  if not keyword_set(QUIET) then print,'NOTE: TRIMALL assumes the result of all specified variables has dimensions from the varfile (with ghosts)'
    variables = 'pc_noghost('+variables+',dim=dim,run2D=run2D)'
  endif
;
  makeobject = "object = "+ $
      "CREATE_STRUCT(name=objectname,['t','x','y','z','dx','dy','dz'" + $
      arraytostring(tags,QUOTE="'") + "],t,"+xyzstring+",dx,dy,dz" + $
      arraytostring(variables) + ")"
;
  if (execute(makeobject) ne 1) then begin
    message, 'ERROR Evaluating variables: '+makeobject
    undefine,object
  endif
;
; If requested print a summary (actually the default - unless being QUIET.)
;
  if (keyword_set(STATS)) or (not (keyword_set(NOSTATS) or keyword_set(QUIET))) then begin
    if (not keyword_set(QUIET)) then print,''
    if (not keyword_set(QUIET)) then print,'VARIABLE SUMMARY:'
    pc_object_stats, object, dim=dim, TRIM=TRIMALL, QUIET=QUIET
  endif
;
end
