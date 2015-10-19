; Create the grid data structure from the grid file.
;
; 27-nov-02/tony: coded
;  2-oct-14/MR: keyword 'down' added to trigger reading of downsampled grid
;
pro pc_read_grid, object=object, dim=dim, param=param, trimxyz=trimxyz, $
    datadir=datadir, proc=proc, print=print, quiet=quiet, help=help, $
    swap_endian=swap_endian, allprocs=allprocs, reduced=reduced, down=down
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, mw, nx, ny, nz, nw, ntmax, date0, time0
  common cdat_limits, l1, l2, m1, m2, n1, n2, nghostx, nghosty, nghostz
  common cdat_grid, dx_1, dy_1, dz_1, dx_tilde, dy_tilde, dz_tilde, lequidist, lperi, ldegenerated
  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
  if (keyword_set(help)) then begin
    print, "Usage:"
    print, ""
    print, "pc_read_grid, t=t, x=x, y=y, z=z, dx=dx, dy=dy, dz=dz, $                      "
    print, "              datadir=datadir, proc=proc, $                                   "
    print, "              /PRINT, /QUIET, /HELP                                           "
    print, ""
    print, "Returns the grid position arrays and grid deltas of a Pencil-Code run. For a  "
    print, "specific processor. Returns zeros and empty in all variables on failure.      "
    print, ""
    print, "  datadir: specify root data directory. Default: './data'          [string]   "
    print, "     proc: specify a processor to get the data from. Default is 0  [integer]  "
    print, " allprocs: specify a collectively written grid file.               [integer]  "
    print, " /reduced: use a reduced dataset grid file.                                   "
    print, ""
    print, "   object: optional structure in which to return all the above     [structure]"
    print, ""
    print, "   /PRINT: instruction to print all variables to standard output              "
    print, "   /QUIET: instruction not to print any 'helpful' information                 "
    print, "    /HELP: display this usage information, and exit                           "
    return
  endif

  ; Default data directory
  if (not keyword_set(datadir)) then datadir = pc_get_datadir()

  ; Default filename
  default, gridfile, 'grid.dat'
  if (keyword_set(down)) then gridfile = 'grid_down.dat'

  ; Check if allprocs is consistent with proc.
  if (keyword_set (allprocs) and (size(proc, /type) ne 0)) then $
      message, "pc_read_grid: 'allprocs' and 'proc' cannot be set both."

  ; Check if reduced keyword is set.
  if (keyword_set(reduced) and (size(proc, /type) ne 0)) then $
      message, "pc_read_grid: /reduced and 'proc' cannot be set both."

  ; Get necessary dimensions.
  if (size(dim, /type) eq 0) then $
      pc_read_dim, object=dim, datadir=datadir, proc=proc, reduced=reduced, QUIET=QUIET, down=down
  if (size(param, /type) eq 0) then $
      pc_read_param, object=param, datadir=datadir, QUIET=QUIET

  if (keyword_set (allprocs) or (size(proc, /type) ne 0) or keyword_set (reduced)) then begin
    multiple_reads = 0
  endif else begin
    multiple_reads = 1
  endelse

  ; Set common blocks for derivative routines
  nx = dim.nx
  ny = dim.ny
  nz = dim.nz
  nw = nx*ny*nz
  mx = dim.mx
  my = dim.my
  mz = dim.mz
  mw = mx*my*mz
  nghostx = dim.nghostx
  nghosty = dim.nghosty
  nghostz = dim.nghostz
  l1 = dim.l1
  l2 = dim.l2
  m1 = dim.m1
  m2 = dim.m2
  n1 = dim.n1
  n2 = dim.n2
  ncpus = dim.nprocx * dim.nprocy * dim.nprocz

  lequidist = (safe_get_tag (param, 'lequidist', default=[1,1,1]) ne 0)
  lperi = (param.lperi ne 0)
  ldegenerated = ([nx,ny,nz] eq 1)

  ; Initialize / set default returns for ALL variables
  t = zero
  x = make_array (mx, type=type_idl)
  y = make_array (my, type=type_idl)
  z = make_array (mz, type=type_idl)
  dx = zero
  dy = zero
  dz = zero
  Lx = zero
  Ly = zero
  Lz = zero
  dx_1 = make_array (mx, type=type_idl)
  dy_1 = make_array (my, type=type_idl)
  dz_1 = make_array (mz, type=type_idl)
  dx_tilde = make_array (mx, type=type_idl)
  dy_tilde = make_array (my, type=type_idl)
  dz_tilde = make_array (mz, type=type_idl)

  for i=0, (ncpus-1)*multiple_reads do begin

    ; Build the full filename and read the data
    if (keyword_set (reduced)) then begin
      procdir = '/reduced/'
      procdim = dim
    end else if (keyword_set (allprocs)) then begin
      procdir = '/allprocs/'
      procdim = dim
    endif else if (size(proc, /type) ne 0) then begin
      procdir = '/proc'+str(proc)+'/'
      pc_read_dim, object=procdim, datadir=datadir, proc=i, QUIET=QUIET, down=down
    endif else begin
      procdir = '/proc'+str(i)+'/'
      pc_read_dim, object=procdim, datadir=datadir, proc=i, QUIET=QUIET, down=down
    endelse

    ; Check for existance and read the data
    filename = datadir + procdir + gridfile
    if (not file_test (filename)) then message, 'ERROR: cannot find file "'+filename+'".'

    ; Set processor box dimensions
    xloc = make_array (procdim.mx, type=type_idl)
    yloc = make_array (procdim.my, type=type_idl)
    zloc = make_array (procdim.mz, type=type_idl)
    dx_1loc = make_array (procdim.mx, type=type_idl)
    dy_1loc = make_array (procdim.my, type=type_idl)
    dz_1loc = make_array (procdim.mz, type=type_idl)
    dx_tildeloc = make_array (procdim.mx, type=type_idl)
    dy_tildeloc = make_array (procdim.my, type=type_idl)
    dz_tildeloc = make_array (procdim.mz, type=type_idl)

    ; Don't overwrite ghost zones of processor to the left (and
    ; accordingly in y and z direction makes a difference on the diagonals)
    if (procdim.ipx eq 0L) then begin
      i0x = 0L
      i1x = i0x+procdim.mx-1L
      i0xloc = 0L
      i1xloc = procdim.mx-1L
    endif else begin
      i0x = procdim.ipx*procdim.nx+procdim.nghostx
      i1x = i0x+procdim.mx-1L-procdim.nghostx
      i0xloc = procdim.nghostx
      i1xloc = procdim.mx-1L
    endelse

    if (procdim.ipy eq 0L) then begin
      i0y = 0L
      i1y = i0y+procdim.my-1L
      i0yloc = 0L
      i1yloc = procdim.my-1L
    endif else begin
      i0y = procdim.ipy*procdim.ny+procdim.nghosty
      i1y = i0y+procdim.my-1L-procdim.nghosty
      i0yloc = procdim.nghosty
      i1yloc = procdim.my-1L
    endelse

    if (procdim.ipz eq 0L) then begin
      i0z = 0L
      i1z = i0z+procdim.mz-1L
      i0zloc = 0L
      i1zloc = procdim.mz-1L
    endif else begin
      i0z = procdim.ipz*procdim.nz+procdim.nghostz
      i1z = i0z+procdim.mz-1L-procdim.nghostz
      i0zloc = procdim.nghostz
      i1zloc = procdim.mz-1L
    endelse

    ; Check for existance and read the data
    filename = datadir + procdir + gridfile
    if (not file_test(filename)) then begin
      message, 'ERROR: cannot find file "'+filename+'".'
    endif

    if (not keyword_set (quiet)) THEN print, 'Reading "' , filename , '"...'

    default, found_Lxyz, 0
    default, found_grid_der, 0
    openr, file, filename, /f77, swap_endian=swap_endian, /get_lun
    if (not multiple_reads) then begin
      readu, file, t, x, y, z
      readu, file, dx, dy, dz
      if (not EOF (file)) then begin
        found_Lxyz = 1
        readu, file, Lx, Ly, Lz
        if (not EOF (file)) then begin
          found_grid_der = 1
          readu, file, dx_1, dy_1, dz_1
          readu, file, dx_tilde, dy_tilde, dz_tilde
        endif
      endif
    endif else begin
      readu, file, t, xloc, yloc, zloc
      readu, file, dx, dy, dz
      x[i0x:i1x] = xloc[i0xloc:i1xloc]
      y[i0y:i1y] = yloc[i0yloc:i1yloc]
      z[i0z:i1z] = zloc[i0zloc:i1zloc]
      if (not EOF (file)) then begin
        found_Lxyz = 1
        readu, file, Lx, Ly, Lz
        if (not EOF (file)) then begin
          found_grid_der = 1
          readu, file, dx_1loc, dy_1loc, dz_1loc
          readu, file, dx_tildeloc, dy_tildeloc, dz_tildeloc
          dx_1[i0x:i1x] = dx_1loc[i0xloc:i1xloc]
          dy_1[i0y:i1y] = dy_1loc[i0yloc:i1yloc]
          dz_1[i0z:i1z] = dz_1loc[i0zloc:i1zloc]
          dx_tilde[i0x:i1x] = dx_tildeloc[i0xloc:i1xloc]
          dy_tilde[i0y:i1y] = dy_tildeloc[i0yloc:i1yloc]
          dz_tilde[i0z:i1z] = dz_tildeloc[i0zloc:i1zloc]
        endif
      endif
    endelse
    close, file
    free_lun, file
  endfor

  ; Trim ghost zones of coordinate arrays.
  if (keyword_set(trimxyz)) then begin
    x = x[l1:l2]
    y = y[m1:m2]
    z = z[n1:n2]
    dx_1 = dx_1[l1:l2]
    dy_1 = dy_1[m1:m2]
    dz_1 = dz_1[n1:n2]
    dx_tilde = dx_tilde[l1:l2]
    dy_tilde = dy_tilde[m1:m2]
    dz_tilde = dz_tilde[n1:n2]
  endif

  ; Build structure of all the variables
  if (found_Lxyz and found_grid_der) then begin
    if (keyword_set(trimxyz)) then begin
      Ox = x[0] - lperi[0] * 0.5 / dx_1[0]
      Oy = y[0] - lperi[1] * 0.5 / dy_1[0]
      Oz = z[0] - lperi[2] * 0.5 / dz_1[0]
    endif else begin
      Ox = x[nghostx] - lperi[0] * 0.5 / dx_1[nghostx]
      Oy = y[nghosty] - lperi[1] * 0.5 / dy_1[nghosty]
      Oz = z[nghostz] - lperi[2] * 0.5 / dz_1[nghostz]
    endelse
    object = create_struct(name="pc_read_grid_" + $
        str((size(x))[1]) + '_' + $
        str((size(y))[1]) + '_' + $
        str((size(z))[1]), $
        ['t','x','y','z','dx','dy','dz','Ox','Oy','Oz','Lx','Ly','Lz', $
         'dx_1','dy_1','dz_1','dx_tilde','dy_tilde','dz_tilde', $
         'lequidist','lperi','ldegenerated'], $
        t,x,y,z,dx,dy,dz,Ox,Oy,Oz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde, $
        lequidist,lperi,ldegenerated)
  endif else if (found_Lxyz) then begin
    object = create_struct(name="pc_read_grid_" + $
        str((size(x))[1]) + '_' + $
        str((size(y))[1]) + '_' + $
        str((size(z))[1]), $
        ['t','x','y','z','dx','dy','dz','Lx','Ly','Lz'], $
        t,x,y,z,dx,dy,dz,Lx,Ly,Lz)
    dx_1 = zero
    dy_1 = zero
    dz_1 = zero
    dx_tilde = zero
    dy_tilde = zero
    dz_tilde = zero
  endif else if (found_grid_der) then begin
    object = create_struct(name="pc_read_grid_" + $
        str((size(x))[1]) + '_' + $
        str((size(y))[1]) + '_' + $
        str((size(z))[1]), $
        ['t','x','y','z','dx','dy','dz','dx_1','dy_1','dz_1', $
         'dx_tilde','dy_tilde','dz_tilde'], $
        t,x,y,z,dx,dy,dz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde)
  endif

  ; If requested print a summary
  fmt = '(A,4G15.6)'
  if (keyword_set(print)) then begin
    if (size(proc, /type) eq 0) then begin
      print, FORMAT='(A,I2,A)', 'For all processors calculation domain:'
    endif else if (keyword_set(reduced)) then begin
      print, FORMAT='(A,I2,A)', 'For reduced calculation domain:'
    endif else begin
      print, FORMAT='(A,I2,A)', 'For processor ',proc,' calculation domain:'
    endelse
    print, '             t = ', t
    print, 'min(x), max(x) = ', min(x), ', ', max(x)
    print, 'min(y), max(y) = ', min(y), ', ', max(y)
    print, 'min(z), max(z) = ', min(z), ', ', max(z)
    print, '    dx, dy, dz = ', dx, ', ', dy, ', ', dz
    print, '   periodicity = ', lperi
    print, '  degeneration = ', ldegenerated
  endif

end

