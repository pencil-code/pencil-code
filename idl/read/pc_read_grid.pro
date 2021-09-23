;;
;; $Id$
;;
;;   Read grid.dat
;;
;;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;;  $Date: 2008-06-12 13:26:15 $
;;  $Revision: 1.17 $
;;
;;  27-nov-02/tony: coded
;;   2-oct-14/MR: keyword parameter down added for use with downsampled data
;;  27-jan-16/MR: added check for FORTRAN consistency of grid data + automatic endian swapping if necessary
;;
pro pc_read_grid, object=object, dim=dim, param=param, trimxyz=trimxyz, $
    datadir=datadir, proc=proc, print=print, quiet=quiet, help=help, $
    swap_endian=swap_endian, allprocs=allprocs, reduced=reduced, down=down, single=single
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat, x, y, z, mx, my, mz, nw, ntmax, date0, time0, nghostx, nghosty, nghostz
  common cdat_limits, l1, l2, m1, m2, n1, n2, nx, ny, nz
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
  common cdat_coords,coord_system
;
; Show some help.
;
  if (keyword_set(help)) then begin
    print, "Usage: "
    print, ""
    print, "pc_read_grid, object=object, dim=dim, param=param, trimxyz=trimxyz, datadir=datadir,"
    print, "              proc=proc, print=print, quiet=quiet, help=help, swap_endian=swap_endian, "
    print, "              allprocs=allprocs, reduced=reduced, down=down, single=single"
    print, "                                                                              "
    print, "Returns the grid position arrays and grid deltas of a Pencil-Code run in object 'object'."
    print, " For a specific processor if requested. Returns zeros and empty in all variables on failure."
    print, "                                                                              "
    print, "  datadir: specify root data directory. Default: './data'          [string]   "
    print, "     proc: specify a processor to get the data from. Default is 0  [integer]  "
    print, " allprocs: specify a collectively written grid file. Default is -1 [integer]  "
    print, "           if allprocs is equal to -1, allprocs will be chosen automatically. "
    print, " /reduced: use a reduced dataset grid file.                                   "
    print, "                                                                              "
    print, "   object: optional structure in which to return all the above     [structure]"
    print, "                                                                              "
    print, " /trimxyz: instruction to remove ghost points from arrays                     "
    print, "   /print: instruction to print all variables to standard output              "
    print, "  /single: instruction to convert the data to single precision on input       "
    print, "    /down: instruction to use downsampled data                                "
    print, "/swap_endian: instruction to swap endian                                      "
    print, "   /quiet: instruction not to print any 'helpful' information                 "
    print, "    /help: display this usage information, and exit                           "
    return
  ENDIF

  default, single, 0
  if is_defined(proc) then $
    if (proc lt 0) then undefine, proc
;
; Default data directory
;
  datadir = pc_get_datadir(datadir)
;
; Get necessary dimensions.
;
  pc_read_dim, object=dim, datadir=datadir, proc=proc, reduced=reduced, QUIET=QUIET, down=down
;
;  Set pc_precision.
;
  pc_set_precision, datadir=datadir, dim=dim, /quiet
;
; Set mx,my,mz in common block for derivative routines
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

  pc_read_param, object=param, datadir=datadir, QUIET=QUIET, single=single
;
; General grid properties.
;
    lequidist = (safe_get_tag (param, 'lequidist', default=[1,1,1]) ne 0)
    lperi = (param.lperi ne 0)
    ldegenerated = ([ nx, ny, nz ] eq 1)
;
;  Set coordinate system.
;
    coord_system=param.coord_system

  NaN = !Values.F_NaN * one
  downprec = (precision eq 'D') and single
;
; Default filename stem.
;
  if keyword_set(down) then $
    gridfile='grid_down' $
  else $
    gridfile='grid'

  filename=datadir+'/'+gridfile+'.h5'

  if file_test(filename) then begin
;
; HDF5 file format - check validity of input object.
;
    if is_valid(object,'GRID',filename) then return

    start = [ 0, 0, 0 ]
    count = [ dim.mx, dim.my, dim.mz ]
    if (is_defined(proc)) then start = [ dim.ipx*dim.nx, dim.ipy*dim.ny, dim.ipz*dim.nz ]

    Lx = pc_read ('grid/Lx', file='grid.h5', datadir=datadir, single=single)
    Ly = pc_read ('grid/Ly', single=single)
    Lz = pc_read ('grid/Lz', single=single)
    if (h5_contains ('grid/Ox')) then begin    ; are overwritten later though
      Ox = pc_read ('grid/Ox', single=single)
      Oy = pc_read ('grid/Oy', single=single)
      Oz = pc_read ('grid/Oz', single=single)
    end else begin
      Ox = param.xyz0[0]
      Oy = param.xyz0[1]
      Oz = param.xyz0[2]
    end

    t = (downprec ? !Values.F_NaN : NaN)
    x = pc_read ('grid/x', start=start[0], count=count[0], single=single)
    y = pc_read ('grid/y', start=start[1], count=count[1], single=single)
    z = pc_read ('grid/z', start=start[2], count=count[2], single=single)
    dx_1 = pc_read ('grid/dx_1', start=start[0], count=count[0], single=single)
    dy_1 = pc_read ('grid/dy_1', start=start[1], count=count[1], single=single)
    dz_1 = pc_read ('grid/dz_1', start=start[2], count=count[2], single=single)
    dx_tilde = pc_read ('grid/dx_tilde', start=start[0], count=count[0], single=single)
    dy_tilde = pc_read ('grid/dy_tilde', start=start[1], count=count[1], single=single)
    dz_tilde = pc_read ('grid/dz_tilde', start=start[2], count=count[2], single=single)
    dx = pc_read ('grid/dx', single=single)
    dy = pc_read ('grid/dy', single=single)
    dz = pc_read ('grid/dz', single=single, /close)
    found_Lxyz = 1
    found_grid_der = 1

  end else begin
;
; Old file format
; Default filename.
;
    gridfile += '.dat'
;
    allprocs_exists = file_test(datadir+'/allprocs/'+gridfile)
    default, swap_endian, 0
;
; Default allprocs.
;
    default, allprocs, -1
    if (allprocs eq -1) then begin
      allprocs=0
      if (allprocs_exists and not is_defined(proc)) then allprocs=1
    end
;
; Check if allprocs is consistent with proc.
;
    if ((allprocs gt 0) and is_defined(proc)) then $
        message, "pc_read_grid: 'allprocs' and 'proc' cannot be set both."
;
; Check if reduced keyword is set.
;
    if (keyword_set(reduced) and is_defined(proc)) then $
        message, "pc_read_grid: /reduced and 'proc' cannot be set both."
;
; Build the full path and filename
;
    ncpus=1 & iterate_cpus=0
    if (keyword_set (reduced)) then begin
      filename=datadir+'/reduced/'+gridfile
    endif else if is_defined(proc) then begin
      filename=datadir+'/proc'+str(proc)+'/'+gridfile
    endif else if ((allprocs gt 0) or allprocs_exists) then begin
    ;;;endif else if ((allprocs gt 0)) then begin
      filename=datadir+'/allprocs/'+gridfile
    endif else begin
      filename=datadir+'/procs/'+gridfile
      ncpus=dim.nprocx*dim.nprocy*dim.nprocz
      iterate_cpus=1
    endelse
;
;  Check validity of input object.
;
    if is_valid(object,'GRID',filename) then return
    filename_loc=filename
    print, 'Reading grid data from ' , filename, ' ...'
;
; Initialize / set default returns for ALL variables, here in precision of the run.
;
    t = NaN
    x = replicate (NaN, mx)
    y = replicate (NaN, my)
    z = replicate (NaN, mz)
    dx = NaN
    dy = NaN
    dz = NaN
    Lx = NaN
    Ly = NaN
    Lz = NaN
    Ox = NaN
    Oy = NaN
    Oz = NaN
    dx_1 = replicate (NaN, mx)
    dy_1 = replicate (NaN, my)
    dz_1 = replicate (NaN, mz)
    dx_tilde = replicate (NaN, mx)
    dy_tilde = replicate (NaN, my)
    dz_tilde = replicate (NaN, mz)
;
; Get a unit number
;
    get_lun, file

    for i=0,ncpus-1 do begin

      if iterate_cpus then begin

        filename_loc=datadir+'/proc'+str(i)+'/'+gridfile
        ; Read processor box dimensions
        undefine, procdim
        pc_read_dim,object=procdim,datadir=datadir,proc=i,QUIET=QUIET, down=down
        xloc = replicate (NaN, procdim.mx)
        yloc = replicate (NaN, procdim.my)
        zloc = replicate (NaN, procdim.mz)
        dx_1loc = replicate (NaN, procdim.mx)
        dy_1loc = replicate (NaN, procdim.my)
        dz_1loc = replicate (NaN, procdim.mz)
        dx_tildeloc = replicate (NaN, procdim.mx)
        dy_tildeloc = replicate (NaN, procdim.my)
        dz_tildeloc = replicate (NaN, procdim.mz)

      endif
      ; Check for existence and read the data
      if (not file_test(filename_loc)) then begin
        FREE_LUN,file
        message, 'ERROR: cannot find file '+ filename_loc
      endif

      check=check_ftn_consistency(filename_loc,swap_endian)
      if check eq -1 then begin
        print, 'File "'+strtrim(filename_loc,2)+'" corrupted!'
        return
      endif else $
        if check eq 1 then $
          print, 'Try to read file "'+strtrim(filename_loc,2)+'" with reversed endian swap!'
     
      if (iterate_cpus and not keyword_set(quiet)) then print, 'Reading ' , filename_loc , ' ...'

      openr,file,filename_loc,/F77,SWAP_ENDIAN=swap_endian

      if (not iterate_cpus) then begin
        readu, file, t, x, y, z
        readu, file, dx, dy, dz
        if (downprec) then begin
          t=float(t) & x=float(x) & y=float(y) & z=float(z) & dx=float(dx) & dy=float(dy) & dz=float(dz)
        endif
        on_ioerror, missing
        readu, file, Lx, Ly, Lz
        if (downprec) then begin
          Lx=float(Lx) & Ly=float(Ly) & Lz=float(Lz)
        endif
        found_Lxyz = 1
        readu, file, dx_1, dy_1, dz_1
        if (downprec) then begin
          dx_1=float(dx_1) & dy_1=float(dy_1) & dz_1=float(dz_1)
        endif
        readu, file, dx_tilde, dy_tilde, dz_tilde
        if (downprec) then begin
          dx_tilde=float(dx_tilde) & dy_tilde=float(dy_tilde) & dz_tilde=float(dz_tilde)
        endif
        found_grid_der = 1
      endif else begin
        if (downprec) then begin
          x=float(x) & y=float(y) & z=float(z)
          dx_1=float(dx_1) & dy_1=float(dy_1) & dz_1=float(dz_1)
          dx_tilde=float(dx_tilde) & dy_tilde=float(dy_tilde) & dz_tilde=float(dz_tilde)
        endif
        readu, file, t, xloc, yloc, zloc
        if ((procdim.ipy eq 0) and (procdim.ipz eq 0)) then begin
        ;
        ;  Don't overwrite ghost zones of processor to the left (and
        ;  accordingly in y and z direction makes a difference on the
        ;  diagonals)
        ;
          if (procdim.ipx eq 0L) then begin
            i0x=0L
            i1x=procdim.mx-1L
            i0xloc=0L
          endif else begin
            i0x=i1x-procdim.nghostx+1L
            i1x=i0x+procdim.mx-1L-procdim.nghostx
            i0xloc=procdim.nghostx
          endelse
          i1xloc=procdim.mx-1L
        ;
          x[i0x:i1x] = xloc[i0xloc:i1xloc]
        endif
        if (procdim.ipx eq 0 and procdim.ipz eq 0) then begin
          if (procdim.ipy eq 0L) then begin
            i0y=0L
            i1y=procdim.my-1L
            i0yloc=0L
          endif else begin
            i0y=i1y-procdim.nghosty+1L
            i1y=i0y+procdim.my-1L-procdim.nghosty
            i0yloc=procdim.nghosty
          endelse
          i1yloc=procdim.my-1L

          y[i0y:i1y] = yloc[i0yloc:i1yloc]
        endif
        if (procdim.ipx eq 0 and procdim.ipy eq 0) then begin
          if (procdim.ipz eq 0L) then begin
            i0z=0L
            i1z=procdim.mz-1L
            i0zloc=0L
          endif else begin
            i0z=i1z-procdim.nghostz+1L
            i1z=i0z+procdim.mz-1L-procdim.nghostz
            i0zloc=procdim.nghostz
          endelse
          i1zloc=procdim.mz-1L

          z[i0z:i1z] = zloc[i0zloc:i1zloc]
        endif

        readu, file, dx, dy, dz
        if (downprec) then begin
          dx=float(dx) & dy=float(dy) & dz=float(dz)
        endif
        on_ioerror, missing
        readu, file, Lx, Ly, Lz
        if (downprec) then begin
          Lx=float(Lx) & Ly=float(Ly) & Lz=float(Lz)
        endif
        found_Lxyz = 1

        readu, file, xloc, yloc, zloc       ; actually d[xyz]_1_loc
        if (procdim.ipy eq 0 and procdim.ipz eq 0) then dx_1[i0x:i1x] = xloc[i0xloc:i1xloc]
        if (procdim.ipx eq 0 and procdim.ipz eq 0) then dy_1[i0y:i1y] = yloc[i0yloc:i1yloc]
        if (procdim.ipx eq 0 and procdim.ipy eq 0) then dz_1[i0z:i1z] = zloc[i0zloc:i1zloc]
        readu, file, xloc, yloc, zloc       ; actually d[xyz]_tilde_loc
        if (procdim.ipy eq 0 and procdim.ipz eq 0) then dx_tilde[i0x:i1x] = xloc[i0xloc:i1xloc]
        if (procdim.ipx eq 0 and procdim.ipz eq 0) then dy_tilde[i0y:i1y] = yloc[i0yloc:i1yloc]
        if (procdim.ipx eq 0 and procdim.ipy eq 0) then dz_tilde[i0z:i1z] = zloc[i0zloc:i1zloc]
        found_grid_der = 1
      endelse

    missing:
      on_ioerror, Null
      close,file

    endfor
;
    free_lun,file
;
;  Construct box origin.
;
    Ox = x[nghostx] - lperi[0] * 0.5 / dx_1[nghostx]
    Oy = y[nghosty] - lperi[1] * 0.5 / dy_1[nghosty]
    Oz = z[nghostz] - lperi[2] * 0.5 / dz_1[nghostz]
;
  endelse
;
;  Add missing parameters for backwards compatibility.
;
  if (not keyword_set (found_Lxyz)) then begin
    Lx = param.xyz0[0]
    Ly = param.xyz0[1]
    Lz = param.xyz0[2]
    found_Lxyz = 1
  end
;
  if (not keyword_set (found_grid_der)) then begin
    if (not any (ldegenerated)) then message, "WARNING: no grid derivatives - please update your run!", /info
;
; replacements for degenerated cases
;
    if (ldegenerated[0]) then begin
      dx_1[*] = 1.0 / dx & dx_tilde[*] = -0.0
    endif
    if (ldegenerated[1]) then begin
      dy_1[*] = 1.0 / dy & dy_tilde[*] = -0.0
    endif
    if (ldegenerated[2]) then  begin
      dz_1[*] = 1.0 / dz & dz_tilde[*] = -0.0
    endif
  end
;
;  Trim ghost zones of coordinate arrays.
;
  if (keyword_set (trimxyz)) then begin
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
;
;  Build structure of all the variables
;
  object = create_struct(name="PC_GRID"+':'+strtrim(filename), $
      ['t','x','y','z','dx','dy','dz','Ox','Oy','Oz','Lx','Ly','Lz', $
       'dx_1','dy_1','dz_1','dx_tilde','dy_tilde','dz_tilde', $
       'lequidist','lperi','ldegenerated'], $
      t,x,y,z,dx,dy,dz,Ox,Oy,Oz,Lx,Ly,Lz,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde, $
      lequidist,lperi,ldegenerated)
;
; If requested print a summary
;
  fmt = '(A,4G15.6)'
  if (keyword_set(print)) then begin
    if (not is_defined(proc)) then begin
      print, FORMAT='(A,I2,A)', 'For all processors calculation domain:'
    endif else if (keyword_set(reduced)) then begin
      print, FORMAT='(A,I2,A)', 'For reduced calculation domain:'
    endif else begin
      print, FORMAT='(A,I2,A)', 'For processor ',proc,' calculation domain:'
    endelse
    print, '             t = ', t
    print, 'min(x), max(x) = ',min(x),', ',max(x)
    print, 'min(y), max(y) = ',min(y),', ',max(y)
    print, 'min(z), max(z) = ',min(z),', ',max(z)
    print, '    dx, dy, dz = ' , dx , ', ' , dy , ', ' , dz
    print, '   periodicity = ' , lperi
    print, '  degeneration = ' , ldegenerated
  endif
;
end
