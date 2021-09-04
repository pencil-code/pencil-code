;
;  $Id$
;
;  Read information from dim.dat
;
;  Author: Tony Mee (A.J.Mee@ncl.ac.uk)
;
;  27-nov-02/tony: coded
;   2-oct-14/MR: keyword parameter down added for use with downsampled data
;  14-Sep-2015/PABourdin: calling 'pc_set_precision' directly after reading
;
pro pc_read_dim, mx=mx, my=my, mz=mz, mw=mw, mvar=mvar, $
    nx=nx, ny=ny, nz=nz, nw=nw, $
    nxgrid=nxgrid, nygrid=nygrid, nzgrid=nzgrid, $
    mxgrid=mxgrid, mygrid=mygrid, mzgrid=mzgrid, $
    precision=precision_new, $
    nghostx=nghostx, nghosty=nghosty, nghostz=nghostz, $
    nprocx=nprocx, nprocy=nprocy, nprocz=nprocz, $
    ipx=ipx,ipy=ipy,ipz=ipz, $
    l1=l1, l2=l2, m1=m1, m2=m2, n1=n1, n2=n2, $
    object=object, datadir=datadir, proc=proc, reduced=reduced, $
    print=print, quiet=quiet, help=help, down=down, ogrid=ogrid
;
COMPILE_OPT IDL2, HIDDEN
;
  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
  if ( keyword_set(HELP) ) then begin
    print, "Usage:"
    print, ""
    print, "pc_read_dim, mx=mx, my=my, mz=mz, mvar=mvar,                                                "
    print, "             nx=nx, ny=ny, nz=nz,                                                           "
    print, "             mxgrid=mxgrid, mygrid=mygrid, mzgrid=mzgrid,                                   "
    print, "             nxgrid=nxgrid, nygrid=nygrid, nzgrid=nzgrid,                                   "
    print, "             precision=precision,                                                           "
    print, "             nghostx=nghostx, nghosty=nghosty, nghostz=nghostz,                             "
    print, "             nprocx=nprocx, nprocy=nprocy, nprocz=nprocz,                                   "
    print, "             object=object,                                                                 "
    print, "             datadir=datadir, proc=proc, /PRINT, /QUIET, /HELP                              "
    print, ""
    print, "Returns the run dimensions of a Pencil-Code run.  Either for the whole calculation,         "
    print, "or if `proc' is defined then for a specific processor.                                      "
    print, "Returns zeros and empty in all variables on failure.                                        "
    print, ""
    print, "  datadir: specify the root data directory. Default is './data'                     [string]"
    print, "     proc: specify a processor to get the data from eg. 0                          [integer]"
    print, "           If unspecified data is read for global calculation.                              "
    print, " /reduced: use a reduced dataset dimension file.                                            "
    print, ""
    print, "       mx: x dimension of processor calculation domain including ghost zones       [integer]"
    print, "       my: y dimension of processor calculation domain including ghost zones       [integer]"
    print, "       mz: z dimension of processor calculation domain including ghost zones       [integer]"
    print, "       mw: defined as mx * my * mz                                                 [integer]"
    print, "       nx: x dimension of processor calculation domain excluding ghost zones       [integer]"
    print, "       ny: y dimension of processor calculation domain excluding ghost zones       [integer]"
    print, "       nz: z dimension of processor calculation domain excluding ghost zones       [integer]"
    print, "       nw: defined as nx * ny * nz                                                 [integer]"
    print, "   nxgrid: x dimension of full calculation domain excluding ghost zones            [integer]"
    print, "   nygrid: y dimension of full calculation domain excluding ghost zones            [integer]"
    print, "   nzgrid: z dimension of full calculation domain excluding ghost zones            [integer]"
    print, "   mxgrid: x dimension of full calculation domain including ghost zones            [integer]"
    print, "   mygrid: y dimension of full calculation domain including ghost zones            [integer]"
    print, "   mzgrid: z dimension of full calculation domain including ghost zones            [integer]"
    print, "     mvar: total number of computed scalar variables (NB. 1 vector = 3 scalars)    [integer]"
    print, "precision: returns 'S' or 'D' for Single or Double precision                     [character]"
    print, "  nghostx: number of points in x direction used for ghost zone at each boundary    [integer]"
    print, "  nghosty: number of points in y direction used for ghost zone at each boundary    [integer]"
    print, "  nghostx: number of points in z direction used for ghost zone at each boundary    [integer]"
    print, "   nprocx: number of communicating processors in the x direction                   [integer]"
    print, "   nprocy: number of communicating processors in the y direction                   [integer]"
    print, "   nprocz: number of communicating processors in the z direction                   [integer]"
    print, "   l1, l2: first & last index of non-ghost-point in x                              [integer]"
    print, "   m1, m2: first & last index of non-ghost-point in y                              [integer]"
    print, "   n1, n2: first & last index of non-ghost-point in z                              [integer]"
    print, ""
    print, "   object: optional structure in which to return all the above as tags           [structure]"
    print, ""
    print, "   /PRINT: instruction to print all variables to standard output                            "
    print, "   /QUIET: instruction not to print any 'helpful' information                               "
    print, "    /HELP: display this usage information, and exit                                         "
    return
  endif
;
;  Set default data directory.
;
  datadir = pc_get_datadir(datadir)
;
;  Initialize / set default returns for ALL variables.
;
  mx = 0L
  my = 0L
  mz = 0L
  mw = 0L
  mvar = 0L
  maux = 0L
  mglobal = 0L
  nx = 0L
  ny = 0L
  nz = 0L
  nw = 0L
  nghostx = 0L
  nghosty = 0L
  nghostz = 0L
  ipx = -1L
  ipy = -1L
  ipz = -1L
  nprocx = 0L
  nprocy = 0L
  nprocz = 0L
;
; Default filename
;
  default, dimfile, 'dim.dat'
  if (keyword_set(down)) then dimfile = 'dim_down.dat'
  if (keyword_set(ogrid)) then dimfile = 'ogdim.dat'
;
;  Build the full path and filename.
;
  if (is_defined(proc)) then begin
    ; Read dimensions on local processor.
    if (keyword_set(reduced)) then $
        message, "pc_read_dim: /reduced and 'proc' cannot be set both."
    filename = datadir+'/proc'+str(proc)+'/'+dimfile
  end else begin
    ; Read global dimensions.
    filename = datadir+'/'+dimfile
    if (keyword_set(reduced)) then filename = datadir+'/reduced/'+dimfile
  end
;
  if (file_test (datadir+'/grid.h5')) then begin
    ; HDF5 format is available
    filename = datadir+'/grid.h5'
    if (not keyword_set(quiet)) then print, 'Reading ' + filename + '...'
    mxgrid = h5_read ('settings/mx', file=filename)
    mygrid = h5_read ('settings/my')
    mzgrid = h5_read ('settings/mz')
    nxgrid = h5_read ('settings/nx')
    nygrid = h5_read ('settings/ny')
    nzgrid = h5_read ('settings/nz')
    mvar = h5_read ('settings/mvar')
    maux = h5_read ('settings/maux')
    mglobal = h5_read ('settings/mglobal')
    nghostx = h5_read ('settings/nghost')
    nghosty = nghostx
    nghostz = nghostx
    nprocx = h5_read ('settings/nprocx')
    nprocy = h5_read ('settings/nprocy')
    nprocz = h5_read ('settings/nprocz')
    precision_new = h5_read ('settings/precision', /close)
    pc_set_precision, precision=precision_new
    if (size (proc, /type) ne 0) then begin
      ipx = proc mod nprocx
      ipy = (proc / nprocx) mod nprocy
      ipz = proc / (nprocx * nprocy)
      mx = nxgrid / nprocx + 2*nghostx
      my = nygrid / nprocy + 2*nghosty
      mz = nzgrid / nprocz + 2*nghostz
    end else begin
      ipx = -1
      ipy = -1
      ipz = -1
      mx = mxgrid
      my = mygrid
      mz = mzgrid
    end
  end else begin
    ; old file format
;
;  Check for existence and read the data.
;
    if (not file_test(filename)) then message, 'ERROR: cannot find file ' + filename
    if (not keyword_set(quiet)) then print, 'Reading ' + filename + '...'
    openr, file, filename, /get_lun
    if (execute('readf,file,mx,my,mz,mvar,maux,mglobal') ne 1) then begin
      ; For backwards compatibility with dim.dat without mglobal.
      print
      print, 'Note: the Input conversion error is of no significance.'
      print, 'Will now read without the mglobal parameter.'
      print
      close, file
      openr, file, filename
      readf, file, mx, my, mz, mvar, maux
      mglobal = 0
    end
    precision_new = ''
    readf, file, precision_new
    pc_set_precision, precision=precision_new
    readf, file, nghostx, nghosty, nghostz
    if (size(proc, /type) ne 0) then begin
      readf, file, ipx, ipy, ipz
    end else begin
      readf, file, nprocx, nprocy, nprocz
    end
    close,file
    free_lun, file
    ;
    if (size(proc, /type) ne 0) then begin
      pc_read_dim, obj=globdim, datadir=datadir, /quiet
      nprocx = globdim.nprocx
      nprocy = globdim.nprocy
      nprocz = globdim.nprocz
      nxgrid = globdim.nxgrid
      nygrid = globdim.nygrid
      nzgrid = globdim.nzgrid
      mxgrid = globdim.mxgrid
      mygrid = globdim.mygrid
      mzgrid = globdim.mzgrid
    end else begin
      mxgrid = mx
      mygrid = my
      mzgrid = mz
      nxgrid = mxgrid - (2 * nghostx)
      nygrid = mygrid - (2 * nghosty)
      nzgrid = mzgrid - (2 * nghostz)
    end
  end
;
;  Calculate any derived quantities
;
  mw = mx * my * mz
  nx = mx - (2 * nghostx)
  ny = my - (2 * nghosty)
  nz = mz - (2 * nghostz)
  nw = nx * ny * nz
  l1 = nghostx & l2 = mx-nghostx-1
  m1 = nghosty & m2 = my-nghosty-1
  n1 = nghostz & n2 = mz-nghostz-1
;
;  Build structure of all the variables.
;
  object = create_struct(name='pc_dim_'+strtrim(filename,2),$
      ['mx','my','mz','mw', $
      'mvar','maux','mglobal', $
      'precision', $
      'nx','ny','nz','nw', $
      'nghostx','nghosty','nghostz', $
      'nxgrid','nygrid','nzgrid', $
      'mxgrid','mygrid','mzgrid', $
      'l1','l2','m1','m2','n1','n2', $
      'ipx','ipy','ipz', $
      'nprocx','nprocy','nprocz'], $
      mx,my,mz,mw, $
      mvar,maux,mglobal, $
      precision, $
      nx,ny,nz,nw, $
      nghostx,nghosty,nghostz, $
      nxgrid, nygrid, nzgrid, $
      mxgrid, mygrid, mzgrid, $
      l1,l2,m1,m2,n1,n2, $
      ipx, ipy, ipz, $
      nprocx,nprocy,nprocz)
;
;  Set status of object to "valid".
;
  setenv, 'PC_VALID_DIM=V'
;
;  Print a summary if requested.
;
  if (keyword_set(print)) then begin
    if (size(proc, /type) ne 0) then begin
      print, 'For processor ',strtrim(proc, 2),' calculation domain:'
    end else if (keyword_set(reduced)) then begin
      print, 'For REDUCED calculation domain:'
    end else begin
      print, 'For GLOBAL calculation domain:'
    end
;
    print, '            (mx,my,mz,mw) = (',mx,',',my,',',mz,',',mw,')'
    print, '    (mvar,maux,precision) = (',mvar,',',maux,',',precision,')'
    print, '            (nx,ny,nz,nw) = (',nx,',',ny,',',nz,',',nw,')'
    print, '      (l1:l2,m1:m2,n1:n2) = (',l1,':',l2,',',m1,':',m1,':',m2,',',n1,':',n2,')'
    print, '(nghostx,nghosty,nghostz) = (',nghostx,',',nghosty,',',nghostz,')'
    print, '   (nxgrid,nygrid,nzgrid) = (',nxgrid,',',nygrid,',',nzgrid,')'
    print, '   (mxgrid,mygrid,mzgrid) = (',mxgrid,',',mygrid,',',mzgrid,')'
    print, '   (nprocx,nprocy,nprocz) = (',nprocx,',',nprocy,',',nprocz,')'
  end
;
end
