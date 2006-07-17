;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   dx_export.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  $Date: 2006-07-17 05:14:59 $
;;;  $Revision: 1.1 $
;;;
;;;  Description:
;;;   Write data into a file such that the DX (IBM Open Data Explorer)
;;;   Import function can import them and write a corresponding header
;;;   file.
;;;
;;;  Usage:
;;;    dx_export, f [, f2[, f3[, f4[, f5[, f6[, f7]]]]]] [,options]
;;;
;;;  Options:
;;;    X=x, Y=y, Z=z -- specify coordinates; assumes equidistant grid
;;;                     at the moment
;;;    T=t           -- document time in a comment
;;;    COMMENT=comm  -- A comment to be included in the .general file
;;;    BASENAME=name -- specify names of data file (<name.dx>) and
;;;                     header file (<name.general>); defaults are
;;;                     `var.dx' and `var.general'
;;;    FILE=file     -- specify name of data file to write (default `var.dx')
;;;    GENFILE=file  -- specify name of header file (default `var.general')
;;;    LABELS=labels -- specify names of the variables
;;;    GHOSTS=ghosts -- number of ghost cells to ignore (scalar ng for
;;;                     all directions, or [ngx,ngy,ngz])
;;;
;;;  Examples:
;;;    dx_export, uu, lam, bb, s
;;;    dx_export, uu, lam, bb, s, X=x, Y=y, Z=z, T=t
;;;    dx_export, uu, bb, FILE='uu-BB.dx', GENFILE='uu-BB.general', $
;;;                       LABELS=['uu','BB']
;;;    dx_export, uu, bb, BASE='uu-BB', LABELS=['uu','BB'], GHOSTS=3
;;;
;;;  Note: data written here are in DX terminology (on Linux and DEC):
;;;    binary; LSF (=little endian); Column order;
;;;    `x0,x1,x2,..,y0,y1,y2,..'; float (scalar or 3-vector)
;;;
;;;  To do:
;;;    - check for byte order (little vs. big endian) and adapt the
;;;      `format' entry accordingly
;;;    - generalise for 2-d data
;;;    - GHOSTS is currently not recorded in the `approximate cmd line'
;;;    - Record the current working directory in a comment
;;;

pro dx_export, f, f2, f3, f4, f5, f6, f7, $
               X=x, Y=y, Z=z, $
               T=t, COMMENT=comment, $
               BASENAME=bname, FILE=fname, GENFILE=genfname, $
               LABELS=labels, $
               GHOSTS=ghosts, $
               HELP=help
  if (keyword_set(help)) then begin
    print, "Usage:"
    print, "  dx_export, f [, f2[, f3[, f4[, f5[, f6[, f7]]]]]] [,options]"
    print, ""
    print, "  Options:"
    print, "  X=x, Y=y, Z=z -- specify coordinates; assumes equidistant grid at the moment"
    print, "  T=t           -- document time in a comment"
    print, "  COMMENT=comm  -- A comment to be included in the .general file
    print, "  BASENAME=name -- specify names of data file (<name.dx>) and header file"
    print, "  (<name.general>); defaults are `var.dx' and `var.general'"
    print, "  FILE=file     -- specify name of data file to write (default `var.dx')"
    print, "  GENFILE=file  -- specify name of header file (default `var.general')"
    print, "  LABELS=labels -- specify names of the variables"
    print, ""
    print, "Examples:"
    print, "  dx_export, uu, lam, bb, s"
    print, "  dx_export, uu, lam, bb, s, COMMENT='Run Magnetic2e'"
    print, "  dx_export, uu, lam, bb, s, X=x, Y=y, Z=z, T=t"
    print, "  dx_export, uu, bb, FILE='uu-BB.dx', GENFILE='uu-BB.general', $"
    print, "                     LABELS=['uu','BB']"
    print, "  dx_export, uu, bb, BASE='uu-BB', LABELS=['uu','BB'], GHOSTS=3"
    print, ""
    goto, done
  endif

  ;; Assumes 3-d data (could be generalised)

  ;; Initialize option strings (for reconstructing the command line)
  bnopt = ""
  fileopt = ""
  genopt = ""
  labopt = ""
  xopt = ""
  yopt = ""
  zopt = ""
  topt = ""
  ghopt = ""
  commopt = ""

  ;; Isolate version number from CVS string
  version = "$Revision: 1.1 $"
  version = strmid(version, 10, strlen(version)-11)
  version = strtrim(version, 2)

  ;; Construct file names
  ibasen = 0                    ; flag for presence of BASENAME option
  if (n_elements(bname) ne 0) then begin
    ibasen = 1
    bnopt = ", BASENAME='" + bname +"'" 
  endif
  ;
  if (n_elements(fname) eq 0) then begin
    if (ibasen) then fname = bname+'.dx' else fname = 'var.dx'
  endif else begin
    fileopt=", FILE='" + fname +"'"
  endelse
  ;
  if (n_elements(genfname) eq 0) then begin
    if (ibasen) then genfname = bname+'.general' else genfname = 'var.general'
  endif else begin
    genopt = ", GENFILE='" + genfname +"'"
  endelse
  ;
  if (n_elements(labels) eq 0) then begin
    labels=['uu','lam','bb','s','var5','var6','var7','var8']
  endif else begin
    labopt = ", LABELS=['" + labels[0] + "'"
    for i=1,n_elements(labels)-1 do begin
      labopt = labopt + ", '" + labels[i] + "'"
    endfor
    labopt = labopt + "]"
  endelse
  ;
  if (n_elements(ghosts) le 0) then ngh=0 else ngh=ghosts
  if (n_elements(ngh) eq 1) then ngh = [1,1,1]*ngh

  ;; Fill the unset arguments by (scalar) zeros (quite unelegant, but
  ;; without pointers there is hardly a better way):
  if (n_elements(f2) eq 0) then f2 = 0
  if (n_elements(f3) eq 0) then f3 = 0
  if (n_elements(f4) eq 0) then f4 = 0
  if (n_elements(f5) eq 0) then f5 = 0
  if (n_elements(f6) eq 0) then f6 = 0
  if (n_elements(f7) eq 0) then f7 = 0
  N_max_args = 7

  if (n_elements(z) le 0) then begin  ; No(t enough) coords specified
    x = [0,1.] & y=x & z=x
  endif else begin
    xopt = ", X=x"
    yopt = ", Y=y"
    zopt = ", Z=z"
  endelse
  dx = x[1]-x[0]
  dy = y[1]-y[0]
  dz = z[1]-z[0]

  if (n_elements(t) gt 0) then topt = ", T=t"
  if (any(ngh gt 0)) then begin
    ghopt = ", GHOSTS="
    if (n_elements(ghosts) eq 1) then begin
      ghopt = ghopt + strtrim(ngh,2)
    endif else begin
      ghopt = ghopt + "[" + strtrim(ngh[0],2) $
                    + "," + strtrim(ngh[1],2) $
                    + "," + strtrim(ngh[0],2) $
                    + "]"
    endelse
  endif
  if (n_elements(comment) gt 0) then commopt = ", COMMENT='" + comment + "'"

  s = size(f)
  if ((s[0] lt 3) or (s[0] gt 4)) then $
      message, 'DX_EXPORT: First argument has wrong dimensions'

  mx = s[1] & my = s[2] & mz = s[3]
  nx = mx-2*ngh[0]              ; size without ghost zones
  ny = my-2*ngh[1]
  nz = mz-2*ngh[2]
  l1 = ngh[0]  &  l2 = l1+nx-1  ; first and last x index
  m1 = ngh[1]  &  m2 = m1+ny-1  ; first and last y index
  n1 = ngh[2]  &  n2 = n1+nz-1  ; first and last z index

  success = 0

  message, /RESET_ERROR_STATE   ; Clean error messages
  on_ioerror, io_err            ; Catch errors
  file = 1                      ; Just to make sure FILE is defined
  openw, file, /GET_LUN, fname
  writeu, file, f[l1:l2,m1:m2,n1:n2,*]
  if (s[0] lt 4) then dims = [1] else dims = [s[4]]
  ;
  nvars = 1                   ; Count number of variables
  for i=2,N_max_args do begin
    case i of                   ; Fill the corresponding argument into var
      2: var = f2
      3: var = f3
      4: var = f4
      5: var = f5
      6: var = f6
      7: var = f7
    endcase
    if (n_elements(var) gt 1) then begin
      writeu, file, var[l1:l2,m1:m2,n1:n2,*]
      s = size(var)
      if (s[0] lt 4) then dims = [dims, 1] else dims = [dims, s[4]]
      nvars = nvars + 1
    endif
  endfor
  ;
  close, file
  success = 1

  io_err:
  on_ioerror, NULL              ; Switch error handling off
  close, file
  free_lun, file
  if (not success) then begin
    print
    print, 'dx_export: Problems ocurred when writing to file ', fname
    print, !error_state.msg
  endif

  ;; Reconstruct the command line (we do not knwo the explicit
  ;; arguments, so we just use the names from LABELS)
  cmdline = "dx_export"
  for i=0,nvars-1 do begin
    cmdline = cmdline + ", " + labels[i]
  endfor
  cmdline = cmdline $
      + xopt + yopt + zopt $
      + topt + commopt $
      + bnopt + fileopt + genopt + labopt

  ;; Get user and host name (Unix-specific, but why bother)
  spawn, "whoami", uid
  spawn, "uname -n", hostn
  spawn, "pwd", cwd

  ;; Write the Import file

  ;; Use this to map DIMS onto `structure' field:
  type = ['unknown', 'scalar', 'unknown', '3-vector', 'unknown']
  on_ioerror, io_err2            ; Catch errors
  openw, file, /GET_LUN, genfname
  ;; Print header with command line
  printf, file, FORMAT='("# Creator: dx_export.pro, version ",A)', version
  printf, file, FORMAT='("# Date: ",A)', systime()
  printf, file, FORMAT='("# For: ",A,"@",A)', $
      strtrim(uid[0],2), strtrim(hostn[0],2)
  printf, file, FORMAT='("# Directory: ",A)', cwd
  printf, file, FORMAT='("# Approximate command line:")'
  printf, file, FORMAT='("#   ",A)', cmdline
  if (n_elements(t) gt 0) then $
      printf, file, FORMAT='("# t = ",A)', strtrim(t,2) 
  if (n_elements(comment) gt 0) then $
      printf, file, FORMAT='("# ",A)', comment

  ;; Print data information
  printf, file, FORMAT='("file = ",A)', fname
  printf, file, FORMAT='("grid = ",I0," x ",I0," x ",I0)', nx, ny, nz
  printf, file, FORMAT='("format = ",A," ",A)', 'lsb', 'ieee'
  printf, file, FORMAT='("interleaving = ",A)', 'record'
  printf, file, FORMAT='("majority = ",A)', 'column'
  printf, file, FORMAT='("field = ",8(A,:,", "))', labels[0:nvars-1]
  printf, file, FORMAT='("structure = ",8(A,:,", "))', type[dims]
  printf, file, FORMAT='("type = ",8(A,:,", "))', replicate('float',nvars)
  printf, file, FORMAT='("dependency = ",8(A,:,", "))', $
      replicate('positions',nvars)
  printf, file, $
      FORMAT='("positions = ",A,", ",A,", ",A,", ",G0,", ",G0,", ",G0,", ",G0,", ",G0,", ",G0)', $
      'regular', 'regular', 'regular', $
      x[ngh[0]], dx, y[ngh[1]], dy, z[ngh[2]], dz
  printf, file, FORMAT='("")'
  printf, file, FORMAT='("end")'



  io_err2:
  on_ioerror, NULL              ; Switch error handling off
  close, file
  free_lun, file
  if (not success) then begin
    print
    print, 'dx_export: Problems ocurred when writing to file ', genfname
    print, !error_state.msg
  endif

  print,'DIMS = ', dims

done:

end

; End of file dx_export.pro
