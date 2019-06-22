;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_read_phiavg.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   11-Aug-2003
;;;
;;;  Description:
;;;    Read phi-averages from file and return them in a structure
;;;  Usage:
;;;    avg = pc_read_phiavg('data/averages/PHIAVG1')
;;;  Keywords:
;;;    TONLY  -- If true, return just the time t (as a scalar)
;;;  Slots of returned structure:
;;;       t        FLOAT              ; time
;;;       rcyl     FLOAT Array[nr]    ; coordinate
;;;       z        FLOAT Array[nz]    ; coordinate (only for root proc)
;;;       nprocz   LONG               ; number of processors in z
;;;       nvars    LONG               ; number of variables
;;;       <var1>   FLOAT Array[nr,nz] ; first averaged variable
;;;       <var2>   FLOAT Array[nr,nz] ; second averaged variable,
;;;       etc.
;;;    The names <var1>, etc are the same as in the file phiaver.in,
;;;    e.g. `b2mphi', etc.
;;;
;;;  File format of PHIAVG files:
;;     3. nr_phiavg, nz_phiavg, nvars, nprocz
;;;    2. t, r_phiavg, z_phiavg, dr, dz
;;;    1. data
;;;    4. len(labels),labels


function parse_labels, line
;
;  Split label line into the individual variable names.
;
  labels = ['']
  ;
  ; strsplit() is not available in IDL prior to 5.3 (and str_sep() was
  ; obsoleted after 5.2..), so we do this manually:
  ;
  while (line ne '') do begin
    sep = strpos(line,',')
    if (sep ge 0) then begin    ; found a comma
      labels = [labels, strtrim(strmid(line,0,sep))]
      line = strmid(line,sep+1)
    endif else begin            ; no comma -> last label
      labels = [labels, strtrim(line)]
      line = ''
    endelse
  endwhile
  ;
  ;  eliminate empty labels
  ;
  good = where(labels ne '')
  return, labels[good]
end
; ---------------------------------------------------------------------- ;
function pc_read_phiavg, file, datadir=datadir, $
                      VARS=vars, TONLY=t_only, $
                      DEBUG=debug, HELP=help

  if (keyword_set(help)) then extract_help, 'pc_read_phiavg'

  default, debug, 0
  default, t_only, 0
  datadir = pc_get_datadir (datadir)

  ; load HDF5 averages, if available
  pos = long (stregex (file, '([0-9]+)$', /extract)) - 1
  h5_file = datadir + '/averages/phi.h5'
  if ((pos ge 0) and file_test (h5_file)) then begin
    last = pc_read ('last', filename='phi.h5', datadir=datadir+'/averages/')
    if (pos gt last) then message, 'pc_read_phiavg: ERROR: "'+h5_file+'" ends after '+str (last+1)+' snapshots!'
    group = str (pos) + '/'
    time = pc_read (group+'time')
    if (keyword_set (t_only)) then begin
      h5_close_file
      return, time
    end
    message, "pc_read_phiavg: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    if (size (vars, /type) ne 7) then vars = h5_content (group)
    found = where (strlowcase (vars) ne 'time', num)
    if (num le 0) then message, 'pc_read_phiavg: ERROR: "'+h5_file+'" contains no known averages!'
    vars = vars[found]
    pc_read_grid, obj=grid, dim=dim, datadir=datadir, /quiet
    r = pc_read ('r', filename='phi.h5', datadir=datadir+'/averages/')
    dr = pc_read ('dr')
    struct = { t:time, last:last, pos:pos, rcyl:r, dr:dr, z:grid.z[dim.n1:dim.n2], nvars:num, labels:vars }
    for pos = 0, num-1 do begin
      struct = create_struct (struct, vars[pos], pc_read (group+vars[pos]))
    end
    h5_close_file
    return, struct
  end

  if (n_elements(vars) gt 0) then begin
    message, /INFO, $
        'VARS keyword (for selecting only certain vars) not yet implemented'
  endif

  get_lun, lun
  close, lun
  openr, lun, file, /F77

  t = 0.
  if (t_only) then begin
    readu, lun
    readu, lun, t
    close, lun
    free_lun, lun
    return, t
  endif

  nr=1L & nz=1L & nvars=1L & nprocz=1L
  readu, lun, nr, nz, nvars, nprocz
  if (debug) then print,'nr,nz,nvars,nprocz=',nr,nz,nvars,nprocz

  rcyl = fltarr(nr)
  z    = fltarr(nz)
  readu, lun, t, rcyl, z
  if (debug) then begin
    print,'t=',t
    print,'rcyl in ', minmax(rcyl)
    print,'z in '   , minmax(z)
  endif

  vars = fltarr(nr,nz,nvars)
  readu, lun, vars
  if (debug) then print, 'vars in ', minmax(vars)

  point_lun, -lun, position       ; save current file position
  llen = 0L                     ; length of labels line
  readu, lun, llen
  if (debug) then print, 'llen=', llen
  if ((llen le 0) or (llen gt 4096)) then $
      message, "Can't believe I'm excpected to read a string of length " $
      + strtrim(llen,2)
  format = '(A' + strtrim(llen,2)+')'
  lline = string('',FORMAT=format) ; predefine string of correct length
  point_lun, lun, position        ; rewind
  readu, lun, llen, lline
  lline = strtrim(lline)
  if (debug) then print, 'lline: <', lline, '>'

  close, lun
  free_lun, lun

  ;; Cycle through labels and construct structure definition
  labels = parse_labels(lline)
  if (n_elements(labels) ne nvars) then $
      message, 'Inconsistency: found ' + strtrim(n_elements(labels),2) + $
      ' labels, but nvars=' + strtrim(nvars,2)
  def = '{t: t, rcyl: rcyl, z: z, nvars: nvars, labels: labels'
  for i=0, nvars-1 do begin
    def = def + ', ' + labels[i] + ': vars[*,*,'+strtrim(i,2)+']'
  endfor
  def = def + '}'
  if (debug) then print, 'def = ', def
  cmd = 'res = ' + def
  if (not execute(cmd)) then $
      message, "Encountered error trying to execute <" + cmd + ">"
  return, res

end
; End of file pc_read_phiavg.pro
