;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   read_phiavg.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   11-Aug-2003
;;;
;;;  Description:
;;;   Read phi-averages from file and return them in a structure
;;;  File format:
;;     3. nr_phiavg, nz_phiavg, nvars, nprocz
;;;    2. t, r_phiavg, z_phiavg, dr, dz
;;;    1. data
;;;    4. labels
;;;  Slots of returned structure:
;;;       t        FLOAT              ; time
;;;       rcyl     FLOAT Array[nr]    ; coordinate
;;;       z        FLOAT Array[nz]    ; coordinate (only for root proc)
;;;       nprocz   LONG               ; number of processors in z
;;;       nvars    LONG               ; number of variables
;;;       <var1>   FLOAT Array[nr,nz] ; first averaged variable
;;;       <var2>   FLOAT Array[nr,nz] ; second averaged variable,
;;;       etc.
;;;  The names <var1>, etc are the same as in the file phiaver.in,
;;;  e.g. `b2mphi', etc.


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
function read_phiavg, file, $
                      VARS=vars, DEBUG=debug

  default, debug, 0

  if (n_elements(vars) gt 0) then begin
    message, /INFO, $
        'VARS keyword (for selecting only certain vars) not yet implemented'
  endif

  close, 1
  openr, 1, file, /F77

  nr=1L & nz=1L & nvars=1L & nprocz=1L
  readu, 1, nr, nz, nvars, nprocz
  if (debug) then print,'nr,nz,nvars,nprocz=',nr,nz,nvars,nprocz

  t = 0.
  rcyl = fltarr(nr)
  z    = fltarr(nz)
  readu, 1, t, rcyl, z
  if (debug) then begin
    print,'t=',t
    print,'rcyl in ', minmax(rcyl)
    print,'z in '   , minmax(z)
  endif

  vars = fltarr(nr,nz*nprocz,nvars)
  readu, 1, vars
  if (debug) then print, 'vars in ', minmax(vars)

  point_lun, -1, position       ; save current file position
  llen = 0L                     ; length of labels line
  readu, 1, llen
  if (debug) then print, 'llen=', llen
  if ((llen le 0) or (llen gt 4096)) then $
      message, "Can't believe I'm excpected to read a string of length " $
      + strtrim(llen,2)
  format = '(A' + strtrim(llen,2)+')'
  lline = string('',FORMAT=format) ; predefine string of correct length
  point_lun, 1, position        ; rewind
  readu, 1, llen, lline
  lline = strtrim(lline)
  if (debug) then print, 'lline: <', lline, '>'

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
; End of file read_phiavg.pro
