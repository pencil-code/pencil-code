;;
;; $Id$
;;
;;   Read phiz-averages from file
;;
pro pc_read_phizaver, object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;;
;;  Default data directory.
;;
datadir = pc_get_datadir(datadir)
default, varfile, 'phizaverages.dat'
default, monotone, 0
default, quiet, 0

  ; load HDF5 averages, if available
  h5_file = datadir + '/averages/phi_z.h5'
  if (file_test (h5_file)) then begin
    last = pc_read ('last', filename='phi_z.h5', datadir=datadir+'/averages/')
    groups = str (lindgen (last + 1))
    times = reform (pc_read (groups+'/time'))
    message, "pc_read_phizaver: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    if (size (vars, /type) ne 7) then vars = h5_content (groups[0])
    found = where (strlowcase (vars) ne 'time', num)
    if (num le 0) then message, 'pc_read_phizaver: ERROR: "'+h5_file+'" contains no known averages!'
    vars = vars[found]
    r = pc_read ('r')
    object = { t:times, last:last, pos:long (groups), rcyl:r, nvars:num, labels:vars }
    for pos = 0, num-1 do begin
      object = create_struct (object, vars[pos], pc_read (groups+'/'+vars[pos]))
    end
    h5_close_file
    return
  end

;;
;;  Get necessary dimensions.
;;
pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
;;
;;  Derived dimensions.
;;
nr=round(dim.nx/2.)
;;
;;  Read variables from phizaver.in
;;
spawn, 'cat '+datadir+'/../phizaver.in', varnames
if (not quiet) then print, 'Preparing to read phiz-averages ', $
    arraytostring(varnames,quote="'",/noleader)
nvar=n_elements(varnames)
;;
;;  Define arrays to put data in.
;;
spawn, 'wc -l '+datadir+'/'+varfile, nlines
nlines=long(nlines[0])
nlin_per_time=1L+ceil(nvar*nr/8.)
nlin_rcyl = ceil(nr/8.)
nit=(nlines-nlin_rcyl)/nlin_per_time
if ((nlines-nlin_rcyl) mod nlin_per_time ne 0) then $
    print, 'Warning: File "'+strtrim(datadir+'/'+varfile,2)+'" corrupted!'

if (not quiet) then print, 'Going to read averages at ', strtrim(nit,2), ' times'

for i=0,nvar-1 do begin
  cmd=varnames[i]+'=fltarr(nr,nit)*one'
  if (execute(cmd,0) ne 1) then message, 'Error defining data arrays'
endfor

rcyl=fltarr(nr)*one
var =fltarr(nr)*one
tt  =fltarr(nit)*one

;;
;;  Prepare for read
;;
filename=datadir+'/'+varfile 
if (not quiet) then print, 'Reading ', filename
if (not file_test(filename)) then begin
  print, 'ERROR: cannot find file '+ filename
  stop
endif
openr, lun, filename, /get_lun

;;
;;  Read phiz-averages and put in arrays.
;;

;; Read radius (first ceil(nr/8) records)
readf, lun, rcyl
for it=0,nit-1 do begin
;; Read time
  readf, lun, t
  tt[it]=t
;; Read data
  for i=0,nvar-1 do begin
    readf, lun, var
    cmd=varnames[i]+'[*,it]=var'
    if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
  endfor
endfor
close, lun
free_lun, lun
;;
;;  Make time monotonous and crop all variables accordingly.
;;  
if (monotone) then begin
  ii=monotone_array(tt)
endif else begin
  ii=lindgen(n_elements(tt))
endelse
;;
;;  Put data in structure.
;;
makeobject="object = CREATE_STRUCT(name=objectname,['t','rcyl'," + $
    arraytostring(varnames,QUOTE="'",/noleader) + "]," + $
    "tt[ii],rcyl[*],"+arraytostring(varnames+'[*,ii]',/noleader) + ")"

if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /INFO
  undefine,object
endif


end
