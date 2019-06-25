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
if nlines-nlin_rcyl mod nlin_per_time ne 0 then $
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
