;;
;; $Id: pc_read_1daver.pro 23239 2015-03-26 20:29:51Z mreinhardt@nordita.org $
;;
;;  Read 1d-averages from file.
;;
pro pc_read_1daver, dir, object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;;
;;  Default data directory.
;;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()

if dir eq 'z' then $
  avdirs='xy' $
else if dir eq 'x' then $
  avdirs='yz' $
else $
  avdirs='xz'

default, varfile, avdirs+'averages.dat'
default, monotone, 0
default, quiet, 0
;;
;;  Get necessary dimensions.
;;
pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
pc_set_precision, dim=dim, quiet=quiet
cmd = 'ndir=dim.n'+dir
dummy = execute(cmd,0)
;;
;;  Read variables from *aver.in
;;
spawn, "echo "+datadir+" | sed -e 's/data\/*$//g'", datatopdir
spawn, 'cat '+datatopdir+'/'+avdirs+'aver.in', varnames

inds = where(varnames ne '')
if inds[0] eq -1 then begin
  print, 'PC_READ_XYAVER: No variables provided!'
  return
endif else $
  varnames = varnames[inds]

if (not quiet) then print, 'Preparing to read '+avdirs+'-averages ', $
    arraytostring(varnames,quote="'",/noleader)
nvar=n_elements(varnames)
;;
;;  Check for existence of data file.
;;
get_lun, file
filename=datadir+'/'+varfile
if (not quiet) then print, 'Reading ', filename
if (not file_test(filename)) then begin
  print, 'ERROR: cannot find file '+ filename
  stop
endif
close, file
openr, file, filename
;;
;;  Define arrays to put data in.
;;
spawn, 'wc -l '+filename, nlines
nlines=long(nlines[0])
nlin_per_time=1L+ceil(nvar*ndir/8.)
nit=nlines/nlin_per_time
if nlines mod nlin_per_time ne 0 then $
  print, 'Warning: File "'+strtrim(filename,2)+'" corrupted!' 
;
if (not quiet) then print, 'Going to read averages at ', strtrim(nit,2), ' times'
;
;  Generate command name. Note that an empty line in the *aver.in
;  file will lead to problems. If this happened, you may want to replace
;  the empty line by a non-empty line rather than nothing, so you can
;  read the data with idl.
;
for i=0,nvar-1 do begin
  cmd=varnames[i]+'=fltarr(ndir,nit)*one'
  if (execute(cmd,0) ne 1) then message, 'Error defining data arrays'
endfor
var=fltarr(ndir*nvar)*one
times =fltarr(nit)*one
;;
;;  Read averages and put in arrays.
;;
for it=0,nit-1 do begin
;; Read time
  readf, file, t
  times[it]=t
;; Read data
  readf, file, var
  for i=0,nvar-1 do begin
    cmd=varnames[i]+'[*,it]=var[i*ndir:(i+1)*ndir-1]'
    if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
  endfor
endfor
;;
;;  Close file.
;;
close, file
free_lun, file
;;
;;  Make time monotonous and crop all variables accordingly.
;;
if (monotone) then begin
  ii=monotone_array(times)
endif else begin
  ii=lindgen(n_elements(times))
endelse
;;
;;  Read grid.
;;
pc_read_grid, obj=grid, /trim, datadir=datadir, /quiet
;;
;;  Put data in structure.
;;
makeobject="object = create_struct(name=objectname,['t','"+dir+"'," + $
    arraytostring(varnames,quote="'",/noleader) + "]," + $
    "times[ii],grid."+dir+","+arraytostring(varnames+'[*,ii]',/noleader) + ")"
;
if (execute(makeobject) ne 1) then begin
  message, 'Error evaluating variables: ' + makeobject, /info
  undefine,object
endif
;
end
