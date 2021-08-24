;
; $Id: pc_read_1d_aver.pro 23239 2015-03-26 20:29:51Z mreinhardt@nordita.org $
;
;  Read 1d-averages from file.
;
pro pc_read_1d_aver, dir, object=object, varfile=varfile, datadir=datadir, $
    monotone=monotone, quiet=quiet, njump=njump, tmin=tmin, dim=dim, grid=grid, single=single
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
;  Get necessary dimensions.
;
if (not is_defined(dim)) then pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
if (not is_defined(grid)) then pc_read_grid, obj=grid, dim=dim, datadir=datadir, quiet=quiet, single=single
;
;  Set pc_precision.
;
  if (not is_defined(precision)) then pc_set_precision, dim=dim, quiet=quiet
;
;  Default data directory.
;
datadir = pc_get_datadir(datadir)

if (dir eq 'z') then begin
  ndir = dim.nz
  avdirs = 'xy'
end else if (dir eq 'y') then begin
  ndir = dim.ny
  avdirs = 'xz'
end else if (dir eq 'x') then begin
  ndir = dim.nx
  avdirs = 'yz'
end else $
  message, 'ERROR: unknown direction "'+dir+'"!'

default, in_file, avdirs+'aver.in'
default, varfile, avdirs+'averages.dat'
default, monotone, 0
default, quiet, 0
default, njump, 1
default, tmin, 0.0
default, single, 0

  ; load HDF5 averages, if available
  if (file_test (datadir+'/averages/'+avdirs+'.h5')) then begin
    message, "pc_read_1d_aver: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info

    last = pc_read ('last', filename=avdirs+'.h5', datadir=datadir+'/averages/')
    groups = str (lindgen (last + 1))
    times = pc_read (groups+'/time', single=single)
    in = (where (times ge tmin))[0:*:njump]
    num = n_elements (in)
    if (num le 0) then message, 'pc_read_1d_aver: ERROR: "'+h5_file+'" no data available after given "tmin"!'
    groups = groups[in]
    times = times[in]
    if (size (vars, /type) ne 7) then vars = h5_content (groups[0])
    found = where (strlowcase (vars) ne 'time', num)
    vars = vars[found]
    object = { t:times, last:last, pos:1+long (groups), num_quantities:num, labels:vars }
    if (dir eq 'z') then begin
      object = create_struct (object, 'z', grid.z[dim.n1:dim.n2])
    end else if (dir eq 'y') then begin
      object = create_struct (object, 'y', grid.y[dim.m1:dim.m2])
    end else begin
      object = create_struct (object, 'x', grid.x[dim.l1:dim.l2])
    end
    for pos = 0, num-1 do begin
      object = create_struct (object, vars[pos], pc_read (groups+'/'+vars[pos],single=single))
    end
    h5_close_file
    return
  end

;
;  Read variables from '*aver.in' file
;
run_dir = stregex ('./'+datadir, '^(.*)data\/', /extract)
varnames = strarr(file_lines(run_dir+in_file))
openr, lun, run_dir+in_file, /get_lun
readf, lun, varnames
close, lun
free_lun, lun
;
; Remove commented and empty elements from allvariables
;
varnames = strtrim (varnames, 2)
inds = where (stregex (varnames, '^[a-zA-Z]', /boolean), nvar)
if (nvar le 0) then message, "ERROR: there are no variables found."
varnames = varnames[inds]
;
;  Check for existence of data file.
;
filename=datadir+'/'+varfile
if (not quiet) then begin
  print, 'Preparing to read '+avdirs+'-averages ', arraytostring(varnames,quote="'",/noleader)
  print, 'Reading ', filename
endif
if (not file_test(filename)) then $
  message, 'ERROR: cannot find file '+ filename
;
;  Define arrays to put data in.
;
nlines=file_lines(filename)
nlin_per_time=1L+ceil(nvar*ndir/8.)
nit=nlines/nlin_per_time/njump
if ((nlines mod nlin_per_time) ne 0) then begin
  print, 'Warning: File "'+strtrim(filename,2)+'" corrupted!'
  corrupt=1
endif else $
  corrupt=0
;
if (not quiet) then print, 'Going to read averages at <=', strtrim(nit,2), ' times'
;
;  Generate command name. Note that an empty line in the *aver.in
;  file will lead to problems. If this happened, you may want to replace
;  the empty line by a non-empty line rather than nothing, so you can
;  read the data with idl.
;
for i=0,nvar-1 do begin
  cmd=varnames[i]+(single ? '=fltarr(ndir,nit)' : '=fltarr(ndir,nit)*one')
  if (execute(cmd,0) ne 1) then message, 'Error defining data arrays'
endfor
var=fltarr(ndir*nvar)*one
times = single ? fltarr(nit) : fltarr(nit)*one
;
;  Read averages and put in arrays.
;
openr, file, filename, /get_lun

nread=0 & dummy=zero & t=zero & line=0
for it=0,nit-1 do begin
  ; Read time
  readf, file, t
  line++
  if corrupt and nread gt 0 then $
    if t lt times[nread-1] then begin
      print, 'Warning: File corrupt before or in line', line
      corrupt=0
    endif
 
  ; Read data
  readf, file, var
  line+=nlin_per_time-1
  if (it mod njump) eq 0 and t ge tmin then begin
    times[nread]=t
    for i=0,nvar-1 do begin
      cmd=varnames[i]+'[*,nread]=var[i*ndir:(i+1)*ndir-1]'
      if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
    endfor
    nread++
  endif
endfor
if nread gt 0 then begin
  times=times[0:nread-1]
  for i=0,nvar-1 do begin
    cmd=varnames[i]+'='+varnames[i]+'[*,0:nread-1]'
    if (execute(cmd,0) ne 1) then message, 'Error trimming data'
  endfor
endif
close, file
free_lun, file
;
;  Make time monotonous and prepare cropping all variables accordingly.
;
if (monotone) then $
  ii=monotone_array(times) $
else $
  ii=lindgen(n_elements(times))
;
;  Trim grid.
;
trimmed_x = grid.x[dim.l1:dim.l2]
trimmed_y = grid.y[dim.m1:dim.m2]
trimmed_z = grid.z[dim.n1:dim.n2]
;
;  Put data in structure.
;
makeobject="object = create_struct(name=objectname,['t','"+dir+"'," + $
    arraytostring(varnames,quote="'",/noleader) + "]," + $
    "times[ii],trimmed_"+dir+","+arraytostring(varnames+'[*,ii]',/noleader) + ")"
;
if (execute(makeobject) ne 1) then begin
  message, 'Error making object: ' + makeobject, /info
  print, 'Variables locally available as times[ii], '+arraytostring(varnames+'[*,ii]',/noleader)+'.'
  undefine,object
  stop
endif
;
end
