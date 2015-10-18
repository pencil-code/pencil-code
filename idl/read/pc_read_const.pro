;;
;; $Id$
;;
;;  Read Pencil Code constants from a file.
;;
;;  Usage:
;;
;;    IDL> pc_read_const, obj=const
;;    IDL> print, const.tt_ion
;;
pro pc_read_const, object=object, varfile=varfile, datadir=datadir, quiet=quiet, specname=specname,specmass=specmass
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
;  Default data directory.
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
default, varfile, 'pc_constants.pro'
default, quiet, 0
fullvarfile=datadir+'/'+varfile
if (not quiet) then $
    print, 'Going to read Pencil Code constants from the file ' $
    + fullvarfile + ' .'
;
;  Find length of data file and define data arrays.
;
nlines = file_lines (fullvarfile)
if (not quiet) then $
    print, 'Data file contains ' + nlines,2) + ' lines.'
array=strarr(nlines)
varnames=strarr(nlines)
values=strarr(nlines)
;
;  Read data file into array.
;
get_lun, file
openr, file, datadir+'/'+varfile
  readf, file, array
close, file
free_lun, file
;
;  Split lines into variable name and associate value.
;
ncom=0
for i=0, nlines-1 do begin
;  Skip comment lines.
  com = strpos(array[i], ';')
  spec = strpos(array[i],'spec')
  if ((com ne 1) and (spec ne 1)) then begin
    pos = strpos(array[i], '=')
    varnames[i] = strmid(array[i], 0, pos)
;  Remove blanks from variable name.
    while (strpos(varnames[i],' ') ne -1) do begin
      j = strpos(varnames[i],' ')
      varnames[i] = strmid(varnames[i], 0, j-1) + $
          strmid(varnames[i], j+1, strlen(varnames[i]))
    endwhile
;  
    values[i] = strmid(array[i], pos+1, strlen(array[i]))
  endif else begin
    if (spec eq 1) then begin
        test=execute(array[i])    
        ncom=ncom+1
    endif else begin
        ncom=ncom+1
    endelse
endelse
endfor
;
;  Put data in structure.
;
makeobject="object = create_struct(name=objectname,[" + $
    arraytostring(varnames[ncom:nlines-1],quote="'",/noleader) + "]," + $
    arraytostring(values[ncom:nlines-1],/noleader) + ")"
;
if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /INFO
  undefine,object
endif
;
end
