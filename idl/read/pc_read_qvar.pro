;
; $Id$
;
;   Read qvar.dat, or other QVAR file
;
pro pc_read_qvar, object=object, varfile=varfile_, datadir=datadir, ivar=ivar, $
    quiet=quiet, qquiet=qquiet, SWAP_ENDIAN=SWAP_ENDIAN
COMPILE_OPT IDL2,HIDDEN
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
;  Defaults.
;
datadir = pc_get_datadir(datadir)
default, quiet, 0
default, qquiet, 0

if n_elements(ivar) eq 1 then begin
  default, varfile_, 'QVAR'
  varfile = varfile_ + strcompress (string (ivar), /remove_all)
  if (file_test (datadir+'/allprocs/'+varfile[0]+'.h5')) then varfile += '.h5'
endif else begin
  default_varfile = 'qvar.dat'
  if (file_test (datadir+'/allprocs/qvar.h5')) then default_varfile = 'qvar.h5'
  default, varfile_, default_varfile
  varfile = varfile_
endelse
;
; Load HDF5 varfile if requested or available.
;
  if (strmid (varfile, strlen(varfile)-3) eq '.h5') then begin
    message, "pc_read_qvar: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    t = pc_read ('time', file=varfile, datadir=datadir)
    object = { t:t, number:pc_read ('number') }
    object = create_struct (object, 'mass', pc_read ('points/mass'))
    object = create_struct (object, 'xx', pc_read ('points/'+['x','y','z']))
    object = create_struct (object, 'vv', pc_read ('points/v'+['x','y','z']))
    return
  end
;
if (qquiet) then quiet=1
;
;  Derived dimensions.
;
;
;  Time Get necessary dimensions.
;
t=0d0
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
pc_read_qdim, obj=qdim, datadir=datadir, /quiet
mqvar =qdim.mqvar
nqpar =qdim.nqpar
;mqpar =0L
;
;  Read qvar indices from qvarname.dat
;
datadir = pc_get_datadir(datadir)
openr, lun, datadir+'/qvarname.dat', /get_lun
while (not eof(lun)) do begin
  line=''
  readf, lun, line, format='(a)'
  found = strsplit (line, /extract)
  if (n_elements (found) ge 2) then begin
    line = strtrim (found[1], 2) + '=' + str (found[0])
    if (execute (line) ne 1) then message, 'There was a problem with "qvarname.dat" ('+line+')', /INF
  end
endwhile
close, lun
free_lun, lun
;
;  Define structure for data
;
varcontent=REPLICATE( $
    {varcontent_all_par, $
    variable   : 'UNKNOWN', $
    idlvar     : 'dummy', $
    idlinit    : 'fltarr(nqpar)*one', $
    skip       : 0}, $
    mqvar+1)

INIT_SCALAR  = 'fltarr(nqpar)*one'
INIT_3VECTOR = 'fltarr(nqpar,3)*one'
;
;  Go through all possible particle variables
;
default, ixq, 0
varcontent[ixq].variable = 'Point mass position (xx)'
varcontent[ixq].idlvar   = 'xx'
varcontent[ixq].idlinit  = INIT_3VECTOR
varcontent[ixq].skip     = 2

default, ivxq, 0
varcontent[ivxq].variable = 'Point mass velocity (vv)'
varcontent[ivxq].idlvar   = 'vv'
varcontent[ivxq].idlinit  = INIT_3VECTOR
varcontent[ivxq].skip     = 2

default, imass, 0
varcontent[imass].variable = 'Particle mass (mass)'
varcontent[imass].idlvar   = 'mass'
varcontent[imass].idlinit  = INIT_SCALAR

varcontent[0].variable    = 'UNKNOWN'
varcontent[0].idlvar      = 'UNKNOWN'
varcontent[0].idlinit     = '0.'
varcontent[0].skip        = 0
;
varcontent = varcontent[1:*]
;
;  Put variable names in array
;
variables = (varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)
;
;  Define arrays from contents of varcontent
;
totalvars = mqvar
for iv=0L,totalvars-1L do begin
  if (varcontent[iv].variable eq 'UNKNOWN') then $
      message, 'Unknown variable at position ' + str(iv) $
      + ' needs declaring in pc_read_qvar.pro', /INF
  if (execute(varcontent[iv].idlvar+'='+varcontent[iv].idlinit,0) ne 1) then $
      message, 'Error initialising ' + varcontent[iv].variable $
      +' - '+ varcontent[iv].idlvar, /INFO
  iv=iv+varcontent[iv].skip
endfor
;
;  Define arrays for temporary storage of data.
;

array=fltarr(nqpar,totalvars)*one
;
if (not keyword_set(quiet)) then $
  print,'Loading ',strtrim(datadir+'/proc0/'+varfile), ')...'

filename=datadir+'/proc0/'+varfile 
;
;  Check if file exists.
;
if (not file_test(filename)) then begin
  print, 'ERROR: cannot find file '+ filename
  stop
endif
;
;  Get a unit number and open file.
;
openr, lun, filename, /F77, /get_lun, SWAP_ENDIAN=SWAN_ENDIAN
;
;  Read the number of particles at the local processor together with their
;  global index numbers.
;
readu, lun, nqpar
;
;  Read particle data (if any).
;
if (nqpar ne 0) then begin
;
;  Read local processor data.
;
  array=fltarr(nqpar,mqvar)*one
  readu, lun, array
  readu, lun, t
  print, 't =', t
;
endif
;
close, lun
free_lun, lun
;
;  Put data into sensibly named arrays.
;
for iv=0L,mqvar-1 do begin
  res=varcontent[iv].idlvar+'=array[*,iv:iv+varcontent[iv].skip]'
  if (execute(res,0) ne 1) then $
      message, 'Error putting data into '+varcontent[iv].idlvar+' array'
  iv=iv+varcontent[iv].skip
endfor
;
;  Put data and parameters in object.
;
makeobject="object = CREATE_STRUCT(name=objectname,['t'," + $
    arraytostring(variables,QUOTE="'",/noleader) + "]," + "t,"+$
    arraytostring(variables,/noleader) + ")"
if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /INFO
  undefine,object
endif

end
