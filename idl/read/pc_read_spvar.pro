;
; $Id$
;
;   Read spvar.dat, or other PVAR file
;
pro pc_read_spvar, object=object, varfile=varfile_, datadir=datadir, ivar=ivar, $
    quiet=quiet, qquiet=qquiet,SWAP_ENDIAN=SWAP_ENDIAN
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
;  Defaults.
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
default, quiet, 0
default, qquiet, 0

if n_elements(ivar) eq 1 then begin
    default,varfile_,'SPVAR'
    varfile=varfile_+strcompress(string(ivar),/remove_all)
endif else begin
    default,varfile_,'spvar.dat'
    varfile=varfile_
endelse

if (qquiet) then quiet=1
;
;  Derived dimensions.
;
;
;  Get necessary dimensions.
;

pc_read_dim, obj=dim, datadir=datadir, /quiet
pc_read_spdim, obj=spdim, datadir=datadir, /quiet
pc_set_precision,dim=dim,/quiet
mspvar=spdim.mspvar
nspar =spdim.nspar
mspar =0L
;
;  Read variable indices from index.pro
;
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
openr, 1, datadir+'/index.pro'
line=''
while ~ eof(1) do begin
  readf, 1, line, format='(a)'
  if (execute(line) ne 1) then $
    message, 'There was a problem with index.pro', /INF
endwhile
close, 1
;
;  Define structure for data
;
varcontent=REPLICATE( $
    {varcontent_all_par, $
    variable   : 'UNKNOWN', $
    idlvar     : 'dummy', $
    idlinit    : 'fltarr(nspar)*one', $
    skip       : 0}, $
    mspvar+1)

INIT_SCALAR  = 'fltarr(nspar)*one'
INIT_3VECTOR = 'fltarr(nspar,3)*one'
;
;  Go through all possible particle variables
;
varcontent[ixp].variable = 'Sink particle position (xxs)'
varcontent[ixp].idlvar   = 'xxs'
varcontent[ixp].idlinit  = INIT_3VECTOR
varcontent[ixp].skip     = 2

varcontent[ivpx].variable = 'Sink particle velocity (vvs)'
varcontent[ivpx].idlvar   = 'vvs'
varcontent[ivpx].idlinit  = INIT_3VECTOR
varcontent[ivpx].skip     = 2

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
totalvars = mspvar
for iv=0L,totalvars-1L do begin
  if (varcontent[iv].variable eq 'UNKNOWN') then $
      message, 'Unknown variable at position ' + str(iv) $
      + ' needs declaring in pc_read_psvar.pro', /INF
  if (execute(varcontent[iv].idlvar+'='+varcontent[iv].idlinit,0) ne 1) then $
      message, 'Error initialising ' + varcontent[iv].variable $
      +' - '+ varcontent[iv].idlvar, /INFO
  iv=iv+varcontent[iv].skip
endfor
;
;  Define arrays for temporary storage of data.
;

array=fltarr(nspar,totalvars)*one
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
get_lun, file
close, file
openr, file, filename, /F77,SWAP_ENDIAN=SWAN_ENDIAN
;
;  Read the number of particles at the local processor together with their
;  global index numbers.
;
readu, file, mspar
;
;  Read particle data (if any).
;
if (mspar ne 0) then begin
;
;  Read local processor data.
;
    array=fltarr(mspar,mspvar)*one
    readu, file, array
;
endif
;
close, file
free_lun, file
;
;  Put data into sensibly named arrays.
;
for iv=1L,mspvar do begin
  res=varcontent[iv].idlvar+'=array[*,iv-1:iv-1+varcontent[iv].skip]'
  if (execute(res,0) ne 1) then $
    message, 'Error putting data into '+varcontent[iv].idlvar+' array'
  iv=iv+varcontent[iv].skip
endfor
;
;  Put data and parameters in object.
;
makeobject="object = CREATE_STRUCT(name=objectname,[" + $
     arraytostring(variables,QUOTE="'",/noleader) + "]," + $
     arraytostring(variables,/noleader) + ")"
if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /INFO
  undefine,object
endif

end
