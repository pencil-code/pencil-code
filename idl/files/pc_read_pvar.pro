; $Id: pc_read_pvar.pro,v 1.7 2005-02-09 15:08:09 ajohan Exp $
;
;   Read pvar.dat, or other PVAR file
;
pro pc_read_pvar, object=object, varfile=varfile, datadir=datadir, $
    QUIET=QUIET
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
;  Default data directory.
;
default, datadir, 'data'
default, varfile, 'pvar.dat'
;
;  Get necessary dimensions.
;
pc_read_dim, obj=dim, datadir=datadir, QUIET=QUIET
pc_read_pdim, obj=pdim, datadir=datadir, QUIET=QUIET
pc_set_precision, dim=dim, QUIET=QUIET
;
;  Derived dimensions.
;
mpvar=pdim.mpvar
npar=pdim.npar
nproc=dim.nprocx*dim.nprocy*dim.nprocz
;
;  Read variable indices from index.pro
;
default,datadir,'data'
cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),     "\n"'+"' "+datadir+'/index.pro'
spawn, cmd, result
res = flatten_strings(result)
if (execute(res) ne 1) then $
    message, 'There was a problem with index.pro', /INF
;
;  Define structure for data
;
varcontent=REPLICATE( $
    {varcontent_all_par, $
    variable   : 'UNKNOWN', $
    idlvar     : 'dummy', $
    idlinit    : 'fltarr(npar)*one', $
    skip       : 0}, $
    mpvar+1)

INIT_SCALAR  = 'fltarr(npar)*one'
INIT_3VECTOR = 'fltarr(npar,3)*one'
;
;  Go through all possible particle variables
;
varcontent[ixxp].variable = 'Particle position (xx)'
varcontent[ixxp].idlvar   = 'xx'
varcontent[ixxp].idlinit  = INIT_3VECTOR
varcontent[ixxp].skip     = 2

varcontent[ivvp].variable = 'Particle velocity (vv)'
varcontent[ivvp].idlvar   = 'vv'
varcontent[ivvp].idlinit  = INIT_3VECTOR
varcontent[ivvp].skip     = 2

varcontent[0].variable    = 'UNKNOWN'
varcontent[0].idlvar      = 'UNKNOWN'
varcontent[0].idlinit     = '0.'
varcontent[0].skip        = 0
;
;  Put variable names in array
;
variables = (varcontent[where((varcontent[*].idlvar ne 'dummy'))].idlvar)[1:*]
;
;  Define arrays from contents of varcontent
;
totalvars = mpvar
for iv=1L,totalvars do begin
  if (varcontent[iv].variable eq 'UNKNOWN') then $
      message, 'Unknown variable at position ' + str(iv) $
      + ' needs declaring in pc_read_pvar.pro', /INF
  if (execute(varcontent[iv].idlvar+'='+varcontent[iv].idlinit,0) ne 1) then $
      message, 'Error initialising ' + varcontent[iv].variable $
      +' - '+ varcontent[iv].idlvar, /INFO
  iv=iv+varcontent[iv].skip
endfor
;
;  Define arrays for temporary storage of data.
;
array=fltarr(npar,totalvars)*one
ipar0=lonarr(npar)
t=zero
npar_loc=0L
;
;  Get a unit number.
;
GET_LUN, file
;
;  Loop over processors.
;
for i=0,nproc-1 do begin

  filename=datadir+'/proc'+strtrim(i,2)+'/'+varfile 
  if (not keyword_set(QUIET)) then print, 'Reading ', filename
;
;  Check if file exists.
;
  dummy=findfile(filename, COUNT=countfile)
  if (not countfile gt 0) then begin
    print, 'ERROR: cannot find file '+ filename
    stop
  endif

  close, file
  openr, file, filename, /F77
;
;  Read the number of particles at the local processor together with their
;  global index numbers.
;
  readu, file, npar_loc 
  ipar_loc=lonarr(npar_loc)
  readu, file, ipar_loc
;
;  Register particle indices for later check if all particles have been read.
;  
  for k=0,npar_loc-1 do begin
    ipar0[ipar_loc[k]-1]=ipar0[ipar_loc[k]-1]+1
  endfor
;
;  Read local processor data.
;
  array_loc=fltarr(npar_loc,mpvar)
  readu, file, array_loc
;
;  Put local processor data into proper place in global data array
;        
  for k=0,npar_loc-1 do begin
    array[ipar_loc[k]-1,*]=array_loc[k,*]
  endfor
;
;  Read time.
;
  readu, file, t

  close, file
  free_lun, file

endfor
;
;  Put data into object structure.
;
for iv=1L,mpvar do begin
  res=varcontent[iv].idlvar+'=array[*,iv-1:iv-1+varcontent[iv].skip]'
  if (execute(res,0) ne 1) then $
    message, 'Error putting data into '+varcontent[iv].idlvar+' array'
  iv=iv+varcontent[iv].skip
endfor
;
;  Check if all particles found exactly once.
;
if ( (max(ipar0) ne 1) or (min(ipar0) ne 1)) then begin
  print, 'Warning: Some particles not found at all or found more'
  print, 'than once in snapshot files.'
  print, 'Particle number---No. of occurences'
  for i=0,npar-1 do begin
    if (ipar0[i] ne 1) then begin
      print, i, ipar0[i]
    endif
  endfor
endif
;
;  Put in to object.
;
makeobject="object = CREATE_STRUCT(name=objectname,[" + $
    arraytostring(variables,QUOTE="'",/noleader) + ",'t']," + $
    arraytostring(variables,/noleader) + ',t' + ")"
if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /INFO
  undefine,object
endif

; If requested print a summary
;if keyword_set(STATS) or (not (keyword_set(NOSTATS) or keyword_set(QUIET))) then begin
;  pc_object_stats,object,dim=dim,QUIET=QUIET
if (not keyword_set(QUIET)) then print,' t = ', t
;endif


end
