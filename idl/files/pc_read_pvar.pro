; $Id: pc_read_pvar.pro,v 1.3 2005-01-12 14:11:30 ajohan Exp $
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
    {varcontent_all, $
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
line=fltarr(totalvars)*one
ipar0=intarr(npar)*1
ipar0_tot=intarr(npar)*1
t=zero
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
;  Read mask containing ones for particles that are present at the
;  current processor and zeros elsewhere.
;
  readu, file, ipar0
;
;  Go through all particles and read data for the ones that are present
;  at the current processor.
;
  for k=0,npar-1 do begin
    if (ipar0[k] eq 1) then begin
      ipar0_tot[k]=ipar0_tot[k]+1
;
;  Read line by line
;      
      readu, file, line
      for iv=1L,mpvar do begin
;
;  Put data from line into proper places of data array
;        
        res=varcontent[iv].idlvar+'[k,*]=line[iv-1:iv-1+varcontent[iv].skip]'
        if (execute(res,0) ne 1) then $
            message, 'Error putting data into '+varcontent[iv].idlvar+' array'
        iv=iv+varcontent[iv].skip
      endfor
    endif
  endfor
;
;  Read time.
;
  readu, file, t

  close, file
  free_lun, file

endfor
;
;  Check if all particles found exactly once.
;
if ( (max(ipar0_tot) ne 1) or (min(ipar0_tot) ne 1)) then begin
  print, 'Warning: Some particles not found at all or found more'
  print, 'than once in snapshot files.'
  print, 'Integrated mask=', ipar0_tot
endif
;
;  Put in to object.
;
makeobject="object = CREATE_STRUCT(name=objectname,[" + $
    arraytostring(variables,QUOTE="'",/noleader) + "]," + $
    arraytostring(variables,/noleader) + ")"
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
