; $Id: pc_read_pvar.pro,v 1.2 2005-01-11 10:36:05 ajohan Exp $
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
;
;  Derived dimensions.
;
mpvar=pdim.mpvar
npar=pdim.npar
nproc=dim.nprocx*dim.nprocy*dim.nprocz
;
;  Define arrays for temporary storage of data.
;
array=fltarr(npar,mpvar)*1.e0
line=fltarr(mpvar)*1.e0
ipar0=intarr(npar)*1
ipar0_tot=intarr(npar)*1
t=0.e0
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
      readu, file, line
      array[k,*]=line
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
  print, 'than once in data file.'
  print, 'Integrated mask=', ipar0_tot
endif
;
;  Put in to object.
;
makeobject="object = CREATE_STRUCT(name=objectname, ['xx','vv'], " + $
    "array[*,0:2], array[*,3:5])"
if (execute(makeobject) ne 1) then begin
  print, 'ERROR Evaluating variables: ' + makeobject
  undefine, object
  stop
endif
;
;  Tidy up memory.
;
undefine, array

; If requested print a summary
;if keyword_set(STATS) or (not (keyword_set(NOSTATS) or keyword_set(QUIET))) then begin
;  pc_object_stats,object,dim=dim,QUIET=QUIET
if (not keyword_set(QUIET)) then print,' t = ', t
;endif


end
