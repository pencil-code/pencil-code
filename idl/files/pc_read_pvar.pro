; $Id: pc_read_pvar.pro,v 1.1 2005-01-09 16:58:18 ajohan Exp $
;
;   Read pvar.dat, or other PVAR file
;
pro pc_read_pvar, object=object, varfile=varfile, datadir=datadir
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one

; Default data directory
default, datadir, 'data'
default, varfile, 'pvar.dat'

; Get necessary dimensions
pc_read_pdim, obj=pdim, datadir=datadir

mpvar=pdim.mpvar
npar=pdim.npar
array=fltarr(npar,mpvar)*1.e0
line=fltarr(mpvar)*1.e0
ipar0=intarr(npar)*1
t=0.e0

; Get a unit number
GET_LUN, file

; Prepare for read
res=''
content=''

filename=datadir+'/proc0/'+varfile 

dummy=findfile(filename, COUNT=countfile)
if (not countfile gt 0) then message, 'ERROR: cannot find file '+ filename

close, file
openr, file, filename, /F77

readu, file, ipar0

for k=0,pdim.npar-1 do begin
  if (ipar0[k] eq k+1) then begin
    readu, file, line
  endif
  array[k,*]=line
endfor

readu, file, t

close, file

free_lun, file

makeobject="object = CREATE_STRUCT(name=objectname, ['xxp','vvp'], " + $
    "array[*,0:2], array[*,3:5])"
if (execute(makeobject) ne 1) then begin
  message, 'ERROR Evaluating variables: ' + makeobject, /INFO
  undefine, object
endif

; If requested print a summary
;if keyword_set(STATS) or (not (keyword_set(NOSTATS) or keyword_set(QUIET))) then begin
;  pc_object_stats,object,dim=dim,QUIET=QUIET
  print,' t = ', t
;endif


end
