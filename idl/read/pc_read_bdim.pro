;;
;; $Id: pc_read_bdim.pro 9839 2008-09-05 07:24:02Z ajohan $
;;
;;  Read block domain decomposition dimension data.
;;
pro pc_read_bdim, object=object, datadir=datadir, quiet=quiet
compile_opt IDL2,HIDDEN
;
; Default data directory.
;
datadir = pc_get_datadir(datadir)
;
; Initialize all variables.
;
nbrickx=0L
nbricky=0L
nbrickz=0L
nblockmax=0L
mxb=0L
myb=0L
mzb=0L
nghostb=0L
nxb=0L
nyb=0L
nzb=0L
l1b=0L
l2b=0L
m1b=0L
m2b=0L
n1b=0L
n2b=0L
;
; Check for existence and read the data.
;
filename=datadir+'/bdim.dat'
if (file_test(filename)) then begin
  if is_valid(object,'BDIM',filename) then return
  if (not keyword_set(quiet)) then print, 'Reading ' + filename + '...'
  get_lun, file
  openr, file, filename
  readf, file, nbrickx, nbricky, nbrickz, nblockmax, $
       mxb, myb, mzb, nghostb, $
       nxb, nyb, nzb, $
       l1b, l2b, m1b, m2b, n1b, n2b
  free_lun,file
endif else $
  message, 'ERROR: cannot find file ' + filename
;
; Build structure of all the variables.
;
object = create_struct(name='PC_BDIM:'+strtrim(filename,2), $
    ['nbrickx','nbricky','nbrickz','nblockmax','mxb','myb','mzb','nghostb', $
    'nxb','nyb','nzb','l1b','l2b','m1b','m2b','n1b','n2b'], $
     nbrickx, nbricky, nbrickz, nblockmax, mxb, myb, mzb, nghostb, $
     nxb, nyb, nzb, l1b, l2b, m1b, m2b, n1b, n2b)
;
end
