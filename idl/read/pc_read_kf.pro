;$Id$
;
;  Read forcing wavenumber kf from the k.dat file
;
pro pc_read_kf,kf,datadir=datadir
;
; Default data directory
;
  if (not keyword_set(datadir)) then datadir='.'
;
;  define
;
nk=0
kf=0.
;
;  read
;
lun=1
openr,lun,datadir+'/k.dat'
readf,lun,nk,kf
close,lun
;
end
