;$Id$
;
;  Read forcing wavenumber kf from the k.dat file
;
pro pc_read_kf,kf
;
;  define
;
nk=0
kf=0.
;
;  read
;
lun=1
openr,lun,'k.dat'
readf,lun,nk,kf
close,lun
;
end
