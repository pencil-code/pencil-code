nxyz=32
lun=41
l=0L & m=0L & n=0L
file='data/proc0/bvec.dat'
close,lun
openr,lun,file,/f77
while not eof(lun) do begin
  readu,lun,l,m,n,bx,by,bz
  ;print,l,m,n
  if l eq 0 then begin
    vecgdv_good,ll-4,mm-4,nn-1,bbx,bby,bbz,indgen(nxyz),indgen(nxyz),indgen(nxyz),len=30
    t=bx
    print,t,n_elements(ll)
    readnew=1
     ;stop
  endif else begin
    if(readnew eq 1) then begin
      ll=l & mm=m & nn=n
      bbx=bx & bby=by & bbz=bz
      readnew=0
    endif else begin
      ll=[ll,l] & mm=[mm,m] & nn=[nn,n]
      bbx=[bbx,bx] & bby=[bby,by] & bbz=[bbz,bz]
    endelse
  endelse
endwhile
close,lun
END
