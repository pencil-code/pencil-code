nxyz=64
lun=41
nread=0
l=0L & m=0L & n=0L
;file='data/proc3/bvec.dat'
file='data/bvec.dat'
close,lun
openr,lun,file,/f77
while not eof(lun) do begin
  readu,lun,l,m,n,bx,by,bz
  if l eq 0 then begin
    t=bx
    print,'nread=',nread
    if nread gt 0 then begin
      vecgdv_good,ll-4,mm-4,nn-1,bbx,bby,bbz,indgen(nxyz),indgen(nxyz),indgen(nxyz),ax=30,az=30,len=1e6
      print,t,nread,n_elements(ll)
    endif
    nread=nread+1
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
