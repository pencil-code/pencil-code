;  $Id: pc_rvector.pro,v 1.7 2003-08-18 15:56:16 brandenb Exp $
;
;  Reads pre-selected vectors and plots in a 3-D box.
;  Data must be preprocessed with read_vectorfiles.x
;
;  18-aug-03/axel: coded
;
pro pc_rvector,nxyz=nxyz,png=png
;
default,nxyz,64*2
print,'Assumed default box size is: nxyz=',nxyz
print,'(If not ok, then set nxyz!)'
;
;  time stamp label
;
siz=2
fo='(f6.1)'
loadct,5
;
;  open MPEG file, if keyword is set
;
dev='x' ;(default)
if keyword_set(png) then begin
  set_plot, 'z'                   ; switch to Z buffer
  device, SET_RESOLUTION=[!d.x_size,!d.y_size] ; set window size
  itpng=0 ;(image counter)
  dev='z'
endif
;
lun=41
nread=0
l=0L & m=0L & n=0L
file='data/bvec.dat'
close,lun
openr,lun,file,/f77
while not eof(lun) do begin
  readu,lun,l,m,n,bx,by,bz
  if l eq 0 then begin
    t=bx
    print,'nread=',nread
    if nread gt 0 then begin
      vecgdv_good,ll-4,mm-4,nn-1,bbx,bby,bbz,$
        indgen(nxyz),indgen(nxyz),indgen(nxyz),$
        ax=30,az=30,len=1e6,back=255
      xyouts,-5,-5,'!8t!6='+string(t,fo=fo),siz=siz,col=1
      wait,.05
      ;
      if keyword_set(png) then begin
        istr2 = strtrim(string(itpng,'(I20.4)'),2) ;(only up to 9999 frames)
        image = tvrd()
        ;
        ;  write png file
        ;
        tvlct, red, green, blue, /GET
        imgname = 'img_'+istr2+'.png'
        write_png, imgname, image, red, green, blue
        itpng=itpng+1 ;(counter)
      endif
      ;
    endif
    nread=0
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
    nread=nread+1
  endelse
endwhile
close,lun
END
