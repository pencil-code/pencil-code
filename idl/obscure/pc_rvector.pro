;  $Id$
;
;  Reads pre-selected vectors and plots in a 3-D box.
;  Data must be preprocessed with read_vectorfiles.x
;
;  18-aug-03/axel: coded
;
pro pc_rvector,nxyz=nxyz,png=png,cltbl=cltbl,backval=backval,dir=dir, $
  infile=infile,len=len,help=help
;
if keyword_set(help) then begin
  print,'pc_rvector,nxyz=nxyz,png=png,cltbl=cltbl,backval=backval,dir=dir,help=help'
  print,'dir is target directory; default is current working directory'
  print,'Typical calling sequence:'
  print,"pc_rvector,nxyz=512,len=15,dir='PNG/',/png"
  return
endif
;
default,dir,''
default,infile,'bvec.dat'
default,nxyz,64*2
default,cltbl,5
default,backval,1
default,len,5
print,'Assumed default box size is: nxyz=',nxyz
print,'(If not ok, then set nxyz!)'
;
;  time stamp label
;
siz=2
fo='(f6.1)'
loadct,cltbl
;
;  open MPEG file, if keyword is set
;
dev='x' ;(default)
if keyword_set(png) then begin
  set_plot, 'z'                   ; switch to Z buffer
  ;device, SET_RESOLUTION=[!d.x_size,!d.y_size] ; set window size
  device, SET_RESOLUTION=[1000,800] ; set window size
  print,'z-buffer; resolution=',!d.x_size,!d.y_size
  itpng=0 ;(image counter)
  dev='z'
endif
;
lun=41
nread=0
l=0L & m=0L & n=0L
file='data/'+infile
close,lun
openr,lun,file,/f77
while not eof(lun) do begin
  readu,lun,l,m,n,bx,by,bz
  ;print,'l,m,n,bx,by,bz=',l,m,n,bx,by,bz
  if l eq 0 then begin
    t=bx
    print,'t=',t
    print,'nread=',nread
    if nread gt 0 then begin
      pc_vectors_selected,ll-4,mm-4,nn-1,bbx,bby,bbz,$
        indgen(nxyz),indgen(nxyz),indgen(nxyz),$
        ax=30,az=30,len=len,back=backval,field=2
      xyouts,-12,-10,'!8t!6='+string(t,fo=fo),siz=siz,col=255-backval
      wait,.05
      ;
      if keyword_set(png) then begin
        istr2 = strtrim(string(itpng,'(I20.4)'),2) ;(only up to 9999 frames)
        image = tvrd()
        ;
        ;  write png file
        ;
        tvlct, red, green, blue, /GET
        imgname = dir+'img_'+istr2+'.png'
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
