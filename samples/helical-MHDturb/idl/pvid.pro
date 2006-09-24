$Id: pvid.pro,v 1.2 2006-09-24 17:09:58 brandenb Exp $
;
;  preprocess run using "src/read_videofiles.x" with "bb"
;
;@param
default,tmin,0.
default,size_label,3.0
default,thlabel,4.
default,xsize,500
default,ysize,420
;
default,bbmax,0.
;
fo='(f5.1)'
rvid_box,'bz',fo=fo,min=-bbmax,max=bbmax,tmin=tmin,$
  xsize=xsize,ysize=ysize,/newwindow,ylabel=1.15,$
  size_label=size_label,thlabel=thlabel,/mono,/png
;
END
