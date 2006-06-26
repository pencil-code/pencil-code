$Id: pvid.pro,v 1.1 2006-06-26 21:38:51 brandenb Exp $
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
default,bbmax,.03
;
fo='(f4.2)'
rvid_box,'bx',fo=fo,min=-bbmax,max=bbmax,tmin=tmin,$
  xsize=xsize,ysize=ysize,/newwindow,ylabel=1.15,$
  size_label=size_label,thlabel=thlabel,/png
;
END
