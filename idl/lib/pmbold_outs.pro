
  pro pmbold_outs, x, y, nstr1, bstr, nstr2, relshift=relshift, _extra=extra

  ; mimics boldface fonts: string bstr is output at point (x,y) of the data coordinate system 
  ;                        in boldface between non-boldface strings nstr1 and nstr2, the latter dropppable.
  ;                        relshift can be used to control the boldness: >1 bolder, <1 more slender.
  ;                        extra may contain valid xyouts keywords except /normal, /device and everything
  ;                        connected with 3D output.

    default, relshift, 1.
    if n_params() eq 4 then nstr2=''

    if nstr1 ne '' then xyouts, x, y, nstr1, width=w1, _extra=extra else w1=0.

    ncoors=convert_coord(x,y,/data,/to_normal)
    xyouts, ncoors[0], ncoors[1], nstr1+bstr, /norm, _extra=extra
    ncoors(0)+=w1+.0012*relshift
    xyouts, ncoors[0], ncoors[1], bstr+nstr2, /norm, _extra=extra

  end
