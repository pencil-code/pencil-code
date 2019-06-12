
  pro pmbold_outs, x, y, nstr1, bstr, nstr2, relshift=relshift, _extra=extra

  ; mimics boldface fonts: string bstr is output at point (x,y) of the data coordinate system 
  ;                        in boldface between non-boldface strings nstr1 and nstr2, the latter dropppable.
  ;                        relshift can be used to control the boldness: >1 bolder, <1 more slender.
  ;                        extra may contain valid xyouts keywords except everything
  ;                        connected with 3D output.

    default, relshift, 1.
    if n_params() eq 4 then nstr2=''

    if nstr1 ne '' then xyouts, x, y, nstr1, width=w1, _extra=extra else w1=0.

    if n_elements(extra) ne 0 then begin
      if is_in(tag_names(extra),'DEVICE',/abbrev) ge 0 then $
        ncoors=convert_coord(x,y,/device,/to_normal) $
      else if is_in(tag_names(extra),'NORMAL',/abbrev) gt 0 then $
        ncoors=[x,y] $
      else $
        ncoors=convert_coord(x,y,/data,/to_normal)
    endif else $
        ncoors=convert_coord(x,y,/data,/to_normal)

    xyouts, ncoors[0], ncoors[1], nstr1+bstr, /norm, _extra=extra
    ncoors(0)+=w1+.0017*relshift
    xyouts, ncoors[0], ncoors[1], bstr+nstr2, /norm, _extra=extra

  end
