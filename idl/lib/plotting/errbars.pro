pro errbars, xvals, yvals, errvals, del, _extra=extra
;
; Plots symmetric errorbars of width 2*errvals[i] about points (xvals[i], yvals[i]), all given in data coordinates, in "double-T" style.
; Width of the T: del in device units.
; Meaningful extra parameters: color, thickness, linestyle, /T3D, z-value.
; 
  for i=0,n_elements(xvals)-1 do begin
    absz=xvals[i] & ord=yvals[i]-errvals[i] & pnt=convert_coord(absz,ord,/data,/to_device) & xrng=pnt[0]+[-del,del] & yrng=[pnt[1],pnt[1]]
    plots, xrng, yrng, /device, _extra=extra 
    ord=yvals[i]+errvals[i] & pnt=convert_coord(absz,ord,/data,/to_device) & xrng=pnt[0]+[-del,del] & yrng=[pnt[1],pnt[1]]
    plots, xrng, yrng, /device, _extra=extra
    oplot, [xvals[i],xvals[i]],[yvals[i]-errvals[i],yvals[i]+errvals[i]],_extra=extra
  endfor

end
