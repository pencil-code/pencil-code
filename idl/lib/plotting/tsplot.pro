;   
;  Just like plot, BUT go through the X values and find places
;  where the X is decreases.  At these points, break the graph
;  and over plot the next bit (up to the next jump!)
;
;  Useful for restarted jobs where there may be some repeatition
;  in the time series. (Shows up restartability nicely!!)
;
;  24-jul-05/tony: coded
;
pro tsplot,x,y, $
         change_symbol=change_symbol, $
         change_line=change_line, $
         change_color=change_color, $
         _extra=_extra

  nsections=0L

  if (n_elements(X) ne n_elements(Y)) then begin
    message,"X and Y arrays differ in size!"
  endif

  startpoint=[0]
  for i=1L,n_elements(X)-1 do begin
    if (X[i] lt X[i-1]) then startpoint=[startpoint,i]
  endfor

  if (n_elements(startpoint) eq 1) then begin
    plot,X,Y,_extra=_extra
    return
  endif

  plot,X[startpoint[0]:startpoint[1]-1],Y[startpoint[0]:startpoint[1]-1],_extra=_extra

  startpoint=[startpoint,n_elements(X)]
  for i=1L,n_elements(startpoint)-2 do begin
    if (keyword_set(change_color)) then begin
     col=128+byte(128*i/(n_elements(startpoint)-2.0))
     oplot,X[startpoint[i]:startpoint[i+1]-1],Y[startpoint[i]:startpoint[i+1]-1],col=col,_extra=_extra
    endif else if (keyword_set(change_line)) then begin
     oplot,X[startpoint[i]:startpoint[i+1]-1],Y[startpoint[i]:startpoint[i+1]-1],lines=i,_extra=_extra
    endif else if (keyword_set(change_symbol)) then begin
     oplot,X[startpoint[i]:startpoint[i+1]-1],Y[startpoint[i]:startpoint[i+1]-1],ps=i,_extra=_extra
    endif else begin
     oplot,X[startpoint[i]:startpoint[i+1]-1],Y[startpoint[i]:startpoint[i+1]-1],_extra=_extra
    endelse
  endfor

end
