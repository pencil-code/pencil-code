;  $Id$
;
;  Summarize data
;  omit ghost zones in the analysis
;
xyz = ['x', 'y', 'z']
fmt = '(A10,A,4G15.6)'
if (quiet le 2) then begin
  print, '  var             minval         maxval          mean           rms'
  ;
  ;
  for iv=0L,totalvars-1L do begin
    if (varcontent[iv].variable eq 'UNKNOWN') then continue
    varname = varcontent[iv].idlvar
    if (varcontent[iv].skip eq 2) then begin
      for j=0,2 do begin
        cmd = "print, FORMAT=fmt,strmid('"+varname+"_'+xyz["+str(j) $
            + "]+'        ',0,10),'=', " $
            + "minmax("+varname+"[l1:l2,m1:m2,n1:n2,"+str(j)+"]), " $
            + "mean("+varname+"[l1:l2,m1:m2,n1:n2,"+str(j)+"],/DOUBLE), " $
            + "rms("+varname+"[l1:l2,m1:m2,n1:n2,"+str(j)+"],/DOUBLE)"
        if (execute(cmd,1) ne 1) then $
            message, 'Error printing stats for ' + varcontent[iv].variable
      endfor
    endif else begin
      cmd = "print, FORMAT=fmt,strmid('"+varname+"        ',0,10),'=', " $
        + "minmax("+varname+"[l1:l2,m1:m2,n1:n2]), " $
        + "mean("+varname+"[l1:l2,m1:m2,n1:n2],/DOUBLE), " $
        + "rms("+varname+"[l1:l2,m1:m2,n1:n2],/DOUBLE)"
      if (execute(cmd,1) ne 1) then $
          message, 'Error printing stats for ' + varcontent[iv].variable
    endelse
    iv=iv+varcontent[iv].skip
  endfor
endif
if (quiet le 3) then begin 
  ;
  print,'t = ',t
  ;
endif
