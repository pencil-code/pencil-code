;
;  Summarize data
;
xyz = ['x', 'y', 'z']
fmt = '(A9,A,4G15.6)'
if (quiet le 2) then begin
  print, '  var             minval         maxval          mean           rms'
  ;
  ;
  for iv=1,totalvars do begin
    varname = varcontent[iv].idlvar
    if (varcontent[iv].skip eq 2) then begin
      for j=0,2 do begin
        cmd = "print, FORMAT=fmt,strmid('"+varname+"_'+xyz["+str(j) $
            + "]+'        ',0,8),'=', " $
            + "minmax("+varname+"(*,*,*,"+str(j)+")), " $
            + "mean("+varname+"(*,*,*,"+str(j)+"),/DOUBLE), " $
            + "rms("+varname+"(*,*,*,"+str(j)+"),/DOUBLE)"
        if (execute(cmd,1) ne 1) then $
            message, 'Error printing stats for ' + varcontent[iv].variable
      endfor
    endif else begin
      cmd = "print, FORMAT=fmt,strmid('"+varname+"        ',0,8),'=', " $
        + "minmax("+varname+"(*,*,*)), " $
        + "mean("+varname+"(*,*,*),/DOUBLE), " $
        + "rms("+varname+"(*,*,*),/DOUBLE)"
      if (execute(cmd,1) ne 1) then $
          message, 'Error printing stats for ' + varcontent[iv].variable
    endelse
    iv=iv+varcontent[iv].skip
  endfor
  ;
  print,'t = ',t
  ;
endif
