xyz = ['x', 'y', 'z']
fmt = '(A9,A,4G15.6)'
if (quiet le 2) then $
    print, '  var             minval         maxval          mean           rms'
;
;
if (quiet le 2) then begin
  for iv=1,totalvars do begin
    if (varcontent[iv].skip eq 2) then begin
      for j=0,2 do begin
        cmd = "print, FORMAT=fmt,strmid('"+varcontent[iv].idlvar+"_'+xyz["+str(j)+"]+'        ',0,8),'=', " $
            + "minmax("+varcontent[iv].idlvar+"(*,*,*,"+str(j)+")), " $
            + "mean("+varcontent[iv].idlvar+"(*,*,*,"+str(j)+"),/DOUBLE), " $
            + "rms("+varcontent[iv].idlvar+"(*,*,*,"+str(j)+"),/DOUBLE)"
        if (execute(cmd,1) ne 1) then $
            message, 'Error printing stats for ' + varcontent[iv].variable
      endfor
    endif else begin
      cmd = "print, FORMAT=fmt,strmid('"+varcontent[iv].idlvar+"        ',0,8),'=', " $
        + "minmax("+varcontent[iv].idlvar+"(*,*,*)), " $
        + "mean("+varcontent[iv].idlvar+"(*,*,*),/DOUBLE), " $
        + "rms("+varcontent[iv].idlvar+"(*,*,*),/DOUBLE)"
      if (execute(cmd,1) ne 1) then $
          message, 'Error printing stats for ' + varcontent[iv].variable         
    endelse
    iv=iv+varcontent[iv].skip
  endfor
  ;
  print,'t = ',t
endif
