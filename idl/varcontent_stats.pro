xyz = ['x', 'y', 'z']
fmt = '(A9,A,4G15.6)'
print, '  var             minval         maxval          mean           rms'
;
;
for i=1,totalvars do begin
  if (varcontent[i].skip eq 2) then begin
      for j=0,2 do begin
          cmd = "print, FORMAT=fmt,strmid('"+varcontent[i].idlvar+"_'+xyz["+str(j)+"]+'        ',0,8),'=', " $
            + "minmax("+varcontent[i].idlvar+"(*,*,*,"+str(j)+")), " $
            + "mean("+varcontent[i].idlvar+"(*,*,*,"+str(j)+"),/DOUBLE), " $
            + "rms("+varcontent[i].idlvar+"(*,*,*,"+str(j)+"),/DOUBLE)"
          if (execute(cmd,1) ne 1) then $
                          message, 'Error printing stats for ' + varcontent[i].variable         
      end
  end else begin
      cmd = "print, FORMAT=fmt,strmid('"+varcontent[i].idlvar+"        ',0,8),'=', " $
        + "minmax("+varcontent[i].idlvar+"(*,*,*)), " $
        + "mean("+varcontent[i].idlvar+"(*,*,*),/DOUBLE), " $
        + "rms("+varcontent[i].idlvar+"(*,*,*),/DOUBLE)"
      if (execute(cmd,1) ne 1) then $
        message, 'Error printing stats for ' + varcontent[i].variable         
  end
  i=i+varcontent[i].skip
end
;
print,'t = ',t
