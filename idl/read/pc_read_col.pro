FUNCTION get_col_lines, filename
  openr, unit, filename, /GET_LUN,/f77
  str = ''  
  count = 0ll  
  WHILE ~ EOF(unit) DO BEGIN  
    READU, unit, str  
    count = count + 1  
  ENDWHILE  
  FREE_LUN, unit  
  RETURN, count  
END  
;
PRO pc_read_col, filename=filename,t=t, coldat=coldat,compdat=compdat
;
;  bands=cs
;  velbands=cv
  default, filename, 'data/collisions.dat'
;
  pc_read_param, obj=par,/param2
  pc_read_pdim, obj=pd
;
  np=float(pd.npar)
  cv=par.colvel & cs=par.colspace & velmult=par.velmult & maxs=par.col_radius
  bands=cs
  velbands=cv
;
  print,'cv=',cv,' cs=',cs,' velmult=',velmult,' maxs=', maxs,' np=', np
  print,''
;
  lines=get_col_lines(filename)
  t=fltarr(lines)
  coldat=fltarr(lines,bands,velbands)
  compdat=coldat
;
  count=0
  tempcoldat=fltarr(bands,velbands)
  tempcompdat=tempcoldat
;
  openr, unit, filename,/GET_LUN,/f77
  WHILE (count lt lines) DO BEGIN
    READU, unit, tempt,tempcoldat,tempcompdat
    t(count)=tempt
    coldat(count,*,*) =tempcoldat
    compdat(count,*,*)=tempcompdat
    count=count+1
  ENDWHILE
  FREE_LUN, unit
END
