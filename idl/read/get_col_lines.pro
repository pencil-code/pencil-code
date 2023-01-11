FUNCTION get_col_lines, filename, formatted=form, comchar=comchar

  default, form, 0

  check = form and is_defined(comchar)
  str = check ? replicate(' ',strlen(comchar)) : ''

  count = 0ll
  openr, unit, filename, /GET_LUN,f77=~form
  WHILE ~ EOF(unit) DO BEGIN
  
    if form then $
      READF, unit, str $
    else $ 
      READU, unit, str

    if check then begin
      if stregex(str,'^ *'+strtrim(comchar,2)) < 0 then count = count + 1
    endif else $
      count = count + 1  
  
  ENDWHILE
  FREE_LUN, unit  
  RETURN, count  
END  

