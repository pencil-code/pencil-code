;$Id$
pro pc_check_math,location=location
;
;  date: 2005/12/14 11:56:20;  author: mee;  state: Exp;
;  Routine to allow regular checking of the floating point error status
;  and optionally give the location where the error was trapped.
;
;Create a string array of error names.  
ERRS = ['Integer divided by zero', $
        'Integer overflow', $
        'UNKNOWN ERROR', $
        'UNKNOWN ERROR', $
        'Divided by zero', $
        'Floating-point underflow', $
        'Floating-point overflow', $
        'Floating-point operand error']
  
;Get math error status.  
math_status = check_math()

;Check to see if an error occurred and print the corresponding  
;error message.  
FOR i = 4, 7 DO IF ISHFT(math_status, -i) AND 1 THEN begin
  if keyword_set(location) then begin
   print,ERRS[i], " occured in ",location,"."
  endif else begin
   print,ERRS[i], " occured."
  endelse
endif  

end
