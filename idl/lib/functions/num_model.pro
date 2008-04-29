;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   num_model.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   11-Jun-2004
;;;
;;;  Description:
;;;   Like F90's epsilon(), huge() and tiny() functions
;;;   Only created one function to keep name space clean
;;;
;;;  Usage:
;;;   epsi = num_model(1.,/epsilon)
;;;
;;;  Keywords:
;;;   EPSILON -- return the smallest psitive number that makes a
;;;              difference when added to 1.
;;;    HUGE   -- reuturn the largest number
;;;    TINY   -- reuturn tha smallest positive number

function num_model, x, $
                    EPSILON=epsilon, HUGE=huge, TINY=tiny, $
                    HELP=help


if (keyword_set(help)) then extract_help, 'num_model'

default, epsilon, 0
default, huge   , 0
default, tiny   , 0

s = size(x)
type = s[s[0]+1]

case type of

  2: begin                      ; int
    epsi = 1
    hug  = 32767
    tiny = 1
  end

  3: begin                      ; long int
    epsi = 1L
    hug  = 2147483647
    tiny = 1L
  end

  4: begin                      ; real
    epsi = 1.1920929e-07
    hug  = 3.4028235e+38
    tiny = 1.1754944e-38
  end

  5: begin                      ; double precision
    epsi = 2.220446049250313d-016
    hug  = 1.7976931348623158d+308
    tiny = 2.225073858507201d-308
  end

  else: begin
    message, "Don't know how to treat type " + strtrim(type,2), /INFO
    return, !values.f_nan
  endelse

endcase

if (epsilon) then begin
  return, epsi
endif else if (huge) then begin
  return, hug
endif else if (tiny) then begin
  return, tiny
endif else begin
  message, "Need to specify one of the /EPSILON, /HUGE, /TINY switches"
endelse

print, 'Here'

end
; End of file epsilon.pro
