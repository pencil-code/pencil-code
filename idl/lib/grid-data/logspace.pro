;;;;;;;;;;;;;;;;;;;;;;;;
;;;   logspace.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   21-Jun-2001
;;;  Version: 0.35 (CVS: $Revision: 1.2 $)
;;;  Description:
;;;     Return a real vector of length N, containing logarithmically
;;;   equidistant values between 10^x1 and 10^x2 inclusively.
;;;   If N is omitted, a value of 100 is adopted.
;;;   Mimics the octave function `logspace'.
;;;     Alllows for the first 2 arguments to be written as a vector,
;;;   so `logspace(minmax(v),10)' works. 
;;;  Keywords:
;;;   GHOST     -- set this to the number of ghost cells before x1 or
;;;                after x2; GHOST can be a 2-element array
;;;                [left,right] or a scalar (applied to both sides)
;;;   UNIQUE    -- flag for returning a list of unique elements.
;;;                This implies that you may get less than N elements
;;;                (in many cases just one).
;;;                Useful if you call
;;;                  contour,..,LEVELS=linspace(minmax(var),N)
;;;                where no repeated levels are allowed
;;;  Note:
;;;   This would best be implemented as a call to linspace:
;;;     return, 10^linspace(x1,x2,n)

function logspace, x1, x2, n, $
                   GHOST=ghost, $
                   UNIQUE=unique
  on_error,2

  if (n_elements(ghost) eq 0) then begin
    nghost = [0,0]
  endif else begin
    nghost = rebin([ghost], 2)  ; replicate if necessary
  endelse
  default, unique, 0

  if (n_elements(x1) ge 2) then begin ; Reshuffle arguments
    xx1 = x1[0]
    xx2 = x1[1]
    if (n_elements(x2) ne 0) then nn = x2
  endif else begin
    xx1 = x1
    xx2 = x2    
    if (n_elements(n) ne 0) then nn = n
  endelse
  default, nn, 100
  n_real = nn - nghost[0] - nghost[1]
  list = 10^(xx1 + (findgen(nn)-nghost[0])*(xx2-xx1)/(n_real-1))

  if (unique) then list = list(uniq(list))

  return, list
end
; End of file logspace.pro
