;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   shocktube.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   27-Jun-2002
;;;
;;;  Description:
;;;   Analytical solution of the shocktube problem
;;;
;;;  Initial state (pressure jump):
;;;
;;;    ------------------+------------------
;;;            pl        |        pr
;;;    ------------------+------------------
;;;                      0
;;;
;;;  Evolved state: membrane has snapped, four domains:
;;;
;;;    -----------------------------------------
;;;          1 (l)   |   2  |   3   |  4 | 5 (r)
;;;    --------------+------+-------+----+------
;;;                  x      x       x    x
;;;                   1      2       3    4
;;;
;;;  Velocity:
;;;                          *************
;;;                        *
;;;                      *                      
;;;                    *
;;;    -**************--------------------***---
;;;
;;;  Nomenclature:
;;;    pl, pr - pressure far left/right from the shock structure
;;;    p2     - pressure in the expansion wave
;;;    p3=p4  - no pressure jump across the contact discontinuity
;;;    p: pressure, T: temperature, rho: density,
;;;    u: flow velocity, cs: sound speed
;;;
;;;    gamma (adiabatic index) is assumed to be constant
;;;
;;;    The borders between the different domains are x1, x2, x3, x4
;;;    and move at velocities ul-csl, u4-c4,  u2, u3, u4, respectively.
;;;
;;;    Warning: This works so far only in the case ul=ur=0

pro shocktube, x, t, p, rho, u, parl, parr, gamma, DEBUG=debug
;
;  X          : coordinate vector (the initial discontinuity is always
;               at x=0, so you may want to call this routine in the form
;                 shocktube, x-x0, t, p, rho, parl, parr, gamma
;  T          : time after membrane snapped
;  P, RHO     : pressure and density at positions X
;  PARL, PARR : parameters [u,p,rho] to left and right of membrane
;  GAMMA      : adiabatic index
;
; Default parameters are for Sod's reference problem

  default, debug, 0


  if (n_elements(parl) le 0) then begin
    message, "No parameters specified, running Sod's reference problem", /INFO
    gamma = 1.4
    ul   = 0. &  ur  = 0.
    pl   = 3. &  pr  = 0.1
    rhol = 1. & rhor = 0.125
  endif else begin
    ul   = parl[0] & ur   = parr[0]
    pl   = parl[1] & pr   = parr[1]
    rhol = parl[2] & rhor = parr[2]
  endelse
;; Warn about imperfections:
  if ((ul ne 0) or (ur ne 0)) then $
      message, "Case ur=0 or ul=0 is not OK yet -- results won't make sense", $
      /INFO

  gamm1=gamma-1.

  csl=sqrt(gamma*pl/rhol)       ; left sound speed

;; declare fields:
  u   = x*0.
  p   = x*0.
  rho = x*0.

;
;  iteratively find p3/pl
;
  p3=pl*(pr/pl)^0.2             ; initial guess
  for i=1,20 do begin
    u3 = csl*2/gamm1*(1-(p3/pl)^(gamm1/2/gamma))
    p3 = pr + (u3-ur)*sqrt(rhor/2*((gamma+1)*p3+gamm1*pr))

    if (debug) then print, p3/pl, u3
  endfor

  rho3 = rhol*(p3/pl)^(1./gamma)
  cs3 = sqrt(gamma*p3/rho3)

  p4 = p3
  u4 = u3
  us = ur + (pr-p4)/(ur-u4)/rhor ; velocity of shock front
  rho4 = -(pr-p4)/(ur-u4)/(u4-us)
  cs4 = sqrt(gamma*p4/rho4)

;; positions of separating faces
  x1 = (ul-csl)*t
  x2 = (u3-cs3)*t
  x3 = u4*t
  x4 = us*t

;; calculate profiles
  left  = where(x le x1)
  reg2  = where((x gt x1) and (x lt x2)) ; expansion region
  reg3  = where((x gt x2) and (x lt x3))
  reg4  = where((x gt x3) and (x lt x4))
  right = where(x gt x4)

  if (left[0] ge 0) then begin
    u[left] = ul
    p[left] = pl
    rho[left] = rhol
  endif

  if (reg2[0] ge 0) then begin  ; expansion region
    u[reg2] = 2/(gamma+1)*(csl+x[reg2]/t+gamm1/2*ur)
    p[reg2] = pl*(1-gamm1/2*u[reg2]/csl)^(2*gamma/gamm1)
    rho[reg2] = rhol*(1-gamm1/2*u[reg2]/csl)^(2/gamm1)
  endif

  if (reg3[0] ge 0) then begin
    u[reg3] = u3
    p[reg3] = p3
    rho[reg3] = rho3
  endif

  if (reg4[0] ge 0) then begin  ; expansion region
    u[reg4] = u4
    p[reg4] = p4
    rho[reg4] = rho4
  endif

  if (right[0] ge 0) then begin
    u[right] = ur
    p[right] = pr
    rho[right] = rhor
  endif

  if (0) then begin
    !p.multi=[0,2,2]
    plot, x, p, /YLOG, YSTYLE=3, YTITLE='!8p!X'
    plot, x, u, YSTYLE=3, YTITLE='!8u!X', ps=-1
    plot, x, rho, YSTYLE=3, /YLOG, YTITLE='!8r!X'
  endif

  if (debug) then begin
    print, 'u3=u4 = ', u3, u4
    print, 'p3=p4 = ', p3, p4
    print, 'rho4 = ', rho4
    print, 'rho3 = ', rho3
    print
    print, 'V1 = ', ul-csl
    print, 'V2 = ', u4-cs3
    print, 'V3 = ', u4
    print, 'V4 = ', us
  endif

end
; End of file shocktube.pro
