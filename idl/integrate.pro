;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   integrate.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   10-Apr-2002
;;;  $Id$
;;;
;;;  Description:
;;;   Calculat a few integral quantities and provide dVol for further
;;;   integrations.
;;;  To do:
;;;   Also construct dSurf for surface integrals.

;; Construct dVol

dVol = make_array(nx,ny,nz,VALUE=dx*dy*dz)

if (nghostx gt 0) then begin
  dVol[0:nghostx-1    ,*,*] = 0.
  dVol[nx-nghostx:nx-1,*,*] = 0.
  dVol[nghostx     ,*,*] = 0.5*dVol[nghostx     ,*,*]
  dVol[nx-nghostx-1,*,*] = 0.5*dVol[nx-nghostx-1,*,*]
endif

if (nghosty gt 0) then begin
  dVol[*,0:nghosty-1    ,*] = 0.
  dVol[*,ny-nghosty:ny-1,*] = 0.
  dVol[*,nghosty     ,*] = 0.5*dVol[*,nghosty     ,*]
  dVol[*,ny-nghosty-1,*] = 0.5*dVol[*,ny-nghosty-1,*]
endif

if (nghostz gt 0) then begin
  dVol[*,*,0:nghostz-1] = 0.
  dVol[*,*,nz-nghostz:nz-1] = 0.
  dVol[*,*,nghostz     ] = 0.5*dVol[*,*,nghostz     ]
  dVol[*,*,nz-nghostz-1] = 0.5*dVol[*,*,nz-nghostz-1]
endif

;; Construct dSurf

dSurfx=(dSurfy=(dSurfz=make_array(nx,ny,nz,VALUE=0.)))

dSurfx[nghostx     ,nghosty:ny-nghosty-1,nghostz:nz-nghostz-1] = -dy*dz
dSurfx[nx-nghostx-1,nghosty:ny-nghosty-1,nghostz:nz-nghostz-1] =  dy*dz
dSurfx[*,nghosty     ,*        ] = 0.5*dSurfx[*,nghosty     ,*        ]
dSurfx[*,ny-nghosty-1,*        ] = 0.5*dSurfx[*,ny-nghosty-1,*        ]
dSurfx[*,*        ,nghostz     ] = 0.5*dSurfx[*,*        ,nghostz     ]
dSurfx[*,*        ,nz-nghostz-1] = 0.5*dSurfx[*,*        ,nz-nghostz-1]

dSurfy[nghostx:nx-nghostx-1,nghosty     ,nghostz:nz-nghostz-1] = -dx*dz
dSurfy[nghostx:nx-nghostx-1,ny-nghosty-1,nghostz:nz-nghostz-1] =  dx*dz
dSurfy[nghostx     ,*,*        ] = 0.5*dSurfy[nghostx     ,*,*        ]
dSurfy[nx-nghostx-1,*,*        ] = 0.5*dSurfy[nx-nghostx-1,*,*        ]
dSurfy[*        ,*,nghostz     ] = 0.5*dSurfy[*        ,*,nghostz     ]
dSurfy[*        ,*,nz-nghostz-1] = 0.5*dSurfy[*        ,*,nz-nghostz-1]

dSurfz[nghostx:nx-nghostx-1,nghosty:ny-nghosty-1,nghostz     ] = -dx*dy
dSurfz[nghostx:nx-nghostx-1,nghosty:ny-nghosty-1,nz-nghostz-1] =  dx*dy
dSurfz[nghostx     ,*        ,*] = 0.5*dSurfz[nghostx     ,*        ,*]
dSurfz[nx-nghostx-1,*        ,*] = 0.5*dSurfz[nx-nghostx-1,*        ,*]
dSurfz[*        ,nghosty     ,*] = 0.5*dSurfz[*        ,nghosty     ,*]
dSurfz[*        ,ny-nghosty-1,*] = 0.5*dSurfz[*        ,ny-nghosty-1,*]

;; Evaluate some integral quantities

Mtot = total(exp(lnrho)*dVol)
Momtot = fltarr(3)
for i=0,2 do Momtot[i] = total(exp(lnrho)*dVol*uu[*,*,*,i])
Ekin = 0.5*total(exp(lnrho)*dot2(uu)*dVol)
Etherm = 1./gamma_m1*total(exp(gamma*(ss+lnrho)))
;
if (n_elements(pot) ne 0) then begin
  Epot = total(exp(lnrho)*pot*dVol)
  Mheavy = total(xder(pot)*dSurfx+yder(pot)*dSurfy+zder(pot)*dSurfz) / (4*!pi)
endif else Epot=(Mheavy=!VALUES.F_NAN)

print, 'Total mass, heavy mass:     ', Mtot, Mheavy
print, 'Momentum: ', Momtot
print, 'E_kin, Etherm:    ', Ekin, Etherm
print, 'Epot :   ', Epot

end
; End of file integrate.pro
