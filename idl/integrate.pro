;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   integrate.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   10-Apr-2002
;;;
;;;  Description:
;;;   Calculat a few integral quantities and provide dVol for further
;;;   integrations.
;;;  To do:
;;;   Also construct dSurf for surface integrals.

;; Construct dVol

dVol = make_array(nx,ny,nz,VALUE=dx*dy*dz)

if (nghx gt 0) then begin
  dVol[0:nghx-1    ,*,*] = 0.
  dVol[nx-nghx:nx-1,*,*] = 0.
  dVol[nghx     ,*,*] = 0.5*dVol[nghx     ,*,*]
  dVol[nx-nghx-1,*,*] = 0.5*dVol[nx-nghx-1,*,*]
endif

if (nghy gt 0) then begin
  dVol[*,0:nghy-1    ,*] = 0.
  dVol[*,ny-nghy:ny-1,*] = 0.
  dVol[*,nghy     ,*] = 0.5*dVol[*,nghy     ,*]
  dVol[*,ny-nghy-1,*] = 0.5*dVol[*,ny-nghy-1,*]
endif

if (nghz gt 0) then begin
  dVol[*,*,0:nghz-1] = 0.
  dVol[*,*,nz-nghz:nz-1] = 0.
  dVol[*,*,nghz     ] = 0.5*dVol[*,*,nghz     ]
  dVol[*,*,nz-nghz-1] = 0.5*dVol[*,*,nz-nghz-1]
endif

;; Construct dSurf

dSurfx=(dSurfy=(dSurfz=make_array(nx,ny,nz,VALUE=0.)))

dSurfx[nghx     ,nghy:ny-nghy-1,nghz:nz-nghz-1] = -dy*dz
dSurfx[nx-nghx-1,nghy:ny-nghy-1,nghz:nz-nghz-1] =  dy*dz
dSurfx[*,nghy     ,*        ] = 0.5*dSurfx[*,nghy     ,*        ]
dSurfx[*,ny-nghy-1,*        ] = 0.5*dSurfx[*,ny-nghy-1,*        ]
dSurfx[*,*        ,nghz     ] = 0.5*dSurfx[*,*        ,nghz     ]
dSurfx[*,*        ,nz-nghz-1] = 0.5*dSurfx[*,*        ,nz-nghz-1]

dSurfy[nghx:nx-nghx-1,nghy     ,nghz:nz-nghz-1] = -dx*dz
dSurfy[nghx:nx-nghx-1,ny-nghy-1,nghz:nz-nghz-1] =  dx*dz
dSurfy[nghx     ,*,*        ] = 0.5*dSurfy[nghx     ,*,*        ]
dSurfy[nx-nghx-1,*,*        ] = 0.5*dSurfy[nx-nghx-1,*,*        ]
dSurfy[*        ,*,nghz     ] = 0.5*dSurfy[*        ,*,nghz     ]
dSurfy[*        ,*,nz-nghz-1] = 0.5*dSurfy[*        ,*,nz-nghz-1]

dSurfz[nghx:nx-nghx-1,nghy:ny-nghy-1,nghz     ] = -dx*dy
dSurfz[nghx:nx-nghx-1,nghy:ny-nghy-1,nz-nghz-1] =  dx*dy
dSurfz[nghx     ,*        ,*] = 0.5*dSurfz[nghx     ,*        ,*]
dSurfz[nx-nghx-1,*        ,*] = 0.5*dSurfz[nx-nghx-1,*        ,*]
dSurfz[*        ,nghy     ,*] = 0.5*dSurfz[*        ,nghy     ,*]
dSurfz[*        ,ny-nghy-1,*] = 0.5*dSurfz[*        ,ny-nghy-1,*]

;; Evaluate some integral quantities

Mtot = total(exp(lam)*dVol)
Momtot = fltarr(3)
for i=0,2 do Momtot[i] = total(exp(lam)*dVol*uu[*,*,*,i])
Ekin = 0.5*total(exp(lam)*dot2(uu)*dVol)
Etherm = 1./gamma1*total(exp(gamma*(ent+lam)))
;
if (n_elements(pot) ne 0) then begin
  Epot = total(exp(lam)*pot*dVol)
  Mheavy = total(xder(pot)*dSurfx+yder(pot)*dSurfy+zder(pot)*dSurfz) / (4*!pi)
endif else Epot=(Mheavy=!VALUES.F_NAN)

print, 'Total mass, heavy mass:     ', Mtot, Mheavy
print, 'Momentum: ', Momtot
print, 'E_kin, Etherm:    ', Ekin, Etherm
print, 'Epot :   ', Epot

end
; End of file integrate.pro
