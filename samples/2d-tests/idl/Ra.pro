;;;;;;;;;;;;;;;;;;;
;;;   Ra.pro   ;;;
;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   13-Jan-2004
;;;
;;;  Description:
;;;   Calculate the Rayleigh number
;;;   Requires thermo.pro to have run before

n_mid = (n1+n2)/2               ; requires odd nz

Tinit = cs20/gamma1*exp(gamma*ssinit+gamma1*(lnrhoinit-lnrho0))

ss_bot   = ssinit[n1]
ss_mid   = ssinit[n_mid]
ss_top   = ssinit[n2]
T_bot   = Tinit[n1]
T_mid   = Tinit[n_mid]
T_top   = Tinit[n2]
rho_bot = exp(lnrhoinit[n1])
rho_mid = exp(lnrhoinit[n_mid])
rho_top = exp(lnrhoinit[n2])

dss   = ss_top-ss_bot
dTemp = T_top-T_bot
drho  = rho_top - rho_bot


ds_dz = (zder(spread(ssinit,[0,1],[1,1])))[0,0,n_mid]


grav = -1.D0
height = 1.D0
Kbot = par2.Kbot

dT_ad = grav

print, 'dss      = ', strtrim(dss,2)
print, 'dT,dT_ad = ', strtrim(dTemp,2), ',  ', strtrim(dT_ad,2)
print, 'drho     = ', strtrim(drho,2)
print, 'T_{bot,mid,top}   = ', $
    strtrim(T_bot,2),', ',strtrim(T_mid,2),', ',strtrim(T_top,2)
print, 'rho_{bot,mid,top} = ', $
     strtrim(rho_bot,2),', ',strtrim(rho_mid,2),', ',strtrim(rho_top,2)

rho_bmt = [rho_bot,rho_mid,rho_top]
T_bmt   = [T_bot  ,T_mid  ,T_top  ]

if (strtrim(par2.ivisc) eq 'rho_nu-const') then begin
  print, 'rho*nu=const'
  mu_bmt = par2.nu*rho0*[1,1,1]
  nu_bmt = mu_bmt/rho_bmt
endif else begin
  print, 'nu=const'
  nu_bmt = par2.nu*[1,1,1]
  mu_bmt = nu_bmt*rho_bmt
endelse
print, 'nu_mid  = ', strtrim(nu_bmt[1],2)
print, 'mu_mid  = ', strtrim(mu_bmt[1],2)
print, 'K_bot   = ', strtrim(Kbot,2)

print,'------------------------------'


Ra_s_bulk = grav*dss*height^3*rho_bmt/nu_bmt/Kbot
print, 'Ra_s,bulk  = ', strtrim(grav*dss*height^3*rho_bmt/nu_bmt/Kbot,2)
print, 'Ra_s,deriv = ', strtrim(grav*ds_dz*height^4*rho_bmt/nu_bmt/Kbot,2)
print, 'Ra_T       = ', strtrim(grav*(dTemp-dT_ad)/T_bmt*height^3*rho_bmt/nu_bmt/Kbot,2)

end
; End of file Ra.pro
