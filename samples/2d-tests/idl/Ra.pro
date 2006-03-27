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
;;;
;;;  Modified 24-3-2006, for use with cp /= 0 runs, ngrs & neamvonk

n_mid = (n1+n2)/2               ; requires odd nz

cp = par2.cp
Tinit = cs20/gamma1/cp*exp(gamma*ssinit/cp + gamma1*(lnrhoinit-lnrho0))  

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
dTemp = T_top-T_bot       ;; nb: this is for the *initial* T profile...
drho  = rho_top - rho_bot

ds_dz = (zder(spread(ssinit,[0,1],[1,1])))[0,0,n_mid]

grav = par2.gravz
height = 1.D0
Kbot = par2.Kbot

R_star=(gamma1/gamma)*cp
beta0=grav/( R_star*(mpoly+1) )
dT_stat = beta0*height   ;; i.e. dT_hydrostatic=dTemp = beta0*height
beta1=grav/cp
dT_ad = beta1*height      ;; i.e. dT_adiabatic = beta1*height, from GMSW
dT_sup=dT_stat-dT_ad      ;; i.e. dT_superadiabatic = beta*height, 
                          ;;       where beta=beta0-beta1

print, 'dss      = ', strtrim(dss,2)
print, 'dT,dT_ad = ', strtrim(dTemp,2), ',  ', strtrim(dT_ad,2)
print, 'drho     = ', strtrim(drho,2)
print, 'T_{bot,mid,top}   = ', $
    strtrim(T_bot,2),', ',strtrim(T_mid,2),', ',strtrim(T_top,2)
print, 'rho_{bot,mid,top} = ', $
     strtrim(rho_bot,2),', ',strtrim(rho_mid,2),', ',strtrim(rho_top,2)

rho_bmt = [rho_bot,rho_mid,rho_top]
T_bmt   = [T_bot  ,T_mid  ,T_top  ]

;; if (strtrim(par2.ivisc) eq 'rho_nu-const') then begin
    print, 'rho*nu=const'
    mu_bmt = par2.nu*rho0*[1,1,1]
    nu_bmt = mu_bmt/rho_bmt
;; endif else begin
;;   print, 'nu=const'
;;   nu_bmt = par2.nu*[1,1,1]
;;   mu_bmt = nu_bmt*rho_bmt
;; endelse
chi_bmt=Kbot/(cp*rho_bmt)
print, 'nu_mid  = ', strtrim(nu_bmt[1],2)
print, 'mu_mid  = ', strtrim(mu_bmt[1],2)
print, 'K_bot   = ', strtrim(Kbot,2)
print, 'chi_bmt   = ', strtrim(chi_bmt,2)

print,'------------------------------'

Ra_s_bulk = ( grav*(dss/cp)*height^3 ) / ( nu_bmt*chi_bmt )
Ra_s_deriv = ( grav*ds_dz/cp*height^4 ) / ( nu_bmt*chi_bmt )
Ra_T = ( grav*dT_sup/T_bmt*height^3 ) / ( nu_bmt*chi_bmt )

print, 'Ra_s,bulk  = ', strtrim(Ra_s_bulk,2)
print, 'Ra_s,deriv = ', strtrim(Ra_s_deriv,2)
print, 'Ra_T       = ', strtrim(Ra_T,2)

print,'------------------------------'

Pr = nu_bmt/chi_bmt
print, 'Pr          = ', strtrim(Pr,2)

print,'------------------------------'

end
; End of file Ra.pro
