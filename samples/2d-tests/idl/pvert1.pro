;;;;;;;;;;;;;;;;;;;;;
;;;   pvert.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   11-Nov-2001
;;;  $Id: pvert1.pro,v 1.2 2005-06-26 17:34:29 eos_merger_tony Exp $
;;;
;;;  Description:
;;;    Plot vertical profiles of uz, lnrho and entropy.
;;;  PostScript output for the manual:
;;;    psa,file='pvert0.eps' & device,XSIZE=14,YSIZE=7
;;;    pvert_layout=[0,4,1]
;;;    .r thermo
;;;    .r pvert
;;;    pse

default, pvert_layout, [0,4,3]
default, nprofs, 50              ; set to N for only N profiles, 0 for all

sym = texsyms()

Temp = cs0^2/gamma1*exp(gamma*ss+gamma1*(lnrho-lnrho0))

nign = 3                        ; Number of close-to-bdry points to ignore

;
; construct vector of vertical pencils to plot
;
Nxmax = mx-2*nign > 1
Nymax = my-2*nign > 1
Nmax = Nxmax*Nymax
if ((nprofs le 0) or (nprofs gt Nmax)) then nprofs = Nmax
ixp=[[mx/2]] & iyp=[[my/2]]     ; case nprofs=1
if (nprofs gt 1) then begin
  Nxp = sqrt(nprofs+1e-5)*Nxmax/Nymax > 1
  Nyp = sqrt(nprofs+1e-5)*Nymax/Nxmax > 1
  Nxyp = Nxp*Nyp
  ;; Correct for case where Nxp or Nyp was <<1 and thus Nxp*Nyp>>nprofs
  if (Nxp eq 1) then begin
    Nyp = Nyp*nprofs/Nxyp > 1
  endif else if (Ny eq 1) then begin
    Nxp = Nxp*nprofs/Nxyp > 1
  endif
  Nxyp = Nxp*Nyp
  ixp = nign + spread( indgen(Nxp)*Nxmax/Nxp, 1, Nyp )
  iyp = nign + spread( indgen(Nyp)*Nymax/Nyp, 0, Nxp )
  ixp = floor(ixp)
  iyp = floor(iyp)
; floor() will have dropped a degenerate trailing dimension (any IDL
; function will, this is crazy..), so let us reform again to be sure:
  ixp = reform(ixp,Nxp,Nyp,/OVERWRITE)
  iyp = reform(iyp,Nxp,Nyp,/OVERWRITE)
endif

save_state
if (!p.charsize eq 0) then !p.charsize=1.3

!p.multi = pvert_layout

if (!d.name eq 'X') then begin
  red = 130/256.*!d.table_size   ; brick red for color table 5
endif else begin
  red = !p.color   ; black for PostScript  
endelse

if (!d.name eq 'X') then begin
  !p.charsize = 1. + (max(!p.multi)-1)*0.3
endif

!y.title = '!8z!X'

for ivar = 0,8 do begin
  xlog = 0                      ; default
  yr = minmax(z)
  case ivar of
    0: begin
      var = exp(lnrho)
      title = '!6Density ' + sym.varrho
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
      if (n_elements(lnrhoinit) gt 0) then xr = minmax([xr,lnrhoinit])
      xlog = 1
    end
    1: begin
      var = ss
      title = '!6Entropy !8s!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
      if (n_elements(ssinit) gt 0) then xr = minmax([xr,ssinit])
    end
    2: begin
      var = Temp
      title = '!6Temperature !8T!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
      if (n_elements(Tinit) gt 0) then xr = minmax([xr,Tinit])
    end
    3: begin
      var = sqrt(cs20*exp(gamma*ss+gamma1*(lnrho-lnrho0)))
      title = '!6Sound speed !8c!6!Ds!N!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
      if (n_elements(ssinit) gt 0) then xr = minmax([xr,ssinit])
    end
    4: begin
      var = cs20*rho0/gamma*exp(gamma*(ss+lnrho-lnrho0))
      title = '!6Pressure !8p!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
      if (n_elements(Tinit) gt 0) then xr = minmax([xr,Tinit])
      xlog = 1
    end
    5: begin
      var = zder(Temp)
      title = '!8dT!6/!8dz!6!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
      yr=[-1.15,-0.95]
    end
    6: begin
      var = uu[*,*,*,0]
      title = '!8u!Dx!N!6!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
    end
    7: begin
      var = zder(uu[*,*,*,0])
      title = '!8du!Dx!N!6/!8dz!6!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
      yr=[-1.15,-0.95]
    end
    8: begin
      var = uu[*,*,*,2]
      title = '!8u!Dz!N!6!X'
      if (ny eq 1) then xr = minmax(var[*,3,n1:n2]) else xr = minmax(var)
    end
    else: begin
      restore_state
      message, 'No such variable'
    end
  endcase

  ;; Draw axes
  plot, z, z, /NODATA, $
      XRANGE=xr, XSTYLE=3, XLOG=xlog, $
      YRANGE=yr, YSTYLE=3,  $
      TITLE=title
  ;; Plot curves along different lines
  for ix=0,(size(ixp))[1]-1 do begin
    for iy=0,(size(ixp))[2]-1 do begin
      oplot, var[ixp[ix,iy],iyp[ix,iy],*], z
    endfor
  endfor
  ;; Overplot averages
  var_xyavg = haver(var[l1:l2,m1:m2,*])
  oplot, var_xyavg, z, COLOR=!p.background, LINE=2, THICK=3
  oplot, var_xyavg, z, COLOR=120, LINE=2, THICK=2
  ;; Horizontal lines for orientation
  ophline, [z0,z1,z2,ztop]
  if (ivar eq 1) then opvline

;; overplot initial profiles
  if (n_elements(Tinit) le 0) then begin
    if (ivar eq 0) then $
        message, 'No Tinit -- you should run thermo.pro', /INFO
  endif else begin
    case ivar of
      0: oplot, lnrhoinit, z, LINE=2, COLOR=red, THICK=2
      1: ;nothing to overplot
      2: oplot, ssinit, z, LINE=2, COLOR=red, THICK=2
      3: oplot, Tinit, z, LINE=2, COLOR=red, THICK=2
      else: ; nothing to overplot
    endcase
  endelse

endfor

;; Add some parameters
if (!d.name eq 'PS') then begin
  !p.charsize = !p.charsize*0.6
endif else begin

endelse

plot, [0,1], [0,1], /NODATA, XSTYLE=4, YSTYLE=4 ; establish coordinate system
dx = 0.6
dy = 0.1

x00 = -0.1
y00 = 1 - dy/2

x0 = x00 & y0 = y00

n12 = (n1+n2)/2

T_bot = mean(Temp[l1:l2,m1:m2,n1])
T_top = mean(Temp[l1:l2,m1:m2,n2])
T_mid = 0.5*(T_bot+T_top)
ss_bot = mean(ss[l1:l2,m1:m2,n1])
ss_top = mean(ss[l1:l2,m1:m2,n2])
rho_bot = exp(mean(lnrho[l1:l2,m1:m2,n1]))
rho_mid = exp(mean(lnrho[l1:l2,m1:m2,n12]))
rho_top = exp(mean(lnrho[l1:l2,m1:m2,n2]))
beta = (T_top-T_bot)/Lz
beta_ad = par.gravz
chi_bot = par2.Kbot/rho_bot
chi_mid = par2.Kbot/rho_mid
chi_top = par2.Kbot/rho_top
delta_T = T_top - T_bot
delta_ss = ss_top-ss_bot
Ra1_mid = (beta-beta_ad)*par.gravz*Lz^4/(T_mid*nu*chi_mid)
Ra2_mid = delta_ss*gravz*Lz^3/(nu*chi_mid)
Ra1_top = (beta-beta_ad)*par.gravz*Lz^4/(T_top*nu*chi_top)
Ra2_top = delta_ss*gravz*Lz^3/(nu*chi_top)

y0=y0    & xyouts, x0+dx/2, y0, '!8t!6 = '     + strtrim(t          ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!8m!6 = '          + strtrim(par.mpoly  ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!8c!6!Ds0!N!6 = '  + strtrim(par.cs0    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!8g!Dz!N!6 = '     + strtrim(par.gravz  ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, sym.gamma + '!6 ='  + strtrim(par.gamma  ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, sym.nu + '!6 ='     + strtrim(par2.nu    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!8K!N=!6const =' + strtrim(par2.Kbot  ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!8F!D!6bot!N ='    + strtrim(par2.Fbot  ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!6d!8T!6/d!8z!6 =' + strtrim(beta       ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!6(d!8T!6/d!8z!6)!Dad!N =' $
                                               + strtrim(beta_ad    ,2) + '!X'

x0 = x0 + dx & y0 = y00

y0=y0-dy & xyouts, x0, y0, sym.delta_c + '!8T!6-' $
                           + sym.delta_c + '!8T!D!6ad!N!6 =' $
                                               + strtrim(delta_T-beta_ad*Lz $
                                                                    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, sym.delta_c + '!8s!6 =' $
                                               + strtrim(delta_ss   ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, sym.chi + '!6!Dbot!N!6 =' $
                                               + strtrim(chi_bot    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, sym.chi + '!6!Dmid!N!6 =' $
                                               + strtrim(chi_mid    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, sym.chi + '!6!Dtop!N!6 =' $
                                               + strtrim(chi_top    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!6Ra!D!8T!6,mid!N!6 = ' $
                                               + strtrim(Ra1_mid    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!6Ra!D!8s!6,mid!N!6 = ' $
                                               + strtrim(Ra2_mid    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!6Ra!D!8T!6,top!N!6 = ' $
                                               + strtrim(Ra1_top    ,2) + '!X'
y0=y0-dy & xyouts, x0, y0, '!6Ra!D!8s!6,top!N!6 = ' $
                                               + strtrim(Ra2_top    ,2) + '!X'


if (all(pvert_layout eq [0,2,2])) then begin
  ; don't know where to place otherwise
  xyouts, 0.45,0.5, '!8t!6=' + strtrim(t,2), /NORMAL
endif

restore_state

end
; End of file pvert.pro


