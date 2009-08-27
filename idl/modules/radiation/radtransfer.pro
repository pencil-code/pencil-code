;;;
;;; radtransfer.pro
;;;
;;; Author: Tobias Heinemann
;;; $Date: 2003-10-02 17:07:15 $
;;; $Revision: 1.3 $
;;;
;;; Description:
;;;   Re-implementation of some radiative-transfer routines in IDL, for
;;;   easier testing.
;;; Usage:
;;;     IDL> f=radtransfer('stratification.dat',/DOUBLE)
;;;   or
;;;     IDL> f=radtransfer('stratification.dat',/DOUBLE,/L2NDORDER)
;;;   or
;;;     IDL> f=radtransfer('stratification.dat',/DOUBLE,/L2NDORDER,/UPWARDS)
;;;   then
;;;     IDL> plot, f.z, f.qrad*f.kaprho
;;;   use
;;;     IDL> help, f, /STRUCT
;;;   for a full list of slots

pro radcalc
;
;  Calculate opacity and source function
;
  common variables,z,lnrho,ss,yH,TT,Srad,kaprho,Qrad0,emtau
  @data/pc_constants.pro

  ioncalc,lnrho,ss,yH,TT

  Srad=Srad0*(TT/TT_ion)^4
  kaprho=.25*exp(2.*lnrho-lnrho_e_)*(TT_ion_/TT)^1.5 $
            *exp(TT_ion_/TT)*yH*(1.-yH)*kappa0

end
;*******************************************************************
pro Qintr,Qrad,nrad,L2NDORDER=l2ndorder,UPWARDS=upwards
;
;  Integration radiation transfer equation along rays
;
  common variables,z,lnrho,ss,yH,TT,Srad,kaprho,Qrad0,emtau
  common startstop,n1,n2,nnstart,nnstop,nsign

  dlength=abs(z[1]-z[0])

  if (nrad gt 0) then begin
    nnstart=n1
    nnstop=n2+1
    nsign=1
  endif

  if (nrad lt 0) then begin
    nnstart=n2
    nnstop=n1-1
    nsign=-1
  endif

  emtau=0*lnrho+1
  Qrad=0*lnrho

  for n=nnstart,nnstop,nsign do begin
    if (keyword_set(l2ndorder)) then begin
      if (keyword_set(upwards)) then begin   ; use points n-1, n, n+1 to calculate point n
        dtau_m=(5*kaprho[n-nrad]+8*kaprho[n]-1*kaprho[n+nrad])*dlength/12
        dtau_p=(5*kaprho[n+nrad]+8*kaprho[n]-1*kaprho[n-nrad])*dlength/12
        dSdtau_m=(Srad[n]-Srad[n-nrad])/dtau_m
        dSdtau_p=(Srad[n+nrad]-Srad[n])/dtau_p
        Srad1st=(dSdtau_p*dtau_m+dSdtau_m*dtau_p)/(dtau_m+dtau_p)
        Srad2nd=2*(dSdtau_p-dSdtau_m)/(dtau_m+dtau_p)
        emdtau=exp(-dtau_m)
        emtau[n]=emtau[n-nrad]*emdtau
        if (dtau_m>1e-5) then begin
          emdtau1=1-emdtau
          emdtau2=emdtau*(1+dtau_m)-1
        endif else begin
          emdtau1=dtau_m-dtau_m^2/2+dtau_m^3/6
          emdtau2=-dtau_m^2/2+dtau_m^3/3-dtau_m^4/8
        endelse
          Qrad[n]=Qrad[n-nrad]*emdtau-Srad1st*emdtau1-Srad2nd*emdtau2 
      endif else begin          ; use points n-2, n-1, n to calculate point n
        dtau01=(5*kaprho[n     ]+8*kaprho[n-  nrad]-kaprho[n-2*nrad])*dlength/12
        ;; It seems to be slighlty more consistent to use points
        ;; n-2,n-1,n here, but I see no discernible difference:
;       dtau12=(5*kaprho[n-nrad]+8*kaprho[n-2*nrad]-kaprho[n-3*nrad])*dlength/12
        dtau12=(5*kaprho[n-2*nrad]+8*kaprho[n-  nrad]-kaprho[n       ])*dlength/12
        dSrad01=Srad[n     ]-Srad[n-  nrad]
        dSrad12=Srad[n-nrad]-Srad[n-2*nrad]
        Srad1st=(dSrad01*dtau12/dtau01+dSrad12*dtau01/dtau12)/(dtau01+dtau12)
        Srad2nd=2*(dSrad01/dtau01-dSrad12/dtau12)/(dtau01+dtau12)
        emdtau=exp(-dtau01)
        emtau[n]=emtau[n-nrad]*emdtau
        if (dtau01>1e-5) then emdtau1=1-emdtau $
                         else emdtau1=dtau01-dtau01^2/2+dtau01^3/6
        if (dtau01>1e-5) then emdtau2=emdtau*(1+dtau01)-1 $
                         else emdtau2=-dtau01^2/2+dtau01^3/3-dtau01^4/8
        Qrad[n]=Qrad[n-nrad]*emdtau $
                - Srad2nd*dtau01 $
                + (Srad2nd-Srad1st)*emdtau1
      endelse
    endif else begin
      dtau=.5*(kaprho[n-nrad]+kaprho[n])*dlength
      emdtau=exp(-dtau)
      emtau[n]=emtau[n-nrad]*emdtau
      if (dtau>1e-5) then tau_term=(1-emdtau)/dtau $
                     else tau_term=1-dtau/2+dtau^2/6
      Qrad[n]=Qrad[n-nrad]*emdtau+tau_term*(Srad[n-nrad]-Srad[n])
    endelse
  endfor
end
;*******************************************************************
pro radboundary,nrad

  common variables,z,lnrho,ss,yH,TT,Srad,kaprho,Qrad0,emtau
  common startstop,n1,n2,nnstart,nnstop

  Qrad0=0*lnrho

  if (nrad gt 0) then begin
    Qrad0[nnstart+nrad]=0
  endif
  if (nrad lt 0) then begin
    Qrad0[nnstart+nrad]=-Srad[nnstart+nrad]
  endif

end

pro Qrev,Qrad,nrad
  ;; Revise Q_rad, adding the (attenuated) boundary value

  common variables,z,lnrho,ss,yH,TT,Srad,kaprho,Qrad0,emtau
  common startstop,n1,n2,nnstart,nnstop,nsign

  for n=nnstart,nnstop,nsign do begin
    Qrad0[n]=Qrad0[n-nrad]
    Qrad[n]=Qrad[n]+Qrad0[n]*emtau[n]
  endfor

end
;*******************************************************************
function radtransfer,datafile,DOUBLE=double,_EXTRA=_extra
  ;; _extra will handle l2ndorder and upwards

  common variables,z,lnrho,ss,yH,TT,Srad,kaprho,Qrad0,emtau
  common startstop,n1,n2,nnstart,nnstop,nsign

  if (keyword_set(double)) then data=input_table(datafile,/double) $
                           else data=input_table(datafile)

  z=reform(data[0,*])
  lnrho=reform(data[1,*])
  ss=reform(data[2,*])

  n1=3 & n2=n_elements(lnrho)-4

  radcalc

  ; downward ray
  nrad=-1
  Qintr,Qrad_down,nrad, _EXTRA=_extra
  radboundary,nrad
  Qrev,Qrad_down,nrad

  ; upward ray
  nrad=+1
  Qintr,Qrad_up,nrad, _EXTRA=_extra
  radboundary,nrad
  Qrev,Qrad_up,nrad
  
  Qrad=.5*(Qrad_up+Qrad_down)

  f={z:         z[n1:n2], $
     lnrho:     lnrho[n1:n2], $
     ss:        ss[n1:n2], $
     yH:        yH[n1:n2], $
     Srad:      Srad[n1:n2], $
     kaprho:    kaprho[n1:n2], $
     emtau:     emtau[n1:n2], $
     Qrad_up:   Qrad_up[n1:n2], $
     Qrad_down: Qrad_down[n1:n2], $
     Qrad:      Qrad[n1:n2], $
     Qrad0:     Qrad0[n1:n2], $
     TT:        TT[n1:n2]}

  return,f

end
