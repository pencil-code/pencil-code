;;;;;;;;;;;;;;;;;;;
;;;   sp.pro   ;;;
;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   13-Jan-2004
;;;
;;;  Description:
;;;   Plot spectrum

spec_uz = abs(fft(uu[l1:l2,3,(n1+n2)/2,2]))
spec_ux = abs(fft(uu[l1:l2,3,(3*n1+n2)/4,0]))

plot, spec_uz, /YLOG, PSYM=-1, $
    XRANGE=[0,10], YRANGE=[1.e-4,1]*max(spec_uz), YSTYLE=3
oplot, spec_ux*max(spec_uz)/max(spec_ux),col=120, PSYM=-1

print, FORMAT='(I3,2F14.8)', $
    [transpose(dindgen(8)), $
     transpose(spec_uz[0:7])/max(spec_uz), $
     transpose(spec_ux[0:7])/max(spec_ux)]

end
; End of file sp.pro
