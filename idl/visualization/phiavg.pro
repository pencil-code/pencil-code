;;;;;;;;;;;;;;;;;;;;;;
;;;   phiavg.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   11-Aug-2003
;;;
;;;  Description:
;;;    Read azimuthal averages and plot
;;;    vname must be set to one of the fields listed in phiaver.in

default, avgdir, datatopdir+'/averages'
default, phiavgfile, 'PHIAVG1'
default, vname, 'b2mphi'

az = pc_read_phiavg(avgdir+'/'+phiavgfile)

pos = aspect_pos((max(az.z)-min(az.z))/(max(az.rcyl)-min(az.rcyl)), $
             MARGIN=0.1)

if (any(tag_names(az) eq strupcase(vname))) then begin
  ;; Extract slot by string vname given above. This is awkward, but
  ;; that's the price for having the values in a structure
  if (execute('var = az.'+vname)) then begin
    contourfill, var, az.rcyl, az.z, $
        POS=pos, XSTYLE=1, YSTYLE=1, /GRID, $
        XTITLE='!8r!X', YTITLE='!8z!X', TITLE=vname
  endif
endif else begin
  print, 'No such variable in phi-averages: vname=<', vname,'>'
endelse

end
; End of file phiavg.pro
