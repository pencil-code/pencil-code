;
; Visualize the temperature
;
!x.style=1
!x.title='x [cm]'
!y.title='T [K]'
!p.charsize=2
rvid_line,'lnTT',/exp,min=500,max=2200,/nocontour,/xgrid
END
