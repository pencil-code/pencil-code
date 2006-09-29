if !d.name eq 'PS' then begin
  device,xsize=30,ysize=8,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  mv idl.ps ~/tex/mhd/alpha_pot/fig/p.ps
;
!x.title='x'
!y.title='z'
!p.charsize=1.6
contour,reform(aa(*,3,*,1)),nlev=20,x,z
oplot,x,x-x-.5,col=122
oplot,x,x-x+.5,col=122
oplot,+!pi+z-z,z,col=122
oplot,-!pi+z-z,z,col=122
END
