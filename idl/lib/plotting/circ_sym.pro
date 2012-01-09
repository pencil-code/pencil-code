;
;       CIRC_SYM
;
;     sets up circle psym symbols
;
pro circ_sym,siz,fill,thick=thick
;
if n_params(0) eq 0 then begin
  print,'circ_sym,siz,fill'
  print,'eg, circ_sym,1,1'
  return
endif
;
th=findgen(46)/45.*2*!pi
xx=cos(th)*siz
yy=sin(th)*siz
if fill eq 1 then usersym,xx,yy,/fill,thick=thick else usersym,xx,yy,thick=thick
return
end
;
;
;
