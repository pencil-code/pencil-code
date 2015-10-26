;
;       CIRC_SYM
;
;     sets up circle psym symbols
;     10/03/26: MR added handing over of all keyword parameters to usersym,
;               if parameter thick missing or 0 value is set from !p.thick
;
pro circ_sym,siz,fill,thick=thick,_extra=extra
;
if (n_params() eq 0) then begin
  print,'circ_sym,siz,fill'
  print,'eg, circ_sym,1,1'
  return
endif
;
th=findgen(46)/45.*2*!pi
xx=cos(th)*siz
yy=sin(th)*siz

if (n_params() eq 2) then $
  if (fill ne 0) then $
    if (exists(extra)) then $
      extra=create_struct(extra, 'fill', 1) $
    else $
      extra=create_struct('fill', 1)

if (keyword_set(extra)) then $
  if (has_tag(extra, 'thick')) then $
    extra=create_struct('thick', !p.thick, extra) $
  else begin
    if (extra.thick eq 0) then $
      extra.thick = !p.thick
  endelse $
else $
  extra=create_struct('thick', !p.thick)

if (fill eq 1) then usersym,xx,yy,/fill,thick=thick else usersym,xx,yy,thick=thick
return
end
;
;
;
