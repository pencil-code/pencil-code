pro oplot_sn,sn=sn,tscale=tscale,_EXTRA=_EXTRA
  default, color, 128
  default, tscale, 1.

  if n_elements(sn) ne 1 then pc_read_sn,obj=sn

  opvline,sn.t*tscale,_EXTRA=_EXTRA
end
