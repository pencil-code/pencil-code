pro oplot_sn,sn=sn,_EXTRA=_EXTRA
  default, color, 128

  if n_elements(sn) ne 1 then pc_read_sn,obj=sn

  opvline,sn.t,_EXTRA=_EXTRA
end
