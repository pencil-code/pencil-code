function pc_proc2box,proc,l,m,n,dim=dim,procdim=procdim,GHOST_INCLUDED=GHOST_INCLUDED

  if n_elements(dim) ne 1 then pc_read_dim,object=dim
  if n_elements(procdim) ne 1 then pc_read_dim,object=procdim,proc=proc

  if keyword_set(GHOST_INCLUDED) then begin
    l_loc=l-3 & m_loc=m-3 & n_loc=n-3
  endif else begin
    l_loc=l   & m_loc=m   & n_loc=n
  endelse
  nproc=dim.nprocx*dim.nprocy*dim.nprocz

  l_grid=procdim.nx*procdim.ipx+l_loc  
  m_grid=procdim.ny*procdim.ipy+m_loc  
  n_grid=procdim.nz*procdim.ipz+n_loc  

  if keyword_set(GHOST_INCLUDED) then begin
    l_grid=l_grid+3 & m_grid=m_grid+3 & n_grid=n_grid+3
  endif
 
  return,[l_grid,m_grid,n_grid]
end

