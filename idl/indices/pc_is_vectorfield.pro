function pc_is_vectorfield,variable,dim=dim,subscripts=subscripts,NOVECTORS=NOVECTORS,_EXTRA=e
  varsize=size(variable)
COMPILE_OPT IDL2,HIDDEN

  if n_elements(dim) ne 1 then pc_read_dim,obj=dim,_EXTRA=e

  if varsize[0] lt 4 then return, 0

  result=0
  if (varsize[1] eq dim.nx) and (varsize[2] eq dim.ny) and (varsize[3] eq dim.nz) then begin
    result=1
    if arg_present(subscripts) then begin
      subscripts=make_array(varsize[0],/STRING,value='*')
    endif
  endif else if (varsize[1] eq dim.mx) and (varsize[2] eq dim.my) and (varsize[3] eq dim.mz) then begin 
    result=1
    if arg_present(subscripts) then begin
      subscripts=make_array(varsize[0],/STRING,value='*')
      subscripts[0]=str(dim.l1)+':'+str(dim.l2) 
      subscripts[1]=str(dim.m1)+':'+str(dim.m2) 
      subscripts[2]=str(dim.n1)+':'+str(dim.n2) 
    endif
  endif

  if result eq 1 then begin
    if (varsize[4] ne 3) or (keyword_set(NOVECTORS)) then result=0
  endif

 

  if result eq 0 then undefine,subscripts
  return, result
end
