pro pc_read_ode, file=file, object=obj

  common pc_precision, zero, one, precision, data_type, data_bytes, type_idl

  datadir = strtrim(pc_get_datadir(datadir),2)
  pc_set_precision, datadir=datadir, quiet=quiet
   
  nvar=0L & num=0L

  if not is_defined(file) then file='ode.dat'

  openr, 11, datadir+'/allprocs/'+strtrim(file), /f77
  readu, 11, nvar, num                 
  names=replicate(strjoin(replicate(' ',30)),num)
  lengs=lonarr(num)
  readu, 11, names
  readu, 11, lengs                     
  f_ode=make_array(nvar, type=type_idl)
  readu, 11, f_ode
  close, 11

  cmd='obj=create_struct('
  iend=lonarr(num+1) & iend[0]=-1
  for i=0,num-1 do begin
    iend(i+1)=iend(i)+lengs(i)
    cmd += '"'+strtrim(names(i))+'",f_ode(iend('+strtrim(string(i),2)+')+1:iend('+strtrim(string(i+1),2)+')),'
  endfor
  cmd += ')'
  cmd = cmd.replace(',)',')')
  
  res = execute(cmd)
end
