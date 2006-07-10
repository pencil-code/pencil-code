;; $Id: pc_read_zaver.pro,v 1.2 2006-07-10 11:24:02 ajohan Exp $
;;
;;   Read z-averages from file.
;;   Default is to only plot the data (with tvscl), not to save it in memory.
;;   The user can get the data returned in an object by specifying nit, the
;;   number of snapshots to save.
;;
pro pc_read_zaver, object=object, varfile=varfile, datadir=datadir, $
    nit=nit, lplot=lplot, iplot=iplot, min=min, max=max, zoom=zoom, $
    quiet=quiet
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;;
;;  Default data directory.
;;
default, datadir, './data'
default, varfile, 'zaverages.dat'
default, nit, 0
default, lplot, 1
default, iplot, 0
default, zoom, 1
default, min, 0.0
default, max, 0.0
default, quiet, 0
;;
;;  Get necessary dimensions.
;;
pc_read_dim, obj=dim, datadir=datadir, quiet=quiet
pc_set_precision, dim=dim, quiet=quiet
;;
;;  Derived dimensions.
;;
nx=dim.nx
ny=dim.ny
;;
;;  Read variables from zaver.in
;;
spawn, "echo "+datadir+" | sed -e 's/data$//g'", datatopdir
spawn, 'cat '+datatopdir+'/zaver.in', varnames
if (not quiet) then print, 'Preparing to read z-averages ', $
    arraytostring(varnames,quote="'",/noleader)
nvar=n_elements(varnames)
;;  Die if attempt to plot variable that does not exist.
if (iplot gt nvar-1) then message, 'iplot must not be greater than nvar-1!'
;;
;;  Define arrays to put data in.
;;  The user must supply the length of the time dimension (nit).
;;
if (nit gt 0) then begin

  if (not quiet) then print, 'Returning averages at ', strtrim(nit,2), ' times'

  tt=fltarr(nit)*one
  for i=0,nvar-1 do begin
    cmd=varnames[i]+'=fltarr(nx,ny,nit)*one'
    if (execute(cmd,0) ne 1) then message, 'Error defining data arrays'
  endfor

endif
;;
;;  Variables to put single time snapshot in.
;;
array=fltarr(nx,ny,nvar)*one
t  =0.0*one
;;
;;  Prepare for read
;;
GET_LUN, file
filename=datadir+'/'+varfile 
if (not quiet) then print, 'Reading ', filename
dummy=findfile(filename, COUNT=countfile)
if (not countfile gt 0) then begin
  print, 'ERROR: cannot find file '+ filename
  stop
endif
close, file
openr, file, filename, /f77
;;
;;  Read z-averages and put in arrays if requested.
;;
it=0
while (not eof(file)) do begin

  readu, file, t
  readu, file, array
;;
;;  Plot requested variable, variable number 0 by default.
;;
  if (lplot) then begin
    array_plot=array[*,*,iplot]
    ii=where(array_plot gt max) & if (ii[0] ne -1) then array_plot[ii]=max
    ii=where(array_plot lt min) & if (ii[0] ne -1) then array_plot[ii]=min
    tvscl, rebin(array_plot,zoom*[nx,ny])
  endif
;;
;;  Diagnostics.
;;
  if (not quiet) then begin
    if (it eq 0 ) then $
        print, '  ------- it ------- ivar -------- t --------- min(var) ------- max(var) -----'
    for ivar=0,nvar-1 do begin
      print, it, ivar, t, min(array[*,*,ivar]), max(array[*,*,ivar])
    endfor
  endif
;;
;;  Split read data into named arrays.
;;
  if ( it le nit-1 ) then begin
    tt[it]=t
    for ivar=0,nvar-1 do begin
      cmd=varnames[ivar]+'[*,*,it]=array[*,*,ivar]'
      if (execute(cmd,0) ne 1) then message, 'Error putting data in array'
    endfor
  endif
;
  it=it+1
endwhile
;;
;;  Put data in structure.
;;
if (nit ne 0) then begin
  makeobject="object = CREATE_STRUCT(name=objectname,['t'," + $
      arraytostring(varnames,QUOTE="'",/noleader) + "]," + $
      "tt,"+arraytostring(varnames,/noleader) + ")"
;
  if (execute(makeobject) ne 1) then begin
    message, 'ERROR Evaluating variables: ' + makeobject, /INFO
    undefine,object
  endif
endif
;
end
