;;
;; $Id: pc_read_pstalk.pro,v 1.2 2007-11-14 07:18:26 ajohan Exp $
;;
;; NAME:
;;      pc_read_pstalk
;;
;; PURPOSE:
;;     Read information about local state of the gas around a
;;     selected group of particles.
;;
;; MODIFICATION HISTORY:
;;     Written by: Anders Johansen (johansen@mpia.de) on 13.07.2007
;;
pro pc_read_pstalk, object=object, datadir=datadir, quiet=quiet
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
; Default values.
;
default, quiet, 0
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
;
; Read dimensions and set precision.
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
pc_read_pdim, obj=pdim, datadir=datadir, /quiet
pc_set_precision, dim=dim, datadir=datadir, /quiet
;
; Read the number of output times from file.
;
tout=zero
nout=0
openr, 1, datadir+'/tstalk.dat'
  readf, 1, tout, nout
close, 1
;
; Read header information from file.
;
header=''
openr, 1, datadir+'/particles_stalker_header.dat'
  readf, 1, header, format='(A)'
close, 1
;
; Extract fields from header in order to know what to read from file.
;
fields=strsplit(header,',',/extract)
nfields=n_elements(fields)
fields_loc=fields+'_loc'
if (not quiet) then begin
  print, 'Going to read the '+strtrim(nfields,2)+' fields: '
  print, '  ', fields
  print, 'at ', strtrim(nout,2), ' times'
endif
;
; Initialize data arrays.
;
t=fltarr(nout)*zero
for ifield=0,nfields-1 do begin
  command=fields[ifield]+'=fltarr(pdim.npar_stalk,nout)*zero' & status=execute(command)
  command=fields_loc[ifield]+'=zero' & status=execute(command)
endfor
;
; Go through all processor directories.
;
for iproc=0,dim.nprocx*dim.nprocy*dim.nprocz-1 do begin
;
; Initialize variables.
;
  it=0
  t_loc=zero
  npar_stalk_loc=0L
  ipar=0L
;
  openr, 1, datadir+'/proc'+strtrim(iproc,2)+'/particles_stalker.dat', /f77
    while (it lt nout and not eof(1)) do begin
      readu, 1, t_loc, npar_stalk_loc

      for k=0,npar_stalk_loc-1 do begin
        command='readu, 1, ipar'
        for ifield=0,nfields-1 do begin
          command=command+', '+fields_loc[ifield]
        endfor
        status=execute(command)
        for ifield=0,nfields-1 do begin
          command=fields[ifield]+'[ipar-1,it]='+fields_loc[ifield]
          status=execute(command)
        endfor
      endfor

      t[it]=t_loc

      it=it+1
    endwhile
  close, 1

endfor
;
; Build structure of all the variables.
;
command="object = create_struct(name=objectname,['t'"+ $
    arraytostring(fields,quote="'")+"],t"+arraytostring(fields)+")"
status=execute(command)
;
end
