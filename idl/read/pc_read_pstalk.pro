;;
;; $Id$
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
pro pc_read_pstalk, object=object, datadir=datadir, it0=it0, it1=it1, $
    swap_endian=swap_endian, quiet=quiet, noutmax=noutmax
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
; Default values.
;
default, swap_endian, 0
default, quiet, 0
default, it1, -1
default, it0, 0
if (not keyword_set(datadir)) then datadir=pc_get_datadir()
;
; Read dimensions and set precision.
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
pc_read_pdim, obj=pdim, datadir=datadir, /quiet
pc_set_precision, dim=dim, datadir=datadir, /quiet
;
; Stop if no particles have been stalked.
;
if (pdim.npar_stalk eq 0) then begin
  print, 'pc_read_pstalk: npar_stalk is zero - set it in cparam.local and rerun'
  return
endif
;
; Read the number of output times from file.
;
tout=zero
default, noutmax, -1
openr, 1, datadir+'/tstalk.dat'
  readf, 1, tout, noutmaxfile
close, 1
if (noutmax eq -1 or noutmax gt noutmaxfile) then nout=noutmaxfile else nout=noutmax
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
  print, 'at ', strtrim(string(nout,format='(i)'),2), ' times'
endif
;
; Initialize data arrays.
;
t=fltarr(nout)*zero
array=fltarr(nfields,pdim.npar_stalk,nout)*zero
;
; Go through all processor directories.
;
for iproc=0,dim.nprocx*dim.nprocy*dim.nprocz-1 do begin

  if (not quiet) then begin
    print, 'Reading data from processor '+ $
        strtrim(iproc,2)+'/'+strtrim(dim.nprocx*dim.nprocy*dim.nprocz-1,2)
    print, '-------- iproc ------ it --------- t ----------- npar ------- '
  endif
;
; Initialize variables.
;
  it=0
  ntread=0
  t_loc=zero
  npar_stalk_loc=0L
  ipar=0L
;
  openr, 1, datadir+'/proc'+strtrim(iproc,2)+'/particles_stalker.dat', /f77, $
      swap_endian=swap_endian
  while (ntread lt nout and not eof(1)) do begin
    readu, 1, t_loc, npar_stalk_loc
;
    if (it ge it0) then begin
      if ( (it1 ne -1) and (it mod it1 eq 0) ) then $
          print, iproc, it, t_loc, npar_stalk_loc
;
      if (npar_stalk_loc ge 1) then begin
        ipar_loc=lonarr(npar_stalk_loc)
        readu, 1, ipar_loc
;
        array_loc=fltarr(nfields,npar_stalk_loc)*zero
        readu, 1, array_loc
        array[*,ipar_loc-1,it-it0]=array_loc
      endif
;
      ntread=ntread+1
;
      t[it-it0]=t_loc
    endif else begin
      if (npar_stalk_loc ge 1) then begin
        dummyinteger=0L & readu, 1, dummyinteger
        dummyreal=zero  & readu, 1, dummyreal
      endif
    endelse
;
    it=it+1
;
  endwhile
  close, 1
;
endfor
;
; Build structure of all the variables.
;
command="object = create_struct(name=objectname,['t'"+ $
    arraytostring(fields,quote="'")+"],t"+ $
    arraytostring('reform(array['+strtrim(indgen(nfields),2)+",*,*])")+")"
status=execute(command)
;
end
