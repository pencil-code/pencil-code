;;
;; $Id$
;+
; NAME:
;      pc_read_pstalk
;
; PURPOSE:
;     Read information about local state of the gas around a
;     selected group of particles.
;-
; MODIFICATION HISTORY:
;     Written by: Anders Johansen (johansen@mpia.de) on 13.07.2007
;
pro pc_read_pstalk, object=object, datadir=datadir, it0=it0, it1=it1, $
    swap_endian=swap_endian, quiet=quiet, noutmax=noutmax, single=single, help=help
COMPILE_OPT IDL2,HIDDEN
;
common pc_precision, zero, one, precision, data_type, data_bytes, type_idl
;
  if (keyword_set(help)) then begin
    doc_library, 'pc_read_pstalk'
    return
  endif
;
; Default values.
;
default, swap_endian, 0
default, quiet, 0
default, it1, -1
default, it0, 0
default, single, 0
datadir = pc_get_datadir(datadir)
;
; Read dimensions.
;
pc_read_dim, obj=dim, datadir=datadir, /quiet
;
; Load HDF5 file if requested or available.
;
  files = file_search (datadir+'/allprocs/PSTALK*.h5')
  if (files[0] ne '') then begin
    message, "pc_read_pstalk: WARNING: please use 'pc_read' to load HDF5 data efficiently!", /info
    num_files = n_elements (files)
    for pos = 0, num_files-1 do begin
      t = pc_read ('time', file='PSTALK'+strtrim (pos, 2)+'.h5', datadir=datadir, single=single)
      if (pos eq 0) then begin
        distribution = pc_read ('proc/distribution')
        num_part = total (distribution)
        quantities = h5_content('stalker')
        num_quantities = n_elements (quantities)
        object = { t:replicate(single ? !Values.F_NaN : !Values.D_NaN*zero, num_files), IPAR:pc_read('stalker/ID') }
        data = make_array(num_part, num_files, type=single ? 4 : type_idl)
        for i = 0, num_quantities-1 do begin
          if (strupcase (quantities[i]) eq 'ID') then continue
          object = create_struct (object, quantities[i], data)
        end
      end
      tags = (tag_names (object))
      for i = 0, num_quantities-1 do begin
        if ((strupcase (tags[i]) eq 'T') or (strupcase (tags[i]) eq 'IPAR')) then continue
        found = where (strupcase (quantities) eq strupcase (tags[i]), num_found)
        if (num_found ne 1) then message, 'pc_read_pstalk: ERROR while reading quantites.'
        tmp = object.(i)
        tmp[*,pos] = pc_read('stalker/'+quantities[found],single=single)
        object.(i) = tmp
      end
      object.t[pos] = t
    end
    h5_close_file
    return
  end
;
pc_read_pdim, obj=pdim, datadir=datadir, /quiet
;
; Read parameters from file.
;
pc_read_param, obj=param, datadir=datadir, /quiet, single=single
;
; Stop if no particles have been stalked.
;
if (pdim.npar_stalk eq 0) then begin
  print, 'pc_read_pstalk: npar_stalk is zero - set it in cparam.local and rerun'
  return
endif
;
; Find out whether or not we stalk sink particles.
;
if (max(tag_names(param) eq 'LSTALK_SINK_PARTICLES') eq 0) then begin
  lstalk_sink_particles=0
endif else begin
  lstalk_sink_particles=param.lstalk_sink_particles
endelse
;
; Read the number of output times from file.
;
tout=zero
default, noutmax, -1
openr, lun, datadir+'/tstalk.dat', /get_lun
  readf, lun, tout, noutmaxfile
close, lun
free_lun, lun
if (noutmax eq -1 or noutmax gt noutmaxfile) then nout=noutmaxfile else nout=noutmax
;
; Read header information from file.
;
header=''
openr, lun, datadir+'/particles_stalker_header.dat', /get_lun
  readf, lun, header, format='(A)'
close, lun
free_lun, lun
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
t=make_array(nout, type=single ? 4 : type_idl)
array=make_array(nfields,pdim.npar_stalk,nout, type=single ? 4 : type_idl)
;
; Sink particles have random particle indices, so we need to keep track of the
; particle index for later sorting.
;
if (lstalk_sink_particles) then begin
  ipar_stalk=lonarr(pdim.npar_stalk,nout)
  ipar_stalk[*,*]=-1
  npar_stalk_read=lonarr(nout)
endif else begin
  ipar_stalk=indgen(pdim.npar_stalk)
endelse
;
; Go through all processor directories.
;
for iproc=0,dim.nprocx*dim.nprocy*dim.nprocz-1 do begin
;
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
  openr, lun, datadir+'/proc'+strtrim(iproc,2)+'/particles_stalker.dat', /f77, /get_lun, swap_endian=swap_endian
  while (ntread lt nout and not eof(lun)) do begin
    readu, lun, t_loc, npar_stalk_loc
;
    if (it ge it0) then begin
      if ( (it1 ne -1) and (it mod it1 eq 0) ) then $
          print, iproc, it, t_loc, npar_stalk_loc
;
      if (npar_stalk_loc ge 1) then begin
        ipar_loc=lonarr(npar_stalk_loc)
        readu, lun, ipar_loc
;
        array_loc=make_array(nfields,npar_stalk_loc, type=type_idl)
        readu, lun, array_loc
;
; Put data from local processor into global data structure.
;
        if (not lstalk_sink_particles) then begin
          array[*,ipar_loc-1,it-it0]=array_loc
        endif else begin ; Sink particles are sorted by index later
          array[*,npar_stalk_read[it-it0]: $
              npar_stalk_read[it-it0]+npar_stalk_loc-1,it-it0]=array_loc
          ipar_stalk[npar_stalk_read[it-it0]:npar_stalk_read[it-it0]+ $
              npar_stalk_loc-1,it-it0]=ipar_loc
          npar_stalk_read[it-it0]=npar_stalk_read[it-it0]+npar_stalk_loc
        endelse
      endif
;
      ntread=ntread+1
;
      t[it-it0]=t_loc
    endif else begin
      if (npar_stalk_loc ge 1) then begin
        dummyinteger=0L
        readu, lun, dummyinteger
        dummyreal=zero
        readu, lun, dummyreal
      endif
    endelse
;
    it=it+1
;
  endwhile
  close, lun
  free_lun, lun
;
endfor
;
; Order sink particles by particle index and trim data array.
;
if (lstalk_sink_particles) then begin
  array2=array
  ipar_stalk2=ipar_stalk
  ipar_stalk=ipar_stalk[sort(ipar_stalk)]
  kuniq=ipar_stalk[uniq(ipar_stalk)]
  if (n_elements(kuniq) gt 1) then begin
    kuniq=kuniq[1:n_elements(kuniq)-1]
    for k=0,n_elements(kuniq)-1 do begin
      for it=1,nout-1 do begin
        kk=where(ipar_stalk2[*,it] eq kuniq[k])
        if (kk[0] ne -1) then array[*,k,it]=array2[*,kk,it] $
        else array[*,k,it]=zero
      endfor
    endfor
    ipar_stalk=kuniq
    array=array[*,0:n_elements(ipar_stalk)-1,*]
  endif else begin
    ipar_stalk=-1
    array=array[*,0,*]
  endelse
endif
;
; Build structure of all the variables.
;
command="object = create_struct(name=objectname,['t','ipar'"+ $
    arraytostring(fields,quote="'")+"],t,ipar_stalk"+ $
    arraytostring('reform(array['+strtrim(indgen(nfields),2)+",*,*])")+")"
status=execute(command)

if (status eq 0) then begin
  message, 'Error: building of object failed, but data locally available as t,ipar_stalk,array', /info
  undefine, object
  stop
endif
;
end
