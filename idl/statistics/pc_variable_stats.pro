;;
;;  $Id$
;;
;;  Summarize data
;;  omit ghost zones in the analysis
;;
pro pc_variable_stats,variable, struct=struct, tag=tag, varname=varname, dim=dim, $
    MINIMUM=MINIMUM_, MAXIMUM=MAXIMUM_, MEAN=MEAN_, RMS=RMS_, $
    TRIM=TRIM, NOHEADER=NOHEADER, NOVECTORS=NOVECTORS, YINYANG=yinyang, _EXTRA=e
COMPILE_OPT IDL2,HIDDEN
;
;  Statistics for scalar and vector variables.
;  Yet missing: multicomponent non-vector variables.
;
if not keyword_set(yinyang) then yinyang=0 
yinyang = logical_true(yinyang)
;
xyz = ['x', 'y', 'z']
;
;;;statvar = (n_params() eq 1 ? "variable" : "struct."+strtrim(tag))

stat_headings=[ $
               '  name            ',$
               '   ',             $
               'minval         ',$
               'maxval          ',$
               'mean         ',$
               'stddev' ]
stats        =[ $
               "' '+strmid(strlowcase(varname)+'          ',0,10)",       $
               "'-->'",         $
               'min(statvar)', $
               'max(statvar)', $
               'mean(statvar,/double)',$
               ;'stddev(statvar)' ]
               ;
               ;  work-around to prevent stddev from exiting on scalars
               ;
               'stddev(statvar*[1,1],/double)' ]
stat_formats =[ $
               'A12',$
               'A',$
               '1G15.6',$
               '1G15.6',$
               '1G15.6',$
               '1G15.6' ]
;
requested_stats=make_array(6,/INT,VALUE=0)
;
if (keyword_set(MINIMUM_) or $
    keyword_set(MAXIMUM_) or $
    keyword_set(MEAN_)    or $
    keyword_set(RMS_) ) then begin
  requested_stats[*]=0
  if keyword_set(MINIMUM_) then requested_stats[2]=1
  if keyword_set(MAXIMUM_) then requested_stats[3]=1
  if keyword_set(MEAN_)    then requested_stats[4]=1
  if keyword_set(RMS_)     then requested_stats[5]=1
endif else begin
  requested_stats[*]=1
endelse
;
if (n_elements(varname) eq 1) then begin
  requested_stats[0:1]=1
endif else begin
  requested_stats[0:1]=0
endelse
;
if (n_elements(dim) ne 1) then pc_read_dim, obj=dim, _EXTRA=e
;
if (not keyword_set(NOHEADER)) then $
    print,arraytostring(stat_headings[where(requested_stats)],LIST='')
;
;  We do not want averages of grid and time parameters.
;
if (varname eq 'X')  or (varname eq 'Y')  or (varname eq 'Z') or $
   (varname eq 'DX') or (varname eq 'DY') or (varname eq 'DZ') then return
  if yinyang then begin
    if strpos(varname,'MERGE') ge 0 then return
    if (varname eq 'YZ') or (varname eq 'TRIANGLES') or (varname eq 'SPHERE_DATA') or $
       (varname eq 'YGHOSTS') or (varname eq 'ZGHOSTS') or (varname eq 'FVAL') then return
  endif
;
if n_params() eq 1 then $
  varsize=size(variable) $
else $
  if keyword_set(struct) and keyword_set(tag) then $
    varsize=sizeoftag(struct,tag) $
  else begin
    print, 'pc_variable_stats: Error - no variable given!!!' & stop
  endelse
;
fmt='('+arraytostring(stat_formats[where(requested_stats)],LIST=',',/NOLEAD)+')'
;if ( pc_is_scalarfield(variable, dim=dim, $
;    subscripts=subscripts, TRIM=TRIM, NOVECTORS=NOVECTORS, YINYANG=yinyang) ) then begin

if ( pc_is_vectorfield(varsize, dim=dim, $
   subscripts=subscripts, TRIM=TRIM, NOVECTORS=NOVECTORS, YINYANG=yinyang) eq 1) then begin
;
;  The array is a vector field...
;
  for j=0,2 do begin
    subscripts[varsize[0]-1-yinyang]=string(j)
    stats[0]="' '+strmid(strlowcase(varname)+'_'+xyz[j]+'          ',0,10)"
;
    if (execute("statvar="+(n_params() eq 1 ? "variable" : "struct."+strtrim(tag))+"[" + $
        arraytostring(subscripts,/NOLEADER)+"]",1) ne 1) then $
        message, 'Error printing stats for ' + varname 
;
    cmd = "print, FORMAT=fmt" + arraytostring(stats[where(requested_stats)])  
    if (execute(cmd,1) ne 1) then message, 'Error printing stats for '+stats[0] 
  endfor
;
endif else begin
;
;  The array is a scalar field....
;
  stats[0]="' '+strmid(strlowcase(varname)+'          ',0,10)"
;
  if (execute("statvar="+(n_params() eq 1 ? "variable" : "struct."+strtrim(tag))+"[" + $
      arraytostring(subscripts,/NOLEADER)+"]",1) ne 1) then $
      message, 'Subscripting array  ' + varname
;
  cmd = "print, FORMAT=fmt" + $
      arraytostring(stats[where(requested_stats)])
  if (execute(cmd,1) ne 1) then message, 'Error printing stats for ' + stats[0] 

endelse
;
end
