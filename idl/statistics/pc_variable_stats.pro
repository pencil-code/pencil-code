;  $Id: pc_variable_stats.pro,v 1.1 2004-06-03 17:23:27 mee Exp $
;
;  Summarize data
;  omit ghost zones in the analysis
;
pro pc_variable_stats,variable, varname=varname, $
                    dim=dim, $
                    MINIMUM=MINIMUM_, $ 
                    MAXIMUM=MAXIMUM_, $
                    MEAN=MEAN_,       $
                    RMS=RMS_,         $
                    NOHEADER=NOHEADER, $
                    NOVECTORS=NOVECTORS, _EXTRA=e
COMPILE_OPT IDL2,HIDDEN

xyz = ['x', 'y', 'z']

stat_headings=[ $
               '  name            ',$
               '   ',             $
               'minval         ',$
               'maxval          ',$
               'mean           ',$
               'rms' ]
stats        =[ $
               "' '+strmid(strlowcase(varname)+'          ',0,10)",       $
               "'-->'",         $
               'min(statvar)', $
               'max(statvar)', $
               'mean(statvar)',$
               'rms(statvar)' ]
stat_formats =[ $
               'A10',$
               'A',$
               '1G15.6',$
               '1G15.6',$
               '1G15.6',$
               '1G15.6' ]
requested_stats=make_array(6,/INT,VALUE=0)

if (keyword_set(MINIMUM_) or $
    keyword_set(MAXIMUM_) or $
    keyword_set(MEAN_)    or $
    keyword_set(RMS_)        $
) then begin
  requested_stats[*]=0
  if keyword_set(MINIMUM_) then requested_stats[2]=1
  if keyword_set(MAXIMUM_) then requested_stats[3]=1
  if keyword_set(MEAN_)    then requested_stats[4]=1
  if keyword_set(RMS_)     then requested_stats[5]=1
endif else begin
  requested_stats[*]=1
endelse

if n_elements(varname) eq 1 then begin
  requested_stats[0:1]=1
endif else begin
  requested_stats[0:1]=0
endelse


if n_elements(dim) ne 1 then $
      pc_read_dim, obj=dim, _EXTRA=e

if not keyword_set(NOHEADER) then $
    print,arraytostring(stat_headings[where(requested_stats)],LIST='')

fmt='('+arraytostring(stat_formats[where(requested_stats)],LIST=',',/NOLEAD)+')'
  ;
  ;
  varsize=size(variable)


  if pc_is_scalarfield(variable, dim=dim,     $
                        subscripts=subscripts, $
                        NOVECTORS=NOVECTORS) then begin ; It is a scalar field
    stats[0]="' '+strmid(strlowcase(varname)+'          ',0,10)"
 
    if (execute("statvar=variable[" + $
                  arraytostring(subscripts,/NOLEADER)+"]",1) ne 1) then $
      message, 'Subscripting array  ' + varname 

    cmd = "print, FORMAT=fmt" + $
                          arraytostring(stats[where(requested_stats)])  
    if (execute(cmd,1) ne 1) then $
      message, 'Error printing stats for ' + stats[0] 
   endif $
  else if pc_is_vectorfield(variable,dim=dim,      $
                             subscripts=subscripts, $
                             NOVECTORS=NOVECTORS) then begin 
                                                     ; It is a vector field
    for j=0,2 do begin
      subscripts[3]=str(j) 
      stats[0]="' '+strmid(strlowcase(varname)+'_'+xyz[j]+'          ',0,10)"

      if (execute("statvar=variable[" + $
                  arraytostring(subscripts,/NOLEADER)+"]",1) ne 1) then $
        message, 'Error printing stats for ' + varname 

      cmd = "print, FORMAT=fmt" + $
                            arraytostring(stats[where(requested_stats)])  
      if (execute(cmd,1) ne 1) then $
        message, 'Error printing stats for ' + stats[0] 
    endfor
  endif

end
