;  $Id: pc_object_stats.pro,v 1.1 2004-06-03 17:23:27 mee Exp $
;
;  Summarize data
;  omit ghost zones in the analysis
;
pro pc_object_stats,object,tags=tags,skiptags=skiptags, $
                    dim=dim,  $
                    NOVECTORS=NOVECTORS, QUIET=QUIET, $
                    _EXTRA=e
COMPILE_OPT IDL2,HIDDEN

  default,tags,tag_names(object)

  if n_elements(dim) ne 1 then $
      pc_read_dim, obj=dim, QUIET=QUIET, _EXTRA=e

  donefirst=0
  ;
  for iv=0L,n_elements(tags)-1 do begin
    res=execute("pc_variable_stats,object." + tags[iv] + $
                ",varname=tags[iv],NOVECTORS=NOVECTORS," + $
                  "dim=dim,NOHEADER=donefirst")
    donefirst=1
  endfor

end
