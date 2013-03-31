;$Id$
;
;  this routine returns the name of the directory
;  without the rest of the full path
;
pro cwd,run,debug=debug
spawn,'pwd',path
if keyword_set(debug) then print,'path="'+path+'"'
b=byte(path)
n=n_elements(b)
slash=max(where(b eq 47))
run=string(b(slash+1:n-1))
END
