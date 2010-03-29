;$Id: cwd.pro,v 1.1 2009-02-25 08:57:44 brandenb Exp $
;
;  this routine returns the name of the directory
;  without the rest of the full path
;
pro cwd,run
spawn,'echo $cwd',path
b=byte(path)
n=n_elements(b)
slash=max(where(b eq 47))
run=string(b(slash+1:n-1))
END
