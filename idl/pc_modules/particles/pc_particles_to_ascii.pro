;
;  $Id: pc_particles_to_ascii.pro,v 1.1 2007-07-27 10:35:27 ajohan Exp $
;
;  Output particles positions in ascii file.
;
;  Author: Anders Johansen
;
pro pc_particles_to_ascii, xxp, filename=filename, npar=npar, $
    lwrite_tauf=lwrite_tauf

default, filename, './particles.dat'
default, npar, n_elements(xxp[*,0])
default, lwrite_tauf, 0

close, 1
openw, 1, filename

  for ipar=0L,npar-1 do begin
    if (lwrite_tauf) then begin
      if (ipar lt 1999999) then begin
        tauf=0.25
      endif else if (ipar lt 3999999) then begin
        tauf=0.50
      endif else if (ipar lt 5999999) then begin
        tauf=0.75
      endif else begin
        tauf=1.00
      endelse
      printf, 1, xxp[ipar,*], tauf, format='(4f9.4)'
    endif else begin
      printf, 1, xxp[ipar,*], format='(3f9.4)'
    endelse
  endfor

close, 1

end
