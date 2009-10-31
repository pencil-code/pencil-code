;$Id$
pro pc_coarsegrain,t,f,nevery,tt,ff,aver=aver,firstindex=firstindex
;
;  Coarsegrain data array in the second (or, optionally the first) index.
;  The keyword /aver allows averaging over the coarsegrained direction.
;
;  Calling example
;    pc_coarsegrain,t,alpij,nevery,tt,alpijc,/aver
;
;  Determine time indices. make sure we go all the way to the end.
;  First downsample time array, 
;
nt=n_elements(t) & ntout=nt/nevery
it2=nt-1 & it1=nt-nevery*ntout
tt=reform((reform(t(it1:it2),nevery,ntout))(nevery-1,*))
print,'coarsegrain: it1,it2,nt=',it1,it2,nt
s=size(f)
;
;  consider option of coarsegraining over first index
;
if keyword_set(firstindex) then begin
  if s[1] ne nt then stop,'wrong f array size'
  nz=s[2]
  if keyword_set(aver) then begin
    ff=total(reform(f(it1:it2,*),nevery,ntout,nz),1)/nevery
  endif else begin
    ff=reform((reform(f(it1:it2,*),nevery,ntout,nz))(nevery-1,*,*))
  endelse
endif else begin
;
;  default: coarsegraining over second index
;  allow fow up to 2 additional indices after the time index
;
  nz=s[1]
  print,'size=',s
  if s[0] eq 2 then begin
    if keyword_set(aver) then begin
      ff=total(reform(f(*,it1:it2),nz,nevery,ntout),2)/nevery
    endif else begin
      ff=reform((reform(f(*,it1:it2),nz,nevery,ntout))(*,nevery-1,*))
    endelse
  endif else if s[0] eq 3 then begin
    i1=s[3]
    if keyword_set(aver) then begin
      ff=total(reform(f(*,it1:it2,*),nz,nevery,ntout,i1),2)/nevery
    endif else begin
      ff=reform((reform(f(*,it1:it2,*),nz,nevery,ntout,i1))(*,nevery-1,*,*))
    endelse
  endif else if s[0] eq 4 then begin
    i1=s[3]
    i2=s[4]
    if keyword_set(aver) then begin
      ff=total(reform(f(*,it1:it2,*,*),nz,nevery,ntout,i1,i2),2)/nevery
    endif else begin
      ff=reform((reform(f(*,it1:it2,*,*),nz,nevery,ntout,i1,i2))(*,nevery-1,*,*,*))
    endelse
  endif
endelse
;
end
