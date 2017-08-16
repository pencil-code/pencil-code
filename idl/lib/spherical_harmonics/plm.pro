FUNCTION plm,th,n,m
;
;  Legendre polynomials
;
      nth=n_elements(th)
      plm=replicate(0.,nth)
      IF (M GT N) then RETURN,plm
      PM0=replicate(0.,nth)
      PM1=replicate(1.,nth)
;
      IF (M ne 0) then begin
        DS=SIN(TH)
        for I=1,M do PM1=PM1*(2*I-1)*DS*(-1.)
      endif
;
      plm=PM1
      IF (M EQ N) then RETURN,plm
      DC=COS(TH)
      M1=M-1
      I1=M+1
      for I=I1,N do begin
        plm=((2*I-1)*DC*PM1-(I+M1)*PM0)/FLOAT(I-M)
        PM0=PM1
        PM1=plm
      endfor
;
return,plm
END
