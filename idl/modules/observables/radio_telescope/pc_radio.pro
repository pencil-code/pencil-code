function pc_radio_beamsmooth,var,beam_width_arcseconds,dx=dx,dy=dy,nsigma=nsigma,xmin=xmin,ymin=ymin,dist_to_obj=dist_to_obj
;FORWARD_FUNCTION filter_image
;
;  Take A and smooth with gaussian of (HPBW) beam_width;  the array B is used to
;  store the smmothed A until all calculations are done.  beam_width should be 
;  provided in angular arcseconds, and then recalculated here from the provided 
;  input / parameters. 
;
; require nx,ny,dx,dy(,xmin,ymin ? - maybe not; just need spacing dx,dy)
;   
    
COMPILE_OPT IDL2,HIDDEN
;

fxln2 = 2.772588722239780 ; = 4*ln(2)
one_arcsec=4.848137e-6


default,nsigma,5.  ; use 5.0\sigma smoothing by default
default,dx,1.
default,dy,1.
default,xmin,0.
default,ymin,0.

;    get beam_width corresponding to sky_length
;
; so providing the units of dx and dist_to_object are identical we are good !
;
  beam_width=beam_width_arcseconds*dist_to_obj*one_arcsec
  return, filter_image( var, FWHM = [beam_width/dx,beam_width/dy], /ALL )

;    units are pc by default; for kpc specify length_unit=1000.0
;
;    bundle the beam/weighting stuff into constants
;      
  gaussian_const= -1.*fxln2/(beam_width^2)
  gaussian_cutoff2= nsigma*nsigma*beam_width^2
  gaussian_weighting=fxln2*dx*dy/(!pi*beam_width^2)
  if ( (beam_width^2) lt (0.25*(dx^2+dy^2)) ) then return,var
;
;   do the smoothing
;  
  print,'Smoothing to HPBW',beam_width_arcseconds,'"'

  nx=(size(var))[1]
  ny=(size(var))[2]
  smoothed=fltarr(nx,ny)  ; Initialize smoothed to be the size of var


  for j=0,ny - 1 do begin
  yy=spread(ymin+(findgen(ny)-j)*dy,0,nx)

  for i=0,nx - 1 do begin
  xx=spread(xmin+(findgen(nx)-i)*dx,1,ny)
  rr2=xx^2+yy^2
  pnts=where(rr2 le gaussian_cutoff2,siz)
  if siz gt 0 then smoothed[i,j]=smoothed[i,j]+total(var[pnts]*exp(rr2[pnts]*gaussian_const))

;    for jj=0,(size(var))[2] - 1 do begin
;      py = dy*(j-jj) 
;      exp1 = exp(py^2*gaussian_const)
;    for ii=0,(size(var))[1] - 1 do begin
;      px = dx*(i-ii)
;      if (px^2 + py^2) lt gaussian_cutoff2 then begin
;        exp2 = exp1*exp(px*px*gaussian_const)
;        smoothed(i,j) = smoothed(i,j) + var(ii,jj) * exp2
;      endif
;    endfor
;    endfor

  endfor
if ((j mod 10) eq 0) then print,j*100./((size(var))[2] - 1.0),'% complete'
  endfor
help,smoothed
  return,smoothed*gaussian_weighting
      
end
 
pro pc_radio,object=object,lambdas=lambdas,dim=dim,param=param,data=data,units=units,test=test,HPBW=HPBW,_extra=e
;  * We require output from simulation: rho, Bx, By, Bz, 
;    ionization frac (optional)
;  * Need to know the size of such arrays at compile time(?), then read the 
;    data from a named file.
;  * Need to initialise the following arrays:    
;    
;     ___2D______________________ 
;
;     Qone, Qtwo, Uone, Utwo
;
;     ___3D______________________
;
;     (Bx, By, Bz) .or. (B(1:3)), rho, epsilon, faraday_depth (perhaps),
;     ionization_frac (optional) 


@data/pc_constants.pro
  default,lambdas,[0.035,0.062]

  default,HPBW,0. ; was 15.
; Err... What are these?
  default,distance_to_object,100.0  ; in parsecs
  default,c1,1.               ; Actually a fn of n_cr 
  default,c2,1.62             ; c2=1.0*2.0*0.81=2*K
  default,gamma,3.
  default,p,.75

if not keyword_set(TEST) then begin
  if (n_elements(dim) ne 1) then pc_read_dim,object=dim,_extra=e
  if (n_elements(param) ne 1) then pc_read_param,object=param,dim=dim,_extra=e
  if (n_elements(units) ne 1) then pc_units,object=units,dim=dim,param=param,_extra=e
  if (n_elements(data) ne 1) then pc_read_var,object=data,variables=['curl(aa)','exp(lnrho)'],tags=['bb','rho'],dim=dim,param=param,_extra=e

  default,distance_to_object,0.100

  m_H=1.6726D-24
  mu=param.mu
  dx=data.dx ;*(units.length/3.086E18)
  dy=data.dy ;*(units.length/3.086E18)
  dz=data.dz ;*(units.length/3.086E18)
;  nx=dim.nx
;  ny=dim.ny
;  nz=dim.nz
;  l1=dim.l1
;  l2=dim.l2
;  m1=dim.m1
;  m2=dim.m2
;  n1=dim.n1
;  n2=dim.n2
  nx=dim.mx
  ny=dim.my
  nz=dim.mz
  l1=0
  l2=nx-1  
  m1=0
  m2=ny-1  
  n1=0
  n2=nz-1  
  xmin=param.xyz0[0]-l1*dx
  ymin=param.xyz0[1]-m1*dy

;  rescale density so that it is n_e
  n_e=data.rho*(units.density/(mu*m_H))
  bb=data.bb*units.magnetic_field*1E+6
;  make n_cr the numerical value of B^2
  n_cr=dot2(bb)

  length_unit=units.length/3.086D18
 ; length_unit=1.
endif else begin


  nx=150 & ny=150 & nz=10
  dx=1./128. & dy=1./128. & dz=1.
  xmin=0. & ymin=0.
  l1=0 & m1=0 & n1=0
  l2=nx-1 & m2=ny-1 & n2=nz-1
  bb=fltarr(nx,ny,nz,3)
  rho=fltarr(nx,ny,nz)
  mu=1.0
  m_H=1.6726
; create a purely azimuthal field
  Br=0.
  Bphi=1.
  for j=m1,m2 do begin
  for i=l1,l2 do begin
    xx=(i-75)*1.01
    yy=(j-75)*1.01
    phi=atan(yy,xx)
    if (phi lt (0.0D0)) then phi=phi+(2.0*!pi)
    bb(i,j,*,0)= Br*cos(phi) - Bphi*sin(phi)
    bb(i,j,*,1)= Br*sin(phi) + Bphi*cos(phi)
    bb(i,j,*,2)=0.
  endfor
  endfor  


  n_e=fltarr(nx,ny,nz)
  n_e[*,*,*]=1.0
;     n_e=1.0
;  make n_cr the numerical value of B^2
  n_cr=dot2(bb)
  length_unit=1.
endelse
;
;
;    make sure c2 is scaled so that it is in pc...etc..
;
  c2=c2*length_unit
;
;
;
;
;
; Prepare some arrays
  eepsilon = fltarr(nx,ny,nz) 
  angle = fltarr(nx,ny,nz)
  Q=fltarr(nx,ny,n_elements(lambdas)) 
  U=fltarr(nx,ny,n_elements(lambdas))
  QQ=fltarr(nx,ny) 
  UU=fltarr(nx,ny)
;
;  loop over 2 values of lambda and obtain Q and U for each
;
  for nl=0,n_elements(lambdas)-1 do begin
    lambda=lambdas[nl]
;    
;   *  Calculate epsilon, the total synchrotron emissivity
;
    eepsilon=(sqrt(bb[l1:l2,m1:m2,n1:n2,0]^2  +  bb[l1:l2,m1:m2,n1:n2,1]^2))^((gamma+1.0)/2.0)
    eepsilon=eepsilon*(lambda^((gamma-1.0)/2.0))
    eepsilon=eepsilon*n_cr[l1:l2,m1:m2,n1:n2]*c1       ; c1 ~ 1.0
;         
;   *  Calculate the Faraday depth   
;   *  Calculate the effective local angle ( 2*intrinsic + F)        
;       

    faraday_depth=fltarr(nx,ny)   ; use a cumulative form of 'depth'
    faraday_depth[*,*]=0D

    for k=nz-1,0,-1 do begin                ; do the sum for 'depth' backwards to save time
      faraday_depth=faraday_depth +                      $
          c2*lambda^2*n_e[l1:l2,m1:m2,n1+k]*bb[l1:l2,m1:m2,n1+k,2]*dz


      
; Angle = 2 * intrinsic polatisation angle + faraday_depth
;      intr_angle=atan(bb[l1:l2,m1:m2,n1+k,1]/bb[l1:l2,m1:m2,n1+k,0]+(1e-9))+0.5*!pi
      intr_angle=atan(bb[l1:l2,m1:m2,n1+k,1],bb[l1:l2,m1:m2,n1+k,0])+0.5*!pi
      pnts=where(intr_angle lt 0.*!pi,siz)
      if siz gt 0 then intr_angle[pnts]=intr_angle[pnts]+!pi
      pnts=where(intr_angle gt !pi,siz)
      if siz gt 0 then intr_angle[pnts]=intr_angle[pnts]-!pi
     
      angle[*,*,k]=2.0*intr_angle  - faraday_depth
    endfor
;
;   *  Calculate the Q and U by summing along the line of sight
;   
    Q[*,*,nl]=0D
    U[*,*,nl]=0D
    for k=0,nz-1 do begin
      Q[*,*,nl]=Q[*,*,nl] + p*eepsilon[*,*,k]*cos(angle[*,*,k])*dz
      U[*,*,nl]=U[*,*,nl] + p*eepsilon[*,*,k]*sin(angle[*,*,k])*dz
    endfor
  endfor

if (HPBW gt 0.) then begin
  for nl=0,n_elements(lambdas)-1 do begin
    Q[*,*,nl]=pc_radio_beamsmooth(Q[*,*,nl],HPBW,dx=dx,dy=dy,xmin=xmin,ymin=ymin,dist_to_obj=distance_to_object/length_unit)
    U[*,*,nl]=pc_radio_beamsmooth(U[*,*,nl],HPBW,dx=dx,dy=dy,xmin=xmin,ymin=ymin,dist_to_obj=distance_to_object/length_unit)
  endfor
endif

;
; Polarization angle
;
  polarization_angle=0.5*atan(U,Q)
  pnts=where(polarization_angle lt 0.,siz)
  if siz gt 0 then polarization_angle[pnts]=polarization_angle[pnts]+!pi
;  polarization_angle=polarization_angle-0.5*!pi ;Don't bother... Will only have to add it back on later.

;   So at this point we have the observed angle of polarisation
;   for the E-vector
;AJWM OR DO WE??? 
;AJWM... Not if we calculate it as above!! 
  polarization_angle=polarization_angle + 0.5*!pi

; Polarized intensity
  P_I=sqrt(Q^2 + U^2)

; Rotation measure - NEEDS CHECKING
  RM=polarization_angle[*,*,1]-polarization_angle[*,*,0]
  pnts=where(RM gt (!pi/2.0),siz)
  if siz gt 0 then RM[pnts]=RM[pnts]-!pi
  pnts=where(RM lt (-!pi/2.0),siz)
  if siz gt 0 then RM[pnts]=RM[pnts]+!pi
  RM=RM/(lambdas[1]^2-lambdas[0]^2)

  bad_pnts0=where(P_I[*,*,0] lt 0.1*mean(P_I[*,*,0]))
  bad_pnts1=where(P_I[*,*,1] lt 0.1*mean(P_I[*,*,1]))

;  RM[bad_pnts0]=!values.f_nan
;  RM[bad_pnts1]=!values.f_nan
; Line of sight averaged B Field
  bb_los=fltarr(nx,ny,3)
  for i=0,nz-1 do bb_los = bb_los + bb[l1:l2,m1:m2,i,*]
  bb_los=bb_los/nz
  
  bb_perp=fltarr(nx,ny,2,2)
  bb_perp[*,*,0,0]=P_I[*,*,0]*cos(polarization_angle[*,*,0])
  bb_perp[*,*,1,0]=P_I[*,*,0]*sin(polarization_angle[*,*,0])
  bb_perp[*,*,0,1]=P_I[*,*,1]*cos(polarization_angle[*,*,1])
  bb_perp[*,*,1,1]=P_I[*,*,1]*sin(polarization_angle[*,*,1])



  bb_parallel=fltarr(nx,ny,1)
  for j=0,ny-1 do begin 
  for i=0,nx-1 do begin 
    bb_parallel[i,j]=RM[i,j]/(nz*dz*total(n_e[i,j,*])*0.81)
  endfor
  endfor

  bb_x=bb_perp[*,*,0,0]
  bb_y=bb_perp[*,*,1,0]

;  bb_x[bad_pnts0]=!values.f_nan
;  bb_x[bad_pnts1]=!values.f_nan
;  bb_y[bad_pnts0]=!values.f_nan
;  bb_y[bad_pnts1]=!values.f_nan

  j_parallel=der(spread(bb_y,2,1),0)-der(spread(bb_x,2,1),1)
  h_c=bb_parallel*j_parallel*dx*dy

;wind=0
;for i=0,(size(lambdas))[1]-1 do begin
;  window,wind
;
;  contour,polarization_angle[*,*,i],nlevels=20,/fill
;  print,minmax(polarization_angle[*,*,i])
;  wind=wind+1
;  window,wind
;  contour,P_I[*,*,i],nlevels=20,/fill
;  print,minmax(P_I[*,*,i])
;  wind=wind+1
;  window,wind
;  contour,Q[*,*,i],nlevels=20,/fill
;  print,minmax(Q[*,*,i])
;  wind=wind+1
;  window,wind
;  print,minmax(U[*,*,i])
;  wind=wind+1
;endfor
;contour,RM,nlevels=64,/fill

object=CREATE_STRUCT(['polarization_angle','polarized_intensity', $
                      'bb_los','bb_perp','bb_parallel','Q','U','RM', $
                      'j_parallel','h_c',                  $
                      'lambda','nlambda'],                 $
                       polarization_angle,P_I,             $
                       bb_los,bb_perp,bb_parallel,Q,U,RM,  $
                       j_parallel,h_c,                     $
                       lambdas,n_elements(lambdas))

end


