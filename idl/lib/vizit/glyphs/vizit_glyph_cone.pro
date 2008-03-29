function vizit_glyph_cone
  default,nphi,3
  default,headsize,0.7
  default,stemradius,0.1

  vertices=fltarr(3,nphi*3+2)  
  triangles=lonarr(nphi*24)  
  ipnt=0

  vertices[0,ipnt]=-0.5 
  vertices[1,ipnt]=0. 
  vertices[2,ipnt]=0. 
  ipnt=ipnt+1 

  for i=0,nphi-1 do begin
   vertices[0,ipnt]=-0.5 
   vertices[1,ipnt]=stemradius*sin(2.*!dpi*i/float(nphi)) 
   vertices[2,ipnt]=stemradius*cos(2.*!dpi*i/float(nphi)) 
   ipnt=ipnt+1 
  endfor

  vertices[0,ipnt]=0.5 
  vertices[1,ipnt]=0. 
  vertices[2,ipnt]=0. 
  ipnt=ipnt+1 

  itri=0
; The base
  for i=0,nphi-1 do begin
   triangles[itri  ]=3
   triangles[itri+1]=0
   triangles[itri+2]=i+1
   triangles[itri+3]=((i+1) mod nphi) + 1
   itri=itri+4 
  endfor


; The head
  for i=0,nphi-1 do begin
   triangles[itri  ]=3
   triangles[itri+1]=i+1
   triangles[itri+2]=nphi+1
   triangles[itri+3]=((i+1) mod nphi)+1
   itri=itri+4 
  endfor
  object=CREATE_STRUCT(['vertices','triangles'],ptr_new(vertices),ptr_new(triangles))
  return,object
end
