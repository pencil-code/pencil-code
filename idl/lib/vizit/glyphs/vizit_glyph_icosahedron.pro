function vizit_glyph_icosahedron
  vertices=fltarr(3,12)
  normals=fltarr(3,12)
  triangles=lonarr(20*4)  

      c0 = 0.0
      c1 = 0.9512
      c2 = 0.8507
      c3 = 0.8090
      c4 = 0.4253
      c5 = 0.2629
      c6 = 0.5
      c7 = 0.6882

   vertices[0:2,0]=[c0,c0,c1]
   vertices[0:2,1]=[c0,c2,c4]
   vertices[0:2,2]=[c3,c5,c4]
   vertices[0:2,3]=[c6,-c7,c4]
   vertices[0:2,4]=[-c6,-c7,c4]
   vertices[0:2,5]=[-c3,c5,c4]
   vertices[0:2,6]=[-c6,c7,-c4]
   vertices[0:2,7]=[c6,c7,-c4]
   vertices[0:2,8]=[c3,-c5,-c4]
   vertices[0:2,9]=[c0,-c2,-c4]
   vertices[0:2,10]=[-c3,-c5,-c4]
   vertices[0:2,11]=[c0,c0,-c1]

  ;if keyword_set(rounded) then begin
    normals=vertices/spread(dot2(transpose(vertices)),0,3)
  ;endif

  triangles = [3,0,1,2, 3,0,2,3, 3,0,3,4, 3,0,4,5, 3,0,5,1, 3,1,7,2, 3,2,8,3, $
               3,3,9,4, 3,4,10,5, 3,5,6,1, 3,1,6,7, 3,2,7,8, 3,3,8,9, $
               3,4,9,10, 3,5,10,6, 3,6,11,7, 3,7,11,8, 3,8,11,9, $
               3,9,11,10, 3,10,11,6 ]


  object=CREATE_STRUCT(['vertices','triangles','normals'],ptr_new(vertices),ptr_new(triangles),ptr_new(normals))
  return,object
end
