function vizit_glyph_cube,rounded=rounded
  vertices=fltarr(3,24)
  normals=fltarr(3,24)
  triangles=lonarr(12)  

  normals[0:2,0:23]=0.

  vertices[0:2,0]=[0,0,0]
  vertices[0:2,1]=[0,0,1]
  vertices[0:2,2]=[0,1,1]
  vertices[0:2,3]=[0,1,0]
  normals[0,0:3]=-1

  vertices[0:2,4]=[1,1,0]
  vertices[0:2,5]=[1,1,1]
  vertices[0:2,6]=[1,0,1]
  vertices[0:2,7]=[1,0,0]
  normals[0,4:7]=1

  vertices[0:2,8]=[0,0,1]
  vertices[0:2,9]=[1,0,1]
  vertices[0:2,10]=[1,1,1]
  vertices[0:2,11]=[0,1,1]
  normals[2,8:11]=1

  vertices[0:2,12]=[0,1,0]
  vertices[0:2,13]=[1,1,0]
  vertices[0:2,14]=[1,0,0]
  vertices[0:2,15]=[0,0,0]
  normals[2,12:15]=-1

  vertices[0:2,16]=[0,1,0]
  vertices[0:2,17]=[0,1,1]
  vertices[0:2,18]=[1,1,1]
  vertices[0:2,19]=[1,1,0]
  normals[1,16:19]=1

  vertices[0:2,20]=[1,0,0]
  vertices[0:2,21]=[1,0,1]
  vertices[0:2,22]=[0,0,1]
  vertices[0:2,23]=[0,0,0]
  normals[1,20:23]=-1

  vertices=vertices-0.5
  if keyword_set(rounded) then begin
    normals=vertices/spread(dot2(transpose(vertices)),0,3)
  endif

  triangles = [5,0,1,2,3,0, 5,4,5,6,7,4, 5,8,9,10,11,8, 5,12,13,14,15,12, $
               5,16,17,18,19,16, 5,20,21,22,23,20 ]

  object=CREATE_STRUCT(['vertices','triangles','normals'],ptr_new(vertices),ptr_new(triangles),ptr_new(normals))
  return,object
end
