function vizit_glyph_sphere,refinement=refinement
  default,refinement,0

  sphere=vizit_glyph_icosahedron()

  for i=1,refinement do begin
    refine_surface,sphere
  endfor
  
  *(sphere.vertices)=*(sphere.vertices)/spread(sqrt(dot2(transpose(*(sphere.vertices)))),0,3)
  *(sphere.normals)=*(sphere.vertices)

  return,sphere
end
