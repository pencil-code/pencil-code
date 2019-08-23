   function tag_exists, struct, tag
;
;  Tests whether a structure has a given tag name.
;
     return, ((where(strupcase (tag_names(struct)) eq strupcase (tag)))[0] ge 0)

   end
