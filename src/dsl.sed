# sed script for generation of field declarations for DSL code
# so far only PDE variables are considered
# multiple register calls for the same field result in multiple declarations -> 
# run-time switches like ldensity_nolog need to be taken into account in addition

/^ *call  *farray_register_pde/ !d
/call  *farray_register_pde/ s/^[^!]*farray_register_pde *( *.\([0-9a-zA-Z_]*\). *,[^0-9]*$/FIELD \U\1/
/call  *farray_register_pde/ s/^[^!]*farray_register_pde *( *.\([0-9a-zA-Z_]*\). *,.*vector *= *\([0-9]*\).*$/FIELD \U\1X,\U\1Y,\U\1Z\n#\Ldefine \U\1 F\Lield\2(\U\1X,\U\1Y,\U\1Z)/
#/call  *farray_register_pde/ s/^[^!]*farray_register_pde *( *.\([0-9a-zA-Z_]*\). *,.*array *= *\([0-9]*\).*$/FIELD \1[\2]/

