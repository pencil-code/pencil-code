/^ *call  *farray_register_pde/ !d
/call  *farray_register_pde/ s/^[^!]*farray_register_pde *( *.\([0-9a-zA-Z_]*\). *,[^0-9]*$/FIELD \U\1/
/call  *farray_register_pde/ s/^[^!]*farray_register_pde *( *.\([0-9a-zA-Z_]*\). *,.*vector *= *\([0-9]*\).*$/FIELD \U\1X,\U\1Y,\U\1Z\n#\Ldefine \U\1 F\Lield\2(\U\1X,\U\1Y,\U\1Z)/
#/call  *farray_register_pde/ s/^[^!]*farray_register_pde *( *.\([0-9a-zA-Z_]*\). *,.*array *= *\([0-9]*\).*$/FIELD \1[\2]/

