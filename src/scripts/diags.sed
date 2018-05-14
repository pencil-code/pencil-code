# export DIAGNOSTICS=$MODULE_PREFIX'diagnostics'$MODULE_INFIX
#
#sed -e '/idiag_[0-9a-zA-Z_]*/ !d' -e's/^ *integer.*idiag_\([0-9a-zA-Z_]*\).*$/if (idiag_\1!=0) \n{\n \/\/Get \1 from \n  min_scal_cuda(d_\1, d_partial_result, d_uu_y); \n  cudaDeviceSynchronize();\n  cudaMemcpy(\&\1, (float*)d_\1, sizeof(float), cudaMemcpyDeviceToHost);\n  cudaDeviceSynchronize(); \n}/' < cdata.f90 > diags
#sed -e's/idiag_\([0-9a-zA-Z_]*\)/$DIAGNOSTICSsave_name(\1,idiag_\1);/'
#
/idiag_[0-9a-zA-Z_]*/ !d
h
s/^ *integer.*idiag_\([0-9a-zA-Z_]*\).*$/if (idiag_\1!=0) \n{\n \/\/Get \1 from \n  min_scal_cuda(d_\1, d_partial_result, d_uu_y); \n  cudaDeviceSynchronize();\n  cudaMemcpy(\&\1, (float*)d_\1, sizeof(float), cudaMemcpyDeviceToHost);\n  cudaDeviceSynchronize();/
p
g
s/^ *integer.*idiag_\([0-9a-zA-Z_]*\).*$/  __diagnostics_MOD_save_name(\1,idiag_\1);\n}/
