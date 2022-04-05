sed -i -e '/module  *procedure/ ! s/^\( *module  *\)[a-zA-Z0-9_]* */\1Special/' \
              -e 's/^\( *endmodule  *\)[a-zA-Z0-9_]* *$/\1Special/' \
              -e "s/[a-zA-Z0-9_]*\.inc/special_dummies.inc/" \
              -e "s/namelist *\/[a-zA-Z0-9_]*[a-zA-Z0-9]_\([ir][a-z]*_pars\)/namelist \/special_\1/" \
              -e "s/NML *= *[a-zA-Z0-9_]*_\([a-z][a-z]*_pars\)/NML=special_\1/" $1
