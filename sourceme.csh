#
#  this should be sourced each time you log in;
#  You may or may not want to put the lins
#    source $HOME/f90/pencil_modular/sourceme.csh
#  into your .cshrc
#
if (! $?_sourceme) then
  set cdpath = ( . ../  ../../ ../../../ ~/ )
  set path=($path ../bin ../../bin ../../../bin ../../../../bin )
  setenv IDL_PATH "+../idl:+../../idl:+../../../idl:tmp:+~/idl:<IDL_DEFAULT>"
  # if ($?TEXINPUTS) then
  #   setenv TEXINPUTS ${TEXINPUTS}:./texinputs
  # else
  #   setenv TEXINPUTS ::./texinputs
  # endif
  #
  #  additional aliases (for axel)
  #
  alias gb 'cd $gt ; set gg=$gb ; set gb=$gt ; set gt=$gg ; echo $gt "->" $gb'
  alias gt 'set gt=$cwd; cd \!^ ; set gb=$cwd ; echo $gt "->" $gb'
  alias d ls -sCF
  alias .. 'set pwd = $cwd ; cd ..'
  alias cd 'cd \!* ; set prompt=`whoami`"@"`hostname`:"$cwd> "; setenv CWD $cwd'
  alias local 'cp -p \!:1 tmp.$$; \rm \!:1; mv tmp.$$ \!:1; chmod u+w \!:1'
endif
set _sourceme = "set"
