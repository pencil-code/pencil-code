#
#  this should be sourced each time you log in;
#  You may or may not want to put the lins
#    source $HOME/f90/pencil_modular/sourceme.csh
#  into your .cshrc
#
set cdpath = ( . ../  ../../ ../../../ ~/ )
set path=(. ../bin ../../bin ../../../bin ../../../../bin $HOME/bin $path)
setenv IDL_PATH "+../idl:+../../idl:+../../../idl:tmp:+~/idl:<IDL_DEFAULT>"
#
#  additional aliases (for axel)
#
alias gb 'cd $gt ; set gg=$gb ; set gb=$gt ; set gt=$gg ; echo $gt "->" $gb'
alias gt 'set gt=$cwd; cd \!^ ; set gb=$cwd ; echo $gt "->" $gb'
alias d ls -sCF
alias .. 'set pwd = $cwd ; cd ..'
alias cd 'cd \!* ; set prompt=`whoami`"@"`hostname`:"$cwd> "; setenv CWD $cwd'
alias local 'cp -p \!:1 tmp.$$; \rm \!:1; mv tmp.$$ \!:1; chmod u+w \!:1'

