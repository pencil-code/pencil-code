#
#  This file tries to set the PENCIL_HOME environment variable if it
#  doesn't exist yet, and then adds stuff to your PATH and IDL_PATH.
#
#  You may or may not want to put the lines
#    setenv PENCIL_HOME [...]
#    source $HOME/f90/pencil_modular/sourceme.csh
#  into your .cshrc
#

#  set cdpath = ( . ../  ../../ ../../../ ~/ )

if (! $?PENCIL_HOME) then
  #
  # Try to identify position of the code's home directory:
  #
  foreach _dir ( . .. ../.. ../../.. ../../../.. pencil pencil_modular \
                f90/pencil f90/pencil_modular)
    if ( (-e $_dir/sourceme.csh) && \
         (-d $_dir/bin)          && \
	 (-d $_dir/doc)          && \
	 (-d $_dir/src)          && \
         (-d $_dir/runs)            \
       ) then
      setenv PENCIL_HOME `cd $_dir; pwd`
      goto found
    endif
  end

  echo "Cannot locate home directory of pencil code."
  echo "Try sourcing me from the home directory itself, or set PENCIL_HOME"
  goto eof
endif
    
found:
echo "PENCIL_HOME = <$PENCIL_HOME>"

if (! $?_sourceme) then		# called for the fist time?
  if (-d $PENCIL_HOME/bin) then
    #  Set shell path
    echo "Adding $PENCIL_HOME/bin to path"
    set path = ( $path $PENCIL_HOME/bin )
    #  Set IDL path
    if ($?IDL_PATH) then
      setenv IDL_PATH "./idl:../idl:+${PENCIL_HOME}/idl:./tmp:$IDL_PATH"
    else
      setenv IDL_PATH "./idl:../idl:+${PENCIL_HOME}/idl:./tmp:<IDL_DEFAULT>"
    endif
    set _sourceme = 'set'
  else
    echo "Not adding $PENCIL_HOME/bin to path: not a directory"
  endif
  #
  #  additional aliases (for axel)
  #
  alias gb 'cd $gt ; set gg=$gb ; set gb=$gt ; set gt=$gg ; echo $gt "->" $gb'
  alias gt 'set gt=$cwd; cd \!^ ; set gb=$cwd ; echo $gt "->" $gb'
  # alias d ls -sCF
  alias .. 'set pwd = $cwd ; cd ..'
  # alias local 'cp -p \!:1 tmp.$$; \rm \!:1; mv tmp.$$ \!:1; chmod u+w \!:1'
endif

#
#  Clean up and exit
#
eof:

unset _dir
