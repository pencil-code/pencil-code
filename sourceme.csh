#$Id$
#
#  This file tries to set the PENCIL_HOME environment variable if it
#  doesn't exist yet, and then adds stuff to your PATH and IDL_PATH.
#  If _sourceme_quiet is set, no output is printed, which enables you to
#  put the lines
#    setenv PENCIL_HOME [...]
#    set _sourceme_quiet; source $PENCIL_HOME/sourceme.csh; unset _sourceme_quiet
#  into you .cshrc file
#

#  set cdpath = ( . ../  ../../ ../../../ ~/ )

if (! $?PENCIL_HOME) then
  unset _sourceme		# tabula rasa without PENCIL_HOME
  #
  # Try to identify position of the code's home directory:
  #
  foreach _dir ( . .. ../.. ../../.. ../../../.. \
                pencil pencil-code \
		f90/pencil f90/pencil-code \
		pencil_modular f90/pencil_modular)
    if ( (-e $_dir/sourceme.csh) && \
         (-d $_dir/bin)          && \
	 (-d $_dir/doc)          && \
	 (-d $_dir/src)          && \
	 (-d $_dir/samples)         \
       ) then
      set back_again = `pwd`     
      cd $_dir; setenv PENCIL_HOME `pwd`; cd $back_again
      goto found
    endif
  end

  echo "sourceme.csh: Cannot locate home directory of pencil code."
  echo "  Try sourcing me from the home directory itself, or set PENCIL_HOME"
  goto eof
endif
    
found:
if (! $?_sourceme_quiet) echo "PENCIL_HOME = <$PENCIL_HOME>"

if (! $?_sourceme) then		# called for the fist time?
  if (-d $PENCIL_HOME/bin) then

    #  Set shell path
    if (! $?_sourceme_quiet) echo "Adding $PENCIL_HOME/{bin,utils{,/axel},scripts} to PATH"
    set path = ( $path $PENCIL_HOME/bin \
                       $PENCIL_HOME/utils \
		       $PENCIL_HOME/utils/axel \
                       $PENCIL_HOME/src/scripts \
		       $PENCIL_HOME/remesh/bin)

		       #  Set path for DX macros
    set _dxpath = "${PENCIL_HOME}/dx/macros:${PENCIL_HOME}/dx/macros/others"
    if ($?DXMACROS) then
      setenv DXMACROS "${_dxpath}:$DXMACROS"
    else
      setenv DXMACROS "${_dxpath}"
    endif
    unset _dxpath

    #  Set IDL path
    if ($?IDL_PATH) then
      setenv IDL_PATH "./idl:../idl:+${PENCIL_HOME}/idl:./data:./tmp:$IDL_PATH"
    else
      setenv IDL_PATH "./idl:../idl:+${PENCIL_HOME}/idl:./data:./tmp:<IDL_DEFAULT>"
    endif

    # Set HDF5 path
    if (! $?HDF5_HOME) setenv HDF5_HOME `which h5fc 2>/dev/null | sed -e 's/\/bin\/h5fc$//'`

    #  Set PYTHON path
    if ($?PYTHONPATH) then
      setenv PYTHONPATH "${PYTHONPATH}:${PENCIL_HOME}/python:${PWD}/python"
    else
      setenv PYTHONPATH "${PENCIL_HOME}/python:${PWD}/python"
    endif
    #  Set library path for linker
    if ($?LD_LIBRARY_PATH) then
      setenv LD_LIBRARY_PATH "${LD_LIBRARY_PATH}:./src:./src/astaroth:./src/astaroth/submodule/build/src/core:./src/astaroth/submodule/build/src/core/kernels:./src/astaroth/submodule/build/src/utils"
    else
      setenv LD_LIBRARY_PATH "./src:./src/astaroth:./src/astaroth/submodule/build/src/core:./src/astaroth/submodule/build/src/core/kernels:./src/astaroth/submodule/build/src/utils"
    endif
#    #  Set Perl module path [no longer needed]
#    set _perl5lib = "${PENCIL_HOME}/lib/perl"
#    if ($?PERL5LIB) then
#      setenv PERL5LIB "${_perl5lib}:$PERL5LIB"
#    else
#      setenv PERL5LIB "${_perl5lib}"
#    endif
#    unset _perl5lib

    # Remember that sourceme has been successfully run
    set _sourceme = 'set'

  else
    echo "Not adding $PENCIL_HOME/bin to PATH: not a directory"
  endif
  #
  #  additional aliases (for Axel)
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
