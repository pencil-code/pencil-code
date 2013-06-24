<?
	include "inc/header.inc";
 ?>
<div class="centcolumnpad">
<h2>FAQ</h2>


<h3 class="sectionHead"><span class="titlemark">1    </span> <a 
 id="x1-10001"></a>Troubleshooting / Frequently Asked Questions</h3>
<a 
 id="x1-1001"></a>
<a 
 id="x1-1002"></a>
<!--l. 390--><p class="noindent" >
<h4 class="subsectionHead"><span class="titlemark">1.1    </span> <a 
 id="x1-20001.1"></a>Download and setup</h4>
<!--l. 392--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.1.1    </span> <a 
 id="x1-30001.1.1"></a>Download forbidden</h5>
<a 
 id="x1-3001"></a>
<a 
 id="x1-3002"></a>
<!--l. 396--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Both Google Code and SourceForge are banned from countries on the United
States Office of Foreign Assets Control sanction list, including Cuba, Iran, Libya,
North Korea, Sudan and Syria; see <a 
href="http://en.wikipedia.org/wiki/Google_Code" class="url" ><span 
class="cmtt-12">http://en.wikipedia.org/wiki/Google_Code</span></a> and
<a 
href="http://en.wikipedia.org/wiki/SourceForge" class="url" ><span 
class="cmtt-12">http://en.wikipedia.org/wiki/SourceForge</span></a>. As a remedy, you might download a tarball from
<a 
href="http://pencil-code.nordita.org/" class="url" ><span 
class="cmtt-12">http://pencil-code.nordita.org/</span></a>; see also Section&#x00A0;<span 
class="pncb7t-x-x-120">??</span>.
<!--l. 410--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.1.2    </span> <a 
 id="x1-40001.1.2"></a>When sourcing the &#8216;<span 
class="cmtt-12">sourceme.sh</span>&#8217;/&#8216;<span 
class="cmtt-12">sourceme.csh</span>&#8217; file or running <span 
class="cmtt-12">pc_setupsrc</span>, I get error
messages from the shell, like &#8216;if: Expression Syntax.&#8217; or &#8216;set: Variable name must begin with a
letter.&#8217;</h5>
<a 
 id="x1-4001"></a>
<!--l. 415--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: This sounds like a buggy shell setup, either by yourself or your system administrator
&#8212; or a shell that is even more idiosyncratic than the ones we have been working
with.
<!--l. 419--><p class="noindent" >To better diagnose the problem, collect the following information before filing a bug report to
us:
<!--l. 422--><p class="noindent" >
      <ol  class="enumerate1" >
      <li 
  class="enumerate" id="x1-4003x1"><span 
class="cmtt-12">uname -a</span>
      </li>
      <li 
  class="enumerate" id="x1-4005x2"><span 
class="cmtt-12">/bin/csh -v</span>
                                                                                            
                                                                                            
      </li>
      <li 
  class="enumerate" id="x1-4007x3"><span 
class="cmtt-12">echo $version</span>
      </li>
      <li 
  class="enumerate" id="x1-4009x4"><span 
class="cmtt-12">echo $SHELL</span>
      </li>
      <li 
  class="enumerate" id="x1-4011x5"><span 
class="cmtt-12">ps -p $$</span>
      </li>
      <li 
  class="enumerate" id="x1-4013x6">If you have problems while sourcing the &#8216;<span 
class="cmtt-12">sourceme</span>&#8217;<a 
 id="x1-4014"></a> script,
           <ol  class="enumerate2" >
           <li 
  class="enumerate" id="x1-4016x1">unset the <span 
class="cmtt-12">PENCIL_HOME</span><a 
 id="x1-4017"></a> variable:
               <dl class="description"><dt class="description">
           <span 
class="pncb7t-x-x-120">for </span><span 
class="pncbo7t-x-x-120">csh</span><a 
 id="x1-4018"></a> <span 
class="pncb7t-x-x-120">and similar:</span>  </dt><dd 
class="description"><span 
class="cmtt-12">unsetenv PENCIL_HOME</span>
               </dd><dt class="description">
           <span 
class="pncb7t-x-x-120">for </span><span 
class="pncbo7t-x-x-120">bash</span><a 
 id="x1-4019"></a> <span 
class="pncb7t-x-x-120">and similar:</span>  </dt><dd 
class="description"><span 
class="cmtt-12">unexport PENCIL_HOME; unset PENCIL_HOME</span></dd></dl>
           </li>
           <li 
  class="enumerate" id="x1-4021x2">switch your shell in verbose mode,
               <dl class="description"><dt class="description">
           <span 
class="pncb7t-x-x-120">for </span><span 
class="pncbo7t-x-x-120">csh</span><a 
 id="x1-4022"></a> <span 
class="pncb7t-x-x-120">and similar:</span>  </dt><dd 
class="description"><span 
class="cmtt-12">set verbose; set echo</span>
               </dd><dt class="description">
           <span 
class="pncb7t-x-x-120">for </span><span 
class="pncbo7t-x-x-120">bash</span><a 
 id="x1-4023"></a> <span 
class="pncb7t-x-x-120">and similar:</span>  </dt><dd 
class="description"><span 
class="cmtt-12">set -v; set -x</span>
               </dd><dt class="description">
            </dt><dd 
class="description">then source again.</dd></dl>
           </li></ol>
      </li>
      <li 
  class="enumerate" id="x1-4025x7">If you have problems with <span 
class="cmtt-12">pc_setupsrc</span><a 
 id="x1-4026"></a>, run it with <span 
class="cmtt-12">csh</span><a 
 id="x1-4027"></a> in verbose mode: <div class="alltt">
      <!--l. 455--><p class="noindent" ><div class="obeylines-v">
      <span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/bin/csh</span><span 
class="cmtt-12">&#x00A0;-v</span><span 
class="cmtt-12">&#x00A0;-x</span><span 
class="cmtt-12">&#x00A0;$PENCIL_HOME/bin/pc_setupsrc</span>
      <br /><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span></div>
      </div></li></ol>
<!--l. 462--><p class="noindent" >
<h4 class="subsectionHead"><span class="titlemark">1.2    </span> <a 
 id="x1-50001.2"></a>Compilation</h4>
<!--l. 465--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.1    </span> <a 
 id="x1-60001.2.1"></a>Linker can&#8217;t find the syscalls functions:</h5>
                                                                                            
                                                                                            
<!--l. 467--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb1"><a 
 id="x1-6002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ld:</span><span 
class="cmtt-12">&#x00A0;0711-317</span><span 
class="cmtt-12">&#x00A0;ERROR:</span><span 
class="cmtt-12">&#x00A0;Undefined</span><span 
class="cmtt-12">&#x00A0;symbol:</span><span 
class="cmtt-12">&#x00A0;.is_nan_c</span>
<br class="fancyvrb" /><a 
 id="x1-6004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ld:</span><span 
class="cmtt-12">&#x00A0;0711-317</span><span 
class="cmtt-12">&#x00A0;ERROR:</span><span 
class="cmtt-12">&#x00A0;Undefined</span><span 
class="cmtt-12">&#x00A0;symbol:</span><span 
class="cmtt-12">&#x00A0;.sizeof_real_c</span>
<br class="fancyvrb" /><a 
 id="x1-6006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ld:</span><span 
class="cmtt-12">&#x00A0;0711-317</span><span 
class="cmtt-12">&#x00A0;ERROR:</span><span 
class="cmtt-12">&#x00A0;Undefined</span><span 
class="cmtt-12">&#x00A0;symbol:</span><span 
class="cmtt-12">&#x00A0;.system_c</span>
<br class="fancyvrb" /><a 
 id="x1-6008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ld:</span><span 
class="cmtt-12">&#x00A0;0711-317</span><span 
class="cmtt-12">&#x00A0;ERROR:</span><span 
class="cmtt-12">&#x00A0;Undefined</span><span 
class="cmtt-12">&#x00A0;symbol:</span><span 
class="cmtt-12">&#x00A0;.get_env_var_c</span>
<br class="fancyvrb" /><a 
 id="x1-6010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ld:</span><span 
class="cmtt-12">&#x00A0;0711-317</span><span 
class="cmtt-12">&#x00A0;ERROR:</span><span 
class="cmtt-12">&#x00A0;Undefined</span><span 
class="cmtt-12">&#x00A0;symbol:</span><span 
class="cmtt-12">&#x00A0;.get_pid_c</span>
<br class="fancyvrb" /><a 
 id="x1-6012r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ld:</span><span 
class="cmtt-12">&#x00A0;0711-317</span><span 
class="cmtt-12">&#x00A0;ERROR:</span><span 
class="cmtt-12">&#x00A0;Undefined</span><span 
class="cmtt-12">&#x00A0;symbol:</span><span 
class="cmtt-12">&#x00A0;.file_size_c</span></div>
<!--l. 478--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The Pencil Code needs a working combination of a Fortran- and a C-compiler. If this is not
correctly set up, usually the linker won&#8217;t find the functions inside the <span 
class="pncro7t-x-x-120">syscalls</span><a 
 id="x1-6013"></a> module. If that
happens, either the combination of C- and Fortran-compiler is inappropriate (e.g. <span 
class="pncro7t-x-x-120">ifort</span><a 
 id="x1-6014"></a>
needs <span 
class="pncro7t-x-x-120">icc</span><a 
 id="x1-6015"></a>), or the compiler needs additional flags, like <span 
class="pncro7t-x-x-120">g95</span><a 
 id="x1-6016"></a> might need the option
&#8216;<span 
class="cmtt-12">-fno-second-underscore</span>&#8217;<a 
 id="x1-6017"></a> and <span 
class="pncro7t-x-x-120">xlf</span><a 
 id="x1-6018"></a> might need the option &#8216;<span 
class="cmtt-12">-qextname</span>&#8217;<a 
 id="x1-6019"></a>. Please refer to Sect.&#x00A0;<span 
class="pncb7t-x-x-120">??</span>,
Table&#x00A0;<span 
class="pncb7t-x-x-120">??</span>.
<!--l. 488--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.2    </span> <a 
 id="x1-70001.2.2"></a>Make gives the following error now:</h5>
<!--l. 489--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb2"><a 
 id="x1-7002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;PGF90-S-0017-Unable</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;open</span><span 
class="cmtt-12">&#x00A0;include</span><span 
class="cmtt-12">&#x00A0;file:</span><span 
class="cmtt-12">&#x00A0;chemistry.h</span><span 
class="cmtt-12">&#x00A0;(nochemistry.f90:</span><span 
class="cmtt-12">&#x00A0;43)</span>
<br class="fancyvrb" /><a 
 id="x1-7004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;0</span><span 
class="cmtt-12">&#x00A0;inform,</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;0</span><span 
class="cmtt-12">&#x00A0;warnings,</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;1</span><span 
class="cmtt-12">&#x00A0;severes,</span><span 
class="cmtt-12">&#x00A0;0</span><span 
class="cmtt-12">&#x00A0;fatal</span><span 
class="cmtt-12">&#x00A0;for</span><span 
class="cmtt-12">&#x00A0;chemistry</span><br class="fancyvrb" /><a 
 id="x1-7006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span>
<br class="fancyvrb" /><a 
 id="x1-7008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Line</span><span 
class="cmtt-12">&#x00A0;43</span><span 
class="cmtt-12">&#x00A0;of</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;nochemistry</span><span 
class="cmtt-12">&#x00A0;routine,</span><span 
class="cmtt-12">&#x00A0;only</span><span 
class="cmtt-12">&#x00A0;has</span><span 
class="cmtt-12">&#x00A0;&#8217;contains&#8217;.</span></div>
<!--l. 498--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: This is because somebody added a new module (together with a corresponding
<span 
class="pncro7t-x-x-120">nomodule.f90</span><a 
 id="x1-7009"></a> and a <span 
class="pncro7t-x-x-120">module.h</span><a 
 id="x1-7010"></a> file (chemistry in this case). These files didn&#8217;t exist before, so
you need to say:
<div class="fancyvrb" id="fancyvrb3"><a 
 id="x1-7012r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;pc_setupsrc</span></div>
<!--l. 504--><p class="noindent" >If this does not help, say first <span 
class="cmtt-12">make clean</span><a 
 id="x1-7013"></a> and then <span 
class="cmtt-12">pc_setupsrc</span><a 
 id="x1-7014"></a>.
<!--l. 508--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.3    </span> <a 
 id="x1-80001.2.3"></a>How do I compile the <span 
class="pncrc7t-x-x-120">P<span 
class="small-caps">e</span><span 
class="small-caps">n</span><span 
class="small-caps">c</span><span 
class="small-caps">i</span><span 
class="small-caps">l</span> C<span 
class="small-caps">o</span><span 
class="small-caps">d</span><span 
class="small-caps">e</span> </span>with the Intel (<span 
class="pncro7t-x-x-120">ifc</span>) compiler under <span 
class="pncro7t-x-x-120">Linux</span>?</h5>
<!--l. 512--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The <span 
class="pncrc7t-x-x-120">P<span 
class="small-caps">e</span><span 
class="small-caps">n</span><span 
class="small-caps">c</span><span 
class="small-caps">i</span><span 
class="small-caps">l</span> C<span 
class="small-caps">o</span><span 
class="small-caps">d</span><span 
class="small-caps">e</span> </span>should compile successfully with <span 
class="pncro7t-x-x-120">ifc</span><a 
 id="x1-8001"></a>&#x00A0;6.x, <span 
class="pncro7t-x-x-120">ifc</span><a 
 id="x1-8002"></a>&#x00A0;7.0, sufficiently recent
versions of <span 
class="pncro7t-x-x-120">ifc</span><a 
 id="x1-8003"></a> 7.1 (you should get the latest version; if yours is too old, you will typically get an
&#8216;internal compiler error&#8217; during compilation of &#8216;<span 
class="cmtt-12">src/hydro.f90</span>&#8217;<a 
 id="x1-8004"></a>), as well as with recent versions
of <span 
class="pncro7t-x-x-120">ifort</span><a 
 id="x1-8005"></a> 8.1 (8.0 may also work).
<!--l. 519--><p class="noindent" >You can find the <span 
class="pncro7t-x-x-120">ifort</span><a 
 id="x1-8006"></a> compiler at <a 
href="ftp://download.intel.com/software/products/compilers/downloads" class="url" ><span 
class="cmtt-12">ftp://download.intel.com/software/products/compilers/downloads</span></a>.
<!--l. 522--><p class="noindent" >On many current (as of November 2003) Linux systems, there is a mismatch between the <span 
class="pncro7t-x-x-120">glibc</span><a 
 id="x1-8007"></a>
versions used by the compiler and the linker. To work around this, use the following flag for
compiling
                                                                                            
                                                                                            
<div class="fancyvrb" id="fancyvrb4"><a 
 id="x1-8009r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;FC=ifc</span><span 
class="cmtt-12">&#x00A0;-i_dynamic</span></div>
<!--l. 528--><p class="noindent" >and set the environment variable
<div class="fancyvrb" id="fancyvrb5"><a 
 id="x1-8011r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;LD_ASSUME_KERNEL=2.4.1;</span><span 
class="cmtt-12">&#x00A0;export</span><span 
class="cmtt-12">&#x00A0;LD_ASSUME_KERNEL</span></div>
<!--l. 532--><p class="noindent" >or
<div class="fancyvrb" id="fancyvrb6"><a 
 id="x1-8013r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;setenv</span><span 
class="cmtt-12">&#x00A0;LD_ASSUME_KERNEL</span><span 
class="cmtt-12">&#x00A0;2.4.1</span></div>
<!--l. 536--><p class="noindent" >This has solved the problems e.g.&#x00A0;on a system with <span 
class="pncri7t-x-x-120">glibc-2.3.2 </span>and kernel <span 
class="pncri7t-x-x-120">2.4.22</span>.
<!--l. 539--><p class="noindent" >Thanks to Leonardo J. Milano (<a 
href="http://udel.edu/~lmilano/" class="url" ><span 
class="cmtt-12">http://udel.edu/</span><span 
class="cmtt-12">~</span><span 
class="cmtt-12">lmilano/</span></a>) for part of this info.
<!--l. 544--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.4    </span> <a 
 id="x1-90001.2.4"></a>I keep getting segmentation faults with &#8216;<span 
class="cmtt-12">start.x</span>&#8217; when compiling with <span 
class="pncro7t-x-x-120">ifort</span>
8.0</h5>
<!--l. 548--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: There was/is a number of issues with <span 
class="pncro7t-x-x-120">ifort</span><a 
 id="x1-9001"></a> 8.0. Make sure you have the latest patches
applied to the compiler. A number of things to consider or try are:
      <ol  class="enumerate1" >
      <li 
  class="enumerate" id="x1-9003x1">Compile with the the &#8216;<span 
class="cmtt-12">-static -nothreads</span>&#8217;<a 
 id="x1-9004"></a> flags.
      </li>
      <li 
  class="enumerate" id="x1-9006x2">Set your stacksize to a large value (but a far too large value may be problematic, too),
      e. g.
      <div class="fancyvrb" id="fancyvrb7"><a 
 id="x1-9008r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;limit</span><span 
class="cmtt-12">&#x00A0;stacksize</span><span 
class="cmtt-12">&#x00A0;256m</span><br class="fancyvrb" /><a 
 id="x1-9010r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ulimit</span><span 
class="cmtt-12">&#x00A0;-s</span><span 
class="cmtt-12">&#x00A0;256000</span></div>
      </li>
      <li 
  class="enumerate" id="x1-9012x3">Set the environment variable KMP&#x02D9;STACKSIZE to a large value (like <span 
class="cmtt-12">100M</span>)</li></ol>
<!--l. 563--><p class="noindent" >See also <a 
href="http://softwareforums.intel.com/ids/board/message?board.id=11&message.id=1375" class="url" ><span 
class="cmtt-12">http://softwareforums.intel.com/ids/board/message?board.id=11&amp;message.id=1375</span></a>
<!--l. 567--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.5    </span> <a 
 id="x1-100001.2.5"></a>When compiling with MPI on a Linux system, the linker complains:</h5>
<!--l. 569--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb8"><a 
 id="x1-10002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;mpicomm.o:</span><span 
class="cmtt-12">&#x00A0;In</span><span 
class="cmtt-12">&#x00A0;function</span><span 
class="cmtt-12">&#x00A0;&#8216;mpicomm_mpicomm_init_&#8217;:</span>
<br class="fancyvrb" /><a 
 id="x1-10004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;mpicomm.o(.text+0x36):</span><span 
class="cmtt-12">&#x00A0;undefined</span><span 
class="cmtt-12">&#x00A0;reference</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;&#8216;mpi_init_&#8217;</span>
<br class="fancyvrb" /><a 
 id="x1-10006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;mpicomm.o(.text+0x55):</span><span 
class="cmtt-12">&#x00A0;undefined</span><span 
class="cmtt-12">&#x00A0;reference</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;&#8216;mpi_comm_size_&#8217;</span>
<br class="fancyvrb" /><a 
 id="x1-10008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;mpicomm.o(.text+0x6f):</span><span 
class="cmtt-12">&#x00A0;undefined</span><span 
class="cmtt-12">&#x00A0;reference</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;&#8216;mpi_comm_rank_&#8217;</span><br class="fancyvrb" /><a 
 id="x1-10010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;[...]</span></div>
<!--l. 579--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: This is the infamous <span 
class="pncro7t-x-x-120">underscore problem</span><a 
 id="x1-10011"></a>. Your <span 
class="pncro7t-x-x-120">MPI</span><a 
 id="x1-10012"></a> libraries have been compiled with <span 
class="pncro7t-x-x-120">G77</span><a 
 id="x1-10013"></a>
without the option &#8216;<span 
class="cmtt-12">-fno-second-underscore</span>&#8217;<a 
 id="x1-10014"></a>, which makes the <span 
class="pncro7t-x-x-120">MPI</span><a 
 id="x1-10015"></a> symbol names
incompatible with other Fortran compilers.
<!--l. 584--><p class="noindent" >As a workaround, use
                                                                                            
                                                                                            
<div class="fancyvrb" id="fancyvrb9"><a 
 id="x1-10017r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;MPICOMM</span><span 
class="cmtt-12">&#x00A0;=</span><span 
class="cmtt-12">&#x00A0;mpicomm_</span></div>
<!--l. 588--><p class="noindent" >in &#8216;<span 
class="cmtt-12">Makefile.local</span>&#8217;<a 
 id="x1-10018"></a>. Or, even better, you can set this globally (for the given computer) by
inserting that line into the file &#8216;<span 
class="cmtt-12">~/.adapt-mkfile.inc</span>&#8217;<a 
 id="x1-10019"></a> (see <span 
class="cmtt-12">perldoc adapt-mkfile </span>for more
details). <a 
 id="x1-10020"></a>
<!--l. 595--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.6    </span> <a 
 id="x1-110001.2.6"></a>Compilation stops with the cryptic error message:</h5>
<!--l. 596--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb10"><a 
 id="x1-11002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;f95</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-O3</span><span 
class="cmtt-12">&#x00A0;-u</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;.f90.f90</span><br class="fancyvrb" /><a 
 id="x1-11004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;:</span><span 
class="cmtt-12">&#x00A0;Could</span><span 
class="cmtt-12">&#x00A0;not</span><span 
class="cmtt-12">&#x00A0;open</span><span 
class="cmtt-12">&#x00A0;sourcefile</span><span 
class="cmtt-12">&#x00A0;.f90.f90</span>
<br class="fancyvrb" /><a 
 id="x1-11006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;compilation</span><span 
class="cmtt-12">&#x00A0;aborted</span><span 
class="cmtt-12">&#x00A0;for</span><span 
class="cmtt-12">&#x00A0;.f90.f90</span><span 
class="cmtt-12">&#x00A0;(code</span><span 
class="cmtt-12">&#x00A0;1)</span><br class="fancyvrb" /><a 
 id="x1-11008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[1]:</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;[.f90.o]</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;1</span></div>
<!--l. 602--><p class="noindent" >What is the problem?
<!--l. 606--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: There are two possibilities:
      <ol  class="enumerate1" >
      <li 
  class="enumerate" id="x1-11010x1">One of the variables for <span 
class="pncro7t-x-x-120">make</span><a 
 id="x1-11011"></a> has not been set, so <span 
class="pncro7t-x-x-120">make</span><a 
 id="x1-11012"></a> expands it to the empty
      string. Most probably you forgot to specify a module in &#8216;<span 
class="cmtt-12">src/Makefile.local</span>&#8217;<a 
 id="x1-11013"></a>. One
      possibility is that you have upgraded from an older version of the code that did not
      have some of the modules the new version has.
      <!--l. 616--><p class="noindent" >Compare your &#8216;<span 
class="cmtt-12">src/Makefile.local</span>&#8217;<a 
 id="x1-11014"></a> to one of the examples that work.
      </li>
      <li 
  class="enumerate" id="x1-11016x2">One of the variables for <span 
class="pncro7t-x-x-120">make</span><a 
 id="x1-11017"></a> has a space appended to it, e. g.&#x00A0;if you use the
      line
           <div class="quote">
           <!--l. 623--><p class="noindent" ><span 
class="cmtt-12">MPICOMM = mpicomm_</span>_</div>
      <!--l. 625--><p class="noindent" >(see <span 
class="cmsy-10x-x-120">§</span>&#x00A0;<a 
href="#x1-100001.2.5">1.2.5<!--tex4ht:ref: Sec-underscore-problem --></a>) with a trailing blank, you will encounter this error message. Remove the
      blank. This problem can also occur if you added a new module (and have an empty
      space after the module name in &#8216;<span 
class="cmtt-12">src/Makefile.src</span>&#8217;<a 
 id="x1-11018"></a>, i.e.&#x00A0;<span 
class="cmtt-12">CHIRAL=nochiral</span>_<a 
 id="x1-11019"></a>),
      in which case the compiler will talk about &#8220;circular dependence&#8221; for the file
      &#8216;<span 
class="cmtt-12">nochiral</span>&#8217;<a 
 id="x1-11020"></a>.
      </li></ol>
<!--l. 637--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.7    </span> <a 
 id="x1-120001.2.7"></a>The code doesn&#8217;t compile,</h5>
<!--l. 638--><p class="noindent" >&hellip;there is a problem with <span 
class="pncro7t-x-x-120">mvar</span><a 
 id="x1-12001"></a>:
<div class="fancyvrb" id="fancyvrb11"><a 
 id="x1-12003r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make</span><span 
class="cmtt-12">&#x00A0;start.x</span><span 
class="cmtt-12">&#x00A0;run.x</span><br class="fancyvrb" /><a 
 id="x1-12005r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;f95</span><span 
class="cmtt-12">&#x00A0;-O4</span><span 
class="cmtt-12">&#x00A0;-u</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;cdata.f90</span><br class="fancyvrb" /><a 
 id="x1-12007r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Error:</span><span 
class="cmtt-12">&#x00A0;cdata.f90,</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;71:</span><span 
class="cmtt-12">&#x00A0;Implicit</span><span 
class="cmtt-12">&#x00A0;type</span><span 
class="cmtt-12">&#x00A0;for</span><span 
class="cmtt-12">&#x00A0;MVAR</span>
<br class="fancyvrb" /><a 
 id="x1-12009r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;detected</span><span 
class="cmtt-12">&#x00A0;at</span><span 
class="cmtt-12">&#x00A0;MVAR@)</span><br class="fancyvrb" /><a 
 id="x1-12011r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;[f95</span><span 
class="cmtt-12">&#x00A0;terminated</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;errors</span><span 
class="cmtt-12">&#x00A0;found</span><span 
class="cmtt-12">&#x00A0;by</span><span 
class="cmtt-12">&#x00A0;pass</span><span 
class="cmtt-12">&#x00A0;1]</span><br class="fancyvrb" /><a 
 id="x1-12013r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[1]:</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;[cdata.o]</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;2</span></div>
                                                                                            
                                                                                            
<!--l. 650--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Check and make sure that &#8216;<span 
class="cmtt-12">mkcparam</span>&#8217;<a 
 id="x1-12014"></a> (directory &#8216;<span 
class="cmtt-12">$PENCIL_HOME/bin</span>&#8217;<a 
 id="x1-12015"></a>) is in your path. If this
doesn&#8217;t help, there may be an <span 
class="pncri7t-x-x-120">empty </span>&#8216;<span 
class="cmtt-12">cparam.inc</span>&#8217;<a 
 id="x1-12016"></a> file in your &#8216;<span 
class="cmtt-12">src</span>&#8217;<a 
 id="x1-12017"></a> directory. Remove
&#8216;<span 
class="cmtt-12">cparam.inc</span>&#8217;<a 
 id="x1-12018"></a> and try again (Note that &#8216;<span 
class="cmtt-12">cparam.inc</span>&#8217;<a 
 id="x1-12019"></a> is automatically generated from the
&#8216;<span 
class="cmtt-12">Makefile</span>&#8217;<a 
 id="x1-12020"></a>).
<!--l. 658--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.8    </span> <a 
 id="x1-130001.2.8"></a>Some samples don&#8217;t even compile,</h5>
<!--l. 659--><p class="noindent" >as you can see on the web, <a 
href="http://www.nordita.org/software/pencil-code/tests.html" class="url" ><span 
class="cmtt-12">http://www.nordita.org/software/pencil-code/tests.html</span></a>.
<div class="fancyvrb" id="fancyvrb12"><a 
 id="x1-13002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;samples/helical-MHDturb:</span><br class="fancyvrb" /><a 
 id="x1-13004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Compiling..</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;not</span><span 
class="cmtt-12">&#x00A0;ok:</span><br class="fancyvrb" /><a 
 id="x1-13006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make</span><span 
class="cmtt-12">&#x00A0;start.x</span><span 
class="cmtt-12">&#x00A0;run.x</span><span 
class="cmtt-12">&#x00A0;read_videofiles.x</span>
<br class="fancyvrb" /><a 
 id="x1-13008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[1]:</span><span 
class="cmtt-12">&#x00A0;Entering</span><span 
class="cmtt-12">&#x00A0;directory</span><span 
class="cmtt-12">&#x00A0;&#8216;/home/dobler/f90/pencil-code/samples/helical-MHDturb/src&#8217;</span>
<br class="fancyvrb" /><a 
 id="x1-13010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/usr/lib/lam/bin/mpif95</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-O3</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;initcond.f90</span><br class="fancyvrb" /><a 
 id="x1-13012r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/usr/lib/lam/bin/mpif95</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-O3</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;density.f90</span>
<br class="fancyvrb" /><a 
 id="x1-13014r7"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;use</span><span 
class="cmtt-12">&#x00A0;Gravity,</span><span 
class="cmtt-12">&#x00A0;only:</span><span 
class="cmtt-12">&#x00A0;gravz,</span><span 
class="cmtt-12">&#x00A0;nu_epicycle</span><br class="fancyvrb" /><a 
 id="x1-13016r8"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;^</span>
<br class="fancyvrb" /><a 
 id="x1-13018r9"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;208</span><span 
class="cmtt-12">&#x00A0;at</span><span 
class="cmtt-12">&#x00A0;(467:density.f90)</span><span 
class="cmtt-12">&#x00A0;:</span><span 
class="cmtt-12">&#x00A0;No</span><span 
class="cmtt-12">&#x00A0;such</span><span 
class="cmtt-12">&#x00A0;entity</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;module</span>
<br class="fancyvrb" /><a 
 id="x1-13020r10"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;355</span><span 
class="cmtt-12">&#x00A0;:</span><span 
class="cmtt-12">&#x00A0;In</span><span 
class="cmtt-12">&#x00A0;procedure</span><span 
class="cmtt-12">&#x00A0;INIT_LNRHO</span><span 
class="cmtt-12">&#x00A0;variable</span><span 
class="cmtt-12">&#x00A0;NU_EPICYCLE</span><span 
class="cmtt-12">&#x00A0;has</span><span 
class="cmtt-12">&#x00A0;not</span><span 
class="cmtt-12">&#x00A0;been</span><span 
class="cmtt-12">&#x00A0;given</span><span 
class="cmtt-12">&#x00A0;a</span><span 
class="cmtt-12">&#x00A0;type</span>
<br class="fancyvrb" /><a 
 id="x1-13022r11"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;355</span><span 
class="cmtt-12">&#x00A0;:</span><span 
class="cmtt-12">&#x00A0;In</span><span 
class="cmtt-12">&#x00A0;procedure</span><span 
class="cmtt-12">&#x00A0;POLYTROPIC_LNRHO_DISC</span><span 
class="cmtt-12">&#x00A0;variable</span><span 
class="cmtt-12">&#x00A0;NU_EPICYCLE</span><span 
class="cmtt-12">&#x00A0;has</span><span 
class="cmtt-12">&#x00A0;not</span><span 
class="cmtt-12">&#x00A0;been</span><span 
class="cmtt-12">&#x00A0;given</span><span 
class="cmtt-12">&#x00A0;a</span><span 
class="cmtt-12">&#x00A0;type</span>
<br class="fancyvrb" /><a 
 id="x1-13024r12"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;3</span><span 
class="cmtt-12">&#x00A0;Errors</span><br class="fancyvrb" /><a 
 id="x1-13026r13"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;compilation</span><span 
class="cmtt-12">&#x00A0;aborted</span><span 
class="cmtt-12">&#x00A0;for</span><span 
class="cmtt-12">&#x00A0;density.f90</span><span 
class="cmtt-12">&#x00A0;(code</span><span 
class="cmtt-12">&#x00A0;1)</span><br class="fancyvrb" /><a 
 id="x1-13028r14"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[1]:</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;[density.o]</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;1</span>
<br class="fancyvrb" /><a 
 id="x1-13030r15"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[1]:</span><span 
class="cmtt-12">&#x00A0;Leaving</span><span 
class="cmtt-12">&#x00A0;directory</span><span 
class="cmtt-12">&#x00A0;&#8216;/home/dobler/f90/pencil-code/samples/helical-MHDturb/src&#8217;</span>
<br class="fancyvrb" /><a 
 id="x1-13032r16"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make:</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;[code]</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;2</span></div>
<!--l. 682--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Somebody may have checked in something without having run auto-test beforehand. The
problem here is that something has been added in one module, but not in the corresponding
no-module. You can of course check with <span 
class="pncro7t-x-x-120">svn</span><a 
 id="x1-13033"></a> who it was&hellip;
<!--l. 687--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.9    </span> <a 
 id="x1-140001.2.9"></a>Internal compiler error with Compaq/Dec F90</h5>
<!--l. 688--><p class="noindent" >The Dec Fortran optimizer has occasional problems with &#8216;<span 
class="cmtt-12">nompicomm.f90</span>&#8217;<a 
 id="x1-14001"></a>:
<div class="fancyvrb" id="fancyvrb13"><a 
 id="x1-14003r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make</span><span 
class="cmtt-12">&#x00A0;start.x</span><span 
class="cmtt-12">&#x00A0;run.x</span><span 
class="cmtt-12">&#x00A0;read_videofiles.x</span><br class="fancyvrb" /><a 
 id="x1-14005r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;f90</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-fast</span><span 
class="cmtt-12">&#x00A0;-O5</span><span 
class="cmtt-12">&#x00A0;-tune</span><span 
class="cmtt-12">&#x00A0;ev6</span><span 
class="cmtt-12">&#x00A0;-arch</span><span 
class="cmtt-12">&#x00A0;ev6</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;cparam.f90</span>
<br class="fancyvrb" /><a 
 id="x1-14007r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;[...]</span><br class="fancyvrb" /><a 
 id="x1-14009r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;f90</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-fast</span><span 
class="cmtt-12">&#x00A0;-O5</span><span 
class="cmtt-12">&#x00A0;-tune</span><span 
class="cmtt-12">&#x00A0;ev6</span><span 
class="cmtt-12">&#x00A0;-arch</span><span 
class="cmtt-12">&#x00A0;ev6</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;nompicomm.f90</span>
<br class="fancyvrb" /><a 
 id="x1-14011r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;otal</span><span 
class="cmtt-12">&#x00A0;vm</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;2755568</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;otal</span><span 
class="cmtt-12">&#x00A0;vm</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;2765296</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;otal</span><span 
class="cmtt-12">&#x00A0;vm</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;2775024</span><br class="fancyvrb" /><a 
 id="x1-14013r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;otal</span><span 
class="cmtt-12">&#x00A0;vm</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;2784752</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;otal...</span>
<br class="fancyvrb" /><a 
 id="x1-14015r7"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Assertion</span><span 
class="cmtt-12">&#x00A0;failure:</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Compiler</span><span 
class="cmtt-12">&#x00A0;internal</span><span 
class="cmtt-12">&#x00A0;error</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;please</span><span 
class="cmtt-12">&#x00A0;submit</span><span 
class="cmtt-12">&#x00A0;problem</span><span 
class="cmtt-12">&#x00A0;r...</span>
<br class="fancyvrb" /><a 
 id="x1-14017r8"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;GEM</span><span 
class="cmtt-12">&#x00A0;ASSERTION,</span><span 
class="cmtt-12">&#x00A0;Compiler</span><span 
class="cmtt-12">&#x00A0;internal</span><span 
class="cmtt-12">&#x00A0;error</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;please</span><span 
class="cmtt-12">&#x00A0;submit</span><span 
class="cmtt-12">&#x00A0;problem</span><span 
class="cmtt-12">&#x00A0;report</span>
<br class="fancyvrb" /><a 
 id="x1-14019r9"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Fatal</span><span 
class="cmtt-12">&#x00A0;error</span><span 
class="cmtt-12">&#x00A0;in:</span><span 
class="cmtt-12">&#x00A0;/usr/lib/cmplrs/fort90_540/decfort90</span><span 
class="cmtt-12">&#x00A0;Terminated</span><br class="fancyvrb" /><a 
 id="x1-14021r10"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;Exit</span><span 
class="cmtt-12">&#x00A0;3</span>
<br class="fancyvrb" /><a 
 id="x1-14023r11"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Stop.</span><br class="fancyvrb" /><a 
 id="x1-14025r12"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;Exit</span><span 
class="cmtt-12">&#x00A0;1</span><br class="fancyvrb" /><a 
 id="x1-14027r13"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Stop.</span></div>
<!--l. 708--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The occurrence of this problem depends upon the grid size; and the problem never seems
to occur with &#8216;<span 
class="cmtt-12">mpicomm.f90</span>&#8217;<a 
 id="x1-14028"></a>, except when <span 
class="cmtt-12">ncpus=1</span>. The problem can be avoided by
switching off the loop transformation optimization (part of the &#8216;<span 
class="cmtt-12">-O5</span>&#8217;<a 
 id="x1-14029"></a> optimization),
via:
<div class="fancyvrb" id="fancyvrb14"><a 
 id="x1-14031r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;#OPTFLAGS=-fast</span><span 
class="cmtt-12">&#x00A0;-O5</span><span 
class="cmtt-12">&#x00A0;-notransform_loops</span></div>
                                                                                            
                                                                                            
<!--l. 722--><p class="noindent" >This is currently the default compiler setting in &#8216;<span 
class="cmtt-12">Makefile</span>&#8217;<a 
 id="x1-14032"></a>, although it has a measurable
performance impact (some 8% slowdown).
<!--l. 725--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.10    </span> <a 
 id="x1-150001.2.10"></a>Assertion failure under SunOS</h5>
<!--l. 726--><p class="noindent" >Under SunOS, I get an error message like
<div class="fancyvrb" id="fancyvrb15"><a 
 id="x1-15002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;user@sun&#x003E;</span><span 
class="cmtt-12">&#x00A0;f90</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;param_io.f90</span><br class="fancyvrb" /><a 
 id="x1-15004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Assertion</span><span 
class="cmtt-12">&#x00A0;failed:</span><span 
class="cmtt-12">&#x00A0;at_handle_table[at_idx].tag</span><span 
class="cmtt-12">&#x00A0;==</span><span 
class="cmtt-12">&#x00A0;VAR_TAG,</span>
<br class="fancyvrb" /><a 
 id="x1-15006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;file</span><span 
class="cmtt-12">&#x00A0;../srcfw/FWcvrt.c,</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;4018</span><br class="fancyvrb" /><a 
 id="x1-15008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;f90:</span><span 
class="cmtt-12">&#x00A0;Fatal</span><span 
class="cmtt-12">&#x00A0;error</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;f90comp:</span><span 
class="cmtt-12">&#x00A0;Abort</span></div>
<!--l. 736--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: This is a compiler bug that we find at least with Sun&#8217;s WorkShop Compiler version &#8216;5.0
00/05/17 FORTRAN 90 2.0 Patch 107356-05&#8217;. Upgrade the compiler version (and possibly also
the operating system): we find that the code compiles and works with version &#8216;Sun WorkShop
6 update 2 Fortran 95 6.2 Patch 111690-05 2002/01/17&#8217; under SunOS version &#8216;5.8
Generic&#x02D9;108528-11&#8217;.
<!--l. 743--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.11    </span> <a 
 id="x1-160001.2.11"></a>After some dirty tricks I got pencil code to compile with MPI, ...</h5>
<!--l. 745--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb16"><a 
 id="x1-16002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span><span 
class="cmtt-12">&#x00A0;Before</span><span 
class="cmtt-12">&#x00A0;that</span><span 
class="cmtt-12">&#x00A0;i</span><span 
class="cmtt-12">&#x00A0;installed</span><span 
class="cmtt-12">&#x00A0;lam-7.1.4</span><span 
class="cmtt-12">&#x00A0;from</span><span 
class="cmtt-12">&#x00A0;source.</span><br class="fancyvrb" /><a 
 id="x1-16004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span>
<br class="fancyvrb" /><a 
 id="x1-16006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Goodness</span><span 
class="cmtt-12">&#x00A0;gracious</span><span 
class="cmtt-12">&#x00A0;me,</span><span 
class="cmtt-12">&#x00A0;you</span><span 
class="cmtt-12">&#x00A0;shouldn&#8217;t</span><span 
class="cmtt-12">&#x00A0;have</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;compile</span><span 
class="cmtt-12">&#x00A0;your</span><span 
class="cmtt-12">&#x00A0;own</span><span 
class="cmtt-12">&#x00A0;MPI</span><span 
class="cmtt-12">&#x00A0;library.</span></div>
<!--l. 753--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Then don&#8217;t use the old LAM-MPI. It is long superseded by open-mpi now. Open-mpi doesn&#8217;t
need a daemon to be running. I am using the version that ships with Ubuntu (e.g.
9.04):
<!--l. 757--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb17"><a 
 id="x1-16008r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;frenesi:~&#x003E;</span><span 
class="cmtt-12">&#x00A0;aptitude</span><span 
class="cmtt-12">&#x00A0;-w</span><span 
class="cmtt-12">&#x00A0;210</span><span 
class="cmtt-12">&#x00A0;search</span><span 
class="cmtt-12">&#x00A0;openmpi</span><span 
class="cmtt-12">&#x00A0;|</span><span 
class="cmtt-12">&#x00A0;grep</span><span 
class="cmtt-12">&#x00A0;&#8217;^i&#8217;</span><br class="fancyvrb" /><a 
 id="x1-16010r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span>
<br class="fancyvrb" /><a 
 id="x1-16012r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;i</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;libopenmpi-dev</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;high</span><span 
class="cmtt-12">&#x00A0;performance</span><span 
class="cmtt-12">&#x00A0;message</span><span 
class="cmtt-12">&#x00A0;passing</span><span 
class="cmtt-12">&#x00A0;library</span><span 
class="cmtt-12">&#x00A0;--</span><span 
class="cmtt-12">&#x00A0;header</span><span 
class="cmtt-12">&#x00A0;files</span>
<br class="fancyvrb" /><a 
 id="x1-16014r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;i</span><span 
class="cmtt-12">&#x00A0;A</span><span 
class="cmtt-12">&#x00A0;libopenmpi1</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;high</span><span 
class="cmtt-12">&#x00A0;performance</span><span 
class="cmtt-12">&#x00A0;message</span><span 
class="cmtt-12">&#x00A0;passing</span><span 
class="cmtt-12">&#x00A0;library</span><span 
class="cmtt-12">&#x00A0;--</span><span 
class="cmtt-12">&#x00A0;shared</span><span 
class="cmtt-12">&#x00A0;library</span>
<br class="fancyvrb" /><a 
 id="x1-16016r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;i</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;openmpi-bin</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;high</span><span 
class="cmtt-12">&#x00A0;performance</span><span 
class="cmtt-12">&#x00A0;message</span><span 
class="cmtt-12">&#x00A0;passing</span><span 
class="cmtt-12">&#x00A0;library</span><span 
class="cmtt-12">&#x00A0;--</span><span 
class="cmtt-12">&#x00A0;binaries</span>
<br class="fancyvrb" /><a 
 id="x1-16018r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;i</span><span 
class="cmtt-12">&#x00A0;A</span><span 
class="cmtt-12">&#x00A0;openmpi-common</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;high</span><span 
class="cmtt-12">&#x00A0;performance</span><span 
class="cmtt-12">&#x00A0;message</span><span 
class="cmtt-12">&#x00A0;passing</span><span 
class="cmtt-12">&#x00A0;library</span><span 
class="cmtt-12">&#x00A0;--</span><span 
class="cmtt-12">&#x00A0;common</span><span 
class="cmtt-12">&#x00A0;files</span>
<br class="fancyvrb" /><a 
 id="x1-16020r7"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;i</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;openmpi-doc</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;-</span><span 
class="cmtt-12">&#x00A0;high</span><span 
class="cmtt-12">&#x00A0;performance</span><span 
class="cmtt-12">&#x00A0;message</span><span 
class="cmtt-12">&#x00A0;passing</span><span 
class="cmtt-12">&#x00A0;library</span><span 
class="cmtt-12">&#x00A0;--</span><span 
class="cmtt-12">&#x00A0;man</span><span 
class="cmtt-12">&#x00A0;pages</span></div>
<!--l. 767--><p class="noindent" >Install that and keep your configuration (Makefile.src and getconf.csh) close to that for &#8216;frenesi&#8217;
or &#8216;norlx50&#8217;. That should work. <a 
 id="x1-16021"></a>
<!--l. 771--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.12    </span> <a 
 id="x1-170001.2.12"></a>Error: Symbol &#8217;mpi&#x02D9;comm&#x02D9;world&#8217; at (1) has no IMPLICIT type</h5>
<!--l. 773--><p class="noindent" >
                                                                                            
                                                                                            
<div class="fancyvrb" id="fancyvrb18"><a 
 id="x1-17002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;I</span><span 
class="cmtt-12">&#x00A0;installed</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;pencil</span><span 
class="cmtt-12">&#x00A0;code</span><span 
class="cmtt-12">&#x00A0;on</span><span 
class="cmtt-12">&#x00A0;Ubuntu</span><span 
class="cmtt-12">&#x00A0;system</span><span 
class="cmtt-12">&#x00A0;and</span><span 
class="cmtt-12">&#x00A0;tested</span><span 
class="cmtt-12">&#x00A0;"run.csh"</span>
<br class="fancyvrb" /><a 
 id="x1-17004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;...\samples\conv-slab.</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Here</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;code</span><span 
class="cmtt-12">&#x00A0;worked</span><span 
class="cmtt-12">&#x00A0;pretty</span><span 
class="cmtt-12">&#x00A0;well.</span>
<br class="fancyvrb" /><a 
 id="x1-17006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Nevertheless,</span><span 
class="cmtt-12">&#x00A0;running</span><span 
class="cmtt-12">&#x00A0;(auto-test),</span><span 
class="cmtt-12">&#x00A0;I</span><span 
class="cmtt-12">&#x00A0;found</span><span 
class="cmtt-12">&#x00A0;there</span><span 
class="cmtt-12">&#x00A0;are</span><span 
class="cmtt-12">&#x00A0;some</span><span 
class="cmtt-12">&#x00A0;errors.</span><br class="fancyvrb" /><a 
 id="x1-17008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><br class="fancyvrb" /><a 
 id="x1-17010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;The</span><span 
class="cmtt-12">&#x00A0;messages</span><span 
class="cmtt-12">&#x00A0;are,</span><br class="fancyvrb" /><a 
 id="x1-17012r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span>
<br class="fancyvrb" /><a 
 id="x1-17014r7"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Error:</span><span 
class="cmtt-12">&#x00A0;Symbol</span><span 
class="cmtt-12">&#x00A0;&#8217;mpi_comm_world&#8217;</span><span 
class="cmtt-12">&#x00A0;at</span><span 
class="cmtt-12">&#x00A0;(1)</span><span 
class="cmtt-12">&#x00A0;has</span><span 
class="cmtt-12">&#x00A0;no</span><span 
class="cmtt-12">&#x00A0;IMPLICIT</span><span 
class="cmtt-12">&#x00A0;type</span><br class="fancyvrb" /><a 
 id="x1-17016r8"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Fatal</span><span 
class="cmtt-12">&#x00A0;Error:</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;count</span><span 
class="cmtt-12">&#x00A0;reached</span><span 
class="cmtt-12">&#x00A0;limit</span><span 
class="cmtt-12">&#x00A0;of</span><span 
class="cmtt-12">&#x00A0;25.</span>
<br class="fancyvrb" /><a 
 id="x1-17018r9"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[2]:</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;[mpicomm_double.o]</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;1</span><br class="fancyvrb" /><a 
 id="x1-17020r10"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[2]:</span><span 
class="cmtt-12">&#x00A0;Leaving</span><span 
class="cmtt-12">&#x00A0;directory</span>
<br class="fancyvrb" /><a 
 id="x1-17022r11"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#8216;/home/pkiwan/Desktop/pencil-code/samples/2d-tests/selfgravitating-shearwave/src&#8217;</span>
<br class="fancyvrb" /><a 
 id="x1-17024r12"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[1]:</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;[code]</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;2</span><br class="fancyvrb" /><a 
 id="x1-17026r13"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make[1]:</span><span 
class="cmtt-12">&#x00A0;Leaving</span><span 
class="cmtt-12">&#x00A0;directory</span>
<br class="fancyvrb" /><a 
 id="x1-17028r14"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#8216;/home/pkiwan/Desktop/pencil-code/samples/2d-tests/selfgravitating-shearwave/src&#8217;</span>
<br class="fancyvrb" /><a 
 id="x1-17030r15"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;make:</span><span 
class="cmtt-12">&#x00A0;***</span><span 
class="cmtt-12">&#x00A0;[default]</span><span 
class="cmtt-12">&#x00A0;Error</span><span 
class="cmtt-12">&#x00A0;2</span><br class="fancyvrb" /><a 
 id="x1-17032r16"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><br class="fancyvrb" /><a 
 id="x1-17034r17"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Finally,</span><span 
class="cmtt-12">&#x00A0;###</span><span 
class="cmtt-12">&#x00A0;auto-test</span><span 
class="cmtt-12">&#x00A0;failed</span><span 
class="cmtt-12">&#x00A0;###</span><br class="fancyvrb" /><a 
 id="x1-17036r18"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><br class="fancyvrb" /><a 
 id="x1-17038r19"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Will</span><span 
class="cmtt-12">&#x00A0;it</span><span 
class="cmtt-12">&#x00A0;be</span><span 
class="cmtt-12">&#x00A0;OK?</span><span 
class="cmtt-12">&#x00A0;Or,</span><span 
class="cmtt-12">&#x00A0;how</span><span 
class="cmtt-12">&#x00A0;can</span><span 
class="cmtt-12">&#x00A0;I</span><span 
class="cmtt-12">&#x00A0;fix</span><span 
class="cmtt-12">&#x00A0;this?</span></div>
<!--l. 797--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Thanks for letting me know about the status, and congratulations on your progress! Those
tests that fail are those that use MPI. If your machine is a dual or multi core machine, you
could run faster by running under MPI. But this is probably not crucial for you at this point. (I
just noticed that there is a ToDo listed in the auto-test command to implement the option not
to run the MPI tests, but this hasn&#8217;t been done yet. So I guess you can start with the science
next.
<!--l. 807--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.2.13    </span> <a 
 id="x1-180001.2.13"></a>Error: Can&#8217;t open included file &#8217;mpif.h&#8217;</h5>
<!--l. 809--><p class="noindent" >It always worked, but now, after some systems upgrade, I get
<div class="fancyvrb" id="fancyvrb19"><a 
 id="x1-18002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;gfortran</span><span 
class="cmtt-12">&#x00A0;-O3</span><span 
class="cmtt-12">&#x00A0;-o</span><span 
class="cmtt-12">&#x00A0;mpicomm.o</span><span 
class="cmtt-12">&#x00A0;-c</span><span 
class="cmtt-12">&#x00A0;mpicomm.f90</span><br class="fancyvrb" /><a 
 id="x1-18004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Error:</span><span 
class="cmtt-12">&#x00A0;Can&#8217;t</span><span 
class="cmtt-12">&#x00A0;open</span><span 
class="cmtt-12">&#x00A0;included</span><span 
class="cmtt-12">&#x00A0;file</span><span 
class="cmtt-12">&#x00A0;&#8217;mpif.h&#8217;</span></div>
<!--l. 815--><p class="noindent" >When I say <span 
class="cmtt-12">locate mpif.h</span><a 
 id="x1-18005"></a> I only get things like
<div class="fancyvrb" id="fancyvrb20"><a 
 id="x1-18007r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/scratch/ntest/1.2.7p1-intel/include/mpif.h</span></div>
<!--l. 819--><p class="noindent" >But since I use <span 
class="cmtt-12">FC=mpif90</span><a 
 id="x1-18008"></a> I thought I don&#8217;t need to worry.
<!--l. 823--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Since you use <span 
class="cmtt-12">FC=mpif90</span><a 
 id="x1-18009"></a> there must definitely be something wrong with their setup. Try
<span 
class="cmtt-12">mpif90 -showme</span><a 
 id="x1-18010"></a> or <span 
class="cmtt-12">mpif90 -show</span><a 
 id="x1-18011"></a>; the &#8216;-I&#8217; option should say where it looks for &#8217;mpif.h&#8217;. If
those directories don&#8217;t exist, it&#8217;s no wonder that it doesn&#8217;t work, and it is time to
complain.
<!--l. 832--><p class="noindent" >
<h4 class="subsectionHead"><span class="titlemark">1.3    </span> <a 
 id="x1-190001.3"></a>Pencil check</h4>
<!--l. 834--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.3.1    </span> <a 
 id="x1-200001.3.1"></a>The pencil check complains for no reason.</h5>
<!--l. 838--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The pencil check only complains for a reason.
                                                                                            
                                                                                            
<!--l. 840--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.3.2    </span> <a 
 id="x1-210001.3.2"></a>The pencil check reports MISSING PENCILS and quits</h5>
<!--l. 844--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: This could point to a serious problem in the code. Check where the missing pencil is used in
the code. Request the right pencils, likely based on input parameters, by adapting one or more
of the <span 
class="cmtt-12">pencil_criteria_MODULE </span>subroutines.
<!--l. 849--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.3.3    </span> <a 
 id="x1-220001.3.3"></a>The pencil check reports unnecessary pencils</h5>
<!--l. 851--><p class="noindent" >The pencil check reports <span 
class="cmtt-12">possible overcalculation... pencil rho ( 43) is requested, but</span>
<span 
class="cmtt-12">does not appear to be required!</span>
<!--l. 857--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Such warnings show that your simulation is possibly running too slowly because it is
calculating pencils that are not actually needed. Check in the code where the unnecessary
pencils are used and adapt one or more of the <span 
class="cmtt-12">pencil_criteria_MODULE </span>subroutines to request
pencils only when they are actually needed.
<!--l. 863--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.3.4    </span> <a 
 id="x1-230001.3.4"></a>The pencil check reports that most or all pencils are missing</h5>
<!--l. 867--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: This is typically a thing that can happen when testing new code development for the
first time. It is usually an indication that the reference <span 
class="cmtt-12">df </span>changes every time you
call <span 
class="cmtt-12">pde</span>. Check whether any newly implemented subroutines or functionality has a
&#8220;memory&#8221;, i.e. if calling the subroutine twice with the same <span 
class="cmtt-12">f </span>gives different output
<span 
class="cmtt-12">df</span>.
<!--l. 873--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.3.5    </span> <a 
 id="x1-240001.3.5"></a>Running the pencil check triggers mathematical errors in the code</h5>
<!--l. 875--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The pencil check puts random numbers in <span 
class="cmtt-12">f </span>before checking the dependence of <span 
class="cmtt-12">df </span>on the
chosen set of pencils. Sometimes these random numbers are inconsistent with the physics and
cause errors. In that case you can set <span 
class="cmtt-12">lrandom_f_pencil_check=F </span>in <span 
class="cmtt-12">&amp;run_pars </span>in &#8216;<span 
class="cmtt-12">run.in</span>&#8217;<a 
 id="x1-24001"></a>. The
initial condition may contain many idealized states (zeros or ones) which then do not trigger
pencil check errors when <span 
class="cmtt-12">lrandom_f_pencil_check=F</span>, even if pencils are missing. But it does
prevent mathematical inconsistencies.
<!--l. 883--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.3.6    </span> <a 
 id="x1-250001.3.6"></a>The pencil check still complains</h5>
                                                                                            
                                                                                            
<!--l. 887--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Then you need to look into the how the code and the pencil check operate. Reduce the
problem in size and dimensions to find the smallest problem that makes the pencil check fail
(e.g. 1x1x8 grid points). At the line of &#8216;<span 
class="cmtt-12">pencil_check.f90</span>&#8217;<a 
 id="x1-25001"></a> when a difference is found between
<span 
class="cmtt-12">df_ref </span>and <span 
class="cmtt-12">df</span>, add some debug lines telling you which variable is inconsistent and in what
place. Often you will be surprised that the pencil check has correctly found a problem in the
simulation.
<!--l. 895--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.3.7    </span> <a 
 id="x1-260001.3.7"></a>The pencil check is annoying so I turned it off</h5>
<!--l. 899--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Then you are taking a major risk. If one or more pencils are not calculated properly, then
the results will be wrong.
<!--l. 904--><p class="noindent" >
<h4 class="subsectionHead"><span class="titlemark">1.4    </span> <a 
 id="x1-270001.4"></a>Running</h4>
<!--l. 908--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.1    </span> <a 
 id="x1-280001.4.1"></a>Why does &#8216;<span 
class="cmtt-12">start.x</span>&#8217; / &#8216;<span 
class="cmtt-12">start.csh</span>&#8217; write data with periodic boundary conditions?</h5>
<!--l. 912--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Because you are setting the boundary conditions in &#8216;<span 
class="cmtt-12">run.in</span>&#8217;<a 
 id="x1-28001"></a>, not in &#8216;<span 
class="cmtt-12">start.in</span>&#8217;<a 
 id="x1-28002"></a>; see Sect.&#x00A0;<span 
class="pncb7t-x-x-120">??</span>.
There is nothing wrong with the initial data &#8212; the ghost-zone values will be re-calculated
during the very first time step.
<!--l. 917--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.2    </span> <a 
 id="x1-290001.4.2"></a>csh problem?</h5>
<!--l. 918--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: On some rare occasions we have problems with csh not being supported on other machines.
(We hope to fix this by contacting the responsible person, but may not be that trivial today!)
Oliver says this is a well known bug of some years ago, etc. But maybe in the long run it would
be good to avoid csh.
<!--l. 927--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: These occasions will become increasingly frequent, and eventually for some architectures,
there may not even be a csh variant that can be installed.
<!--l. 931--><p class="noindent" >We never pushed people to use <span 
class="cmtt-12">pc_run</span><a 
 id="x1-29001"></a> and friends (and to report corresponding bugs and get
them fixed), but if we don&#8217;t spend a bit of effort (or annoy users) now, we create a future
emergency, where someone needs to run on some machine, but there is no csh and he or she
just gets stuck.
                                                                                            
                                                                                            
<!--l. 937--><p class="noindent" >We don&#8217;t have that many csh files, and for years now it should be possible to compile run
without csh (using <span 
class="cmtt-12">bin/pc_run</span><a 
 id="x1-29002"></a>) &#8212; except that people still fall back on the old way of doing
things. This is both cause and consequence of the &#8216;new&#8217; way not being tested that much, at
least for the corner cases like &#8216;<span 
class="cmtt-12">RERUN</span>&#8217;<a 
 id="x1-29003"></a>, &#8216;<span 
class="cmtt-12">NEWDIR</span>&#8217;<a 
 id="x1-29004"></a>, &#8216;<span 
class="cmtt-12">SCRATCH_DIR</span>&#8217;<a 
 id="x1-29005"></a>.
<!--l. 944--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.3    </span> <a 
 id="x1-300001.4.3"></a>&#8216;<span 
class="cmtt-12">run.csh</span>&#8217; doesn&#8217;t work:</h5>
<!--l. 945--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb21"><a 
 id="x1-30002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Invalid</span><span 
class="cmtt-12">&#x00A0;character</span><span 
class="cmtt-12">&#x00A0;&#8217;&#8217;&#8217;</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;NAMELIST</span><span 
class="cmtt-12">&#x00A0;input</span><br class="fancyvrb" /><a 
 id="x1-30004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Program</span><span 
class="cmtt-12">&#x00A0;terminated</span><span 
class="cmtt-12">&#x00A0;by</span><span 
class="cmtt-12">&#x00A0;fatal</span><span 
class="cmtt-12">&#x00A0;I/O</span><span 
class="cmtt-12">&#x00A0;error</span><br class="fancyvrb" /><a 
 id="x1-30006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Abort</span></div>
<!--l. 951--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The string array for the boundary condition, e.g.&#x00A0;<span 
class="pncro7t-x-x-120">bcx</span><a 
 id="x1-30007"></a> or <span 
class="pncro7t-x-x-120">bcz</span><a 
 id="x1-30008"></a> is too long. Make sure it has
exactly as many elements as <span 
class="pncro7t-x-x-120">nvar</span><a 
 id="x1-30009"></a> is big.
<!--l. 955--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.4    </span> <a 
 id="x1-310001.4.4"></a>Namelist problem under IRIX</h5>
<!--l. 956--><p class="noindent" >Under IRIX, I get
<div class="fancyvrb" id="fancyvrb22"><a 
 id="x1-31002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;lib-4001</span><span 
class="cmtt-12">&#x00A0;:</span><span 
class="cmtt-12">&#x00A0;UNRECOVERABLE</span><span 
class="cmtt-12">&#x00A0;library</span><span 
class="cmtt-12">&#x00A0;error</span><br class="fancyvrb" /><a 
 id="x1-31004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><br class="fancyvrb" /><a 
 id="x1-31006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Encountered</span><span 
class="cmtt-12">&#x00A0;during</span><span 
class="cmtt-12">&#x00A0;a</span><span 
class="cmtt-12">&#x00A0;namelist</span><span 
class="cmtt-12">&#x00A0;READ</span><span 
class="cmtt-12">&#x00A0;from</span><span 
class="cmtt-12">&#x00A0;unit</span><span 
class="cmtt-12">&#x00A0;1</span>
<br class="fancyvrb" /><a 
 id="x1-31008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Fortran</span><span 
class="cmtt-12">&#x00A0;unit</span><span 
class="cmtt-12">&#x00A0;1</span><span 
class="cmtt-12">&#x00A0;is</span><span 
class="cmtt-12">&#x00A0;connected</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;a</span><span 
class="cmtt-12">&#x00A0;sequential</span><span 
class="cmtt-12">&#x00A0;formatted</span><span 
class="cmtt-12">&#x00A0;text</span><span 
class="cmtt-12">&#x00A0;file:</span><span 
class="cmtt-12">&#x00A0;"run.in"</span><br class="fancyvrb" /><a 
 id="x1-31010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;IOT</span><span 
class="cmtt-12">&#x00A0;Trap</span><br class="fancyvrb" /><a 
 id="x1-31012r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Abort</span></div>
<!--l. 968--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: This is a compiler bug that has been found at least with the MIPSpro F90 compiler version
7.3.1.3m. The problem seems to have been fixed in version 7.4.20m.
<!--l. 975--><p class="noindent" >The error comes and goes, depending on the configuration (and possibly even the input
parameters) you are using. Until SGI fix their compiler, you can experiment with adding new
variables to the module <span 
class="pncro7t-x-x-120">Param&#x02D9;IO</span><a 
 id="x1-31013"></a>; this has solved the problem once for us. If this trick does
not help, you will need to turn your namelist input (at least &#8216;<span 
class="cmtt-12">run.in</span>&#8217;<a 
 id="x1-31014"></a> into Fortran statements,
include them into a replacement version of &#8216;<span 
class="cmtt-12">param_io.f90</span>&#8217;<a 
 id="x1-31015"></a>, and recompile each time you make
changes.
<!--l. 983--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.5    </span> <a 
 id="x1-320001.4.5"></a>Code crashes after restarting</h5>
<!--l. 985--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb23"><a 
 id="x1-32002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span><span 
class="cmtt-12">&#x00A0;removing</span><span 
class="cmtt-12">&#x00A0;mu_r</span><span 
class="cmtt-12">&#x00A0;from</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;namelist</span><span 
class="cmtt-12">&#x00A0;just</span><span 
class="cmtt-12">&#x00A0;&#8216;like</span><span 
class="cmtt-12">&#x00A0;that&#8217;</span><span 
class="cmtt-12">&#x00A0;makes</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;code</span><br class="fancyvrb" /><a 
 id="x1-32004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span><span 
class="cmtt-12">&#x00A0;backwards</span><span 
class="cmtt-12">&#x00A0;incompatible.</span><br class="fancyvrb" /><a 
 id="x1-32006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span>
<br class="fancyvrb" /><a 
 id="x1-32008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span><span 
class="cmtt-12">&#x00A0;That</span><span 
class="cmtt-12">&#x00A0;means</span><span 
class="cmtt-12">&#x00A0;that</span><span 
class="cmtt-12">&#x00A0;we</span><span 
class="cmtt-12">&#x00A0;can</span><span 
class="cmtt-12">&#x00A0;never</span><span 
class="cmtt-12">&#x00A0;get</span><span 
class="cmtt-12">&#x00A0;rid</span><span 
class="cmtt-12">&#x00A0;of</span><span 
class="cmtt-12">&#x00A0;a</span><span 
class="cmtt-12">&#x00A0;parameter</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;start.in</span><span 
class="cmtt-12">&#x00A0;once</span><span 
class="cmtt-12">&#x00A0;we</span><br class="fancyvrb" /><a 
 id="x1-32010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#x003E;</span><span 
class="cmtt-12">&#x00A0;have</span><span 
class="cmtt-12">&#x00A0;introduced</span><span 
class="cmtt-12">&#x00A0;it,</span><span 
class="cmtt-12">&#x00A0;right?</span></div>
<!--l. 995--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: In the current implementation, without a corresponding cleaning procedure, unfortunately
yes.
<!--l. 998--><p class="noindent" >Of course, this does not affect users&#8217; private changes outside the central svn tree.
                                                                                            
                                                                                            
<!--l. 1032--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.6    </span> <a 
 id="x1-330001.4.6"></a>auto-test gone mad...?</h5>
<!--l. 1034--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: Have you ever seen this before:
<div class="fancyvrb" id="fancyvrb24"><a 
 id="x1-33002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;giga01:/home/pg/n7026413/cvs-src/pencil-code/samples/conv-slab&#x003E;</span><span 
class="cmtt-12">&#x00A0;auto-test</span><br class="fancyvrb" /><a 
 id="x1-33004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;.</span><br class="fancyvrb" /><a 
 id="x1-33006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span>
<br class="fancyvrb" /><a 
 id="x1-33008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/home/pg/n7026413/cvs-src/pencil-code/samples/conv-slab:</span><br class="fancyvrb" /><a 
 id="x1-33010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Compiling..</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ok</span>
<br class="fancyvrb" /><a 
 id="x1-33012r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;No</span><span 
class="cmtt-12">&#x00A0;data</span><span 
class="cmtt-12">&#x00A0;directory;</span><span 
class="cmtt-12">&#x00A0;generating</span><span 
class="cmtt-12">&#x00A0;data</span><span 
class="cmtt-12">&#x00A0;-&#x003E;</span><span 
class="cmtt-12">&#x00A0;/var/tmp/pencil-tmp-25318</span>
<br class="fancyvrb" /><a 
 id="x1-33014r7"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Starting..</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ok</span><br class="fancyvrb" /><a 
 id="x1-33016r8"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Running..</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;ok</span>
<br class="fancyvrb" /><a 
 id="x1-33018r9"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Validating</span><span 
class="cmtt-12">&#x00A0;results..Malformed</span><span 
class="cmtt-12">&#x00A0;UTF-8</span><span 
class="cmtt-12">&#x00A0;character</span><span 
class="cmtt-12">&#x00A0;(unexpected</span><span 
class="cmtt-12">&#x00A0;continuation</span>
<br class="fancyvrb" /><a 
 id="x1-33020r10"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;byte</span><span 
class="cmtt-12">&#x00A0;0x80,</span><span 
class="cmtt-12">&#x00A0;with</span><span 
class="cmtt-12">&#x00A0;no</span><span 
class="cmtt-12">&#x00A0;preceding</span><span 
class="cmtt-12">&#x00A0;start</span><span 
class="cmtt-12">&#x00A0;byte)</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;split</span><span 
class="cmtt-12">&#x00A0;at</span>
<br class="fancyvrb" /><a 
 id="x1-33022r11"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/home/pg/n7026413/cvs-src/pencil-code/bin/auto-test</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;263.</span>
<br class="fancyvrb" /><a 
 id="x1-33024r12"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Malformed</span><span 
class="cmtt-12">&#x00A0;UTF-8</span><span 
class="cmtt-12">&#x00A0;character</span><span 
class="cmtt-12">&#x00A0;(unexpected</span><span 
class="cmtt-12">&#x00A0;continuation</span><span 
class="cmtt-12">&#x00A0;byte</span><span 
class="cmtt-12">&#x00A0;0x80,</span><span 
class="cmtt-12">&#x00A0;with</span><span 
class="cmtt-12">&#x00A0;no</span>
<br class="fancyvrb" /><a 
 id="x1-33026r13"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;preceding</span><span 
class="cmtt-12">&#x00A0;start</span><span 
class="cmtt-12">&#x00A0;byte)</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;split</span><span 
class="cmtt-12">&#x00A0;at</span><br class="fancyvrb" /><a 
 id="x1-33028r14"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/home/pg/n7026413/cvs-src/pencil-code/bin/auto-test</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;263.</span></div>
<!--l. 1052--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: You are running on a RedHat 8 or 9 system, right?
<!--l. 1055--><p class="noindent" >Set <span 
class="cmtt-12">LANG=POSIX</span><a 
 id="x1-33029"></a> in your shell&#8217;s startup script and life will be much better.
<!--l. 1058--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.7    </span> <a 
 id="x1-340001.4.7"></a>Can I restart with a different number of cpus?</h5>
<!--l. 1060--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: I am running a simulation of nonhelical turbulence on the cluster using MPI. Suppose if I
am running a <span 
class="cmr-12">128</span><sup><span 
class="cmr-8">3</span></sup> simulation on 32 cpus/cores i.e.
<div class="fancyvrb" id="fancyvrb25"><a 
 id="x1-34002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;integer,</span><span 
class="cmtt-12">&#x00A0;parameter</span><span 
class="cmtt-12">&#x00A0;::</span><span 
class="cmtt-12">&#x00A0;ncpus=32,nprocy=2,nprocz=ncpus/nprocy,nprocx=1</span>
<br class="fancyvrb" /><a 
 id="x1-34004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;integer,</span><span 
class="cmtt-12">&#x00A0;parameter</span><span 
class="cmtt-12">&#x00A0;::</span><span 
class="cmtt-12">&#x00A0;nxgrid=128,nygrid=nxgrid,nzgrid=nxgrid</span></div>
<!--l. 1067--><p class="noindent" >And I stop the run after a bit. Is there a way to resume this run with different number of cpus
like this :
<div class="fancyvrb" id="fancyvrb26"><a 
 id="x1-34006r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;integer,</span><span 
class="cmtt-12">&#x00A0;parameter</span><span 
class="cmtt-12">&#x00A0;::</span><span 
class="cmtt-12">&#x00A0;ncpus=16,nprocy=2,nprocz=ncpus/nprocy,nprocx=1</span>
<br class="fancyvrb" /><a 
 id="x1-34008r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;integer,</span><span 
class="cmtt-12">&#x00A0;parameter</span><span 
class="cmtt-12">&#x00A0;::</span><span 
class="cmtt-12">&#x00A0;nxgrid=128,nygrid=nxgrid,nzgrid=nxgrid</span></div>
<!--l. 1073--><p class="noindent" >I understand it has to be so in a new directory but making sure that the run starts from where
I left it off in the previous directory.
<!--l. 1076--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The answer is no, if you use the standard distributed io. There is also parallel io, but I never
used it. That would write the data in a single file, and then you could use the data for restart
in another processor layout.
<!--l. 1082--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.8    </span> <a 
 id="x1-350001.4.8"></a>Can I restart with a different number of cpus?</h5>
                                                                                            
                                                                                            
<!--l. 1084--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: Is it right that once the simulation is resumed, pencil-code takes the last data from var.dat
(which is the current snapshot of the fields)? If that is true, then, is it not possible to give that
as the initial condition for the run in the second directory (with changed &#8221;ncpus&#8221;)? Is there a
mechanism already in place for that?
<!--l. 1091--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Yes, the code restarts from the last var.dat. It is written after a successful completion of the
run, but it crashes or you hit a time-out, there will be a var.dat that is overwritten every isave
timesteps. If the system stops during writing, some var.dat files may be corrupt or have the
wrong time. In that case you could restart from a good VAR file, if you have one, using,
e.g.,
<div class="fancyvrb" id="fancyvrb27"><a 
 id="x1-35002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;restart-new-dir-VAR</span><span 
class="cmtt-12">&#x00A0;.</span><span 
class="cmtt-12">&#x00A0;46</span></div>
<!--l. 1101--><p class="noindent" >where 46 is the number of your VAR file, i.e., VAR46 im this case. To restart in another
directory, you say, from the old run directory,
<div class="fancyvrb" id="fancyvrb28"><a 
 id="x1-35004r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;restart-new-dir</span><span 
class="cmtt-12">&#x00A0;../another_directory</span></div>
<!--l. 1106--><p class="noindent" >Hope this helps. Look into pencil-code/bin/restart-new-dir to see what it is doing.
<!--l. 1109--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.9    </span> <a 
 id="x1-360001.4.9"></a>fft_xyz_parallel_3D: nygrid needs to be an integer multiple...</h5>
<!--l. 1111--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: I just got an:
<div class="fancyvrb" id="fancyvrb29"><a 
 id="x1-36002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;fft_xyz_parallel_3D:</span><span 
class="cmtt-12">&#x00A0;nygrid</span><span 
class="cmtt-12">&#x00A0;needs</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;be</span><span 
class="cmtt-12">&#x00A0;an</span><span 
class="cmtt-12">&#x00A0;integer</span><span 
class="cmtt-12">&#x00A0;multiple</span><span 
class="cmtt-12">&#x00A0;of</span><span 
class="cmtt-12">&#x00A0;nprocy*nprocz</span></div>
<!--l. 1116--><p class="noindent" >In my case, nygrid=2048, nprocy=32, and nprocz=128, so nprocy*nprocz=4096. In other words,
2048 needs to be a multiple of 4096. But isn&#8217;t this the case then?
<!--l. 1120--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: No, because 2048 = 0.5 * 4096 and 0.5 is not an integer. Maybe try either setting nprocz=64
or nprocy=64. You could compensate the change of ncpus with the <span 
class="cmmi-12">x</span>-direction. For <span 
class="cmr-12">2048</span><sup><span 
class="cmr-8">3</span></sup>
simulations, nprocy=32 and nprocz=64 would be good.
<!--l. 1126--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.4.10    </span> <a 
 id="x1-370001.4.10"></a>Unit-agnostic calculations?</h5>
<!--l. 1128--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: The manual speaks about unit-agnostic calculations, stating that one may choose to
interpret the results in any (consistent) units, depending on the problem that is solved at
hand. So, for example, if I chose to run the &#8216;<span 
class="cmtt-12">2d-tests/battery_term</span>&#8217;<a 
 id="x1-37001"></a> simulation for an arbitrary
number of time-steps and then choose to examine the diagnostics, am I correct in assuming the
following:
                                                                                            
                                                                                            
<div class="verbatim" id="verbatim-1">
1)&#x00A0;[Brms]&#x00A0;=&#x00A0;Gauss&#x00A0;(as&#x00A0;output&#x00A0;by&#x00A0;unit_magnetic,&#x00A0;before&#x00A0;the&#x00A0;run&#x00A0;begins)
&#x00A0;<br />2)&#x00A0;[t]&#x00A0;=&#x00A0;s&#x00A0;(since&#x00A0;the&#x00A0;default&#x00A0;unit&#x00A0;system&#x00A0;is&#x00A0;left&#x00A0;as&#x00A0;CGS)
&#x00A0;<br />3)&#x00A0;[urms]&#x00A0;=&#x00A0;cm/s&#x00A0;(again,&#x00A0;as&#x00A0;output&#x00A0;by&#x00A0;unit_velocity,&#x00A0;before&#x00A0;the&#x00A0;run&#x00A0;begins)
&#x00A0;<br />4)&#x00A0;and&#x00A0;etc.&#x00A0;for&#x00A0;the&#x00A0;units&#x00A0;of&#x00A0;the&#x00A0;other&#x00A0;diagnostics</div>
<!--l. 1141--><p class="nopar" >
<!--l. 1143--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Detailed correspondence on this item can be found on:
<a 
href="https://groups.google.com/forum/?fromgroups#!topic/pencil-code-discuss/zek-uYNbgXI" class="url" ><span 
class="cmtt-12">https://groups.google.com/forum/?fromgroups#!topic/pencil-code-discuss/zek-uYNbgXI</span></a> There is also working
material on unit systems under <a 
href="http://www.nordita.org/~brandenb/teach/PencilCode/MixedTopics.html" class="url" ><span 
class="cmtt-12">http://www.nordita.org/</span><span 
class="cmtt-12">~</span><span 
class="cmtt-12">brandenb/teach/PencilCode/MixedTopics.html</span></a>
with a link to <a 
href="http://www.nordita.org/~brandenb/teach/PencilCode/material/AlfvenWave_SIunits/" class="url" ><span 
class="cmtt-12">http://www.nordita.org/</span><span 
class="cmtt-12">~</span><span 
class="cmtt-12">brandenb/teach/PencilCode/material/AlfvenWave_SIunits/</span></a>
Below is a pedagogical response from Wlad Lyra:
                                                                                            
                                                                                            
<div class="verbatim" id="verbatim-2">
In&#x00A0;the&#x00A0;sample&#x00A0;battery-term,&#x00A0;the&#x00A0;sound&#x00A0;speed&#x00A0;cs0=1&#x00A0;sets&#x00A0;the&#x00A0;unit&#x00A0;of
&#x00A0;<br />velocity.&#x00A0;Together&#x00A0;with&#x00A0;the&#x00A0;unit&#x00A0;of&#x00A0;length,&#x00A0;that&#x00A0;sets&#x00A0;your&#x00A0;unit&#x00A0;of
&#x00A0;<br />time.&#x00A0;The&#x00A0;unit&#x00A0;of&#x00A0;magnetic&#x00A0;field&#x00A0;follows&#x00A0;from&#x00A0;the&#x00A0;unit&#x00A0;of&#x00A0;velocity,
&#x00A0;<br />density,&#x00A0;and&#x00A0;your&#x00A0;choice&#x00A0;of&#x00A0;magnetic&#x00A0;permittivity,&#x00A0;according&#x00A0;to&#x00A0;the
&#x00A0;<br />definition&#x00A0;of&#x00A0;the&#x00A0;Alfven&#x00A0;velocity.</div>
<!--l. 1157--><p class="nopar" >
                                                                                            
                                                                                            
<div class="verbatim" id="verbatim-3">
If&#x00A0;you&#x00A0;are&#x00A0;assuming&#x00A0;cgs,&#x00A0;you&#x00A0;are&#x00A0;saying&#x00A0;that&#x00A0;your&#x00A0;sound&#x00A0;speed&#x00A0;cs0=1
&#x00A0;<br />actually&#x00A0;means&#x00A0;[U]=1&#x00A0;cm/s.&#x00A0;Your&#x00A0;unit&#x00A0;of&#x00A0;length&#x00A0;is&#x00A0;equivalently&#x00A0;1&#x00A0;cm,
&#x00A0;<br />and&#x00A0;therefore&#x00A0;the&#x00A0;unit&#x00A0;of&#x00A0;time&#x00A0;is&#x00A0;[t]&#x00A0;=&#x00A0;[L]/[U]=1&#x00A0;s.&#x00A0;The&#x00A0;unit&#x00A0;of
&#x00A0;<br />density&#x00A0;is&#x00A0;[rho]&#x00A0;=&#x00A0;1&#x00A0;g/cm^3.&#x00A0;Since&#x00A0;in&#x00A0;cgs&#x00A0;vA=B/sqrt(4*pi&#x00A0;*&#x00A0;rho),&#x00A0;your
&#x00A0;<br />unit&#x00A0;of&#x00A0;magnetic&#x00A0;field&#x00A0;is&#x00A0;[B]&#x00A0;=&#x00A0;[U]&#x00A0;*&#x00A0;sqrt([rho]&#x00A0;*&#x00A0;4*pi)&#x00A0;~=&#x00A0;3.5
&#x00A0;<br />sqrt(g/cm)&#x00A0;/&#x00A0;s&#x00A0;=&#x00A0;3.5&#x00A0;Gauss.</div>
<!--l. 1166--><p class="nopar" >
                                                                                            
                                                                                            
<div class="verbatim" id="verbatim-4">
If&#x00A0;instead&#x00A0;you&#x00A0;are&#x00A0;assuming&#x00A0;SI,&#x00A0;you&#x00A0;have&#x00A0;cs0=1&#x00A0;assuming&#x00A0;that&#x00A0;means
&#x00A0;<br />[U]=1&#x00A0;m/s&#x00A0;and&#x00A0;rho0=1&#x00A0;assuming&#x00A0;that&#x00A0;to&#x00A0;mean&#x00A0;[rho]=1&#x00A0;kg/m^3.&#x00A0;Using&#x00A0;[L]=1
&#x00A0;<br />m,&#x00A0;you&#x00A0;have&#x00A0;still&#x00A0;[t]=1&#x00A0;s,&#x00A0;but&#x00A0;now&#x00A0;what&#x00A0;appears&#x00A0;as&#x00A0;B=1&#x00A0;in&#x00A0;your&#x00A0;output
&#x00A0;<br />is&#x00A0;actually&#x00A0;[B]&#x00A0;=&#x00A0;[U]&#x00A0;*&#x00A0;sqrt&#x00A0;(mu&#x00A0;*&#x00A0;[rho])&#x00A0;=&#x00A0;1&#x00A0;m/s&#x00A0;*&#x00A0;sqrt(4*pi&#x00A0;*&#x00A0;1e-7
&#x00A0;<br />N*A-2&#x00A0;&#x00A0;1&#x00A0;kg/m^3)&#x00A0;&#x00A0;~=&#x00A0;0.0011210&#x00A0;kg/(s^2*A)&#x00A0;~&#x00A0;11&#x00A0;Gauss.</div>
<!--l. 1174--><p class="nopar" >
                                                                                            
                                                                                            
<div class="verbatim" id="verbatim-5">
You&#x00A0;can&#x00A0;make&#x00A0;it&#x00A0;more&#x00A0;interesting&#x00A0;and&#x00A0;use&#x00A0;units&#x00A0;relevant&#x00A0;to&#x00A0;the
&#x00A0;<br />problem.&#x00A0;Say&#x00A0;you&#x00A0;are&#x00A0;at&#x00A0;the&#x00A0;photosphere&#x00A0;of&#x00A0;the&#x00A0;Sun.&#x00A0;You&#x00A0;may&#x00A0;want&#x00A0;to
&#x00A0;<br />use&#x00A0;dimensionless&#x00A0;cs0=1&#x00A0;meaning&#x00A0;a&#x00A0;sound&#x00A0;speed&#x00A0;of&#x00A0;10&#x00A0;km/s.&#x00A0;&#x00A0;Your
&#x00A0;<br />appropriate&#x00A0;length&#x00A0;can&#x00A0;be&#x00A0;a&#x00A0;megameter.&#x00A0;Now&#x00A0;your&#x00A0;time&#x00A0;unit&#x00A0;is
&#x00A0;<br />[t]=[L]/[U]&#x00A0;=&#x00A0;1e3&#x00A0;km/&#x00A0;10&#x00A0;km/s&#x00A0;=&#x00A0;10^2&#x00A0;s,&#x00A0;i.e.,&#x00A0;roughly&#x00A0;1.5&#x00A0;minute.&#x00A0;For
&#x00A0;<br />density,&#x00A0;assume&#x00A0;rho=2x10-4&#x00A0;kg/m^3,&#x00A0;typical&#x00A0;of&#x00A0;the&#x00A0;solar&#x00A0;photosphere.
&#x00A0;<br />Your&#x00A0;unit&#x00A0;of&#x00A0;magnetic&#x00A0;field&#x00A0;is&#x00A0;therefore&#x00A0;[B]&#x00A0;=&#x00A0;[U]&#x00A0;*&#x00A0;sqrt([rho]&#x00A0;*
&#x00A0;<br />4*pi)&#x00A0;=&#x00A0;1e6&#x00A0;cm/s&#x00A0;*&#x00A0;sqrt(4*pi&#x00A0;*&#x00A0;2e-7&#x00A0;g/cm^3)&#x00A0;~&#x00A0;1585.33&#x00A0;Gauss.</div>
<!--l. 1185--><p class="nopar" >
                                                                                            
                                                                                            
<div class="verbatim" id="verbatim-6">
Notice&#x00A0;that&#x00A0;for&#x00A0;mu0=1&#x00A0;and&#x00A0;rho0=1&#x00A0;you&#x00A0;simply&#x00A0;have&#x00A0;vA=B.&#x00A0;Then&#x00A0;you&#x00A0;can
&#x00A0;<br />conveniently&#x00A0;set&#x00A0;the&#x00A0;field&#x00A0;strength&#x00A0;by&#x00A0;your&#x00A0;choice&#x00A0;of&#x00A0;plasma&#x00A0;beta&#x00A0;(=
&#x00A0;<br />2*cs^2/vA^2).&#x00A0;There&#8217;s&#x00A0;a&#x00A0;reason&#x00A0;why&#x00A0;we&#x00A0;like&#x00A0;dimensionless&#x00A0;quantities!</div>
<!--l. 1191--><p class="nopar" >
<!--l. 1195--><p class="noindent" >
<h4 class="subsectionHead"><span class="titlemark">1.5    </span> <a 
 id="x1-380001.5"></a>Visualization</h4>
<!--l. 1198--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.5.1    </span> <a 
 id="x1-390001.5.1"></a>&#8216;<span 
class="cmtt-12">start.pro</span>&#8217; doesn&#8217;t work:</h5>
<!--l. 1199--><p class="noindent" >
<div class="fancyvrb" id="fancyvrb30"><a 
 id="x1-39002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Reading</span><span 
class="cmtt-12">&#x00A0;grid.dat..</span><br class="fancyvrb" /><a 
 id="x1-39004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Reading</span><span 
class="cmtt-12">&#x00A0;param.nml..</span><br class="fancyvrb" /><a 
 id="x1-39006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;\%</span><span 
class="cmtt-12">&#x00A0;Expression</span><span 
class="cmtt-12">&#x00A0;must</span><span 
class="cmtt-12">&#x00A0;be</span><span 
class="cmtt-12">&#x00A0;a</span><span 
class="cmtt-12">&#x00A0;structure</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;this</span><span 
class="cmtt-12">&#x00A0;context:</span><span 
class="cmtt-12">&#x00A0;PAR.</span>
<br class="fancyvrb" /><a 
 id="x1-39008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;\%</span><span 
class="cmtt-12">&#x00A0;Execution</span><span 
class="cmtt-12">&#x00A0;halted</span><span 
class="cmtt-12">&#x00A0;at:</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;\$MAIN\$</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;104</span>
<br class="fancyvrb" /><a 
 id="x1-39010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/home/brandenb/pencil-code/runs/forced/hel1/../../../idl/start.pro</span></div>
<!--l. 1209--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: You don&#8217;t have the subdirectory &#8216;<span 
class="cmtt-12">data</span>&#8217;<a 
 id="x1-39011"></a> in your IDL variable <span 
class="pncro7t-x-x-120">!path</span><a 
 id="x1-39012"></a>. Make sure you source
&#8216;<span 
class="cmtt-12">sourceme.csh</span>&#8217;<a 
 id="x1-39013"></a>/&#8216;<span 
class="cmtt-12">sourceme.sh</span>&#8217;<a 
 id="x1-39014"></a> or set a sufficient IDL path otherwise.
<!--l. 1214--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.5.2    </span> <a 
 id="x1-400001.5.2"></a>&#8216;<span 
class="cmtt-12">start.pro</span>&#8217; doesn&#8217;t work:</h5>
<!--l. 1216--><p class="noindent" >Isn&#8217;t there some clever (or even trivial) way that one can avoid the annoying error
messages that one gets, when running e.g. &#8221;.r rall&#8221; after a new variable has been
introduced in &#8221;idl/varcontent.pro&#8221;? Ever so often there&#8217;s a new variable that can&#8217;t be
found in my param2.nml &#8211; this time it was IECR, IGG, and ILNTT that I had to
circumvent&hellip;
<a 
 id="x1-40001"></a>
<!--l. 1225--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The simplest solution is to invoke &#8216;<span 
class="cmtt-12">NOERASE</span>&#8217;<a 
 id="x1-40002"></a>, i.e.&#x00A0;say
<div class="fancyvrb" id="fancyvrb31"><a 
 id="x1-40004r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;touch</span><span 
class="cmtt-12">&#x00A0;NOERASE</span><br class="fancyvrb" /><a 
 id="x1-40006r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;start.csh</span></div>
<!--l. 1230--><p class="noindent" >or, alternatively, <span 
class="cmtt-12">start_run.csh</span>. What it does is that it reruns <span 
class="cmtt-12">src/start.x </span>with a new version
of the code; this then produces all the necessary auxiliary files, but it doesn&#8217;t overwrite or erase
the &#8216;<span 
class="cmtt-12">var.dat</span>&#8217;<a 
 id="x1-40007"></a> and other &#8216;<span 
class="cmtt-12">VAR</span>&#8217;<a 
 id="x1-40008"></a> and &#8216;<span 
class="cmtt-12">slice</span>&#8217;<a 
 id="x1-40009"></a> files.
                                                                                            
                                                                                            
<!--l. 1236--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.5.3    </span> <a 
 id="x1-410001.5.3"></a>Something about tag name undefined:</h5>
<!--l. 1238--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: In one of my older run directories I can&#8217;t read the data with idl anymore. What should I do?
Is says something like
<div class="fancyvrb" id="fancyvrb32"><a 
 id="x1-41002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Reading</span><span 
class="cmtt-12">&#x00A0;param.nml..</span><br class="fancyvrb" /><a 
 id="x1-41004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;%</span><span 
class="cmtt-12">&#x00A0;Tag</span><span 
class="cmtt-12">&#x00A0;name</span><span 
class="cmtt-12">&#x00A0;LEQUIDIST</span><span 
class="cmtt-12">&#x00A0;is</span><span 
class="cmtt-12">&#x00A0;undefined</span><span 
class="cmtt-12">&#x00A0;for</span><span 
class="cmtt-12">&#x00A0;structure</span><span 
class="cmtt-12">&#x00A0;&#x003C;Anonymous&#x003E;.</span>
<br class="fancyvrb" /><a 
 id="x1-41006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;%</span><span 
class="cmtt-12">&#x00A0;Execution</span><span 
class="cmtt-12">&#x00A0;halted</span><span 
class="cmtt-12">&#x00A0;at:</span><span 
class="cmtt-12">&#x00A0;$MAIN$</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;182</span><br class="fancyvrb" /><a 
 id="x1-41008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;/people/disk2/brandenb/pencil-code/idl/start.pro</span></div>
<!--l. 1250--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Go into &#8216;<span 
class="cmtt-12">data/param.nml</span>&#8217;<a 
 id="x1-41009"></a> and add <span 
class="cmtt-12">, LEQUIDIST=T </span>anywhere in the file (but before the last
slash).
<!--l. 1253--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.5.4    </span> <a 
 id="x1-420001.5.4"></a>Something INC in start.pro</h5>
<!--l. 1255--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: start doesn&#8217;t even work:
<div class="fancyvrb" id="fancyvrb33"><a 
 id="x1-42002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;%</span><span 
class="cmtt-12">&#x00A0;Compiled</span><span 
class="cmtt-12">&#x00A0;module:</span><span 
class="cmtt-12">&#x00A0;$MAIN$.</span><br class="fancyvrb" /><a 
 id="x1-42004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;nname=</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;11</span><br class="fancyvrb" /><a 
 id="x1-42006r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Reading</span><span 
class="cmtt-12">&#x00A0;grid.dat..</span><br class="fancyvrb" /><a 
 id="x1-42008r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Reading</span><span 
class="cmtt-12">&#x00A0;param.nml..</span>
<br class="fancyvrb" /><a 
 id="x1-42010r5"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;Can&#8217;t</span><span 
class="cmtt-12">&#x00A0;locate</span><span 
class="cmtt-12">&#x00A0;Namelist.pm</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;INC</span><span 
class="cmtt-12">&#x00A0;(INC</span><span 
class="cmtt-12">&#x00A0;contains:</span><span 
class="cmtt-12">&#x00A0;/etc/perl</span><span 
class="cmtt-12">&#x00A0;/usr/local/lib/perl/5.8.4</span><span 
class="cmtt-12">&#x00A0;/usr/local/share/perl/5.8.4</span><span 
class="cmtt-12">&#x00A0;/usr/lib/perl5</span><span 
class="cmtt-12">&#x00A0;/usr/share/perl5</span><span 
class="cmtt-12">&#x00A0;/usr/lib/perl/5.8</span><span 
class="cmtt-12">&#x00A0;/usr/share/perl/5.8</span><span 
class="cmtt-12">&#x00A0;/usr/local/lib/site_perl</span><span 
class="cmtt-12">&#x00A0;.)</span><span 
class="cmtt-12">&#x00A0;at</span><span 
class="cmtt-12">&#x00A0;/home/brandenb/pencil-code/bin/nl2idl</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;49.</span>
<br class="fancyvrb" /><a 
 id="x1-42012r6"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;BEGIN</span><span 
class="cmtt-12">&#x00A0;failed--compilation</span><span 
class="cmtt-12">&#x00A0;aborted</span><span 
class="cmtt-12">&#x00A0;at</span><span 
class="cmtt-12">&#x00A0;/home/brandenb/pencil-code/bin/nl2idl</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;49.</span></div>
<!--l. 1268--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Go into &#8216;<span 
class="cmtt-12">$PENCIL_HOME</span>&#8217;<a 
 id="x1-42013"></a> and say <span 
class="cmtt-12">svn up sourceme.csh</span><a 
 id="x1-42014"></a> and/or <span 
class="cmtt-12">svn up sourceme.sh</span><a 
 id="x1-42015"></a>. (They
were just out of date.)
<!--l. 1272--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.5.5    </span> <a 
 id="x1-430001.5.5"></a>nl2idl problem when reading param2.nml</h5>
<!--l. 1274--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: Does anybody encounter a backward problem with nl2idl? The file param*.nml files are
checked in under &#8216;<span 
class="cmtt-12">pencil-code/axel/couette/SStrat128a_mu0.20_g2</span>&#8217;<a 
 id="x1-43001"></a> and the problem is
below.
<div class="fancyvrb" id="fancyvrb34"><a 
 id="x1-43003r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;at</span><span 
class="cmtt-12">&#x00A0;/people/disk2/brandenb/pencil-code/bin/nl2idl</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;120</span>
<br class="fancyvrb" /><a 
 id="x1-43005r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;HCOND0=</span><span 
class="cmtt-12">&#x00A0;0.0,HCOND1=</span><span 
class="cmtt-12">&#x00A0;1.000000,HCOND2=</span><span 
class="cmtt-12">&#x00A0;1.000000,WIDTHSS=</span><span 
class="cmtt-12">&#x00A0;1.192093E-06,MPOLY0=</span>
<br class="fancyvrb" /><a 
 id="x1-43007r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;^------</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;HERE</span><br class="fancyvrb" /><a 
 id="x1-43009r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;at</span><span 
class="cmtt-12">&#x00A0;/people/disk2/brandenb/pencil-code/bin/nl2idl</span><span 
class="cmtt-12">&#x00A0;line</span><span 
class="cmtt-12">&#x00A0;120</span></div>
<!--l. 1288--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: The problem is the stupid ifc compiler writing the following into the namelist
file:
<div class="fancyvrb" id="fancyvrb35"><a 
 id="x1-43011r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;COOLING_PROFILE=&#8217;gaussian</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#8217;,COOLTYPE=&#8217;Temp</span>
<br class="fancyvrb" /><a 
 id="x1-43013r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#8217;COOL=</span><span 
class="cmtt-12">&#x00A0;0.0,CS2COOL=</span><span 
class="cmtt-12">&#x00A0;0.0,RCOOL=</span><span 
class="cmtt-12">&#x00A0;1.000000,WCOOL=</span><span 
class="cmtt-12">&#x00A0;0.1000000,FBOT=</span><span 
class="cmtt-12">&#x00A0;0.0,CHI_T=</span><span 
class="cmtt-12">&#x00A0;0.0</span></div>
<!--l. 1295--><p class="noindent" >If you add a comma after the closing quote:
<div class="fancyvrb" id="fancyvrb36"><a 
 id="x1-43015r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;COOLING_PROFILE=&#8217;gaussian</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#8217;,COOLTYPE=&#8217;Temp</span>
<br class="fancyvrb" /><a 
 id="x1-43017r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;&#8217;,COOL=</span><span 
class="cmtt-12">&#x00A0;0.0,CS2COOL=</span><span 
class="cmtt-12">&#x00A0;0.0,RCOOL=</span><span 
class="cmtt-12">&#x00A0;1.000000,WCOOL=</span><span 
class="cmtt-12">&#x00A0;0.1000000,FBOT=</span><span 
class="cmtt-12">&#x00A0;0.0,CHI_T=</span><span 
class="cmtt-12">&#x00A0;0.0</span></div>
                                                                                            
                                                                                            
<!--l. 1300--><p class="noindent" >things will work.
<!--l. 1302--><p class="noindent" >Note that ifc cannot even itself read what it is writing here, so if this happened to occur in
param.nml, the code would require manual intervention after each start.csh.
<!--l. 1307--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.5.6    </span> <a 
 id="x1-440001.5.6"></a>Spurious dots in the time series file</h5>
<!--l. 1309--><p class="noindent" ><span 
class="pncb7t-x-x-120">Q</span>: Wolfgang, you explained it to me once, but I forget. How can one remove spurious dots after
the timestep number if the time format overflows?
<!--l. 1316--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: I don&#8217;t know whether it exists anywhere, but it&#8217;s easy. In Perl you&#8217;d say
<div class="fancyvrb" id="fancyvrb37"><a 
 id="x1-44002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;perl</span><span 
class="cmtt-12">&#x00A0;-pe</span><span 
class="cmtt-12">&#x00A0;&#8217;s/^(\s*[-0-9]+)\.([-0-9eEdD])/$1</span><span 
class="cmtt-12">&#x00A0;$2/g&#8217;</span></div>
<!--l. 1322--><p class="noindent" >and in sed (but that&#8217;s harder to read)
<div class="fancyvrb" id="fancyvrb38"><a 
 id="x1-44004r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;sed</span><span 
class="cmtt-12">&#x00A0;&#8217;s/^\(</span><span 
class="cmtt-12">&#x00A0;*[-0-9]\+\)\.\([-0-9eEdD]\)/\1</span><span 
class="cmtt-12">&#x00A0;\2/g&#8217;</span></div>
<!--l. 1330--><p class="noindent" >
<h4 class="subsectionHead"><span class="titlemark">1.6    </span> <a 
 id="x1-450001.6"></a>General questions</h4>
<!--l. 1333--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.6.1    </span> <a 
 id="x1-460001.6.1"></a>&#8220;Installation&#8221; procedure</h5>
<!--l. 1334--><p class="noindent" >Why don&#8217;t you use GNU <span 
class="pncro7t-x-x-120">autoconf/automake</span><a 
 id="x1-46001"></a> for installation of the <span 
class="pncrc7t-x-x-120">P<span 
class="small-caps">e</span><span 
class="small-caps">n</span><span 
class="small-caps">c</span><span 
class="small-caps">i</span><span 
class="small-caps">l</span> C<span 
class="small-caps">o</span><span 
class="small-caps">d</span><span 
class="small-caps">e</span></span>?
<!--l. 1339--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: What do you mean by &#8220;installation&#8221;? Unlike the applications that normally use <span 
class="pncro7t-x-x-120">autoconf</span><a 
 id="x1-46002"></a>,
the <span 
class="pncro7t-x-x-120">Pencil Code</span><a 
 id="x1-46003"></a> is neither a binary executable, nor a library that you compile once and then
dump somewhere in the system tree. <span 
class="pncro7t-x-x-120">Autoconf</span><a 
 id="x1-46004"></a> is the right tool for these applications, but
not for numerical codes, where the typical compilation and usage pattern is very
different:
<!--l. 1346--><p class="noindent" >You have different directories with different &#8216;<span 
class="cmtt-12">Makefile.local</span>&#8217;<a 
 id="x1-46005"></a> settings, recompile after
introducing that shiny new term in your equations, etc. Moreover, you want to sometimes
switch to a different compiler (but just for that run directory) or another <span 
class="pncro7t-x-x-120">MPI</span><a 
 id="x1-46006"></a> implementation.
Our <span 
class="cmtt-12">adapt-mkfile</span><a 
 id="x1-46007"></a> approach gives you this flexibility in a reasonably convenient way, while
doing the same thing with <span 
class="pncro7t-x-x-120">autoconf</span><a 
 id="x1-46008"></a> would be using that system against most of its design
principles.
<!--l. 1354--><p class="noindent" >Besides, it would really get on my (WD&#8217;s) nerves if I had to wait two minutes for <span 
class="pncro7t-x-x-120">autoconf</span><a 
 id="x1-46009"></a> to
finish before I can start compiling (or maybe 5&#8211;10 minutes if I worked on a NEC
machine&hellip;).
<!--l. 1358--><p class="noindent" >Finally, if you have ever tried to figure out what a &#8216;<span 
class="cmtt-12">configure</span>&#8217;<a 
 id="x1-46010"></a> script does, you will appreciate a
comprehensible configuration system.
                                                                                            
                                                                                            
<!--l. 1362--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.6.2    </span> <a 
 id="x1-470001.6.2"></a>Small numbers in the code</h5>
<!--l. 1363--><p class="noindent" >What is actually the difference between epsi, tini and tiny?
<!--l. 1367--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>:
                                                                                            
                                                                                            
<div class="verbatim" id="verbatim-7">
F90&#x00A0;has&#x00A0;two&#x00A0;functions&#x00A0;epsilon()&#x00A0;and&#x00A0;tiny(),&#x00A0;with
&#x00A0;<br />
&#x00A0;<br />&#x00A0;&#x00A0;epsilon(x)&#x00A0;=&#x00A0;1.1920929e-07
&#x00A0;<br />&#x00A0;&#x00A0;tiny(x)&#x00A0;&#x00A0;&#x00A0;&#x00A0;=&#x00A0;1.1754944e-38
&#x00A0;<br />(and&#x00A0;then&#x00A0;there&#x00A0;is&#x00A0;huge(x)&#x00A0;=&#x00A0;3.4028235e+38)
&#x00A0;<br />for&#x00A0;a&#x00A0;single-precision&#x00A0;number&#x00A0;x.
&#x00A0;<br />
&#x00A0;<br />epsilon(x)&#x00A0;is&#x00A0;the&#x00A0;smallest&#x00A0;number&#x00A0;that&#x00A0;satisfies
&#x00A0;<br />&#x00A0;&#x00A0;1+epsilon(1.)&#x00A0;/=&#x00A0;1&#x00A0;,
&#x00A0;<br />while&#x00A0;tiny(x)&#x00A0;is&#x00A0;the&#x00A0;smallest&#x00A0;number&#x00A0;that&#x00A0;can&#x00A0;be&#x00A0;represented&#x00A0;without
&#x00A0;<br />precision&#x00A0;loss.
&#x00A0;<br />
&#x00A0;<br />In&#x00A0;the&#x00A0;code&#x00A0;we&#x00A0;have&#x00A0;variants&#x00A0;hereof,
&#x00A0;<br />&#x00A0;&#x00A0;&#x00A0;epsi=5*epsilon(1.0)
&#x00A0;<br />&#x00A0;&#x00A0;&#x00A0;tini=5*tiny(1.0)
&#x00A0;<br />&#x00A0;&#x00A0;&#x00A0;huge1=0.2*huge(1.0)
&#x00A0;<br />that&#x00A0;have&#x00A0;added&#x00A0;safety&#x00A0;margins,&#x00A0;so&#x00A0;we&#x00A0;don&#8217;t&#x00A0;have&#x00A0;to&#x00A0;think&#x00A0;about&#x00A0;doing
&#x00A0;<br />things&#x00A0;like&#x00A0;1/tini.
&#x00A0;<br />
&#x00A0;<br />So&#x00A0;in&#x00A0;sub.f90,
&#x00A0;<br />&#x00A0;&#x00A0;-&#x00A0;&#x00A0;&#x00A0;&#x00A0;&#x00A0;&#x00A0;evr&#x00A0;=&#x00A0;evr&#x00A0;/&#x00A0;spread(r_mn+epsi,2,3)
&#x00A0;<br />did&#x00A0;(minimally)&#x00A0;affect&#x00A0;the&#x00A0;result&#x00A0;for&#x00A0;r_mn=O(1),&#x00A0;while&#x00A0;the&#x00A0;correct&#x00A0;version
&#x00A0;<br />&#x00A0;&#x00A0;+&#x00A0;&#x00A0;&#x00A0;&#x00A0;&#x00A0;&#x00A0;evr&#x00A0;=&#x00A0;evr&#x00A0;/&#x00A0;spread(r_mn+tini,2,3)
&#x00A0;<br />only&#x00A0;avoids&#x00A0;overflow.</div>
<!--l. 1393--><p class="nopar" >
<!--l. 1395--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.6.3    </span> <a 
 id="x1-480001.6.3"></a>Why do we need a <span 
class="cmtt-12">/lphysics/ </span>namelist in the first place?</h5>
<!--l. 1397--><p class="noindent" >Wolfgang answered on 29 July 2010: &#8220;&#8216;<span 
class="cmtt-12">cdata.f90</span>&#8217;<a 
 id="x1-48001"></a> has the explanation&#8221;
<div class="fancyvrb" id="fancyvrb39"><a 
 id="x1-48003r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;!</span><span 
class="cmtt-12">&#x00A0;Constant</span><span 
class="cmtt-12">&#x00A0;&#8217;parameters&#8217;</span><span 
class="cmtt-12">&#x00A0;cannot</span><span 
class="cmtt-12">&#x00A0;occur</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;namelists,</span><span 
class="cmtt-12">&#x00A0;so</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;order</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;get</span><span 
class="cmtt-12">&#x00A0;the</span>
<br class="fancyvrb" /><a 
 id="x1-48005r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;!</span><span 
class="cmtt-12">&#x00A0;now</span><span 
class="cmtt-12">&#x00A0;constant</span><span 
class="cmtt-12">&#x00A0;module</span><span 
class="cmtt-12">&#x00A0;logicals</span><span 
class="cmtt-12">&#x00A0;into</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;lphysics</span><span 
class="cmtt-12">&#x00A0;name</span><span 
class="cmtt-12">&#x00A0;list...</span>
<br class="fancyvrb" /><a 
 id="x1-48007r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;!</span><span 
class="cmtt-12">&#x00A0;We</span><span 
class="cmtt-12">&#x00A0;have</span><span 
class="cmtt-12">&#x00A0;some</span><span 
class="cmtt-12">&#x00A0;proxies</span><span 
class="cmtt-12">&#x00A0;that</span><span 
class="cmtt-12">&#x00A0;are</span><span 
class="cmtt-12">&#x00A0;used</span><span 
class="cmtt-12">&#x00A0;to</span><span 
class="cmtt-12">&#x00A0;initialize</span><span 
class="cmtt-12">&#x00A0;private</span><span 
class="cmtt-12">&#x00A0;local</span><span 
class="cmtt-12">&#x00A0;variables</span>
<br class="fancyvrb" /><a 
 id="x1-48009r4"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;!</span><span 
class="cmtt-12">&#x00A0;called</span><span 
class="cmtt-12">&#x00A0;lhydro</span><span 
class="cmtt-12">&#x00A0;etc,</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;lphysics</span><span 
class="cmtt-12">&#x00A0;namelist!</span></div>
<!--l. 1404--><p class="noindent" >So the situation is this: we want to write parameters like ldensity to param.nml
so IDL (and potentially octave, python, etc.) can know whether density was on or
not. To avoid confusion, we want them to have exactly their original names. But
we cannot assemble the original ldensity etc. constants in a namelist, so we have
to define a local ldensity variable. And to provide it with the value of the original
cdata.ldensity, we need to transfer the value via <span 
class="pncro7t-x-x-120">ldensity&#x02D9;var</span><a 
 id="x1-48010"></a>. That&#8217;s pretty scary, although
it seems to work fine. I can track the code back to the big <span 
class="pncro7t-x-x-120">eos&#x02D9;merger</span><a 
 id="x1-48011"></a> commit, so
it may originate from that branch. One obvious problem is that you have to add
                                                                                            
                                                                                            
code in a number of places (the <span 
class="cmtt-12">ldensity</span><a 
 id="x1-48012"></a> <span 
class="cmsy-10x-x-120">&rarr; </span><span 
class="pncro7t-x-x-120">ldensity&#x02D9;var</span><a 
 id="x1-48013"></a> assignment and the local
definition of ldensity) to really get what you need. And when adding a new boolean of
that sort to &#8216;<span 
class="cmtt-12">cdata.f90</span>&#8217;<a 
 id="x1-48014"></a>, you may not even have a clue that you need all the other
voodoo.
<!--l. 1421--><p class="noindent" >There may be a cleaner solution involving generated code. Maybe something like
<div class="fancyvrb" id="fancyvrb40"><a 
 id="x1-48016r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;logical</span><span 
class="cmtt-12">&#x00A0;::</span><span 
class="cmtt-12">&#x00A0;ldensity</span><span 
class="cmtt-12">&#x00A0;!</span><span 
class="cmtt-12">&#x00A0;INCLUDE_IN_LPHYSICS</span></div>
<!--l. 1425--><p class="noindent" >could later generate code (in some param&#x02D9;io&#x02D9;extra.inc file) that looks like this:
<div class="fancyvrb" id="fancyvrb41"><a 
 id="x1-48018r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;write(unit,</span><span 
class="cmtt-12">&#x00A0;*)</span><span 
class="cmtt-12">&#x00A0;&#8217;ldensity</span><span 
class="cmtt-12">&#x00A0;=</span><span 
class="cmtt-12">&#x00A0;&#8217;,</span><span 
class="cmtt-12">&#x00A0;ldensity</span></div>
<!--l. 1430--><p class="noindent" >i.e.&#x00A0;we can manually write in namelist format. But maybe there are even simpler
solutions?
<!--l. 1433--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.6.4    </span> <a 
 id="x1-490001.6.4"></a>Can I run the code on a Mac?</h5>
<!--l. 1437--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: Macs work well for Linux stuff, except that the file structure is slightly different. Problems
when following Linux installs can usually be traced to the PATH. For general reference, if you
need to set an environment variable for an entire OS-X login session, google environment.plist.
That won&#8217;t be needed here.
<!--l. 1443--><p class="noindent" >For a Mac install, the following should work:
      <ul class="itemize1">
      <li class="itemize">Install Dev Tools (an optional install on the MacOS install disks). Unfortunately,
      last time I checked the svn version that comes with DevTools is obsolete. So:
      </li>
      <li class="itemize">Install   MacPorts   (download   from   web).   Note   that   MacPorts   installs   to   a
      non-standard location, and will need to be sourced. The installation normally drops
      an appropriate line in .profile. If it does so, make sure that that line gets sourced.
      Otherwise
      <!--l. 1456--><p class="noindent" >
      <div class="fancyvrb" id="fancyvrb42"><a 
 id="x1-49002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;export</span><span 
class="cmtt-12">&#x00A0;PATH=/opt/local/bin:/opt/local/sbin:$PATH</span><br class="fancyvrb" /><a 
 id="x1-49004r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;export</span><span 
class="cmtt-12">&#x00A0;MANPATH=/opt/local/share/man:$MANPATH</span></div>
      </li>
      <li class="itemize">Install g95 (download from web). Make sure it is linked in /bin.
      </li>
      <li class="itemize">execute macports svn install
      </li>
      <li class="itemize">download the pencil-code and enjoy.</li></ul>
<!--l. 1469--><p class="noindent" >Note: the above way to get svn works. It takes a while however, so there are certainly faster
ways out there. If you already have a non-obsolete svn version, use that instead.
                                                                                            
                                                                                            
<!--l. 1473--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.6.5    </span> <a 
 id="x1-500001.6.5"></a>Pencil Code discussion forum</h5>
<!--l. 1475--><p class="noindent" >Do I just need to send an email somewhere to subscribe or what?
<!--l. 1479--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>&#8221; The answer is yes; just go to:
<div class="fancyvrb" id="fancyvrb43"><a 
 id="x1-50002r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;http://groups.google.com/group/pencil-code-discuss</span></div>
<!--l. 1484--><p class="noindent" >
<h5 class="subsubsectionHead"><span class="titlemark">1.6.6    </span> <a 
 id="x1-510001.6.6"></a>The manual</h5>
<a 
 id="x1-51001"></a>
<!--l. 1487--><p class="noindent" >It would be a good idea to add this useful information in the manual, no?
<!--l. 1491--><p class="noindent" ><span 
class="pncb7t-x-x-120">A</span>: When you have added new stuff to the code, don&#8217;t forget to mention this in the
&#8216;<span 
class="cmtt-12">pencil-code/doc/manual.tex</span>&#8217;<a 
 id="x1-51002"></a> file.
<!--l. 1494--><p class="noindent" >Again, the answer is yes; just go to:
<div class="fancyvrb" id="fancyvrb44"><a 
 id="x1-51004r1"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;cd</span><span 
class="cmtt-12">&#x00A0;pencil-code/doc/</span><br class="fancyvrb" /><a 
 id="x1-51006r2"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;vi</span><span 
class="cmtt-12">&#x00A0;manual.tex</span><br class="fancyvrb" /><a 
 id="x1-51008r3"></a><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;</span><span 
class="cmtt-12">&#x00A0;svn</span><span 
class="cmtt-12">&#x00A0;ci</span><span 
class="cmtt-12">&#x00A0;-m</span><span 
class="cmtt-12">&#x00A0;"explanations</span><span 
class="cmtt-12">&#x00A0;about</span><span 
class="cmtt-12">&#x00A0;a</span><span 
class="cmtt-12">&#x00A0;new</span><span 
class="cmtt-12">&#x00A0;module</span><span 
class="cmtt-12">&#x00A0;in</span><span 
class="cmtt-12">&#x00A0;the</span><span 
class="cmtt-12">&#x00A0;code"</span></div>
                                                                                            
</div>
<?
	include "inc/footer.inc";
 ?>
