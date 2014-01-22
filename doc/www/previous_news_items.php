<!-- $Id$ -->
<?php
	 include "inc/header.inc";
 ?>
<div class="centcolumnpad">
<h2>News archive</h2>

<ul>

<li class="separated">21 Sep 2008:<br>
Now the code is available under
<a href="http://pencil-code.googlecode.com/">http://pencil-code.googlecode.com/</a>.
Please consult the <a href="http://code.google.com/p/pencil-code/wiki/SvnTransition">following notes</a>
to do the transition.</li>

<li class="separated">28 July 2006:<br>
We have incorporated all the feature of the eos branch into the main trunk.
Until 28 July the public repository was not updated (to prevent
users seeing broken code), but meanwhile we have been working
on production runs ourselves without problems, so we deem the code
fit for public consumption.
There are some minor changes that will be necessary in the .in-files,
such as adding eos_init_pars (possibly with gamma=1.) and
viscosity_run_pars (with the value of nu). It will be best to check
out the samples directory for how to do it.
Alternatively, wait until Wolfgang will have prepared his script
that converts from the old files to the new ones.</li>

</ul>

</div>
<?php
	include "inc/footer.inc";
 ?>

