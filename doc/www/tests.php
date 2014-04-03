<!-- $Id$ -->
<?php
	 include "inc/header.inc";
 ?>
<div class="centcolumnpad">
<h2>Automatic tests</h2>

<p>To ensure reproducability, the
<a href="http://pencil-code.nordita.org/">Pencil Code</a>
is tested daily for a number of
sample applications. This is important for us in order to make sure certain
improvements in some parts of the code do not affect the functionality
of other parts. For other users who suspect that a new problem has emerged
it could be useful to first see whether this problem also shows up in
our own tests. The latest test results for a can be seen online:</p>

<ul STYLE="font-size:13px;">
<!--
<li><a href="http://norlx50.albanova.se/~wdobler/pencil-code/tests/g95_debug.html">norlx50 (Linux, g95, by Wolfgang)</a>
-->
<!-- Anders Johansen test -->
<li><a href="http://www.astro.ku.dk/~ajohan/pencil-test.html">gfortran with open-mpi, quad-core with Intel Xeon 2.40 GHz cores, by Anders Johansen)</a>
<!-- Philippe Bourdin test -->
<li><a href="http://www.pab-software.de/Pencil/pc_auto-test.txt">GNU Fortran (Ubuntu 4.4.3-4ubuntu5.1) 4.4.3-4 (by Philippe Bourdin)</a>
<!-- Boris Dintrans test -->
<li><a href="http://www.ast.obs-mip.fr/users/dintrans/tmp/test_runs.html">Copernic (Linux/CentOS5 on 4 x Hexacore Intel Xeon E7540@2.00GHz, ifort 64 bits v12.0.1.107, by Boris Dintrans, regular level 2 test)</a>
<!-- Boris Dintrans test (2) -->
<li><a href="http://www.ast.obs-mip.fr/users/dintrans/tmp/our_tests.html">Copernic (Linux/CentOS5 on 4 x Hexacore Intel Xeon E7540@2.00GHz, ifort 64 bits v12.0.1.107, by Boris Dintrans, 16 separate tests)</a>
<!-- Weekly (big) test -->
<li><a href="http://www.svenbingert.de/auto-test.html">Linux/Ubuntu13.4 on Intel Core 2 Quad Q9000@2.00GHz, ifort 64bit v13.0.1 (Sven Bingert, standard + personal tests)</a>
<li><a href="http://norlx51.albanova.se/~brandenb/pencil-code/tests/g95_debug.html">Nordita Big Test (norlx51, gfortran, openmpi, by Wolfgang/Axel)</a>
<!-- Weekly (big) test -->
<li><a href="http://norlx51.albanova.se/~brandenb/pencil-code/tests/gfortran_hourly.html">Nordita Hourly Test (norlx51, gfortran, openmpi, by Wolfgang/Axel)</a>
<!--
<li><a href="http://www.nordita.org/~brandenb/pencil-code/normac.html">Nordita Mac Mini (os10, g95, lammpi, by Axel)</a>
[<a href="http://www.nordita.org/~brandenb/pencil-code/normac_previous.html">previous</a>]
<li><a href="http://www.nordita.org/~brandenb/pencil-code/nor52.html">Nordita PowerMac (os10, g95, ompi, by Axel)</a>
[<a href="http://www.nordita.org/~brandenb/pencil-code/nor52_previous.html">previous</a>]
-->
<!-- Nordita PowerMac -->
<li><a href="http://norlx51.albanova.se/~brandenb/pencil-code/tests/nor52.html">Nordita PowerMac (os10, g95, ompi, by Axel)</a>
[<a href="http://norlx51.albanova.se/~brandenb/pencil-code/tests/nor52_previous.html">previous</a>]
<!--
<li><a href="http://www.nordita.org/~brandenb/pencil-code/qmul.html">Queen Mary Cluster (London, by Axel)</a>
[<a href="http://www.nordita.org/~brandenb/pencil-code/qmul_previous.html">previous</a>]
<li><a href="http://www.nordita.org/~brandenb/pencil-code/fend.html">DCSC cluster in Copenhagen (pgf90 -fastsse -tp k8-64e, by Axel)</a>
[<a href="http://www.nordita.org/~brandenb/pencil-code/fend_previous.html">previous</a>]
<li><a href="http://www.capca.ucalgary.ca/~theine/pencil-code/dcsc.html">Horseshoe (Linux cluster, ifc 6.0 compiler, by Tobi)</a>
<li><a href="http://www.astro.ku.dk/~ajohan/pencil-code/dcsc.html">Horseshoe (Linux cluster, ifort compiler, by Anders)</a>
<li><a href="http://www.astro.ku.dk/~ajohan/pencil-code/aida25.html">aida25 (Linux cluster, ifort compiler with MPICH, by Anders)</a>
<li><a href="http://bohr.phys.ntnu.no/~nilshau/pencil-code/gridur.html">Gridur (SGI machine in Trondheim, by Nils)</a>
<li><a href="http://www.tac.dk/~brandenb/pencil-code/tac.html">tacsg2 (SGI machine, always some problems...)</a>
-->
</ul>

<p>Note: before checking in your own changes, you should at least do
the very minimal auto-test:</p>

<div class="codescroll"><code>pc_auto-test --level=0 --no-pencil-check -C</code></div>
</div>

<a name="samples"></a>
<div class="centerdivider"></div>
<div class="centcolumnpad">
<h2>Results from tests</h2>
We have prepared a quick overview on the results of some of our automatic tests:
<?php
	include "inc/samples_overview.inc";
?>
</div>
<?php
	include "inc/footer.inc";
 ?>

