<!-- $Id$ -->
<?
	 include "header.inc";
 ?>
<div class="centcolumnpad">
<h2>Automatic test results</h2>

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
<li><a href="http://www.pab-software.de/Pencil/pc_auto-test.txt">GNU Fortran (Ubuntu 4.4.1-4ubuntu9) 4.4.1 (by Philippe Bourdin)</a>
<!-- Boris Dintrans test -->
<li><a href="http://www.ast.obs-mip.fr/users/dintrans/tmp/test_runs.html">Copernic (Linux/CentOS5 on 4 x Hexacore Intel Xeon E7540@2.00GHz, ifort 64 bits v12.0.1.107, by Boris Dintrans, regular level 2 test)</a>
<!-- Boris Dintrans test (2) -->
<li><a href="http://www.ast.obs-mip.fr/users/dintrans/tmp/our_tests.html">Copernic (Linux/CentOS5 on 4 x Hexacore Intel Xeon E7540@2.00GHz, ifort 64 bits v12.0.1.107, by Boris Dintrans, 16 separate tests)</a>
<!-- Weekly (big) test -->
<li><a href="http://www.svenbingert.de/auto-test.html">Linux/Ubuntu10.4 on Intel Core 2 Quad Q9000@2.00GHz, ifort 64bit v11.1 (Sven Bingert, standard + personal tests)</a>
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

<p><code>pc_auto-test --level=0 --no-pencil-check -C</code></p>

<h2>Results from tests</h2>

<h3><tt>samples/1d-tests</tt></h3>
<ul>
<li><i>Sod shock tube tests</i> (checked in under <b>samples/1d-tests/sod_10</b>
to <b>sod_1000</b>). Initial condition is a smoothed (width=0.03) isothermal
pressure jump ranging from 10:1 to 1000:1.
<div align="center">
<table border=0 cellspacing=15>
  <tr>
    <th> pressure jump 10:1 </th>
    <th> pressure jump 100:1 </th>
    <th> pressure jump 1000:1 </th>
  </tr>
  <tr>
    <td align="center"><A href="/samples/1d-tests/sod_10.png"><img src="/samples/1d-tests/sod_10_thumb.jpg" width="100" height="80" border="0" alt="Sod shock tube 10:1"></a></td>
    <td align="center"><A href="/samples/1d-tests/sod_100.png"><img src="/samples/1d-tests/sod_100_thumb.jpg" width="100" height="80" border="0" alt="Sod shock tube 100:1"></a></td>
    <td align="center"><A href="/samples/1d-tests/sod_1000.png"><img src="/samples/1d-tests/sod_1000_thumb.jpg" width="100" height="80" border="0" alt="Sod shock tube 1000:1"></a></td>
  </tr>
  <tr>
    <td> &nu;=.02, &chi;=.0005, <i>t</i>=2.7: </td>
    <td> &nu;=.04, &chi;=.0005, <i>t</i>=1.9: </td>
    <td> &nu;=.08, &chi;=.0005, <i>t</i>=1.5: </td>
  </tr>
</table>
</div>
The values of viscosity are chosen rather conservatively;
for weak shocks one can get away with less viscosity
(<A href="/samples/1d-tests/sod_10_nu0.014.png">&nu;=.014 for 10:1</a> and
<A href="/samples/1d-tests/sod_100_nu0.028.png">&nu;=.028 for 100:1</a>).
For strong shocks (pressure jumps of 1000:1 and above) the
discrepancy compared with the inviscid analytic solution
becomes quite noticeable.</li>

<li><i>Rarefaction shocks</i> (checked in under samples/1d-tests/expans, expans_bfield and riemann_bfield).
<div align="center">
<table border=0 cellspacing=15>
  <tr>
    <th>no B-field</th>
    <th>with B-field</th>
    <th>shock with B-field</th>
  </tr>
  <tr>
    <td align="center"><A href="/samples/1d-tests/expans.png"><img src="/samples/1d-tests/expans_thumb.jpg" width="100" height="80" border="0" alt="no B-field"></a></td>
    <td align="center"><A href="/samples/1d-tests/expans_bfield.png"><img src="/samples/1d-tests/expans_bfield_thumb.jpg" width="100" height="80" border="0" alt="with B-field"></a></td>
    <td align="center"><A href="/samples/1d-tests/riemann_bfield.png"><img src="/samples/1d-tests/riemann_bfield_thumb.jpg" width="100" height="80" border="0" alt="shock with B-field"></a></td>
  </tr>
  <tr>
    <td>cf. Fig. 1 of <A href="http://xxx.lanl.gov/abs/astro-ph/0207419">Falle (2002)</a></td>
    <td>cf. Fig. 2 of <A href="http://xxx.lanl.gov/abs/astro-ph/0207419">Falle (2002)</a></td>
    <td>cf. Fig. 6 of <A href="http://xxx.lanl.gov/abs/astro-ph/0207419">Falle (2002)</a></td>
  </tr>
</table>
</div></li>

<li>Conditions for non-magnetic rarefaction shock.<br>
<i>Left state</i>: &rho;=1, p=10, ux=-3.<br>
<i>Right state</i>: &rho;=.87469, p=8, ux=-2.46537.<br>
This corresponds to s/cp=1.68805 on both sides. 800 points.<br>
&nu;=0.05, &chi;=0.0002.</li>

<li>Conditions for magnetic rarefaction shock.<br>
<i>Left state</i>: &rho;=1, p=.2327, ux=-4.6985, uy=-1.085146, Bx=-0.7, By=1.9680.<br>
<i>Right state</i>: &rho;=.7270, p=.1368, ux=-4.0577, uy=-0.8349, Bx=-0.7, By=1.355.<br>
This corresponds to s/cp=-0.5682 on both sides. 800 points.<br>
&nu;=&chi;=&eta;=0.07.</li>

<li>Conditions for magnetic Riemann problem.<br>
<i>Left state</i>: &rho;=0.5, p=10, ux=0, uy=2, Bx=2, By=2.5.<br>
<i>Right state</i>: &rho;=0.1, p=0.1, ux=-10, uy=0, Bx=2, By=2.<br>
This corresponds to s/cp=2.38119 on the left and 1.22753 on the right. 800 points.<br>
&nu;=&chi;=&eta;=2.</li>

<li><i>Note</i>: in the tests above, uniform viscosity is used. The
viscosity has to be chosen such as to cope with the strongest
compression in the domain. By using
a <a href="/samples/1d-tests/visc_shock"> nonuniform (artificial)
viscosity</a>, both compression and expansion shocks can be made as
sharp as possible.</li>
</ul>


<h3><tt>samples/conv-slab-2d</tt></h3>
Vertical cross-section at <i>t</i>=920: K=0.008, &nu;=0.004.
<ul>
<li><i>2-D convection</i> (checked in under
samples/2d-tests/conv-slab-2d and conv-slab-2d2).
<div align="center">
<table border=0 cellspacing=15>
  <tr>
    <th>no B-field</th>
    <th>with B-field</th>
  </tr>
  <tr>
    <td align="center"><A href="/samples/2d-tests/conv-slab-2d.gif"><img src="/samples/2d-tests/conv-slab-2d_thumb.jpg" width="75" height="130" border="0" alt="no B-field"></a></td>
    <td align="center"><A href="/samples/2d-tests/conv-slab-2d2.gif"><img src="/samples/2d-tests/conv-slab-2d2_thumb.jpg" width="75" height="130" border="0" alt="with B-field"></a></td>
  </tr>
  <tr>
    <td>(<i>t</i>=920)</td>
    <td>(<i>t</i>=320)</td>
  </tr>
</table>
</div>
</li>

<li><i>2-D convection</i> (checked in under samples/2d-tests/A3+chi11+Ra1e5).
<div align="center">
<table border=0 cellspacing=15>
  <tr>
    <th>aspect ratio 3, density ratio 11, Ra=10<sup>5</sup></th>
  </tr>
  <tr>
    <td align="center"><A href="/samples/2d-tests/A3+chi11+Ra1e5.png"><img src="/samples/2d-tests/A3+chi11+Ra1e5_thumb.jpg" width="265" height="90" border="0" alt="large Ra"></a></td>
  </tr>
  <tr>
    <td>(<i>t</i>=530, K=&nu;=0.0011, 150x51 points)</td>
  </tr>
</table>
</div>
</li>
</ul>


<h3><tt>samples/turbulence/helical-MHDturb32-4procs</tt></h3>
Low resolution helical MHD turbulence run, 2x2 processors, initial field: random.
<div align="center">
<table border=0 cellspacing=15>
  <tr>
    <th>t=200-1700</th>
    <th>t=200-800</th>
    <th>initial field: Beltrami</th>
  </tr>
  <tr>
  <td align="center"><A href="/samples/turbulence/helical-MHDturb32-4procs/bz_till1700.mpg"><img src="/samples/turbulence/helical-MHDturb32-4procs/bz_thumb.jpg" width="100" height="80" border="0" alt="t=200-1700"></a></td>
  <td align="center"><A href="/samples/turbulence/helical-MHDturb32-4procs/bz_till800.mpg"><img src="/samples/turbulence/helical-MHDturb32-4procs/bz_thumb.jpg" width="100" height="80" border="0" alt="t=200-800"></a></td>
  <td align="center"><A href="/samples/turbulence/rot512_Om0a/I8_0-37.mpg"><img src="/samples/turbulence/rot512_Om0a/img_0000_thumb.jpg" width="160" height="120" border="0" alt="Beltrami initial field"></a></td>
  </tr>
  <tr>
    <td align="center">full sequence(7.3Mb)</td>
    <td align="center">every 2nd frame (1.4Mb)</td>
    <td align="center">hires run: 512<sup>3</sup> (1.1Mb)<br>rot512_Om0a</td>
  </tr>
</table>
</div>

<h3><tt>samples/turbulence/hydro512f</tt></h3>
A non-helical hydro turbulence run, 512<sup>3</sup> meshpoints, 128 processors.
<div align="center">
<table border=0 cellspacing=15>
  <tr>
    <th></th>
    <th></th>
    <th></th>
  </tr>
  <tr>
  <td align="center"><A href="/samples/turbulence/hydro512f/uz_0-131.mpg"><img src="/samples/turbulence/hydro512f/uz_thumb.jpg" width="134" height="112" border="0" alt="v_z early"></a></td>
  <td align="center"><A href="/samples/turbulence/hydro512f/uz_0-148.mpg"><img src="/samples/turbulence/hydro512f/uz_thumb2.jpg" width="160" height="120" border="0" alt="v_z late"></a></td>
  <td align="center"><A href="/samples/turbulence/hydro512f/lnrho_0-148.mpg"><img src="/samples/turbulence/hydro512f/lnrho_thumb.jpg" width="144" height="108" border="0" alt="lnrho"></a></td>
  </tr>
  <tr>
    <td align="center">z-comp of velocity (4.7Mb)<br>t=500-565</td>
    <td align="center">z-comp of velocity (6.4Mb)<br>t=566-640</td>
    <td align="center">log of density: lnrho (4.1Mb)<br>t=566-640</td>
  </tr>
</table>
</div>


<h3><tt>samples/shearing-box/BH256_3D_mean_Bz=0b1</tt></h3>

Shearing box simulation, 256<sup>3</sup> meshpoints, 32 processors.
<div align="center">
<table border=0 cellspacing=15>
  <tr>
    <th></th>
  </tr>
  <tr>
  <td align="center"><A href="/samples/shearing-box/BH256_3D_mean_Bz=0b1.mpg"><img src="/samples/shearing-box/BH256_3D_mean_Bz=0b1_thumb.jpg" width="144" height="108" border="0" alt="v_z"></a></td>
  </tr>
  <tr>
    <td align="center">z-comp of velocity (2.8Mb)<br>t=296-393</td>
  </tr>
</table>
</div>


</div>
<?
	include "footer.inc";
 ?>

