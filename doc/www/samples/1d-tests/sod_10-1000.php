<!-- $Id$ -->
<?
	 include "inc/header.inc";
 ?>
<div class="centcolumnpad">
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
</div>
<?
	include "inc/footer.inc";
 ?>

