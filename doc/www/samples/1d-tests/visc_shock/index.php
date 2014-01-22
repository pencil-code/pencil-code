<!-- $Id$ -->
<?php
	include "inc/header.inc";
 ?>
<div class="centcolumnpad">
<h2>Shock viscosity</h2>
<h3><tt>samples/1d-tests/sod_10s (to sod_1000s)</tt></h3>
In Makefile.local, set "VISCOSITY=visc_shock" to invoke shock viscosity.
<div align="center">
<table border="0" cellpadding="0" cellspacing="15">
  <tr>
    <th>pressure jump 10:1</th>
    <th>pressure jump 100:1</th>
    <th>pressure jump 1000:1</th>
  </tr>
  <tr>
    <td align="center"><A href="sod_10s.png"><img src="sod_10s_thumb.jpg"></a></td>
    <td align="center"><A href="sod_100s.png"><img src="sod_100s_thumb.jpg"></a></td>
    <td align="center"><A href="sod_1000s.png"><img src="sod_1000s_thumb.jpg"></a></td>
  </tr>
  <tr>
    <td>&nu;<sub>sh</sub>=4, &nu;=2 10<sup>-3</sup></td>
    <td>&nu;<sub>sh</sub>=3, &nu;=2 10<sup>-3</sup></td>
    <td>&nu;<sub>sh</sub>=3, &nu;=5 10<sup>-4</sup></td>
  </tr>
  <tr>
    <td>&chi;=5 10<sup>-4</sup>, <i>t</i>=2.4:</td>
    <td>&chi;=3 10<sup>-4</sup>, <i>t</i>=1.9:</td>
    <td>&chi;=5 10<sup>-5</sup>, <i>t</i>=1.5:</td>
  </tr>
</table>
</div>
</div>
<?php
	include "inc/footer.inc";
 ?>
