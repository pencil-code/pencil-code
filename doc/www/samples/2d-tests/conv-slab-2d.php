<!-- $Id$ -->
<?php
	 include "inc/header.inc";
 ?>
<div class="centcolumnpad">
<h2>Results from tests</h2>
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
    <td align="center"><A href="/samples/2d-tests/conv-slab-2d.png"><img src="/samples/2d-tests/conv-slab-2d_thumb.jpg" width="75" height="130" border="0" alt="no B-field"></a></td>
    <td align="center"><A href="/samples/2d-tests/conv-slab-2d2.png"><img src="/samples/2d-tests/conv-slab-2d2_thumb.jpg" width="75" height="130" border="0" alt="with B-field"></a></td>
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
</div>
<?php
	include "inc/footer.inc";
 ?>

