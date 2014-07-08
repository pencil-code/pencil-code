<!-- $Id: index.php 21445 2014-01-22 16:34:30Z Iouri.Belokopytov@gmail.com $ -->
<?php
	 include "inc/header.inc";
 ?>
<div class="centcolumnpad">
  <h4 id="toc0"><span>Supernova driven turbulence in the interstellar medium</span></h4>
  <p>Movies of supernova (SN) driven turbulence in the interstellar medium (ISM).
  The movie first from the left (in blue) depicts the density of the ISM, in
  which SN remnants appear as dark hollows, where the explosion has forced out the
  plasma leaving a very diffuse remnant interior. 
  The next movie (in red) depicts the temperature. 
  In this case the SN remants appear as light flashes where the explosions heat 
  the ISM to an excess of 10 Million degrees Kelvin. 
  Each face shows a 2D slice through the centre of the box. 
  The midplane of the galaxy is in the centre and our domain extend 
  1 kiloparsec (kpc) either side of the mid plane and covers an area 1 square kpc
  representing the radial and azimuthal directions relative to the galactic 
  centre. 
  We model conditions typical of the solar neighbourhood. (1 kpc = 30.86 thousand million million km)</p>
  <div align="center">
    <table border=0 cellspacing=15>
      <tr>
        <td align="center"><a href="http://fagent.wdfiles.com/local--files/astro/rho.avi"><img src="http://fagent.wdfiles.com/local--files/astro/rho.jpg" width="240" height="200" border="0" alt="plasma density"></a></td>
        <td align="center"><a href="http://fagent.wdfiles.com/local--files/astro/tt.avi"><img src="http://fagent.wdfiles.com/local--files/astro/tt.jpg" width="240" height="200" border="0" alt="plasma temperature"></a></td>
        <td align="center"><a href="http://fagent.wdfiles.com/local--files/astro/e_b.avi"><img src="http://fagent.wdfiles.com/local--files/astro/b2.jpg" width="240" height="200" border="0" alt="magnetic field strength"></a></td>
      </tr>
      <tr>
        <td align="center">plasma density</td>
        <td align="center">plasma temperature</td>
        <td align="center">magnetic field strength</td>
      </tr>
    </table>
  </div>
  <p>The latter movie depicts the magnetic field and shows the dynamo action
  of the ISM from a small seed field. 
  This shows the magnetic energy distribution.
  A small seed field of 0.001 micro Gauss is amplified by the turbulence and 
  differential rotation of the galaxy to about 2 micro Gauss.</p>
  <p>The model successfully simulates a galactic dynamo over a period of 1+ 
  Gyrs, yielding a significantly ordered magnetic 
  field aligned with the rotation of the galaxy.</p>
  <div align="center">
    <table border=0 cellspacing=30>
      <tr>
        <td align="center"><img src="/samples/ism/pavB.png" width="700" height="350" border="0" alt="mean B"></td>
      </tr>
    </table>
  </div>

</div>
<?php
	include "inc/footer.inc";
 ?>

