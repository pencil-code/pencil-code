<?
	 include "header.inc";
 ?>
<div class="centcolumnpad">
<h1>How to do benchmarks and timings.</h1>

<h2>pencil-code/samples/helical-MHDturb</h2>

<p>Most of the timings reported in the <a href="doc/manual.pdf">manual</a>
of the <a href="http://www.nordita.dk/software/pencil-code/">Pencil Code</a>
(currently page 68) are based on the run checked in under
<i>pencil-code/samples/helical-MHDturb</i>.</p>

<p>To run MPI, you want to edit <i>src/Makefile.local</i> and change
<i>MPICOMM=nompicomm</i> into <i>MPICOMM=mpicomm</i>.
You also want to edit the file <i>src/cparam.local</i> and change
<i>ncpus=1,nprocy=1</i> for example to
<i>ncpus=256,nprocy=4</i>, which means that you will be using 256 procs
in total and 4 in the y direction (so 64 in the z direction).
Note that the layout may be very important (see the manual).</p>

<p>Of course, on a bigger computer you also want to run a bigger mesh,
so change <i>nxgrid=32</i> for example to <i>nxgrid=512</i>, which means
512<sup>3</sup> meshpoints.</p>

<p>To run for more time steps, just change in the <i>run.in</i> file the
entries <i>nt=10, it1=2</i> for example to
<i>nt=300, it1=10</i>, which means you are running 300 time steps and
output on the console come every 10 time steps.
(If <i>ialive=1</i> is set, you can check the progress of the code
conveniently by monitoring the content of the file <i>data/proc0/alive.info</i>,
which contains the time step number on processor 0 (or any other processor
if you look under another <i>data/proc?</i> directory.)</p>

<p>The timings that are reported in the manual and shown in the figure below
(for an Origin 3000 and a Linux cluster with clock frequencies between 2.0 and
3.2 MHz) are obtained by looking at the output at the end, where it says for example
<i>Wall clock time/timestep/meshpoint [microsec] =  5.4</i>, which means it
ran 5.4 microseconds per meshpoint and per time step. (Make sure the code
doesn't produce NaNs, of course...)</p>

<p align="center"><img src="pics/ptimings.png"></p>
</div>
<?
	include "footer.inc";
 ?>

