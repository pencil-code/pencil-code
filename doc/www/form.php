<!-- $Id$ -->
<?
	include "inc/header.inc";
 ?>
<div class="centcolumnpad">
<h1>
    <img src="/pics/pencils_65x30.png">
    Obtaining the <em>Pencil Code</em></h1>

<!-- Form -->
<form action="/cgi-bin/people/brandenb/form-eval.cgi" method="get">

<table width="90%" cellspacing=3 cellpadding=15 align=center
  summary="Table for layout (wider margins)">
<tr><td>

<p>
<tr><td>
    To get an idea about the user base of the Pencil Code, we need to
    know a few things about you.
    Please take the time and answer the following questions (questions
    with <font color=red>*</font> are mandatory).

<!-- Identification -->
    <p>
    <table bgcolor="#d8d8f0" width="90%" border=1 cellpadding=20 align=center
      summary="Table for layout (background color)">
    <tr><td>
	<table summary="Table for layout (arrange input fields)">
	  <tr><td>Name: <font color=red>*</font></td>
	      <td><input name="name" size=40 maxlength=60></td>
	  </tr>

	  <tr><td>Email: <font color=red>*</font></td>
	      <td><input name="email" size=40 maxlength=60></td>
	  </tr>

	  <tr><td>Country: <font color=red>*</font></td>
	      <td><input name="country" size=40 maxlength=60></td>
	  </tr>

	  <tr><td>Institution:</td>
	      <td><textarea name="institute" rows=3 cols=35
		  wrap=virtual></textarea></td>
	  </tr>

<!--  	  <tr><td>Type of usage: <font color=red>*</font></td> -->
<!--                <td> -->
<!--  		  <select size=1 name="usage_type"> -->
<!--  		    <option                  > Academic -->
<!--  		    <option value="Testing"  > Just testing -->
<!--  		    <option value="Upgrading"> Upgrading from older version -->
<!--  		    <option selected         > Other -->
<!--  		  </select> -->
<!--  	      </td> -->
<!--  	  </tr> -->
        </table>

<!-- Research interests -->
    <tr><td>Research interests:
<!--  	<td> -->
<!--  	<select name="res_interest" multiple> -->
<!--  	  <option>                   Turbulence -->
<!--  	  <option>                   MHD -->
<!--  	  <option value="Dynamos">   Dynamo theory -->
<!--  	  <option>                   Convection -->
<!--  	  <option value="Accretion"> Accretion discs -->
<!--  	  <option>                   Other -->
<!--  	</select> -->
        <table cellspacing=10
	  summary="Table for layout (three-column arrangement)">
	  <tr>
	    <td><input type=checkbox name="research"
                       value="Turbulence">Turbulence</td>
	    <td><input type=checkbox name="research"
	               value="MHD">MHD</td>
	    <td><input type=checkbox name="research"
	               value="Dynamos">Dynamo theory</td>
	  </tr>
	  <tr>
	    <td><input type=checkbox name="research"
                       value="Convection">Convection</td>
	    <td><input type=checkbox name="research"
		       value="Accretion">Accretion discs</td>
	    <td>Other: <input name="research_other" size=20 maxlength=60></td>
	  </tr>
        </table>
      </td>
    </tr>
    <tr>
      <td>
	<table cellspacing=10
	  summary="Table for layout (three-column arrangement)">
	  <tr><td>Platform you plan to run the code on:</td>
	      <td><input name="platform" size=20 maxlength=60></td>
	  </tr>
        </table>
      </td>
    </tr>

<!-- Submit/Reset button-->
    <tr><td>
	<table summary="Table for layout (arrangement of submit/reset buttons)">
	  <tr>
	    <td width="20%">
	    <input type=submit value="Submit">
	    </td>
	    <td width="20%">
	    <input type=reset value="Clear">
	    </td>
	  </tr>
        </table>
        </td>
    </tr>
    </table>

</tr>

</table>
<p>

</form>
</div>
<?
	include "inc/footer.inc";
 ?>

