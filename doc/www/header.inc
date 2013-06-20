<?
	if (!isset ($right_column)) { $right_column = ""; }
 ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
   "http://www.w3.org/TR/html4/loose.dtd">

<!-- $Id: index.html,v 1.1 2013-06-18 21:36:54 illa Exp $ -->

<!-- <html xmlns="http://www.w3.org/1999/xhtml"> -->
<html>
<head>
  <meta name="language" content="en">
  <meta http-equiv="content-type" content="text/html; charset=ISO-8859-1" >
  <title>The Pencil Code</title>
  <link rel="shortcut icon" href="pics/pencils_16x16.png">
  <!link rel="stylesheet" href="pencil.css" type="text/css">
  <link rel="stylesheet" href="styles.css" type="text/css">
</head>



<body>
<!-- Logo and Title -->
<!--<h1>
  <img src="pics/pencils_65x30.png" alt="Pencilqqq logo" />
  &nbsp;&nbsp;&nbsp;
  The Pencil Code
  &nbsp;&nbsp;&nbsp;
  <img src="pics/pencils_65x30.png">
</h1>
-->

<div id="bodywatermark">
<div id="mcontainerouter">
<!--<table cellpadding="0" cellspacing="10" border="0">
<tr><td align="center" valign="top" width="210">
<!-- Navigation Bar -->

<!--</td>-->
<!--<td align="left" valign="top">-->
<div id="mcontainer">

        <div id="header"><!--<div id="headerwatermark">-->
        
                <div id="headerpad">
                						
                        <h1>
                        <img src="pics/pencils_65x30.png" alt="Pencil logo" />
  															&nbsp;&nbsp;&nbsp;
  																The Pencil Code
  																&nbsp;&nbsp;&nbsp;
  																<img src="pics/pencils_65x30.png">
                        </h1>

                </div><!--End of Header Watermark -->
        </div><!--End of Header -->


<!-- left column stuff-->

<div id="leftcolumn">
         <? include "menu.inc"; ?>      
</div>


  
<!-- right column stuff-->  

<?
	if ($right_column) { include $right_column.".inc"; }
 ?>
