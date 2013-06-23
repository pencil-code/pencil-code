<?
	# $Id$

	require "utils.inc";

	if (!isset ($right_column)) { $right_column = ""; }
 ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html lang=EN>
<head>
	<meta name="language" content="en">
	<meta http-equiv="content-type" content="text/html; charset=UTF-8" >
	<title>The Pencil Code</title>
	<link rel="icon" type="image/png" href="/pics/pencils_32x32.png">
	<link rel="stylesheet" href="/styles.css" type="text/css">
</head>

<body>
<div id="bodywatermark">
<div id="mcontainerouter">
<div id="mcontainer">

<!-- Logo and Title -->
<div id="header">
<div id="headerpad">
<img src="/pics/pencils_100x45.png" alt="Pencil Code logo" width="100" height="45" border="0" vspace="5" hspace="30" align="left" STYLE="position:fixed">
<img src="/pics/pencils_100x45.png" alt="Pencil Code logo" width="100" height="45" border="0" vspace="5" hspace="20" align="right">
<h1>The&nbsp;Pencil&nbsp;Code</h1>
<h4>a high-order finite-difference code for compressible MHD</h4>
</div>
</div>

<!-- left column stuff -->
<div id="leftcolumn">
<!-- Navigation Bar -->
<?
	include "menu.inc";
 ?>
</div>

<!-- right column stuff -->
<?
	if ($right_column) { include $right_column.".inc"; }
 ?>

<div class="centcolumn">

