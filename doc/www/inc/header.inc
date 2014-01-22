<?php
	# $Id$

	date_default_timezone_set ("Europe/Stockholm");
	require "utils.inc";

	if (!isset ($right_column)) { $right_column = ""; }
 ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html lang=EN>
<head>
	<meta name="language" content="en">
	<meta http-equiv="content-type" content="text/html; charset=UTF-8" >
	<title>The Pencil Code</title>
	<meta name="Description" content="The Pencil Code is a high-order finite-difference code for compressible fluid dynamics (CFD) with magnetic fields (MHD).">
	<meta name="Keywords" content="Pencil, Code, MHD">
	<meta name="Copyright" content="&copy; Copyright 2001-<?php echo strftime ('%Y'); ?> by the Pencil Code development team. The source code and all content of this website are unter GPL.">
	<meta name="author" content="Philippe-A. Bourdin, Illa Rivero Losada, and others">
	<meta name="Robots" content="FOLLOW, All, Index">
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
<?php
	include "inc/menu.inc";
 ?>
</div>

<!-- right column stuff -->
<?php
	if ($right_column) { include $right_column.".inc"; }
 ?>

<div class="centcolumn">

