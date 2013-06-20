<?
	$menu = array(
		"Home"           => "/",
		"News"           => "/news.php",
		"Documentation"  => "/doc.php",
		"Autotests"      => "/tests.php",
		"Download"       => "/download.php",
		"Meetings"       => "/meetings.php",
		"Latest changes" => "http://code.google.com/p/pencil-code/source/list",
		"References"     => "/references.php",
		"Contact"        => "/contact.php",
	);

	$own = "/".basename ($_SERVER['SCRIPT_NAME']);
?>
<div class="mmenu" style="position:fixed;">
<ul>
<?
	foreach ($menu as $text => $link) {
 ?>
	<li<? if ($own == $link) { echo ' CLASS="mmenuselect"'; } ?>><a href="<? echo $link; ?>"><? echo $text; ?></a></li>
<?
	}
 ?>
</ul>
<div class="clear"></div>
</div>

