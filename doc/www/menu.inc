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

	$own = parse_url ($_SERVER['SCRIPT_NAME'], PHP_URL_PATH);
	if ($own == "/index.php") { $own = "/"; }
?>
<div class="mmenu" style="position:fixed;">
<ul>
<?
	foreach ($menu as $text => $link) {
 ?>
	<li><a href="<? echo $link; ?>"<? if ($own == $link) { echo ' class="mmenuselect"'; } ?>><? echo $text; ?></a></li>
<?
	}
 ?>
</ul>
<div class="clear"></div>
</div>

