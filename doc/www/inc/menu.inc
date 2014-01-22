<?php
	$menu = array(
		"Home"                => "/",
		"News"                => "/news.php",
		"Documentation"       => "/doc.php",
		"Highlights"          => "/highlights/",
		"Autotests"           => "/tests.php",
		"Download"            => "/download.php",
		"Meetings"            => "/meetings.php",
		"References"          => "/references.php",
		"Contact"             => "/contact.php",
		"Latest changes ..."  => "http://code.google.com/p/pencil-code/source/list",
	);

	$own = parse_url ($_SERVER['SCRIPT_NAME'], PHP_URL_PATH);
	$own = preg_replace ('/\/index\.php$/i', "/", $own);
?>
<div class="mmenu">
<ul>
<?php
	foreach ($menu as $text => $link) {
 ?>
	<li><a href="<?php echo $link; ?>"<?php if ($own == $link) { echo ' class="mmenuselect"'; } ?>><?php echo $text; ?></a></li>
<?php
	}
 ?>
</ul>
<div class="clear"></div>
</div>

