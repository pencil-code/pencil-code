<!-- $Id$ -->
<?
	include "header.inc";

	$meetings = array (
		array ( title => "10th meeting", date => "summer 2014", link => "", notes => "", city => "?", venue => "", country => "" ),
		array ( title => "9th meeting", date => "17-20 Jun, 2013", link => "http://www.astro.lu.se/~michiel/PC2013/", notes => "", city => "Lund", venue => "Lund Observatory", country => "Sweden" ),
		array ( title => "8th meeting", date => "18-21 Jun, 2012", link => "http://agenda.albanova.se/conferenceDisplay.py?confId=3128", notes => "/meeting_2012.php", city => "Helsinki", venue => "Physics Department", country => "Finland" ),
		array ( title => "7th meeting", date => "24-28 Oct, 2011", link => "http://norlx51.albanova.se/~dintrans/Meeting2011/", notes => "/meeting_2011.php", city => "Toulouse", venue => "Observatoire Midi-Pyr&eacute;n&eacute;es", country => "France" ),
		array ( title => "6th meeting", date => "26-30 Jul, 2010", link => "http://www.research.amnh.org/astrophysics/pencil/pencil.html", notes => "/meeting_2010.php", city => "New York", venue => "American Museum of National History", country => "USA" ),
		array ( title => "5th meeting", date => "24-28 Aug, 2009", link => "http://agenda.albanova.se/conferenceDisplay.py?confId=1190", notes => "/meeting_2009.php", city => "Heidelberg", venue => "Max Planck Institute for Astronomy", country => "Germany" ),
		array ( title => "4th meeting", date => "19-22 Aug, 2008", link => "http://www.strw.leidenuniv.nl/~ajohan/pencil2008/", notes => "/meeting_2008.php", city => "Leiden", venue => "Leiden Observatory", country => "Netherlands" ),
		array ( title => "3rd meeting", date => "14-17 Aug, 2007", link => "http://agenda.albanova.se/conferenceDisplay.py?confId=185", notes => "/meeting_2007.php", city => "Stockholm", venue => "Nordita", country => "Sweden" ),
		array ( title => "2nd meeting", date => "13-15 Jul, 2006", link => "http://www.nordita.dk/conference/PencilCode06/", videos => "http://videos.albanova.se/conference/PencilCode06/", city => "Copenhagen", venue => "Nordita", country => "Denmark" ),
		array ( title => "1st meeting", date => "26-28 Jun, 2005", link => "http://www.nordita.org/~brandenb/get-together/meetings/pencil-workshop05a.html", notes => "http://www.nordita.dk/~brandenb/get-together/meetings/pencil-workshop05a.html", city => "Copenhagen", venue => "Nordita", country => "Denmark" ),
	);
 ?>
<div class="centcolumnpad">
<h2>Meetings.</h2>

<p align="center">
<img src="/pics/2006b.jpg" width="240" height="180" border="0" alt="2006"><img src="/pics/2011a.jpg" width="240" height="180" border="0" alt="2011"><img src="/pics/2011b.jpg" width="240" height="180" border="0" alt="2011">
</p>

<table border="0" cellspacing="10" cellpadding="0">
<?
	foreach ($meetings as $meeting) {
 ?>
<tr>
<td align="right" STYLE="white-space:nowrap;"><? ifecho ("", $meeting['date'], ":"); ?></td>
<td STYLE="white-space:nowrap;"><? echolink ($meeting['link'], $meeting['title']); ?></td>
<td STYLE="white-space:nowrap;"><? iflink  ($meeting['notes'], "[notes]"); ?><? iflink ($meeting['videos'], "[videos]"); ?></td>
<td STYLE="font-size:12px;"><? ifecho ("in ", $meeting['city'], ""); ?><? ifecho (", ", $meeting['venue'], ""); ?><? ifecho (" (", $meeting['country'], ")"); ?>.</td>
</tr>
<?
	}
 ?>
</table>

<p align="center">
<img src="/pics/2005a.jpg" width="240" height="180" border="0" alt="2006"><img src="/pics/2005b.jpg" width="240" height="180" border="0" alt="2006"><img src="/pics/2006a.jpg" width="240" height="180" border="0" alt="2006">
</p>
</div>
<?
	include "footer.inc";
 ?>

