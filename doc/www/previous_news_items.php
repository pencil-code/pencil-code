<?
	 include "header.inc";
 ?>
<div class="centcolumn">
<div class="centcolumnpad">
<h2>Previous News.</h2>

<li>13 Jan 2012<br>
New website for the Pencil Code User Meeting 2012 on 18-21 June
in Helsinki.  For more information and registering, see:
<a href=http://agenda.albanova.se/conferenceDisplay.py?confId=3128>http://agenda.albanova.se/conferenceDisplay.py?confId=3128</a>.</li>

<li>13 Jan 2011<br>
New website for the Pencil Code User Meeting 2011 on 24-28 October
in Toulouse.  For more information and registering, see:
<a href=http://norlx51.albanova.se/~dintrans/Meeting2011/index.html>http://norlx51.albanova.se/~dintrans/Meeting2011/index.html</a>.</li>

<li>2 Oct 2010<br>
The Pencil Code User Meeting 2011 will be 24-28 October in Toulouse.
For more information and registering, see: http://www.ast.obs-mip.fr/users/dintrans/PencilCode2011/ (now outdated address).</li>

<li>6 Aug 2010<br>
<a href=meeting2010.php>Notes from the Pencil Code User Meeting 2010</a> in New York.</li>

<li>22 Jan 2010<br>
The <a href=http://www.research.amnh.org/astrophysics/pencil/pencil.html>Pencil Code User Meeting 2010</a> will be 26-30 July in New York.</li>

<li>22 Mar 2009<br>
The Pencil Code User Meeting 2009 will be 24-28 August in Heidelberg.
Please register on: <a href=http://agenda.albanova.se/conferenceDisplay.py?confId=1190>http://agenda.albanova.se/conferenceDisplay.py?confId=1190</a>.</li>

<li>21 Sep 2008<br>
Now the code is available under
<a href=http://pencil-code.googlecode.com/>http://pencil-code.googlecode.com/</a>.
Please consult the <a href="http://code.google.com/p/pencil-code/wiki/SvnTransition">following notes</a>
to do the transition.</li>

<li>12 Sep 2008<br>
Following our discussions in Leiden, Tobi continued negotiating
with the google people, who have now fixed the bug in getting
the google svn repository in sync with ours. Now it is working.
To test things first, Tobi has set up a test account on:
<a href=http://code.google.com/p/pencil-code-playground/>http://code.google.com/p/pencil-code-playground/</a>.
Please sign in and give it a try.
Tobi can help if there are problems.</li>

<li>19 Aug 2008<br>
The Pencil Code is now being served via SVN.
Users may want to do the <a href="svn-transition.html">following change</a>.
CVS is still up, and only the f90/pencil-code is now deactivated
and replaced by
<a href=https://anonymous@svn.nordita.org:/svn/pencil-code/trunk>https://brandenb@svn.nordita.org:/svn/pencil-code/trunk</a>.
The minutes of the sessions are in this
<a href=http://pcum08.blogspot.com/>blog</a>.</li>

<li>22 Jun 2008<br>
As you may have noticed, CVS is now being served from Stockholm, but we are
still waiting with the transfer to SVN. To ease the transition for everyone,
and to provide help, we decided to do it during the
<a href="http://www.strw.leidenuniv.nl/~ajohan/pencil2008/index.php
">Pencil Code User Meeting 2008</a>.
It will be held in Leiden and you are encouraged to sign up!</li>
 
<li>9 Mar 2008<br>
The <a href="http://www.strw.leidenuniv.nl/~ajohan/pencil2008/index.php
">Pencil Code User Meeting 2008</a> will be held in Leiden; please sign up!</li>

<li>2 May 2008<br>
Change to SVN now imminent! As discussed during the last
<a href=http://agenda.albanova.se/conferenceDisplay.py?confId=185>
Pencil Code User Meeting</a> there are only good reasons to
change to SVN.
It is now anticipated that this will happen some time in May.
We'll be in touch with a number of users on the regular
CVS notification email list to make sure it works.
The background is that the server is now in place and will
be installed starting tomorrow.
Note that the
<a href=http://www.nordita.org/~brandenb/Meetings/2007/Pencil/>
temporary SVN server</a> used last year was only a temporary test and
has nothing to do with the new one.  Our current experimental svn
server is
<a href=https://svn.nordita.org/svn/pencil-code/trunk>https://svn.nordita.org/svn/pencil-code/trunk</a>,
but this too is only a temporary name.</li>

<li>29 Aug 2007<br>
A detailed <a href="http://www.nordita.org/~brandenb/Meetings/2007/Pencil/
">Report of the Pencil Code User Meeting 2007</a> is now available.
An important item concerns the planned migration of the Pencil Code from
CVS to Subversion (SVN). Please report potential issues and problems to
one of us.</li>

<li>3 Jul 2007<br>
The CVS server name has changed. For details see the
<a href="http://wiki.nordita.org/index.php/Breaking_News"
>Nordita Computing News</a>.</li>

<li>14 Mar 2007<br>
Unfortunately there was a problem with the registration
page between 27 February and 12 Mach 2007. If you tried to register during
that time and got an error message, please try again. Our apologies for that!</li>
  
<li>29 Dec 2006<br>
The next <a href="http://indico.albanova.se/conferenceDisplay.py?confId=185"
>Pencil Code User Meeting</a> will be held in Stockholm 14-17 August.</li>

<li>4 May 2006<br>
There will be a <a href="http://www.nordita.dk/conference/PencilCode06/"
>Pencil Code User Meeting</a> during 13-15 July. When a detailed
program becomes available, you will see which parts may be of interest
to you and which are not. Everybody is welcome to attend either for the
whole time or just to pop in for some of the sessions.</li>

<li>28 July 2006<br>
Those who remember the news since
<a href="previous_news_items.html">18 June 2005</a> will remember that
the cvs repository had an extra branch (the eos branch).
This has now been closed and during the last
<a href="http://www.nordita.dk/~brandenb/get-together/meetings/pencil-workshop05a.html">Pencil Code Development</a> Meeting
we have incorporated all the feature of the eos branch into the main trunk.
Until 28 July the public repository was not updated (to prevent
users seeing broken code), but meanwhile we have been working
on production runs ourselves without problems, so we deem the code
fit for public consumption.
There are some minor changes that will be necessary in the .in-files,
such as adding eos_init_pars (possibly with gamma=1.) and
viscosity_run_pars (with the value of nu). It will be best to check
out the samples directory for how to do it.
Alternatively, wait until Wolfgang will have prepared his script
that converts from the old files to the new ones.</li>

<li>18 June 2005<br>
The cvs repository has an extra branch (cvs up -r eos). Around 26 June
2005 we plan to incorporate many of the feature of the eos branch into
the main trunk.
This is done in connection with the
<a href="http://www.nordita.dk/~brandenb/get-together/meetings/pencil-workshop05a.html">Pencil Code Development</a>
workshop.
During that time the public repository will not be updated until we are sure
it is fit for public consumption.</li>
</ul>

<hr/>


 </div>
</div>
  
      

<? include "footer.inc" ?>

