# CVS: $Id: tst.csh,v 1.2 2004-03-22 12:55:34 brandenb Exp $
(echo "\$Id: tst.csh,v 1.2 2004-03-22 12:55:34 brandenb Exp $";echo `whoami`"@"`hostname`":"`pwd`;echo "";start.csh;echo "";run.csh) >>.TEST.`timestr`
