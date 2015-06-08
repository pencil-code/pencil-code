#!/usr/bin/perl

# Authors:
#
# Philippe Bourdin
# http://www.Bourdin.ch/Philippe/

$defaultpattern = "*.f90";
$defaultdirectory = ".";

$found = 0;
$nothing = 0;
$testonly = 0;
$recursive = 0;
$dir_depth = 0;
$dir_depth_str = '../';

$pattern = $defaultpattern;
$startdirectory = $defaultdirectory;

while ($ARGV[0] =~ /^-/) {

	$argument = uc (shift ());
	if ($argument eq "-R") { $recursive = 1; }
	elsif ($argument eq "-T") { $testonly = 1; }
	elsif ($argument eq "-P") {
		$pattern = shift();
		if (!$pattern) { &usage ("No Pattern given ..."); }
		$pattern =~ s/\#\?/\*/g;
		$pattern =~ s/\*\*+/\*/g;
	}
	elsif ($argument eq "-D") {
		$startdirectory = shift();
		if (!$startdirectory) { &usage ("No Directory given ..."); }
		$startdirectory =~ s/\*+//g;
		$startdirectory =~ s/\/+$//;
	}
	else {
		&usage ("Option \"".$argument."\" unkown ...");
	}
}

if (shift()) { &usage ("Wrong amount of arguments ..."); }

@affected_files = ();

@directorylist = split (/\|/, $startdirectory);
foreach $directoryitem (@directorylist) {
	if ($directoryitem) {
		if ($directoryitem eq ".") { $directoryitem = ""; }
		elsif ($directoryitem !~ /\/$/) { $directoryitem .= "/"; }
	}

	if ($directoryitem) { print "in directory ".$directoryitem."\n"; }
	&subdir($directoryitem, @affected_files);
}

if (($found == 0) && ($nothing == 0)) {
	print "Sorry, there were no files that matched the pattern \"".$startdirectory.$pattern."\" !\n";
}
elsif ($found == 0) {
	print "Sorry, $nothing files matched the pattern, but the search-string hasn't been found.\n"
}
elsif ($nothing == 0) {
	print "All ".$found." files that matched the pattern contained the search-string.\n";
}
else {
	print "Finished. ".$found." of ".($found + $nothing)." matching files contained the search-string.\n";
}

exit;


sub subdir {

	my $directory = $_[0];
	my $filename = "";
	my $out;
	my $line;
	my $count;
	my $save;
	my $parent_dir;

	if ($directory) {
		if ($directory eq ".") { $directory = ""; }
		elsif ($directory !~ /\/$/) { $directory .= "/"; }
	}
	
	my @dirlist = glob ($directory."*");
	if (($recursive == 1) && ($#dirlist > 0)) {
		$dir_depth++;
		foreach $filename (@dirlist) {
			if (-d $filename) {
				&subdir ($filename);
			}
		}
		$dir_depth--;
	}

	$parent_dir = "";
	if ($dir_depth >= 1) {
		for ($pos = 1; $pos <= $dir_depth; $pos++) {
			$parent_dir .= $dir_depth_str;
		}
	}
	else {
		return;
	}

	@patternlist = split (/\|/, $pattern);
	foreach $patternitem (@patternlist) {
		my @filelist = glob ($directory.$patternitem);
		foreach $filename (@filelist) {

			$out = "";
			open (FILE, $filename);
			while (defined ($line = <FILE>)) { $out .= $line; }
			close (FILE);

			$count = 0;
			while ($out =~ s/\n *(\!+ *\*+\s*\n\s*subroutine +read_[^\s]+_pars.*?\n +include ")(parallel_unit.h" *\n.*?\n +end *subroutine +write_[^\s]+_pars *\n\s*\!+ *\*+) *\n/\n${1}${parent_dir}${2}\n/is) {
				$code = $1.$2;
				$count++;
			}

			if ($count > 0) {
				$found++;
				print "$filename\tfound $count times";
				$out_filename  = $filename;
				if ($testonly) { $out_filename .= ".new"; }
				print " ... ";
				open (FILE, "> ".$out_filename);
				print FILE $out;
				close (FILE);
				print "done.\n";
				push (@affected_files, $filename);
			}
			else {
				$nothing++;
			}
		}
	}
	return;
}


sub usage {

	my $message = $_[0];

	print "\nSorry, ".$message."\n\n";
	print "Usage:\tperl mass-changer.pl [-r] [-t] [-d Directory] [-p Pattern]\n\n";
	print "\t-r = Recursive - descent into directories\n";
	print "\t-t = Test only - create a new file, don't replace old file\n";
	print "\t-d = Directory to start with - default is ".$defaultdirectory."\n";
	print "\t     A list like \"dirA|dirB|dirE|dirF\" is also allowed here.\n";
	print "\t-p = Pattern to match - default is \"".$defaultpattern."\"\n\n";
	print "Example:\n";
	print "\tmass-changer.pl -r -d pages -p \"in*.(f??|h)\"\n";
	print "This recursively collects all files in the directory \"pages\" starting with\n";
	print "\"in\" and ending with \".f\",\".f90\",\".f2k\" or \".h\" and then searches all\n";
	print "strings starting with \"abc\" and ending with \"efg\" by ignoring the case.\n\n";
	exit;
}

