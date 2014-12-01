#!/usr/bin/perl -w

# Description:
#   Run all tests

use strict;
use File::Find;
use TAP::Harness;

my %args = ();

my $harness = TAP::Harness->new(\%args);
(my $t_dir = $0) =~ s{/[^/]*$}{};

my @tests;
sub is_test {
    push @tests, $File::Find::name
        if $_ =~ /\.t$/ && -x $_;
}
find({
      wanted => \&is_test,
      follow => 1,              # follow symlinks
     },
     $t_dir
    );

$harness->runtests(sort @tests);

