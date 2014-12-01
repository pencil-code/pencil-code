#!/usr/bin/perl -w

# Description:
#   Unit tests for Tests::NumericFileComparator
#
#   These tests compare two files and verify that the accuracies are
#   interpreted as expected.

use strict;

use Test::More;
use feature qw/say/;

my $pc_home;

BEGIN {
    # Make sure ${PENCIL_HOME}/lib/perl is in the Perl path
    if (-d "$ENV{PENCIL_HOME}/lib/perl") {
        $pc_home = $ENV{PENCIL_HOME};
    } else {
        if ($0 =~ m!(.*[/\\])!) {
            $pc_home = "$1../../../../../../..";
        }
    }
    unshift @INC, "$pc_home/lib/perl" if $pc_home;
}

BEGIN {
    use_ok('Test::NumericFileComparator');
}


(my $test_dir = $0) =~ s{/[^/]*$}{};


compare_files_2_1();
compare_files_3_1();

done_testing();


sub compare_files_2_1 {
    my $file = "$test_dir/data/time-series-1.dat";
    my $ref_file = "$test_dir/data/time-series-2.dat";
    my $comparator = Test::NumericFileComparator->new($ref_file);
    my @message = $comparator->compare($file);
    ok(! @message, 'time-series-{1,2} equal')
      or diag(join("\n  ", "Files $ref_file, $file differ:", @message));
}


sub compare_files_3_1 {
    my $file = "$test_dir/data/time-series-1.dat";
    my $ref_file = "$test_dir/data/time-series-3.dat";
    my $comparator = Test::NumericFileComparator->new($ref_file);
    my @message = $comparator->compare($file);
    ok(@message, 'time-series-{1,3} differ')
      or diag("Files $ref_file, $file should differ, but are equal");
}

