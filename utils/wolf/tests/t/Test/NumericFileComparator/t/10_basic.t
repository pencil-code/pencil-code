#!/usr/bin/perl -w

# Description:
#   Unit tests for Tests::NumericFileComparator

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


read_column_data();
read_line_data1();
read_line_data2();

done_testing();


sub read_column_data {
    my $dir = "$pc_home/tests";
    my $ts_file = "$dir/input/serial-1/time-series.dat";
    my $comp2 = Test::NumericFileComparator->new($ts_file);
    ok($comp2, "Parse $ts_file");

    my @vars = @{$comp2->{VARS}};
    ok(@vars == 7, "$ts_file: column count");
    my @cols = (qw/it t dt urms rhom ecrm ecrmax/);
    for (my $i = 0; $i < @cols; $i++) {
        ok(
           $vars[$i] eq $cols[$i],
           "$ts_file: column $i, expected name '$cols[$i]', got $vars[$i]"
          );
    }

    ok($vars[1] eq 't', "$ts_file: column 1 name, expected 't', got($vars[1]");

    my @it = @{$comp2->{VALUES}->{'it'}};
    my @t = @{$comp2->{VALUES}->{'t'}};
    ok($it[3] == 150, "ts_file: it[3] value, expected 150, got $it[3]");
    ok($t[3] == 1.480, "ts_file: t[3] value, expected 1.480, got $t[3]");

    my $t_acc = $comp2->{ACCS}->{'t'};
    ok($t_acc == 0.0015, "ts_file: t accuracy, expected 0.0015, got $t_acc");
}


sub read_line_data1 {
    my $dir = "$pc_home/tests";

    my $ref_file = "$dir/python/read-time-series/test1.ref";
    my $comp1 = Test::NumericFileComparator->new($ref_file);
    ok($comp1, "Parse $ref_file");
}


sub read_line_data2 {
    my $dir = "$pc_home/samples/helical-MHDturb/tests";
    my $file = "$dir/read_data.ref";
    my $comp1 = Test::NumericFileComparator->new($file);
    ok($comp1, "Parse $file");
}
