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


(my $test_dir = $0) =~ s{/[^/]*$}{};
my $pressure_values = [
    1.0000000,
    0.5000000,
    0.8333333,
    0.5833333,
    0.2833333,
    0.1066667,
    0.04759523,
    0.01345238,
    0.003124355,
    ];
my $temperature_values = [
    1.0000,
    0.6667,
    0.8667,
    0.7238,
    0.8349,
    0.7440,
    0.8209,
    0.7942,
    0.7734,
    ];


read_column_data_bare();
read_column_data_fancy();
read_column_data_timeseries();

read_line_data_simple();
read_line_data_fancy();

read_line_data1();
read_line_data2();

compare_files_2_1();
compare_files_3_1();

done_testing();


sub read_column_data_bare {
    my $file = "$test_dir/data/column_bare.dat";
    my $comp = Test::NumericFileComparator->new($file);
    ok($comp, "Parse $file");

    _compare_string_lists(
        ['variable_0', 'variable_1'],
        $comp->{VARS},
        $file
        );

    _compare_number_lists(
        $pressure_values,
        $comp->{VALUES}->{'variable_0'},
        "$file:variable_0"
        );
    _compare_number_lists(
        $temperature_values,
        $comp->{VALUES}->{'variable_1'},
        "$file:variable_1"
        );

    _compare_number_lists(
        [1.5e-7, 0],
        $comp->{ACCS}->{'variable_0'},
        'accuracy[variable_0]'
        );
    _compare_number_lists(
        [1.5e-4, 0],
        $comp->{ACCS}->{'variable_1'},
        'accuracy[variable_1]'
        );
}


sub read_column_data_fancy {
    my $file = "$test_dir/data/column_fancy.dat";
    my $comp = Test::NumericFileComparator->new($file);
    ok($comp, "Parse $file");

    _compare_string_lists(
        ['Pressure', 'Temperature'],
        $comp->{VARS},
        $file
        );

    _compare_number_lists(
        $pressure_values,
        $comp->{VALUES}->{'Pressure'},
        "$file:Pressure"
        );
    _compare_number_lists(
        $temperature_values,
        $comp->{VALUES}->{'Temperature'},
        "$file:Temperature"
        );

    _compare_number_lists(
        [0, 5.0e-6],
        $comp->{ACCS}->{'Pressure'},
        'accuracy[Pressure]'
        );
    _compare_number_lists(
        [8.0e-3, 0],
        $comp->{ACCS}->{'Temperature'},
        'accuracy[Temperature]'
        );
}


sub read_column_data_timeseries {
    my $file = "$test_dir/data/column_timeseries.dat";
    my $comp = Test::NumericFileComparator->new($file);
    ok($comp, "Parse $file");

    _compare_string_lists(
        $comp->{VARS},
        [qw/it t urms umax rhom ssm dtc dtu/],
        $file
        );
    my @it = @{$comp->{VALUES}->{'it'}};
    my @t = @{$comp->{VALUES}->{'t'}};
    my @ssm = @{$comp->{VALUES}->{'ssm'}};
    _compare_numbers(30, $it[3], "$file: it[3] value,");
    _compare_numbers(0.20, $t[3], "$file: t[3] value,");
    _compare_numbers(-0.4471, $ssm[3], "$file: ssm[3] value,");

    _compare_number_lists(
        [0.015, 0],
        $comp->{ACCS}->{'t'},
        "$file: t accuracy,"
        );
    _compare_number_lists(
        [0.00015, 0],
        $comp->{ACCS}->{'ssm'},
        "$file: ssm accuracy,",
        1e-12
        );
}


sub read_line_data_simple {
    my $file = "$test_dir/data/line_simple.dat";
    my $comp = Test::NumericFileComparator->new($file);
    ok($comp, "Parse $file");

    _compare_string_lists(
        ['Pressure', 'Temperature'],
        $comp->{VARS},
        $file
        );

    _compare_number_lists(
        [0.6345238],
        $comp->{VALUES}->{'Pressure'},
        "$file:Pressure"
        );
    _compare_number_lists(
        [0.7543],
        $comp->{VALUES}->{'Temperature'},
        "$file:Temperature"
        );

    _compare_number_lists(
        [1.5e-7, 0],
        $comp->{ACCS}->{'Pressure'},
        'accuracy[Pressure]'
        );
    _compare_number_lists(
        [1.5e-4, 0],
        $comp->{ACCS}->{'Temperature'},
        'accuracy[Temperature]'
        );
}


sub read_line_data_fancy {
    my $file = "$test_dir/data/line_fancy.dat";
    my $comp = Test::NumericFileComparator->new($file);
    ok($comp, "Parse $file");

    _compare_string_lists(
        ['Pressure', 'Temp(5,7,3:5)'],
        $comp->{VARS},
        $file
        );

    _compare_number_lists(
        [0.6345238, 0.7345238, 0.655238],
        $comp->{VALUES}->{'Pressure'},
        "$file:Pressure"
        );
    _compare_number_lists(
        [0.7543, 0.7411, 0.7123],
        $comp->{VALUES}->{'Temp(5,7,3:5)'},
        "$file:Temperature"
        );

    _compare_number_lists(
        [0, 1.5e-7],
        $comp->{ACCS}->{'Pressure'},
        'accuracy[Pressure]'
        );
    _compare_number_lists(
        [8.0e-3, 0],
        $comp->{ACCS}->{'Temp(5,7,3:5)'},
        'accuracy[Temperature]'
        );
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


## Helper functions ##

sub _compare_number_lists {
    my ($exp, $got, $label, $abs_acc) = @_;

    my @vars = @$got;
    my @exp_vars = @$exp;

    my $count = @vars + 0;
    my $exp_count = @exp_vars + 0;
    _compare_numbers($exp_count, $count, "$label: value count --");
    for (my $i = 0; $i < @exp_vars; $i++) {
        _compare_numbers(
            $exp_vars[$i],
            $vars[$i],
            "$label: element $i,",
            $abs_acc
            );
    }
}


sub _compare_string_lists {
    my ($exp, $got, $label) = @_;

    my @vars = @$got;
    my @exp_vars = @$exp;

    my $count = @vars + 0;
    my $exp_count = @exp_vars + 0;
    _compare_numbers($exp_count, $count, "$label: value count --");
    for (my $i = 0; $i < @exp_vars; $i++) {
        _compare_strings($exp_vars[$i], $vars[$i], "$label: column $i,");
    }
}


sub _compare_numbers {
    my ($exp, $got, $label, $abs_acc) = @_;

    if (defined $abs_acc) {
        ok(
            abs($got - $exp) <= $abs_acc,
            "$label expected $exp, got $got, differs by > $abs_acc"
            );
    } else {
        ok($exp == $got, "$label expected $exp, got $got");
    }
}


sub _compare_strings {
    my ($exp, $got, $label) = @_;

    ok($exp eq $got, "$label expected '$exp', got '$got'");
}
