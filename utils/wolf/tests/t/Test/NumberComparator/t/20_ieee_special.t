#!/usr/bin/perl -w

# Description:
#   Unit tests for Test::NumberComparator
#
#   Test NaN and ±Inf comparisons

use strict;

use Test::More;

BEGIN {
    # Make sure ${PENCIL_HOME}/lib/perl is in the Perl path
    if (-d "$ENV{PENCIL_HOME}/lib/perl") {
        unshift @INC, "$ENV{PENCIL_HOME}/lib/perl";
    } else {
        if ($0 =~ m!(.*[/\\])!) {
            unshift @INC, "$1../../../../../../../lib/perl";
        }
    }
}

BEGIN {
    use_ok('Test::NumberComparator');
}


my $abs_accuracy = 1.0e-4;
my $rel_accuracy = 1.0e-2;
my $comparator = Test::NumberComparator->new($abs_accuracy, $rel_accuracy);

my ($a, $b);

for my $special_value (qw/NaN Inf +Inf -Inf/) {
    _check_special_number(
        $special_value, $special_value, [0],
        '$special_value equal to itself', 'equality'
        );

    my $real_1 = 1.23;
    _check_special_number(
        $special_value, $real_1, [-1, 1],
        "$special_value different from real number $real_1", 'difference'
        );

    my $real_2 = -3.45e-6;
    _check_special_number(
        $special_value, $real_2, [-1, 1],
        "$special_value different from real number $real_2", 'difference'
        );
}

for my $special_value (qw/+Inf -Inf/) {
    _check_special_number(
        'Inf', $special_value, [0],
        '$special_value equal to Inf', 'equality'
        );
}

for my $special_value (qw/Inf +Inf -Inf/) {
    _check_special_number(
        'NaN', $special_value, [-1, 1],
        "$special_value different from NaN", 'difference'
        );
}

_check_special_number(
    '-Inf', '+Inf', [-1, 1],
    "-Inf different from +Inf", 'difference'
    );


sub _check_special_number {
    my ($a, $b, $expected, $name_fragment, $diag_fragment) = @_;

    _check_forward($a, $b, $expected, $name_fragment, $diag_fragment);
    _check_forward($b, $a, $expected, $name_fragment, $diag_fragment);
}


sub _check_forward {
    my ($a, $b, $expected, $name_fragment, $diag_fragment) = @_;

    my $name = "special number $name_fragment";
    my $got = $comparator->compare($a, $b);
    foreach my $exp (@$expected) {
        if ($exp == $got) {
            ok(1, $name);
            return;
        }
    }
    ok(0, $name);
    diag(
         "Expected $diag_fragment, got ",
         $comparator->format_comparison($a, $b)
         );
}


done_testing();

