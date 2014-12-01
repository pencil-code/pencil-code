#!/usr/bin/perl -w

# Description:
#   Unit tests for Test::NumberComparator

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

($a, $b) = (1.0, 0.995);
ok($comparator->compare($a, $b) == 0, "equal within abs. and rel. accuracy")
    or diag("Expected equality, got ", $comparator->format_comparison($a, $b));

($a, $b) = (100., 99.5);
ok($comparator->compare($a, $b) == 0, "equal within relative accuracy")
    or diag("Expected equality, got ", $comparator->format_comparison($a, $b));

($a, $b) = (1.0e-3, 0.99e-3);
ok($comparator->compare($a, $b) == 0, "equal within relative accuracy")
    or diag("Expected equality, got ", $comparator->format_comparison($a, $b));

($a, $b) = (1.0e-3, 0.9e-3);
ok($comparator->compare($a, $b) == +1, "greater than")
    or diag("Expected $a > $b, got ", $comparator->format_comparison($a, $b));

($a, $b) = (1.0e-3, 1.1e-3);
ok($comparator->compare($a, $b) == -1, "less than")
    or diag("Expected $a < $b, got ", $comparator->format_comparison($a, $b));

done_testing();

