#
#                         ScriptTester.pm
#                         ---------------
#
# Description:
#   Run script tests inside a set of directories.
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Test::ScriptTester;

use warnings;
use strict;
use Carp;
use vars qw($VERSION);

use feature 'say';

##use critic

$VERSION = '0.1';


=head1 NAME

Test::ScriptTester - Run script tests in a set of directories

=head1 SYNOPSIS

  use Test::ScriptTests;

  $tester = Test::ScriptTester->new(
      [dir1, dir2, ...],
      {'python' => 'python', 'idl' => '/usr/bin/gdl'}
      );

  $tester->run();

  my @tests = $tester->list_tests();

  $tester->run_tests(@tests[0..2);

  my %default_interpreters = Test::ScriptTests::get_default_interpreters()

  my $idl_interpreter = Test::ScriptTests::find_interpreter_for('idl')


=head1 DESCRIPTION

Scan the given directories for subdirectories named 'tests/'; in each
tests directory, run all script tests.

A I<script test> consists of a test file that

=over 4

=item *

is named <test_name>.<suff>, where <suff> is a known suffix <suff>
(currently supported: py, pro),

=item *

is executable and can be run from the tests directory (thus, many of the
test scripts read data from '../data/'),

=item *

when run, writes a file <test_name>.out in the same directory.

=item *

There exist a file <test_name>.ref in the same directory that defines the
reference data and possibly accuracy.

=back

Each test script is run (if we find an appropriate interpreter)
and the <test_name>.out file is compared to the reference data
<test_name>.ref .


=head2 Methods

=over 4

=cut


=item B<Test::ScriptTester-E<gt>new>($dirs, $interpreters)

Create a new object that searches the given directories

  $dirs = [dir1, dir2, ...]

and uses the given interpreters, e.g.

  {'python' => 'python', 'idl' => '/usr/bin/gdl'} )

Only test types appearing as keys in the interpreters map are run.
To get the default map of interpreters, use
I<Test::ScriptTests::get_default_interpreters>().

=cut

sub new {
#
#   Test::ScriptTester->new($dirs_ref, $interpreters_ref);
#   $tester->new($dirs_ref, $interpreters_ref);
#
    my $proto = shift;          # either classref or object ref or string
    my @argv  = @_;

    my $self = {};
    my $class;
    my $parent = {};
    if (ref($proto)) {
        $class = ref($proto);
    } else {
        $class = $proto;
    }

    if (@argv == 2) {
        ($self->{ABS_ACC}, $self->{REL_ACC}) = @argv;
    } else {
        croak('Usage: Test::NumberComparator->new($abs_acc, $rel_acc)');
    }

    $self->{DEBUG} = 0;

    bless $self, $class;
    return $self;
}


=item B<$tester-E<gt>run>()

Run all tests supported by this $tester.

=cut

sub run {
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    if (_is_special_ieee($a)) {
        return _compare_special_ieee($a, $b);
    }
    if (_is_special_ieee($b)) {
        return - _compare_special_ieee($b, $a);
    }

    if ($self->equal($a, $b)) {
        return 0;
    } else {
        return $a <=> $b;
    }
}


=item B<$tester-E<gt>list_tests>()

List all tests supported by this $tester.
Returns a list of list refs

  ([$tests_dir1, $script1], ...)

=cut

sub list_tests {
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    return $self->_equal_abs($a, $b) || $self->_equal_rel($a, $b);
}


=item B<$tester-E<gt>run_tests>(@tests)

Run the given tests.

Each test is a list ref [$tests_dir, $script], where $script is the path
to an executable script file, either relative to $tests_dir, or absolute.

=cut

sub run_tests {
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    my $result = $self->compare($a, $b);

    if ($result == 0) {
        return sprintf(
            "$a == $b up to both accuracies (abs: %g, rel: %g)",
            $self->{ABS_ACC}, $self->{REL_ACC}
        );
    }

    my @violated_accuracies = ();
    if (! $self->_equal_abs($a, $b)) {
        push @violated_accuracies,
          sprintf("absolute accuracy %g", $self->{ABS_ACC})
            if $self->{ABS_ACC} > 0;
    }
    if (! $self->_equal_rel($a, $b)) {
        push @violated_accuracies,
          sprintf("relative accuracy %g", $self->{REL_ACC})
            if $self->{REL_ACC} > 0;
    }
    @violated_accuracies = ("strict comparison") unless @violated_accuracies;
    my $reason = join(', ', @violated_accuracies);

    if ($result == -1) {
        return sprintf("$a < $b according to %s", $reason);
    } elsif ($result == 1) {
        return sprintf("$a > $b according to %s", $reason);
    } else {
        croak "Unexpected: \$self->compare($a, $b) returned nonsense.";
    }

}


=item B<get_default_interpreters>()

Return a hash ref

  { type1 => interpreter1,
    type2 => interpreter1, ...
  }

representing the default interpreters (hash values) for all known test
types (hash keys).

=cut

sub get_default_interpreters {


=item B<find_interpreter_for>($test_type)

Return the default interpreter (path of an executable) for the given
$test_type, or undef.

=cut

sub find_interpreter_for {
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    my $deviation = abs($a - $b);
    return $deviation <= $self->{ABS_ACC};

}


# ---------------------------------------------------------------------- #

=back

=head1 EXAMPLES

  use Test::ScriptTester;

  my $tester = Test::ScriptTester->new(
      ["$env{PENCIL_HOME}/tests"],
      \Test::ScriptTester::get_default_interpreters()
      );
  $tester->run();  # Run all test scripts under $PENCIL_HOME/tests/

=cut


1;

__END__



# End of file
