#
#                         NumberComparator.pm
#                         --------------------
#
# Description:
#   Compare numbers up to some specified accuracy.
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Test::NumberComparator;

use warnings;
use strict;
use Carp;
use vars qw($VERSION);

##use critic

$VERSION = '0.1';


=head1 NAME

Test::NumberComparator - Compare numbers up to some specified accuracy.


=head1 SYNOPSIS

  use Test::NumberComparator;

  my $abs_accuracy = 1.0e-4;
  my $rel_accuracy = 1.0e-2;
  my $comparator = Test::NumberComparator->new($abs_accuracy, $rel_accuracy);

  my ($a, $b) = (1.0, 0.995);  # equal within absolute and relative accuracy
  $comparator->compare($a, $b) == 0 \
      or warn "Expected equality, got ", $comparator->format_comparison($a, $b);


=head1 DESCRIPTION

Two numbers are considered equal if they are close enough by either
absolute or relative accuracy.

Equality is symmetric, but cannot be transitive.


=head2 IEEE special numbers

This module is designed for verifying consistency of new date with respect
to reference data and hence does not honour IEEE semantics for e.g. NaN:
If we got NaN before and still get it, that is consistent, even though
IEEE declares that NaN == NaN is false.

We also consider -Inf equal to Inf and +Inf equal to Inf, but -Inf is not
equal to +Inf.
In the light of -0 being numerically equal to +0 (up to any accuracy), we
may want to revise this decision, but for most applications, no NaN or Inf
values are expected anyway.


=head2 Methods

=over 4

=cut


=item B<Test::NumberComparator-E<gt>new>(I<$abs_acc>, I<$rel_acc>)

Create a new object with the given absolute and relative accuracies.

=cut

sub new {
#
#   Test::NumberComparator->new($abs_acc, $rel_acc);
#   $comparator->new($abs_acc, $rel_acc);
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


=item B<$comparator-E<gt>compare>($a, $b)

Compare $a and $b. If the modulus of the difference is smaller than either
I<abs_acc>, or I<rel_acc> * Max(|$a|, |$b|), return 0.
Otherwise return -1 if $a < $b, or +1 if $a > $b.

The B<compare()> method is (anti)symmetric in the sense that

  compare($a, $b) == - compare($b, $a)

It cannot be transitive.

One could call this method 'a sloppy version of $a <=> $b'.

Comparison to NaN or +/-Inf is essentially lexical and ignores accuracy.
The following table shows equality for combinations of IEEE numbers:

  |      | real | Inf | +Inf | -Inf | NaN |
  |------+------+-----+------+------+-----|
  | real |      | !=  | !=   | !=   | !=  |
  | Inf  |      | ==  | ==   | ==   | !=  |
  | +Inf |      |     | ==   | !=   | !=  |
  | -Inf |      |     |      | ==   | !=  |
  | NaN  |      |     |      |      | ==  |

=cut

sub compare {
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


sub _is_special_ieee {
    my ($x) = @_;

    return $x =~ m{ nan | inf }ix;
}


sub _compare_special_ieee {
    my ($a, $b) = @_;

    if ($a =~ m/NaN/i) {
        if ($b =~ m/NaN/i) {
            return 0;           # NaN == NaN for us
        } else {
            return +1;          # NaN > other (abritrary, but must be != 0)
        }
    } elsif ($a =~ m{ (?<sign_a> [-+])? Inf }ix) {
        my $sign_a = $+{sign_a};
        if ($b =~ m{ (?<sign_b> [-+])? Inf }ix) {
            # Comparing +/-Inf to +/-Inf or Inf
            my $sign_b = $+{sign_b};
            if (defined $sign_a && defined $sign_b) {
                return $sign_a cmp $sign_b;
            } else {
                return 0;       # Ind = +Inf or Inf = -Inf
            }
        } else {
            # Comparing Inf to something else
            if (defined $sign_a && $sign_a eq '-') {
                return -1;      # -Inf < $b
            } else {
                return 1;       # Inf = +Ibf > $b
            }
        }
    }
}


=item B<$comparator-E<gt>equal>($a, $b)

Compare $a and $b. If the modulus of the difference is smaller than either
I<abs_acc>, or I<rel_acc> * Max(|$a|, |$b|), return 1
(true), otherwise return '' (false).

The B<equal()> method is symmetric in the sense that

  equal($a, $b) == equal($b, $a)

It cannot be transitive.

=cut

sub equal {
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    return $self->_equal_abs($a, $b) || $self->_equal_rel($a, $b);
}


=item B<$comparator-E<gt>format_comparison>($a, $b)

Compare $a and $b and return a string explaining how the comparison turned
out.

=cut


sub format_comparison {
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


sub _equal_abs {
#
# If $a and $b are equal up to our absolute accuracy, return true, else
# return false.
#
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    my $deviation = abs($a - $b);
    return $deviation <= $self->{ABS_ACC};

}


sub _equal_rel {
#
# If $a and $b are equal up to our relative accuracy, return true, else
# return false.
#
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    my $deviation = abs($a - $b);
    return $deviation <= $self->{REL_ACC} * $self->_max(abs($a), abs($b));

}


sub _max {
#
# Return the maximum of two numbers
#
    my $self = shift;
    my ($a, $b) = @_[0, 1];

    if ($a >= $b) {
        return $a;
    } else {
        return $b;
    }

}

# ---------------------------------------------------------------------- #

=back

=head1 EXAMPLES

  use Test::NumberComparator;

  my $abs_accuracy = 1.0e-4;
  my $rel_accuracy = 1.0e-2;
  my $comparator = Test::NumberComparator->new($abs_accuracy, $rel_accuracy);
  my $abs_comparator = Test::NumberComparator->new($abs_accuracy, 0.0);
  my $rel_comparator = Test::NumberComparator->new(0.0, $rel_accuracy);

  my ($a, $b);

  ($a, $b) = (1.0, 0.995);  # equal within absolute and relative accuracy
  $comparator->compare($a, $b) == 0 \
      or warn "Expected equality, got ", $comparator->format_comparison($a, $b);

  ($a, $b) = (100., 99.5);  # equal within relative accuracy
  $comparator->compare($a, $b) == 0 \
      or warn "Expected equality, got ", $comparator->format_comparison($a, $b);

  ($a, $b) = (1.0e-3, 0.99e-3);  # equal within relative accuracy
  $comparator->compare($a, $b) == 0 \
      or warn "Expected equality, got ", $comparator->format_comparison($a, $b);

  ($a, $b) = (1.0e-3, 0.9e-3);  # not equal
  $comparator->compare($a, $b) == +1 \
      or warn "Expected $a > $b, got ", $comparator->format_comparison($a, $b);

  ($a, $b) = (1.0e-3, 1.1e-3);  # not equal
  $comparator->compare($a, $b) == -1 \
      or warn "Expected $a < $b, got ", $comparator->format_comparison($a, $b);

=cut


1;

__END__



# End of file
