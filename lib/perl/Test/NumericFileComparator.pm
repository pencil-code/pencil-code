#
#                       NumericFileComparator.pm
#                       ------------------------
#
# Description:
#   Compare numbers in two files testing the relative and absolue errors.
#   The required accuracy can be defined in the data files, specified by
#   the caller, or inferred from the data.
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#


package Test::NumericFileComparator;

use warnings;
use strict;
use Carp;
use Math::Complex;
use Test::NumberComparator;

use vars qw($VERSION);

use feature 'say';

##use critic

$VERSION = '0.1';


=head1 NAME

Test::NumericFileComparator - Compare numbers in files.

=head1 SYNOPSIS

  use Test::NumericFileComparator;

  # Infer required accuracy from the reference file
  my $comparator = Test::NumericFileComparator->new('reference.out');

  # Specify accuracy
  my $comparator = Test::NumericFileComparator->new({
      file => 'reference.out',
      relative
  );

  # Compare file to reference data
  if ($comparator->check_file('actual.out')) {
      say 'ok';
  } else {
      say "not ok: %s", $comparator->message();
  }

=head1 DESCRIPTION

Test::NumericFileComparator compares numerical data from files.
The input file can be either in column format:

  # [Maybe
  #  some comment lines]
  [#:names: var1 var2 var3 [...]]
  [#--------var1-var2-var3-[...]]
  [#:accuracies: acc1 acc2 acc3 [...]]
  <var1_val1> <var2_val1> <var3_val1> [...]
  <var1_val2> <var2_val2> <var3_val2> [...]
  [...]

e.g.

  1.0000000  1.0000
  0.5000000  0.6667
  0.8333333  0.8667
  0.5833333  0.7238
  0.7833333  0.8349
  0.6166667  0.7440
  0.7595238  0.8209
  0.6345238  0.7543

or

  # Reference data for pressure and temperature.
  # We print 4 decimals for temperature, but due to the phase of Jupiter's
  # moons, we can only expect an accuracy of about 1e-2.
  #
  #:name:     Pressure   Temperature
  #:accuracy: 1.5e-7     0.8e-2
  # --------------------------------
              1.0000000  1.0000
              0.5000000  0.6667
              0.8333333  0.8667
              0.5833333  0.7238
              0.7833333  0.8349
              0.6166667  0.7440
              0.7595238  0.8209
              0.6345238  0.7543

or

  #--it-----t----urms----umax----rhom----ssm----dtc---dtu--
      0    0.00 0.0063  0.0956 14.4708 -0.4460 0.978 0.025
     10    0.07 0.0056  0.0723 14.4708 -0.4464 0.978 0.019
     20    0.14 0.0053  0.0471 14.4709 -0.4467 0.978 0.019
     30    0.20 0.0056  0.0413 14.4708 -0.4471 0.978 0.017
     40    0.27 0.0058  0.0449 14.4708 -0.4475 0.978 0.013


Alternativly, the input file can be in line format:

  # [Maybe
  #  some comment lines]
  var1: value1 [acc1] [# comment]
  var2: value2 [acc2] [# comment]
  [...]

e.g.

  Pressure:    0.6345238
  Temperature: 0.7543

or

  # Reference data for pressure and temperature.
  # We print 4 decimals for temperature, but due to the phase of Jupiter's
  # moons, we can only expect an accuracy of less than 1e-2.
  Pressure: 0.6345238  1.5e-7  # should be 1.3e-7, but that doesn't work yet
  Temperature: 0.7543  0.8e-2

Note that the accuracies acc1, etc. are absolute accuracies.


=head2 Methods

=over 4

=cut

# ---------------------------------------------------------------------- #
##
## Object constructor
##

=item B<Test::NumericFileComparator-E<gt>new>('I<reference_file>')

Create a new object, reading the reference data from I<reference_file>.
The reference data define the variable list, the accuracies, and the
count of data values.

=cut

sub new {
#
#   Test::NumericFileComparator->new(file);
#   $reader->new(file);
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

    if (@argv == 1) {
        $self->{FILE}  = pop(@argv);
    } else {
        croak("Usage: Test::NumericFileComparator->new(file)");
    }

    $self->{DEBUG} = 0;

    ($self->{VARS}, $self->{VALUES}, $self->{ACCS}) = _parse($self->{FILE});

    bless($self, $class);
    return($self);
}

# ====================================================================== #

##
## Methods
##

=item B<$reader-E<gt>compare>($file)

Compare the given file to the reference data.

The data file may define extra columns or values compared to the reference
data, but the comparison fails if any of the expected columns or values is
missing.

Return an array of descriptive error messages with one entry for each
variable that differs.

=cut

sub compare {
    my $self = shift();
    my ($file) = @_;

    my %problems;
    foreach my $var (@{$self->{VARS}}) {
        $problems{$var} = [];
    }

    my (undef, $val_ref, undef) = _parse($file);
    my %reference = %{$self->{VALUES}};
    my %actual = %$val_ref;
    VAR: foreach my $var (@{$self->{VARS}}) {
        my @column_ref = @{$reference{$var}};
        my @column_act = @{$actual{$var}};
        if (@column_act < @column_ref) {
            push @{$problems{$var}},
                sprintf("Expected %d values, got %d", );
            next VAR;
        }

        my $abs_acc = $self->{ACCS}->{$var};
        my $rel_acc = 0.0;
        my $comparator = Test::NumberComparator->new($abs_acc, $rel_acc);
        for (my $i = 0; $i < @column_ref; $i++) {
            my $ref = $column_ref[$i];
            my $act = $column_act[$i];
            my $comparison = $comparator->compare($ref, $act);
            if ($comparison != 0) {
                my $diagnostic = $comparator->format_comparison($ref, $act);
                if (@column_ref > 1) {
                    $diagnostic = sprintf("Row %d: %s", ($i + 1), $diagnostic);
                }
                push @{$problems{$var}}, $diagnostic;
            }
        }
    }

    my @summary;
    foreach my $var (@{$self->{VARS}}) {
        my $prob = $problems{$var};
        if (@$prob) {
            my $message;
            if ($self->{DEBUG}) {
                $message = join("\n  ", @$prob);  # all diffs for this $var
            } else {
                $message = $prob->[0];  # only the first diff
            }
            push @summary, "$var: $message";
        }
    }

    return @summary;
}

# ====================================================================== #

##
## Private utility subroutines:
##

my $cfloat = '([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?'; # regexp
                                                                   # for C float
my $ieee_float = "(?:$cfloat|[+-]?(?:NaN|Inf))"; # C float | ±NaN | ±Inf


sub _parse {
#
# Parse the given data file and return a list
#   (\@variables, \%values, \%accuracies)
# where @variables lists the variable names in the order in which they are
# defined in the file, %values is a map var_name => [values], and
# %accuracies is a map var_name => accuracy.
#
    my ($file) = @_;

    croak "No such file: $file\n" unless -e $file;
    croak "Cannot open file $file\n" unless -r $file;

    return _read_file_in_column_format($file);
}


sub _read_file_in_column_format {
#
# Read and parse a file in column format (e.g. some time_series.dat file).
# Return
#   (\@variables, \%values, \%accuracies)
#
    my ($file) = @_;

    my $fh;
    unless (open($fh, '<', $file)) {
        croak "Cannot open $file for reading: $!\n";
    }

    my (@variables, %values, %accuracies);

    my @columns;

    while (defined(my $line = <$fh>)) {
        next if $line =~ /^\s*$/;  # empty line
        next if $line =~ /^\s*#/;  # comment line: skip for now
        chomp($line);
        if ($line =~ /^\s*((?:$ieee_float(\s|$)+)+)$/) {
            my @items = split(/\s+/, $1);
            @variables = _create_var_names(@items) unless @variables;
            if ($#items != $#variables) {
                croak(
                    "Inconsistent column count:",
                    " expected ", 0 + @variables,
                    ", found ", 0 + @items,
                    " in line <$line>\n"
                );
            }
            for (my $i = 0; $i < @items; $i++) {
                my $var = $variables[$i];
                $columns[$i] = [] unless defined $columns[$i];
                push $columns[$i], $items[$i];
            }
        } else {
            croak "File $file: Unexpected line in column format: <$line>\n";
        }
    }
    close $fh;

    %accuracies = _infer_accuracies(\@variables, \@columns) unless (%accuracies);
    %values = _numerical_values(\@variables, \@columns);

    return (\@variables, \%values, \%accuracies);
}


sub _create_var_names {
#
# Create variable names by enumeration
#
    my (@items) = @_;

    my @var_names;
    for (my $i = 0; $i < @items; $i++) {
        push @var_names, "variable_$i";
    }
    return @var_names;
}


sub _infer_accuracies {
#
# For each data column, look at value strings and extract an accuracy as
# follows:
# - Select the maximl value (by modulus)
# - Extract the mantissa part
# - Choose 1.5 digits in the last decimal as (absolute) accuracy
#
# This accuracy allows a difference of one in the last decimal, provided
# we the precision of double arithmetics (in Perl) has at least one bit
# more than corresponds to the numbers in the column.
#
# Note that the current way of determining the accuracy from the largest
# value will only work well for static formats (e, d, f); for the dynamic
# g format, we will often underestimate the accuracy. This could be
# improved by analyzing all value strings for a given column.
#
    my ($var_ref, $col_ref) = @_;
    my @variables = @$var_ref;
    my @columns = @$col_ref;

    my %accuracies;
    for (my $i = 0; $i < @variables; $i++) {
        my $var = $variables[$i];
        my $max_string = _max_abs_as_string($columns[$i]);
        $accuracies{$var} = _accuracy_from_num_string($max_string) * 1.5;
    }

    return %accuracies;
}


sub _max_abs_as_string {
#
# Given an array ref of numeric strings, return the element corresponding
# to the largest modulus value.
#
    my ($arr_ref) = @_;
    my @strings = @$arr_ref;

    my $max_val = -$Math::Complex::Inf;
    my $max_string = '<UNDEFINED>';
    foreach my $string (@strings) {
        my $abs_val = abs(0 + $string);
        if ($abs_val > $max_val) {
            $max_string = $string;
            $max_val = $abs_val;
        }
    }
    return $max_string;
}


sub _accuracy_from_num_string {
#
# Return the accuracy of the given numeric string, i.e. a step of 1 in the
# last decimal.
#
# Note that trailing zeros in the mantissa are not treated sepcial, i.e. a
# value of 100 is considered equal to 1.00e2, not 1.e-2.
#
# NaN or Inf values are not handled yet.
#
    my ($num_string) = @_;

    $num_string =~ /^
                    \s *
                    ( ?<mantissa> [-+.0-9] + )
                    ( [eEdD]
                      ( ?<exponent>[-+]?[0-9] * )
                    )?
                    \s *
                    $
                   /x;
    my $exponent = $+{exponent};
    $exponent = 0 unless defined($exponent);
    my $mantissa = $+{mantissa};
    $mantissa =~ /^
                  ( ?<int> [-+0-9] * )
                  ( \.
                    ( ?<frac> [0-9] * )
                  ) ?
                  $
                 /x;
    my $frac = $+{frac};
    $frac = '' unless defined $frac;

    return 10**(-length($frac) + $exponent);
}


sub _numerical_values {
#
# Given a list ref of variable names and a 2d list ref of value strings,
# return the map of numerical values
#   {
#       'var1' => [ var1_val1, var1_val2, ... ],
#       'var2' => [ var2_val1, var2_val2, ... ],
#       ...
#   };
#
    my ($var_ref, $col_ref) = @_;
    my @variables = @$var_ref;
    my @columns = @$col_ref;

    my %values;
    for (my $i = 0; $i < @variables; $i++) {
        my $var = $variables[$i];
        my @col = @{$columns[$i]};
        $values{$var} = [];
        foreach my $val_string (@col) {
            push @{$values{$var}}, 0 + $val_string;
        }
    }

    return %values
}


# ---------------------------------------------------------------------- #

1;

__END__


=back


=head1 BUGS AND LIMITATIONS

=over 4

=item *

None worth mentioning (so far).

=back


=head1 AUTHOR

Wolfgang Dobler <wdobler#cpan:org =~ tr/:#/.@/>


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2007, Wolfgang Dobler.
All rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same conditions as Perl or under the GNU General Public
License, version 3 or later.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=cut


# End of file
