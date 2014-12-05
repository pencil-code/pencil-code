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

# TODO:
# - We currently assume that either no or all accuracies are specified.
#   Need to implement the mixed case.

package Test::NumericFileComparator;

use warnings;
use strict;
use Carp;
use Math::Complex;
use Test::NumberComparator;

use vars qw($VERSION);

use 5.010;                      # for the ~~ operator
use feature 'say';

##use critic

$VERSION = '0.1';


=head1 NAME

Test::NumericFileComparator - Compare numbers in files.

=head1 SYNOPSIS

  use Test::NumericFileComparator;

  # Infer required accuracy from the reference file
  my $comparator = Test::NumericFileComparator->new('reference.out');

  # Compare file to reference data
  my @message = $comparator->compare('actual.out');
  if (@message) {
      say "not ok: @message";
  } else {
      say 'ok';
  }

=head1 DESCRIPTION

Test::NumericFileComparator compares numerical data from files.
The input file can be either in column format:

  [# Maybe]
  [#  some comment lines]
  [#:name: var1 var2 var3 [...]]
  [#--------var1-var2-var3-[...]]
  [#:accuracy: acc1 acc2 acc3 [...]]
  <var1_val1> <var2_val1> <var3_val1> [...]
  <var1_val2> <var2_val2> <var3_val2> [...]
  [...]

e.g.

  1.0000000    1.0000
  0.5000000    0.6667
  0.8333333    0.8667
  0.5833333    0.7238
  0.2833333    0.8349
  0.1066667    0.7440
  0.04759523   0.8209
  0.01345238   0.7942
  0.003124355  0.7734

or

  # Reference data for pressure and temperature.
  # We print 4 decimals for temperature, but due to the phase of Jupiter's
  # moons, we can only expect an accuracy of about 1e-2.
  #
  #:name:     Pressure     Temperature
  #:accuracy: 5.0e-6:r     0.8e-2
  # ----------------------------------
              1.0000000    1.0000
              0.5000000    0.6667
              0.8333333    0.8667
              0.5833333    0.7238
              0.2833333    0.8349
              0.1066667    0.7440
              0.04759523   0.8209
              0.01345238   0.7942
              0.003124355  0.7734

or

  #--it-----t----urms----umax----rhom----ssm----dtc---dtu--
      0    0.00 0.0063  0.0956 14.4708 -0.4460 0.978 0.025
     10    0.07 0.0056  0.0723 14.4708 -0.4464 0.978 0.019
     20    0.14 0.0053  0.0471 14.4709 -0.4467 0.978 0.019
     30    0.20 0.0056  0.0413 14.4708 -0.4471 0.978 0.017
     40    0.27 0.0058  0.0449 14.4708 -0.4475 0.978 0.013


Alternativly, the input file can be in line format:

  [# Maybe
  [# some comment lines]
  var1 : [acc1 :] <var1_val1> [<var1_val2> ... ] [# comment]
  var2 : [acc2 :] <var2_val1> [<var2_val2> ... ] [# comment]
  [...]

e.g.

  Pressure:    0.6345238
  Temperature: 0.7543

or

  # Reference data for pressure and temperature.
  # We print 4 decimals for temperature, but due to the phase of Jupiter's
  # moons, we can only expect an accuracy of about 1e-2.
  Pressure      : 1.5e-7:r : 0.6345238 0.7345238 0.655238  # slice (5,7,3:5)
  Temp(5,7,3:5) : 0.8e-2   : 0.7543 0.7411 0.7123

Note that the colon (:) as separator must be followed by at least one
space. Thus,

  Temp(5,7,3:5): 0.8e-2: 0.7543 0.7411 0.7123

is parsed correctly, but

  Temp(5,7,3:5) :0.8e-2 :0.7543 0.7411 0.7123

is not, nor is

  Temp(5,7,3: 5) : 0.8e-2 : 0.7543 0.7411 0.7123


For both formats, the accuracies acc1, etc. are absolute accuracies if
they are plain numbers or followed by ':a' without a space (e.g. '1.2e-2'
or '1.2e-3:a'). Numbers followed by ':r' indicate relative precision (e.g.
'3.5e-2:r'). A combination of both is possible ('1.2e-3:a|3.5e-2:r'); in
this case, two numbers are considered sufficiently equal if they are close
enough by either absolute or relative accuracy.

If no accuracy is specified, for each column / line, an absolute accuracy
of 1.5 times the last digit of the largest number (by modulus) is used.


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

Compare the given file to the reference data, using the accuracies defined
by the reference data (i.e. any accuracies specified in $file are ignored).

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

        my ($abs_acc, $rel_acc) = @{$self->{ACCS}->{$var}};
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
my $ieee_float = "(?:$cfloat|[+-]?(?:NaN|Inf))"; # C float | +/-NaN | +/-Inf


sub _parse {
#
# Parse the given data file and return a list
#   (\@variables, \%values, \%accuracies)
# with
#   @variables  -- list of variables names
#   %values     -- map variable_name => [value1, value2, ...]
#   %accuracies -- map variable_name => [abs_acc, rel_acc]
#
    my ($file) = @_;

    croak "No such file: $file\n" unless -e $file;
    croak "Cannot open file $file\n" unless -r $file;

    if (_has_line_format($file)) {
        return _read_file_in_line_format($file);
    } else {
        return _read_file_in_column_format($file);
    }
}


sub _has_line_format {
#
# Does the given file look like it is in line format?
# We simply look for a ':' character in a non-comment line.
#
    my ($file) = @_;

    my $fh;
    unless (open($fh, '<', $file)) {
        croak "Cannot open $file for reading: $!\n";
    }

    my $line_format = 0;
    while (defined(my $line = <$fh>)) {
        next if $line =~ /^\s*$/;  # empty line
        next if $line =~ /^\s*#/;  # comment line
        if ($line =~ /:/) {
            $line_format = 1;
            last;
        }
    }
    close $fh;
    return $line_format;
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

    my (@variables, %values);

    my (@columns, @accuracies);

    my $linenum = 0;
    while (defined(my $line = <$fh>)) {
        $linenum++;
        chomp($line);
        next if $line =~ /^ \s * $/x;  # empty line
        if ($line =~ /^
                      \s *
                      \#
                      (?:
                          - +
                          [^ - ] +
                      ) +
                      - *
                      $
                     /x) {
            @variables = _split_line($line, qr'^#[-]+', qr'[-]+');
            next;
        }

        if ($line =~ /^
                      \s *
                      \#:name:
                      \s +
                     /x) {
            @variables = _split_line($line, qr'^#:name:\s+', qr'\s+');
            next;
        }

        if ($line =~ /^
                      \s *
                      \#:accuracy:
                      \s +
                     /x) {
            @accuracies = _extract_accuracies($line);
            next;
        }

        next if $line =~ /^\s*#/;  # comment line: skip for now
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
                push @{$columns[$i]}, $items[$i];
            }
        } else {
            croak "File $file: Unexpected line $linenum in column format: <$line>\n";
        }
    }
    close $fh;

    my %accuracy_map;
    for (my $i=0; $i < @accuracies; $i++) {
        $accuracy_map{$variables[$i]} = $accuracies[$i];
    }
    return _extract_data(\@variables, \@columns, \%accuracy_map, $file);
}


sub _read_file_in_line_format {
#
# Read and parse a file in line format (one variable per line).
# Return
#   (\@variables, \%values, \%accuracies)
#
    my ($file) = @_;

    my $fh;
    unless (open($fh, '<', $file)) {
        croak "Cannot open $file for reading: $!\n";
    }

    my (@variables, %values, %accuracies);
    my %known_vars;

    my @columns;                # in fact rows, but stick to names used above

    my $linenum = 0;
    while (defined(my $line = <$fh>)) {
        $linenum++;
        next if $line =~ /^\s*$/;  # empty line
        next if $line =~ /^\s*#/;  # comment line
        $line =~ s{#.*}{};  # strip trailing comments
        chomp($line);
        if ($line =~ /^
                      \s *
                      ( . +? )
                      \s *
                      : \s +

                      (?:
                          ( . +? )
                          \s *
                          : \s +
                      )?

                      (
                          (?:
                              $ieee_float
                              (?: \s +? | $ )
                          ) +
                      )
                      $
                     /x) {
            my ($var, $acc) = ($1, $2);
            my @vals = split('\s+', $3);
            if (defined $known_vars{$var}) {
                croak "Duplicate variable '$var' in file $file:$linenum\n";
            }

            push @variables, $var;
            $known_vars{$var}++;

            push @columns, \@vals;

            if (defined $acc) {
                $accuracies{$var} = _parse_accuracy($acc);
            }
        } else {
            croak "File $file: Unexpected line $linenum in line format: <$line>\n";
        }
    }
    close $fh;

    return _extract_data(\@variables, \@columns, \%accuracies, $file);
}


sub _extract_data {
    my ($vars_ref, $cols_ref, $accs_ref, $file) = @_;

    croak("No data found in $file") unless @$cols_ref;

    my %accuracies;
    if (%$accs_ref) {
        %accuracies = %$accs_ref;
    } else {
        %accuracies = _infer_accuracies($vars_ref, $cols_ref);
    }
    my %values = _numerical_values($vars_ref, $cols_ref);

    return ($vars_ref, \%values, \%accuracies);
}


sub _extract_accuracies {
#
# Extract accuracies from a header line like
#   #:accuracy:  1  1e-3  1e-3:r  1e-4  1e-6:a|1e-3:r  1e-3:r
# Returns a list of pairs [abs_acc, rel_acc], e.g.
#   ([1, 0], [1e-3, 0], [0, 1e-3], [1e-4, 0], [1e-6, 1e-3], [0, 1e-3])
#
    my ($line) = @_;

    my @strings = _split_line($line, qr'^#:accuracy:\s+', qr'\s+');
    return map { _parse_accuracy($_) } @strings;
}


sub _split_line {
#
# Remove from a line the start sequence and extract the individual
# elements separated by $separator.
# Return the list of elements (as strings).
#
    my ($line, $start, $separator) = @_;

    chomp $line;
    if ($line !~ m{^$start}) {
        croak("Line does not start with $start: $line");
    }
    $line =~ s{$start (?: $separator) ?}{}x;
    $line =~ s{(?: $separator) $}{}x;
    return split($separator, $line);
}


sub _parse_accuracy {
#
# Parse an accuracy string like '1e-3', '1e-3:r', or '1e-6:a|1e-3:r' to a
# pair of numbers like (), (), or ()
    my ($string) = @_;

    my ($abs, $rel) = (0, 0);
    foreach my $acc (split(/\|/, $string)) {
        if ($acc =~ m{^ ($ieee_float) : r $}x) {
            $rel = $1;
        } elsif ($acc =~ m{^ ($ieee_float) (?: :a)? $ }x) {
            $abs = $1;
        } else {
            croak("Cannot parse accuracy string '$acc'");
        }
    }

    return [$abs, $rel];
}


sub _create_var_names {
#
# Create variable names by enumeration
#
    my @items = @_;

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
        $accuracies{$var} = [
            _accuracy_from_num_string($max_string) * 1.5,  # absolute
            0                                              # relative
            ];
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
# Note that trailing zeros in the mantissa are not treated special, i.e. a
# value of 100 is considered equal to 1.00e2, not 1.e-2.
#
# For values NaN or Inf, return 0.
#
    my ($num_string) = @_;

    if ($num_string =~ m{(inf|nan)}i) {
        return 0.0;
    }

    $num_string =~ /^
                    \s *
                    ( [-+.0-9] + )
                    (?: [eEdD]
                        ( [-+]?[0-9] * )
                    )?
                    \s *
                    $
                   /x;
    my $exponent = $2;
    $exponent = 0 unless defined($exponent);
    my $mantissa = $1;
    $mantissa =~ /^
                  ( [-+0-9] * )
                  (?: \.
                      ( [0-9] * )
                  ) ?
                  $
                 /x;
    # my $int = $1;  # not needed
    my $frac = $2;
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



# End of file
