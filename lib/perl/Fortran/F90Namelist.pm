#
#                            F90Namelist.pm
#                            --------------
#
# Description:
#   Parse F90 namelist into a hash and export in different formats.

package Fortran::F90Namelist;

=head1 NAME

Fortran::F90Namelist - Parse F90 namelists into hash and export in different formats

=head1 SYNOPSIS

  use Fortran::F90Namelist;
  my $nl = Fortran::F90Namelist->new() or die "Couldn't get object\n";

  $nl->parse("&runpars\nx=2,y=3\nvec1=1,2,3\nvec2=3*1.3\n/");

  # Operate on each namelist in $text (only works with [mutable]
  # strings, not with files)
  my $text = "&spars\nx=2,y=3\n/\n&runpars\nvec1=1,2,3\nvec2=3*1.3\n/";
  while ($nl->parse($text)) {
      print $nl->name(), "\n";
  }

Dump in arbitrary order:

  use Data::Dumper;
  print "F90Namelist ", $nl->name(), " has ", $nl->nslots(), " slots:\n";
  print Dumper($nl->hash());

Retain original order:

  print "&",$nl->name(),"\n";
  my $nl_hash = $nl->hash();
  foreach my $var (@{$nl->slots()}) {
    print "  $var: ", Dumper($nl_hash->{$var});
  }
  print "/\n";

Read from file:

  # Read one namelist from file `one_list.nml'
  $nl->parse(file => 't/files/one_list.nml');

  # Read one namelist from file handle
  open(my $fh , "< t/files/one_list.nml") or die "Couldn't get file handle\n";
  $nl->parse(file => $fh);
  # or
  open(NAMELIST , "< t/files/one_list.nml") or die "Couldn't open file\n";
  $nl->parse(file => \*NAMELIST);

Read all namelists from file `start.in' and merge into one namelist
called `nlist'

  $nl->parse(file     => 't/files/start.in',
             all      => 1,
             namelist => 'nlist');
  print "Merged namelist ", $nl->name, " contains:\n",
      join(",  ", @{$nl->slots}), "\n";

Merge two namelists

  my $nl2 = Fortran::F90Namelist->new() or die "Couldn't get object\n";
  $nl2->parse(file => 't/files/one_list.nml');
  $nl->merge($nl2,
             { dups_ok => 1 } );
  print $nl->name, " now has ", $nl->nslots, " slots\n";


Write namelist:

  # Write namelist in F90 namelist format
  print "F90 format:\n", $nl->output();

  # Write namelist as Python class
  print "Python format:\n", $nl->output(format => 'python', name => 'par2');

  # Write namelist as IDL structure
  print "IDL format:\n", $nl->output(format => 'idl', name => 'par2');


=head1 DESCRIPTION

Fortran::F90Namelist is a module for parsing Fortran90 namelists into hashs and
re-exporting these hashs in different formats. Currently, the following
data types are supported:

=over 4

=item *

integer

=item *

float/double

=item *

complex numbers

=item *

strings [character(LEN=*)], possibly containing all sorts of quotation
marks

=item *

logical

=back

The following backends exist for re-exporting (or importing into other
languages):

=over 4

=item *

F90 namelist

=item *

Python class

=item *

IDL struct

=back

This module is used with the I<Pencil Code>
(L<http://www.nordita.dk/software/pencil-code/>) to import the values of
all available input parameters into Python, GDL/IDL or other visualization
software.

=head2 Methods

=over 4


=item B<$nl-E<gt>new>()

Create a new namelist object


=item B<$nl-E<gt>parse>(I<string>)

=item B<$nl-E<gt>parse>(text => I<string>)

=item B<$nl-E<gt>parse>(file =>(I<fname>|I<FHANDLE>))

=item B<$nl-E<gt>parse>(file => (I<fname>|I<FHANDLE>) [, I<options> ])

=item B<$nl-E<gt>parse>(\%options)

Parse I<string> or the file represented by I<fname> or I<FHANDLE> (a file
handle), returns the name of the namelist parsed, or undef if parsing
failed.

When reading from a mutable text string $text, the string is modified and
contains everything following the namelist just parsed.

This allows C<while> loops like

  while ($nl->parse($text)) {
      print $nl->name(), "\n";
  }

to work.
This does however not work for files or immutable strings, so

=for test ignore

  while ($nl->parse(file => "t/files/start.in")) {
      print $nl->name(), "\n";
  }

=for test

and

=for test ignore

  while ($nl->parse("&nl1\nx=5.\n/\n&nl2\n/")) {
      print $nl->name(), "\n";
  }

=for test

will fail.

Generally speaking, L<Fortran::F90Namelist::Group|Fortran::F90Namelist::Group>
is the more appropriate tool for handling several namelists in one file or
string.


Additional I<options> are:

=over 8

=item B<merge>

If true, merge data from namelist with any data that may already be
stored in the object.
See L<Fortran::F90Namelist::Group|Fortran::F90Namelist::Group> for a more
flexible framework for dealing with groups of namelists.

=item B<all>

If true, parse all namelists from string or file and merge them into one
namelist object.

=item B<name>

Set name of resulting namelist (default: name of first namelist read).

=item B<dups_ok>

With B<merge>, don't warn if new slots have same names, but different
values as existing slots.

=item B<broken>

Try to parse broken namelists as produced by ifc 7.x, where you can get
something like

=for test ignore

   COOLING_PROFILE='gaussian              ',COOLTYPE='Temp    
   'COOL= 0.0,CS2COOL= 0.0,RCOOL= 1.000000

=for test

if the closing quote for a string (`Temp    ') would end up in column 81.

All options can be passed in a hash(-ref):

  my %options = ( file   => 't/files/one_list.nml',
                  name   => 'broken_nlist',
                  broken => 1 );
  $nl->parse(%options);
  $nl->parse(\%options);  # the same

=back


=item B<$nl-E<gt>merge>($nl2 [, I<options>])

Merge namelist object $nl2 into $nl.

I<Options> are:

=over 8

=item B<name>

Set name of resulting namelist (default: name of $nl).

=item B<dups_ok>

With B<merge>, don't warn if new slots have same names, but different
values as existing slots.

=back


=item B<$nl-E<gt>name>()

=item B<$nl-E<gt>name>($newname)

Return or set name of namelist.


=item B<$nl-E<gt>nslots>()

Return number of slots in namelist


=item B<$nl-E<gt>slots>()

Return ref to list of variable (slot) names in original order


=item B<$nl-E<gt>hash>()

Return namelists as Perl hashref.
See L<HASH FORMAT|/"HASH FORMAT"> below for details of the hash format.


=item B<$nl-E<gt>output>([options])

Write namelist in given I<format>.

Options are

=over 8

=item B<format>=I<format>

Set the output format.
Currently supported formats are `f90' (default), 'python', and `idl'.

=item B<name>=I<name>

Set the name of the namelist (default: C<$nl-E<gt>name>()).

=item B<trim>

Trim all trailing whitespace from strings.

=item B<double>

Write all floating point numbers as double precision numbers.

=item B<oneline>

Print whole namelist in one line (if compatible with the output format).

=item B<maxslots>=I<N>

Print only N slots per line.
Useful for programs like IDL that have restrictions on the length of lines
read from a pipe, so B<oneline> is dangerous.

=back

=back


=head1 HASH FORMAT

The B<hash> method returns a hash reference of the following structure:

=for test ignore

    { 'name of var1' => { 'value' => [ value1, value2, ..],
                          'type'  => numerical_type,
                          'stype' => "type string"
                        },
      'name of var2' => { 'value' => [ value1, value2, ..],
                          'type'  => numerical_type
                          'stype' => "type string"
                        },
      ...
    }

=for test

Here I<numerical_type> is a number identifying each data type, while
I<stype> is a textual description of the given data type.

E.g.

=for test ignore

    { 'xyz0' => { 'value' => [ 0., -3.141593, 0. ],
                  'type'  => 6,
                  'stype' => 'single precision float'
                },
      'nt'   => { 'value' => [ '1000' ],
                  'type'  => 4,
                  'stype' => 'integer'
                }
    }

=for test

Note: This is currently just the internal format used to represent
namelists and can thus change in the future.
In particular the C<type> numbers should not considered to be stable
between releases.


=head1 TO DO

=over 4

=item 1.

new(), parse(), output(), etc. should check for unknown args and complain,
not silently ignore them as is currently the case.

=item 2.

More output methods:

=over 8

=item *

Octave/matlab , C structs, YAML, XML(?), ...

=back

=back


=head1 BUGS AND LIMITATIONS

=over 4

=item *

No user-defined types (records) are supported, so if you have these LaTeX
comment characters in your namelist data, you are out of luck.

=back


=head1 AUTHOR

Wolfgang Dobler <wdobler [at] cpan.org>


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2014, Wolfgang Dobler <wdobler [at] cpan.org>.
All rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same conditions as Perl or under the GNU General Public
License, version 2 or later.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=head1 SEE ALSO

L<Fortran::Namelist> by Victor Marcello Santillan.
That module has a more limited scope (reading a namelist group from file,
inserting namelists, and writing the resulting group to another file [my
interpretation]), but is way faster on large files.

=cut


use warnings;
use strict;
use Carp;
use vars qw($VERSION);

# Cannot use Perl5.8's constant { x => 1, y=>2 , ..} because 5.6
# is very popular still
#
# Possible states of parser [used at all?]
##no critic ProhibitConstantPragma
use constant  UNDEF   => -1;
use constant  START   =>  0;    # initial state of parser
use constant  VAR     =>  1;    # at beginning of variable name
use constant  VALUE   =>  2;    # at beginning of value
use constant  SQUOTE  =>  3;    # in string after opening single quote
use constant  DQUOTE  =>  4;    # in string after opeing double quote
use constant  BRACKET =>  5;    # after opening bracket (e.g. complex number)
use constant  COMMENT =>  6;    # after exclamation mark (F90 comment)
use constant  NL_END  =>  7;    # after closing `/'
#
# F90 data types
use constant  UNKNOWN   => 0;
use constant  SQ_STRING => 1;
use constant  DQ_STRING => 2;
use constant  LOGICAL   => 3;
use constant  INTEGER   => 4;
use constant  FLOAT     => 5;   # a float here can be single or double
use constant  SINGLE    => 6;
use constant  DOUBLE    => 7;
use constant  COMPLEX   => 8;
use constant  DCOMPLEX  => 9;
use constant  MULTIPLE  => 20;
use constant  EMPTY     => 21;
#
use constant  ID        => 100; # variable name (_not_ a data type)

##use critic

$VERSION = '0.6.1';

## Regexps for integer and floating-point numbers
# general float:
my $numeric   = qr/(?:[+-]?)(?=\d|\.\d)\d*(?:\.\d*)?(?:[EeDd](?:[+-]?\d+))?/;
# float:
my $numeric_e = qr/(?:[+-]?)(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee](?:[+-]?\d+))?/;
# double:
my $numeric_d = qr/(?:[+-]?)(?=\d|\.\d)\d*(?:\.\d*)?(?:[Dd](?:[+-]?\d+))?/;
# float with decimal point, but w/o exponential part:
my $float     = qr/(?:[-+]?(?:\d+\.\d*|\d*\.\d+))/;

## Extend floating-point numeric tpes by one- or two-point
## compactification of real numbers (any mathematicians here?), aka IEEE
## denormalized numbers (for engineers):
my $NaN = qr/NaN/;
my $Inf = qr/(?:[-+]?)Inf/;
my $ieee_denorm = qr/(?:$NaN|$Inf)/;
#$numeric_e = qr/(?:$numeric_e|$ieee_denorm)/;
#$numeric_d = qr/(?:$numeric_d|$ieee_denorm)/;
$numeric   = qr/(?:$numeric|$ieee_denorm)/;
$float     = qr/(?:$float|$ieee_denorm)/;

## Regexps for the different data type values. Make sure all brackets are
## marked grouping-but-non-capturing (?:...), or else the parsing
## algorithm will fail.
my @regexp;
$regexp[SQ_STRING] = qr/'(?:[^']|'')*'/; # even covers 'toto''s quote'
$regexp[DQ_STRING] = qr/"(?:[^"]|"")*"/; # ditto for double quotes
$regexp[DCOMPLEX]  = qr/\(\s*$numeric_d\s*,\s*$numeric_d\s*\)/;
$regexp[COMPLEX]   = qr/\(\s*$numeric\s*,\s*$numeric\s*\)/;
$regexp[LOGICAL]   = qr/(?:T|F|\.(?:true|TRUE|false|FALSE)\.)/;
$regexp[MULTIPLE]  = qr/[0-9]+\*/; # also need special treatment...
$regexp[INTEGER]   = qr/[+-]?[0-9]+/;
$regexp[DOUBLE]    = qr/$numeric_d/;
$regexp[SINGLE]    = qr/$numeric_e/;
$regexp[FLOAT]     = qr/$float/;
$regexp[ID]        = qr/[a-zA-Z](?:[a-zA-Z0-9_%])*/; # allowed namelist/var. names

## Corresponding regexp for compatible type class (numeric, complex, ...)
my @regexp2 = @regexp;          # same regexp by default
$regexp2[DCOMPLEX]  = qr/\(\s*$numeric\s*,\s*$numeric\s*\)/;
$regexp2[COMPLEX]   = qr/\(\s*$numeric\s*,\s*$numeric\s*\)/;
$regexp2[INTEGER]   = qr/$numeric/;
$regexp2[DOUBLE]    = qr/$numeric/;
$regexp2[SINGLE]    = qr/$numeric/;
$regexp2[FLOAT]     = qr/$numeric/;

# Hash for looking up symbolic names for type constants. The constants are
# only expanded as numbers if adding 0.
my %stypes = ( UNKNOWN   + 0 => 'unknown',
               SQ_STRING + 0 => 'single-quote string',
               DQ_STRING + 0 => 'double-quote string',
               LOGICAL   + 0 => 'logical',
               INTEGER   + 0 => 'integer',
               FLOAT     + 0 => 'unspecified float',
               SINGLE    + 0 => 'single precision float',
               DOUBLE    + 0 => 'double precision float',
               COMPLEX   + 0 => 'complex number',
               DCOMPLEX  + 0 => 'double precision complex number',
               MULTIPLE  + 0 => 'multiple data (array)',
               EMPTY     + 0 => 'empty array data',
             );

# Global variables related to output() method:
my ($cmplx_pref,$cmplx_suff) = ('', ''); # default delimiters for complex nums

# ---------------------------------------------------------------------- #
##
## Object constructor
##
## Internal structure of Namlist objects (update me):
##   DATA    -- variable names, values, and types (hashref, see below)
##   SLOTS   -- ordered list of variable names (array ref)
##   NSLOTS  -- number of slots
##   NAME    -- name of namelist
##   PARSED_ -- flag indicating that argument has been parsed
##   DEBUG_  -- debug flag
##
## Structure of DATA slot: Note: One namelist object holds only one
## namelist -- use {$nl1, $nl2, ..} to group them.
##
##   $self->{DATA} = data_hash;
##   data_hash = { 'name of var1' => { 'value' => [ value1, value2, ..],
##                                     'type'  => numerical_type,
##                                     'stype' => "type string"
##                                   }
##                 'name of var2' => { 'value' => [ value1, value2, ..],
##                                     'type'  => numerical_type
##                                     'stype' => "type string"
##                                   }
##                 ...
##               };
##
sub new {
# Return new F90Namelist object. By default, the object is unparsed and has
# no name. Use `empty => 1' and `name => $name' to create a complete empty
# namelist with name $name'.
#
# my $nl = Fortran::F90Namelist::new();
# my $nl = Fortran::F90Namelist::new(name => 'toto', empty => 1);
# my $nl = Fortran::F90Namelist::new({name => 'toto', empty => 1});
# IMPLEMENT US:
# my $nl = Fortran::F90Namelist::new({file => $filename);
# my $nl = Fortran::F90Namelist::new({text => $text, debug => 1});
#
    my $proto = shift;          # either classref or object ref or string
    my @argv  = @_;
    my $class = ref($proto) || $proto;
    my $self = {};

    my %data   = ();
    my @slots  = ();
    my $nslots = undef;
    my $parsed = 0;

    my $short_usage =
        "Usage:\n" .
        "  Fortran::F90Namelist::new()\n" .
        "  Fortran::F90Namelist::new(name => \$name)\n" .
        "  Fortran::F90Namelist::new({name => \$name})\n" ;

    # Parse argument(s) (name => <nlname>); may be list or hashref
    my %args;
    if (@argv) {
        if (ref($argv[0]) eq 'HASH') { # parse($hashref)
            %args = %{$argv[0]};
        } else {                # parse(%hash) or parse(@list)
            %args = @argv;
        }
    }
    #
    my $name  = ($args{name}  || '' );
    my $empty = ($args{empty} || '' );
    my $debug = ($args{debug} || 0  );

    if ($empty) {               # Valid but empty namelist
        $nslots = 0;
        $parsed = 1;
    }

    ##
    ## Finish the object
    ##
    # public data of object
    $self->{DATA}   = \%data;
    $self->{SLOTS}  = \@slots;
    $self->{NSLOTS} = $nslots;
    $self->{NAME}   = $name;

    # internal data
    $self->{PARSED_} = $parsed;
    $self->{DEBUG_}  = $debug;

    bless($self,$class);
    return($self);
}

# ====================================================================== #

##
##  Methods
##

sub parse {
#
#   $obj->parse($text)
#   $obj->parse(file => $filename)
#   $obj->parse(file => $filehandle)
#   $obj->parse(file => $filename|$filehandle, merge => 1[, name => $name])
#   $obj->parse({file => $filename|$filehandle, merge => 1[, name => $name]})
#   $obj->parse({text => $textstring, merge => 1[, name => $name]})
#
# IMPLEMENT ME:
#   $obj->parse({text => $textstring, debug => 1})
#
# Parse text or file containing F90 namelist(s)
#
    my $self = shift;
    my @args = @_;              # can't use shift() since we change value
                                # of $text

    my $state = START;
    my $debug = $self->{DEBUG_};

    my %args;
    my $text;
    my $textarg = 0;

    # Parse arguments (file => <filename>, etc.); may be single string,
    # list, hash or hashref
    if (ref($args[0]) eq 'HASH') { # parse($hashref)
        %args = %{$args[0]};
    } else {
        if (@_ == 1) {          # parse($string)
            $textarg = 1;
            $text = $args[0];
        } else {                # parse(%hash) or parse(@list)
            %args = @args;
        }
    }
    my $file    = ($args{file}    || '' );
    my $merge   = ($args{merge}   || 0  );
    my $all     = ($args{all}     || 0  );
    my $name    = ($args{name}    || '' );
    my $dups_ok = ($args{dups_ok} || '' );
    my $broken  = ($args{broken}  || 0  );
    my $format  = ($args{format}  || 'f90' );
    my $double  = ($args{double}  || 0 );

    if (lc($format) eq 'idl' || lc($format) eq 'python' ) {
        if ($double) {
            if (lc($format) eq 'python') {
              $cmplx_pref = "double complex"; # double complex number prefix
            } else {
              $cmplx_pref = "dcomplex"; # double complex number prefix
            }
        } else {
            $cmplx_pref = "complex"; # complex number prefix
        }
    }

    # Get text from file if necessary
    $text      = ($args{text}  || $text );
    if (!defined($text)) {
        croak "\$nl->parse(): need text or file argument\n"
            unless ($file ne '');
        local $/ = undef;
        if (ref($file) eq 'GLOB') { # file handle
            $text = <$file>;
        } else {                    # file name
            ##no critic ProhibitBarewordFileHandles ProhibitTwoArgOpen
            open(FH, "< $file") or croak "Cannot open file <$file> for reading";
            ##use critic
            $text = <FH>;
            close(FH);
        }
    }

    if ($merge) {
        $name ||= $self->{NAME};  # default to previous name
    } else {                      # clear/reset all data
        $self->{DATA}   = {};
        $self->{SLOTS}  = [];
        $self->{NSLOTS} = 0;
        $self->{NAME}   = $name;
    }

    my $done = 0;

    do {
        my ($name1, $nslots1, @slots1);
        my $href = parse_namelist(\$text, \$name1, \$nslots1, \@slots1,
                                  $broken, $debug, $format);
        if (defined($href)
            && defined($name1)
            && $name1 ne ''
           ) { # really read a namelist

            $name ||= $name1; # choose first name if not set yet

            # Call merge method to do the actual work
            $self->merge([\@slots1, $href],
                         { name    => $name,
                           dups_ok => $dups_ok }
                        );
            $done = 1 unless ($all);
        } else {    # read nothing useful --> get out of this loop
            $done = 1;
        }
    } until ($done);

    # Is there a way to find out whether $text is mutable (i.e. no
    # constant)? Until I find one, just use brute force:
    eval { $_[0] = $text };              # Return remaining $text
    $@ = '';                             # We won't use this anyway..

    # Don't mimic success if we haven't read anything useful
    return if ($name eq '');

    if ($debug) {
        print STDERR
          "Fortran::F90Namelist->parse: Successfully read namelist <$name>\n";
        print STDERR "=================================\n";
    }

    $self->{PARSED_} = 1;

    return $self->{NAME};
}

# ---------------------------------------------------------------------- #

sub merge {
#
#   $obj->merge($obj2);
#   $obj->merge($obj2,   name => $name, $dups_ok => 1  );
#   $obj->merge($obj2, { name => $name, $dups_ok => 1 });
#
# Merge another namelist into this one
#
    my $self = shift();
    my $nl2  = shift();
    my @args = @_;              # remaining argument(s)

    # Arg $nl2 can be a namelist or just a data hashref
    my (@slots2,$hashref2);
    if (ref($nl2) eq 'Fortran::F90Namelist') {
        @slots2   = @{$nl2->{SLOTS}};
        $hashref2 = $nl2->{DATA};
    } elsif (ref($nl2) eq 'ARRAY') {
        @slots2   = @{$$nl2[0]};
        $hashref2 = $$nl2[1];
    } else {
        croak "Fortran::F90Namelist->merge(): "
          . "expected Fortran::F90Namelist object or hashref\n";
    }

    # Parse arguments (name => <name>, etc.); may be hash or hashref
    my %args;
    if (ref($args[0]) eq 'HASH') { # parse($hashref)
        %args = %{$args[0]};
    } else {
        %args = @args;
    }
    my $name    = ($args{name}    || $self->{NAME} );
    my $dups_ok = ($args{dups_ok} || ''            );

    my $nslots  = $self->{NSLOTS};
    my @slots   = @{$self->{SLOTS}};
    my $hashref = $self->{DATA};
    my $debug   = $self->{DEBUG_};

    if ($debug) {
        print STDERR
          "Fortran::F90Namelist->merge: "
            , "Merging ", @slots2 + 0,
              "-slots namelist into $nslots-slots namelist\n";
    }

    # Eliminate repeated slots and warn if values don't agree
  slot: foreach my $slot (@slots2) {
        if (defined($$hashref{$slot})) { # slot already known
            my @val1=@{$$hashref{$slot}{'value'}};
            my @val2=@{$$hashref2{$slot}{'value'}};
            while (@val1 and @val2) {
                my $v1 = pop(@val1);
                my $v2 = pop(@val2);
                if (($v1 ne $v2) && ! $dups_ok) {
                    carp "WARNING: Conflicting slots" .
                      " $slot = [@{$$hashref{$slot}{'value'}}]" .
                        " vs. [@{$$hashref2{$slot}{'value'}}]\n";
                    next slot;
                }
            }
        } else {        # new slot
            push @slots, $slot;
            $$hashref{$slot} = $$hashref2{$slot};
            $nslots++;
        }
    }

    # Wrap it up
    $self->{NAME}   = $name;
    $self->{NSLOTS} = $nslots;
    $self->{SLOTS}  = \@slots;
    $self->{DATA}   = $hashref;

    if ($debug) {
        print STDERR
          "Fortran::F90Namelist->merge: "
          . "Successfully merged into namelist <$name>\n";
        print STDERR "=================================\n";
    }

    $self->{PARSED_} = 1;

    return 1;                   # success
}

# ---------------------------------------------------------------------- #

sub name {
# Get or set name of parsed namelists
    my $self = shift();

    if (@_) { $self->{NAME} = shift };
    return $self->{NAME};
}

# ---------------------------------------------------------------------- #

sub nslots {
# Get number of slots in namelist
    my $self = shift();
    return $self->{NSLOTS}
}

# ---------------------------------------------------------------------- #

sub slots { # FIXME
# Return array ref of variable names in slots
    my $self = shift();
    return $self->{SLOTS}
}

# ---------------------------------------------------------------------- #

sub hash {
# Return hash with parsed namelist contents
    my $self = shift;
    return $self->{DATA};
}

# ---------------------------------------------------------------------- #

sub output {
# Write namelist in specified format (defaults to 'f90')
    my $self = shift();

    # Optional arguments:
    #   format   => format   ('f90' [default], 'python', or 'idl')
    #   name     => nl_name  (name of nlist/struct [default: get from nlist])
    #   trim     => 0/1      (trim trailing whitespace off strings)
    #   double   => 0/1      (mark all floats as double precision)
    #   oneline  => 0/1      (write all in one line? [only for some formats])
    #   maxslots => N        (similar to oneline, but split every N slots)
    my @argv = @_;

    # Parse arguments (file => <filename>, etc.); may be list, hash or hashref
    my %args;
    if (ref($argv[0]) eq 'HASH') {
        %args = %{$argv[0]};
    } else {
        %args = @argv;
    }
    my $format   = ($args{format}   || 'f90');
    my $name     = ($args{name}     || $self->name() || '');
    my $trim     = ($args{trim}     || 1);
    my $double   = ($args{double}   || 0);
    my $oneline  = ($args{oneline}  || 0);
    my $maxslots = ($args{maxslots} || 0);
    $oneline = 0 if ($maxslots);

    # Sanity check
    unless ($self->{PARSED_}) {
        croak "Called method output() on unparsed namelist";
    }

    # Get name of namelist(s?)
    my ($name1,$hashref) = (%{$self->{DATA}}); # hash (name=>valhash) ->
                                               # 2-element array; should
                                               # possibly be a loop over all
                                               # name=>hash pairs?
    # Format-dependent settings
    # We are printing the following:
    # $head_pref
    #   <header>
    # $head_suff
    # $slot_pref
    #   <slot1>
    # $slot_join
    #   <slot2>
    # [...]
    #   <slot_maxslots>
    # $slot_suff
    # [...]
    #   <slotN>
    # $last_suff
    # $foot_pref
    #   <footer>
    # $foot_suff

    my ($header,$footer);
    my ($head_pref,$head_suff);
    my ($slot_pref,$slot_join,$slot_suff,$last_suff);
    my ($foot_pref,$foot_suff);

    my ($newline,$indent);      # to play tricks with $oneline
    if ($oneline) {
        $newline = " ";
        $indent  = "";
    } else {
        $newline = "\n";
        $indent  = "    ";
    }

    if      (lc($format) eq 'f90') {
        $header     = "\&$name";
        $footer     = "/";
        #
        $head_pref  = "";
        $head_suff  = "$newline";
        $slot_pref  = "$indent";
        $slot_join  = ", ";
        $slot_suff  = ",$newline";
        $last_suff  = "$newline";
        $foot_pref  = "";
        $foot_suff  = "\n";
    } elsif (lc($format) eq 'python') {
        $header = "class $name:";
        $footer = "";
        #
        $head_pref = "";
        $head_suff = "$newline";
        $slot_pref = "$indent";
        $slot_join = ", ";
        $slot_suff = "$newline";
        $last_suff = "$newline";
        $foot_pref = "";
        $foot_suff = "";
        if ($double) {
            $cmplx_pref = "double complex"; # double complex number prefix
        } else {
            $cmplx_pref = "complex"; # complex number prefix
        }
    } elsif (lc($format) eq 'idl') {
        $header = "$name = {";
        $footer = "}";
        #
        $head_pref = "";
        $head_suff = " \$$newline";
        $slot_pref = "$indent";
        $slot_join = ", ";
        $slot_suff = ", \$$newline";
        $last_suff = " \$$newline";
        $foot_pref = "";
        $foot_suff = "\n";
        #
        if ($oneline) {
            $head_suff = "$newline";
            $slot_suff = ",$newline";
            $last_suff = "$newline";
        }
        #
        if ($double) {
            $cmplx_pref = "dcomplex"; # double complex number prefix
        } else {
            $cmplx_pref = "complex"; # complex number prefix
        }
    } else                         {
        croak "output(): Format <$format> unknown";
    }

    my @slots = format_slots($self,$format,$double,$trim);

    # Take care of $maxslots
    @slots = aggregate_slots(\@slots,$maxslots,$slot_join);

    # Now construct output string
    my $string;
    $string .= $head_pref
               . $header
               . $head_suff;
    if (@slots) {
        $string .= $slot_pref;
        $string .= join($slot_suff . $slot_pref, @slots);
        $string .= $last_suff;
    }
    $string .= $foot_pref
               . $footer
               . $foot_suff;

    return $string;
}

sub debug {
#
#   $obj->debug(1)     # debugging on
#   $obj->debug(0)     # debugging off
#
# Undocumented: Set/get debug flag
    my $self = shift();
    if (@_) { $self->{DEBUG_} = shift };
    return $self->{DEBUG_}
}


# ====================================================================== #

## Private utility subroutines:

sub parse_namelist {
#
# Parse first F90 namelist from text string; return reference to hash
#
#   parse_namelist(\$text,\$name,\$nslots,\@slots,$broken,$debug,$format);
#

    my $textref   = shift;
    my $nameref   = shift;
    my $nslotsref = shift;
    my $slotsref  = shift;
    my $broken    = shift;
    my $debug     = shift;
    my $format    = shift;

    my %hash;
    my $nslots = 0;
    my $state  = START;
    my $id = $regexp[ID];       # allowed namelist/variable names

    my ($status,$var,@values,$type);

    ## Reset to reasonable default values
    $$nslotsref = 0;
    @$slotsref  = ();
    $$nameref   = '';

    ## Get name of nl
    $$nameref = extract_nl_name($textref,$debug) or return;

    $status = VAR;
    ## Extract variable slots

    my $text = $$textref;

    ## Apply fix for brokenness
    if ($broken) {
        $text =~ s{\n'}{',}g;
    }

    while ($text ne '') {
        print STDERR "--------------------\nTop of while loop..\n" if ($debug);
        strip_space_and_comment($text);
        if ($text =~ s/^($id)(\([0-9, \t]+\))?\s*=\s*//s) {
            # string starts with <var=...> or <var(idx,idy,...)=...>
            $var = lc($1);
            $var =~ s/%/_/;
            # any array indices following the variable name?
            if (defined($2)) {
                my $indices = $2;
                $indices =~ s/\s+//g; # squeeze out whitespace
                $var = $var . $indices;
            }
            $status = VALUE;
            if ($debug) {
                print STDERR "parse_namelist 1: \$var=<$var>\n";
                print STDERR "parse_namelist 1: \$text=<",
                  printable_substring($text,50), ">\n";
            }

            # Get values and check
            @values = get_value(\$text,\$type,$var,$debug,$format); # drop $debug here..
            if (@values or $type == EMPTY) {
                $nslots++;
                push @$slotsref, $var;
            } else {
                show_error("Couldn't read value", "", $text, 1);
                return;
            }

        } elsif ($text =~ s{\s*(/|\$end)\s*}{}) { # string is </> or <$end>
            $status = NL_END;
            last;               # end of namelist

        } else {
            show_error("Expected var=[...] not found ", "", $text, 1);
            return;
        }

        print STDERR "[",join(',',@values), "] -> \$hash{$var}\n" if ($debug);
        my $stype = ($stypes{$type} || 'Type inconsistency!');
        $hash{$var} = { type  => $type,
                        stype => $stype,
                        value => [@values]
                      };
    }

    unless ($status == NL_END) {
        carp "Aborted parsing at <",
             printable_substring($text,50),">\n",
             "trying to read slot `$var'\n";
        return;
    }

    print STDERR
      "parse_namelist: Namelist <$$nameref> parsed succesfully\n"
        if ($debug);

    $$textref   = $text;        # propagate remainder of $text back
    $$nslotsref = $nslots;      # propagate number of slots back

    return \%hash;
}

# ---------------------------------------------------------------------- #
sub extract_nl_name {
# Extract namelist name (the part starting with `&' or `$')

    my $textref = shift;
    my $debug   = shift;

    my $text = $$textref;
    my $name;
    my $id = $regexp[ID];       # allowed namelist/variable names

    print STDERR "extract_nl_name 1: \$text=<",
                 printable_substring($text,50),">\n" if ($debug);
    strip_space_and_comment($text);

    print STDERR "extract_nl_name 2: \$text=<",
                 printable_substring($text,50), ">\n" if ($debug);

    if ($text =~ s/^(?:&|\$)($id)//) {
        $name = lc($1);
    } else {                    # empty (comment/whitespace) or erroneous
        if ($text eq '') {
            print STDERR "Empty text (at most some comments)" if ($debug);
            $$textref = $text; # propagate remainder of $text back
            return;
        } else {
            show_error("Namelist does not start with &\n","",$text,1);
            return;             # never got here..
        }
    }
    strip_space_and_comment($text);

    if ($debug) {
        print STDERR "extract_nl_name 3: \$name=<$name>\n";
        print STDERR "extract_nl_name 3: \$text=<",
                     printable_substring($text,50), ">\n";
    }

    $$textref = $text; # propagate remainder of $text back

    return $name;
}

# ---------------------------------------------------------------------- #

sub strip_space_and_comment {
# Strip leading space and anything from possible leading exclamation mark
# till end of line.
    $_[0] =~ s/^(\s*(![^\n]*)?)*//s;
    return;
}

# ---------------------------------------------------------------------- #

sub get_value {
# Extract one or several values from string that starts off immediately
# after the equal sign of a slot assignment
    my $txtptr  = shift;
    my $typeptr = shift;
    my $varname = shift;
    my $debug   = shift;    # Need to somewhow get rid of this argument...
    my $format  = shift;

    my $text = $$txtptr;

    # Sshortcut for special case of empty array, e.g. 'BCX='
    my $id = $regexp[ID];       # allowed namelist/variable names
    if ($text =~ m{^($id)(\([0-9, \t]+\))?\s*=\s*}s
        # string starts with <var=...> or <var(idx,idy,...)=...>
        or
        $text =~ m{^\s*(/|\$end)\s*}
        # string is </> or <$end>
       ) {
        $$typeptr = EMPTY;
        return ();
    }

    my @values;

    strip_space_and_comment($text); # (are comments really allowed here?)
    my $type = infer_data_type($text);
    if ($debug) {               # pretty-printing of type
        print STDERR
          "Found data of type $type (",
          elucidate_type($type),
          ") in <",
          printable_substring($text,40), ">\n";
    }

    if ($type == UNKNOWN) {
        show_error("Cannot identify data type","$varname=","$text");
        croak();
    }

    # Extract data
    my $multiregexp = $regexp[MULTIPLE]; # qr// wouldn't expand the CONSTANT...
    my $re_type     = qr/$regexp2[$type]/;

    while ($text =~ s/^
                      ($multiregexp)?($re_type)
                      \s*
                      (
                          (?: ,? \s* ! [^\n]* | , | \s+ )
                      |
                          (?=\/|\$end)
                      )
                      \s*
                     //sx) {

        my $mul = ( $1 || 1);
        my ($val,$rest) = ($2,$3);
        $mul =~ s/\*//;
        if ($debug) {
            print STDERR "\$mul=<$mul>, \$val=<$val>\n";
            print STDERR "\$rest=<", printable_substring($rest,2),
              ">, \$text=<", printable_substring($text,40), ">\n";
        }

        # `Widen' data type if necessary (e.g. integer -> double for
        # `var = 1, 2.D0')
        my $valtype = infer_data_type($val);
        $type = $valtype if ($valtype > $type);
        if ($debug) {           # pretty-printing of type
            print STDERR
              "Data type is now ($valtype) $type (",
              elucidate_type($type),
              ")\n";
        }

        # Remove quotes around (and doubled in) strings
        if ($type == SQ_STRING) {
            #$val =~ s/''/'/gs;
            $val =~ s/^''(.*)''$/\'$1\'/s;
        }
        if ($type == DQ_STRING) {
            $val =~ s/^"(.*)"$/\'$1\'/s;
            #$val =~ s/""/"/gs;
        }
        # Remove trailing whitespace embedded in strings
        if ($type == SQ_STRING || $type == DQ_STRING) {
           $val =~ s/\s+(["'])\s*$/$1/;
        }

        # Remove embedded newlines from strings (Anders' strange namelist
        # samples from Pencil Code with dust density)
        if (($type == SQ_STRING) || ($type == DQ_STRING)) {
            $val =~ s/\n//g;
        }

        $val = encapsulate_value($val,$format,$type); # preprocess values
        
        if ($mul le '1') {
           push @values, $val
        } else {
          if ($format eq 'idl') {
            push @values, "replicate($val,$mul)"
          } else {
            if ($format eq 'python') {
              #push @values, ($val) x $mul
              push @values, "numpy.repeat($val,$mul)";
            }
          }
        }
        $text =~ s/.*\n// if ($rest eq '!'); # comment
        print STDERR "<<", ($mul||'1'), "x>><<$val>> <<",
          printable_substring($text), ">>\n" if ($debug);
    }

    $$txtptr = $text;           # return remaining unparsed string
    $$typeptr = $type;          # return type

    return @values;
}

# ---------------------------------------------------------------------- #

sub elucidate_type {
# Expand numerical constants into understandable type names
    my $type = shift;

    my @tp;
    $tp[UNKNOWN]   = 'UNKNOWN';
    $tp[SQ_STRING] = 'SQ_STRING';
    $tp[DQ_STRING] = 'DQ_STRING';
    $tp[LOGICAL  ] = 'LOGICAL';
    $tp[INTEGER  ] = 'INTEGER';
    $tp[FLOAT    ] = 'FLOAT';
    $tp[SINGLE   ] = 'SINGLE';
    $tp[DOUBLE   ] = 'DOUBLE';
    $tp[COMPLEX  ] = 'COMPLEX';
    $tp[DCOMPLEX ] = 'DCOMPLEX';
    $tp[MULTIPLE ] = 'MULTIPLE';
    $tp[EMPTY    ] = 'EMPTY';

    return $tp[$type];
}

# ---------------------------------------------------------------------- #

sub infer_data_type {
# Determine the F90 data type of first item in string, skipping multiplier
# if present
    my $text = shift;
    $text =~ s/^\s*[0-9]+\*//;  # ignore multiplier for finding data type
    if      ($text =~ /^\s*'/)                     { return SQ_STRING;
    } elsif ($text =~ /^\s*"/)                     { return DQ_STRING;
    } elsif ($text =~ /^\s*\(\s*$numeric_e\s*,/)   { return COMPLEX;
    } elsif ($text =~ /^\s*\(/)                    { return DCOMPLEX;
    } elsif ($text =~ /^\s*(T|F|.(true|false).)/i) { return LOGICAL;
    } elsif ($text =~ /^\s*[+-]?[0-9]+(\s|,|!|$)/) { return INTEGER;
    } elsif ($text =~ /^\s*$float(\s|,|!|$)/)      { return FLOAT;
    } elsif ($text =~ /^\s*$numeric_e(\s|,|!|$)/)  { return SINGLE;
    } elsif ($text =~ /^\s*$numeric_d(\s|,|!|$)/)  { return DOUBLE;
    } else                                         { return UNKNOWN;
    }
}

# ---------------------------------------------------------------------- #

sub show_error {
# Print error message and beginning of string with marker of current
# position [Slightly ridiculous, as the marker will always point to
# beginning of line]
    my $errmsg = shift;
    my $prefix = shift;
    my $text   = shift;
    my $die    = (shift || 0);

    chomp($errmsg);

    # Escape newlines and only print 75 chars:
    my $subtext = $prefix . $text;
    $subtext =~ s/\n/\\n/g;
    $subtext = substr($subtext,0,75) . "\n";

    # Splice in marker line
    my $marker = " " x length($prefix) . "^------  HERE\n";
    $subtext =~ s/\n/\n$marker/;

    # Prefix error message:
    $subtext = "\e[01m$errmsg:\e[00m\n" . $subtext;

    # Now die
    if ($die) {
        croak "$subtext";       # die
    } else {
        carp "$subtext";        # warn
    }

    return;
}

# ---------------------------------------------------------------------- #

sub printable_substring {
# Extract substring and quote newlines for diagnostic printing
    my $string = shift;
    my $length = shift || 40;

    $string =~ s/\n/\\n/g;
    my $oldlen = length($string);
    $string = substr($string,0,$length);
    substr($string,-3,3) = '...' if ($length<$oldlen);

    return $string;
}

# ---------------------------------------------------------------------- #

sub assign_slot_val {
#
# Assignment of value to slot variable for output (format-dependent: `='
# for f90, `:' for IDL records, etc.)
#
    my $var    = shift;
    my @vals   = @{shift()};
    my $format = shift;
    my $type   = shift;

    my $assmnt = "$var";

    if ($format eq 'f90') {
        $assmnt .= "=";
    } elsif ($format eq 'idl') {
        $assmnt .= ": ";        # structure syntax
    } elsif ($format eq 'python') {
        $assmnt .= " = ";        # structure syntax
    } else {
        croak "assign_slot_val: Unknown format <$format>\n";
    }

    #encapsulate_values(\@vals,$format,$type); # preprocess values
    if ($format eq 'idl' and @vals == 0) {
        # IDL not only lacks empty arrays, but even a usable 'undef'
        # concept (we cannot assign '[]' or '!NULL' to a structure slot).
        # So we follow
        #   http://objectmix.com/idl-pvwave/785581-empty-arrays.html
        # and use an empty pointer instead.
        $assmnt .= 'ptr_new(/allocate_heap)';
    } elsif (@vals == 1) {  # scalar
        $assmnt .= $vals[0];
    } else {           # array
        if ($format eq 'python') {
          my $tolist=0;
          for (@vals) {
            if ($_ =~ /^ *numpy\.repeat.*$/) {
              $tolist=1;
              last;
            }
          }
          if ($tolist==1) {
            for (@vals) {
              if ($_ =~ /^ *numpy\.repeat.*$/) {
                $_ = "list($_)"
              } else {
                $_ = "[$_]"
              }
            }
            $assmnt .= add_array_bracket(join("+", @vals), 'pylist');
          } else {  
            $assmnt .= add_array_bracket(join(",", @vals), $format);
          }
        } else {
          $assmnt .= add_array_bracket(join(",", @vals), $format);
        }
    }

    return $assmnt;
}

# ---------------------------------------------------------------------- #

sub encapsulate_values {
#
# Format-specific preprocessing of data values, e.g. quoting strings,
# mapping logicals to integers for IDL, etc.
#
    my $valref = shift;
    my $format = shift;
    my $type   = shift;
    my @vals = @$valref;

    ## Actions for all formats
    if ($type==COMPLEX or $type==DCOMPLEX) {
        @vals = map { "${cmplx_pref}$_${cmplx_suff}" } @vals;
    }

    ## Actions specific for some formats
    if ($type==SQ_STRING or $type==DQ_STRING) {
        @vals = map { quote_string($_, $format) } @vals;
    }
    if ($type==LOGICAL) {
        @vals = map { encaps_logical($_, $format) } @vals;
    }

    @$valref = @vals;

    return;
}

# ---------------------------------------------------------------------- #

sub encapsulate_value {
#
# Format-specific preprocessing of data values, e.g. quoting strings,
# mapping logicals to integers for IDL, etc.
#
    my $val    = shift;
    my $format = shift;
    my $type   = shift;

    #my $val = $$valref;

    ## Actions for all formats
    if ($type==COMPLEX or $type==DCOMPLEX) {
        $val = ${cmplx_pref}.$val.${cmplx_suff};
          #print STDERR $val, "\n";
    }

    ## Actions specific for some formats
    if ($type==SQ_STRING or $type==DQ_STRING) {
        #$valref = quote_string($valref, $format);
    }
    if ($type==LOGICAL) {
        $val = encaps_logical($val, $format);
    }

    return $val;
}

# ---------------------------------------------------------------------- #

sub format_slots {
#
# Format all slots for printing and collect in a list
#
    my $obj    = shift;
    my $format = (shift || 0);
    my $double = (shift || 0);
    my $trim   = (shift || 1);

    return () unless ($obj->{NSLOTS});

    my @slots;
    my $slot;
    foreach my $var (@{$obj->{SLOTS}}) {
        my $valhash = $obj->{DATA}->{$var};
        my @vals = @{$$valhash{'value'}};
        my $type = $$valhash{'type'};

        # Trim trailing whitespace
        if ($trim) {
            for (@vals) { s/\s+$//; }
        }

        # Replace E[0-9]+ by, or append `D0' where necessary
        if ($double && $format ne 'python') {
            if (($type == FLOAT)   ||
                ($type == SINGLE)  ||
                ($type == DOUBLE)  ||
                ($type == COMPLEX) ||
                ($type == DCOMPLEX))  {
                for (@vals) {
                    s/([0-9\.])[eE]/$1D/g;
                    s/(^|\s|replicate\s*\(|,)($float)($|\s|,)/$1$2D0$3/g;
                }
            }
        }

        # Trim unnecessary zeros and spaces in float and double
        if ($trim) {
            if (($type == FLOAT)   ||
                ($type == SINGLE)  ||
                ($type == DOUBLE)  ||
                ($type == COMPLEX) ||
                ($type == DCOMPLEX))  {
                for (@vals) {
                    s/(\.[0-9]+?)0+(?=\D|$)/$1/g;
                    s/([0-9\.][eEdD])[\-\+]?0+(?=\D|$)/${1}0/g;
                    s/([0-9\.][eEdD][\-\+]?)0+(?=[1-9])/$1/g;
                    s/\s+,/,/g;
                    s/,\s+/,/g;
                }
            }
            if (($type == COMPLEX) ||
                ($type == DCOMPLEX))  {
                for (@vals) {
                    s/\(\s+/\(/g;
                    s/\s+\)/\)/g;
                }
            }
        }

        $slot = assign_slot_val($var,\@vals,$format,$type);
        push @slots, $slot;
    }

    return @slots;
}

# ---------------------------------------------------------------------- #
sub aggregate_slots {
#
# Take list of formatted slot strings, group every $maxslots strings
# together into one string
#
    my $slotsref  = shift;
    my $maxslots  = shift;
    my $slot_join = shift;

    my @slots = @$slotsref;

    # Short-circuit if nothing to do
    return @slots unless ($maxslots>0);

    my @new_slots;
    while (@slots) {
        # Use a loop here, as @slots[1..$maxslots] would generate trailing
        # undefs
        my @group_us;
        foreach my $i (1..$maxslots) {
            push @group_us, shift @slots if (@slots);
        }
        my $aggregated_slot = join($slot_join, @group_us);
        push @new_slots, $aggregated_slot;
    }

    return @new_slots;
}
# ---------------------------------------------------------------------- #

sub add_array_bracket {
# Add format-specific array delimiters around string
    my $string = shift;
    my $format = shift;

    if      ($format eq 'f90') {
        # No delimiters
    } elsif ($format eq 'pylist') {
        $string = "numpy.array($string)";
    } elsif ($format eq 'idl' || $format eq 'python') {
        $string = "[$string]";
    } else {
        croak "add_array_bracket: Unknown format <$format>\n";
    }

    return $string;
}

# ---------------------------------------------------------------------- #

sub encaps_logical {
# Represent logical
    my ($val, $format) = @_;
    
    my ($false, $true);
    if ($format eq 'f90') {
        ($false, $true) = ('F', 'T');
    } elsif ($format eq 'python') {
        ($false, $true) = ('False', 'True');
    } elsif ($format eq 'idl') {
        ($false, $true) = ('0L', '-1L');
    } else {
        croak "encaps_logical: Unknown format <$format>\n";
    }

    $val =~ s{ (\.false\. | F ) }{$false}xi;
    $val =~ s{ (\.true\.  | T ) }{$true}xi;

    return $val;
}

# ---------------------------------------------------------------------- #

sub quote_string {
# Enclose string by quotation marks, doubling any existing quotation marks
# for Fortran and IDL
    my ($val, $format) = @_;

    if ($format eq 'f90' or $format eq 'idl') {
        $val =~ s/'/''/g;
    } elsif ($format eq 'python') {
        $val =~ s/'/\\'/g;
    } else {
        croak "quote_strings: Unknown format <$format>\n";
    }

    return "'$val'";
}

# ---------------------------------------------------------------------- #

## Done.

1;

# End of file
