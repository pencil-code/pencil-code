#
#                         F90Namelist/Group.pm
#                         --------------------
#
# Description:
#   Parse a group of F90 namelists into F90Namelist.pm objects and export in
# different formats.
# Author: wd (wdobler [at] cpan.org)
# $Date: 2007-08-13 10:12:21 $
# $Revision: 1.1 $
# [Date and CVS revision are now pretty irrelevant, as I keep the code
#  under Darcs now]

package Fortran::F90Namelist::Group;

=head1 NAME

Fortran::F90Namelist::Group - Parse F90 namelist groups and export in different formats

=head1 SYNOPSIS

  use Fortran::F90Namelist::Group;
  my $nlgrp = Fortran::F90Namelist::Group->new() or die "Couldn't get object\n";

  $nlgrp->parse(<<'  HERE');
    &runpars
      x=2,y=3
      vec1=1,2,3
      vec2=3*1.3
    /
    &more_runpars
      z=7
      vec1=0,1,2,3
    /
  HERE

Read from file:

  # Read namelist group from file `some_lists.nml':
  $nlgrp->parse(file => 't/files/some_lists.nml');

  # Read namelist group from file handle
  open(my $fh , "< t/files/some_lists.nml") or die "Couldn't get file handle\n";
  $nlgrp->parse(file => $fh);
  # or
  open(NLGROUP , "< t/files/some_lists.nml") or die "Couldn't get file handle\n";
  $nlgrp->parse(file => \*NLGROUP);

  # Print names of all namelists in file `start.in'
  $nlgrp->parse(file => 't/files/start.in') or die "Couldn't parse\n";
  print join(" ", $nlgrp->names), "\n";

Extract or merge namelists from group and return I<Fortran::F90Namelist> object:

  my $nl_1   = $nlgrp->first();   # Extract first namelist from group
  my $nl_3   = $nlgrp->nth(3);    # Extract 4th namelist from group
  my $nl_all = $nlgrp->flatten(); # Flatten all namelists into one

Write namelist group:

  # Write namelist group in F90 namelist format
  print $nlgrp->output();

  # Write namelist as IDL structure
  print $nlgrp->output(format => 'idl');


=head1 DESCRIPTION

I<Fortran::F90Namelist::Group> is a module for parsing Fortran90 namelist
groups into an internal format based on
L<Fortran::F90Namelist|Fortran::F90Namelist>, and re-exporting in different
formats.
Parsing is done by L<Fortran::F90Namelist|Fortran::F90Namelist>, see the
documentation of that module for more details.


=head2 Methods

=over 4


=item B<$nlgrp-E<gt>new>()

Create a new namelist group object

=item B<$nlgrp-E<gt>parse>(I<string>)

=item B<$nlgrp-E<gt>parse>(text => I<string>)

=item B<$nlgrp-E<gt>parse>(file => I<fname>|I<fhandle>)

=item B<$nlgrp-E<gt>parse>(file => I<fname>|I<fhandle> [, <options> ])

Parse I<string> or the file represented by I<fname> or i<fhandle> (a file
handle or L<File::Handle|File::Handle> object [not yet implemeted]);
returns number of parsed namelists, or undef if parsing failed.

Additional I<options> are:

=over 8

=item B<append>

If true, append newly parsed namelists to already existing data in the
object.

=back


=item B<$nlgrp-E<gt>nlists>()

Return number of namelists in group.


=item B<$nlgrp-E<gt>names>()

Return list of namelist names in group (in original order).


=item B<$nlgrp-E<gt>insert>(nl [, pos])

Insert namelist into namelist group at position POS (defaults to appending
nl at end of group.
Returns 1 if successfull, 0 or undef otherwise.


=item B<$nlgrp-E<gt>delete>(nl)

=item B<$nlgrp-E<gt>delete>(name)

=item B<$nlgrp-E<gt>delete>(num)

Delete namelist (identified by namelist object, name, or position in
$nlgrp->names) from namelist group.
Returns 1 if successfull, 0 otherwise.


=item B<$nlgrp-E<gt>first>()

Return the first namelist in the group as
L<Fortran::F90Namelist|Fortran::F90Namelist> object.


=item B<$nlgrp-E<gt>nth>(n)

Return the namelist with index n from the group as
L<Fortran::F90Namelist|Fortran::F90Namelist> object.
Indices count from 0, so this returns the (n+1)st namelist.


=item B<$nlgrp-E<gt>pop>(n)

Return the first namelist in the group as
L<Fortran::F90Namelist|Fortran::F90Namelist> object and remove it from the
group.
Returns C<undef> for an empty group.
This allows to write:

  while (my $nl = $nlgrp->pop()) {
      print $nl->name(), " has ", $nl->nslots(), "slots\n";
  }


=item B<$nlgrp-E<gt>flatten([options])>

Merge all namelist data in the group into one
L<Fortran::F90Namelist|Fortran::F90Namelist> object.
Thus,

  $nlgrp->parse(file => 't/files/some_lists.nml');
  my $nl = $nlgrp->flatten();

is another way of doing

  my $nl = Fortran::F90Namelist->new();
  $nl->parse(file => 't/files/some_lists.nml',
             all  => 1    );

I<Options> are:

=over 8

=item B<name>

Set name of resulting namelist (default: name of first namelist read).

=item B<dups_ok>

Don't warn if new slots have same names as existing slots.

=back


=item B<$nlgrp-E<gt>hash>()

Return namelist group as Perl hash.
See L<HASH FORMAT|/"HASH FORMAT"> below for details of the hash format.


=item B<$nlgrp-E<gt>output>(format => I<format>)

Write namelist group in given I<format>.
Currently supported formats are `f90' (default), and `idl'


=back


=head1 HASH FORMAT

The B<hash> method returns a hash reference of the following structure:

=for test ignore

    { namelist1 => { var1 => { 'value' => [ value1, value2, ..],
                               'type'  => numerical_type,
                               'stype' => "type string"
                             },
                     var2 => { 'value' => [ value1, value2, ..],
                               'type'  => numerical_type
                               'stype' => "type string"
                             },
                     ...
                   },
      namelist2 => { var1 => { 'value' => [ value1, value2, ..],
                               'type'  => numerical_type,
                               'stype' => "type string"
                             },
                     ...
                   },
      ...
    }

=for test

Here I<numerical_type> is a number identifying each data type, while
I<stype> is a textual description of the given data type.

E.g.

=for test ignore

    { 'hydro_init_pars'   => { 'u0' => { 'value' => [ 0., -3.141593, 0. ],
                                         'type'  => 6,
                                         'stype' => 'single precision float'
                                       },
                               'u1' => { 'value' => [ 0., 0., 0.],
                                         'type'  => 6,
                                         'stype' => 'single precision float'
                                       },
                             },
      'density_init_pars' => { 'rho0'   => { 'value' => [ -2.78 ],
                                             'type'  => 6,
                                             'stype' => 'single precision float'
                                           },
                               'ilnrho' => { 'value' => [ 3 ],
                                             'type'  => 4,
                                             'stype' => 'integer'
                                           },
                            },
    }

=for test

Note: This is currently just the internal format used to represent
namelists and can thus change in the future.
In particular the C<type> numbers should not be considered to be stable
between releases.


=head1 TO DO

=over 4

=item 1.

More output methods:

=over 8

=item *

Octave/Matlab , C struct, YAML, XML(?), ...

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

Copyright (c) 2007, Wolfgang Dobler <wdobler [at] cpan.org>.
All rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same conditions as Perl or under the GNU General Public
License, version 2 or later.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=cut


use warnings;
use strict;
use Carp;
use vars qw($VERSION);
use Fortran::F90Namelist;

$VERSION = '0.5.2';

# ---------------------------------------------------------------------- #
##
## Object constructor
##
## Internal structure of Namlist::Group objects:
##   DATA    -- variable names, values, and types (hashref, see below)
##   NAMES   -- ordered list of variable names (array ref)
##   NLISTS  -- number of namelists
##   PARSED_ -- flag indicating that argument has been parsed
##   DEBUG_  -- debug flag
##
##   {
##     'NAMES'  => [ 'nlname1', 'nlname2', ..., 'nlnameN' ],
##     'DATA'   => { 'nlname1' => $nlobj1,
##                   'nlname2' => $nlobj2,
##                   ...
##                   'nlnameN' => $nlobjN },
##     'NLISTS' => N,
##     'PARSED_' => 0/1,
##     'DEBUG_'  => 0/1
##   }
##
sub new {
    my $proto = shift;          # either classref or object ref or string
    my @argv  = @_;
    my $class = ref($proto) || $proto;
    my $self = {};

    my %data   = ();
    my @names  = ();
    my $nlists = undef;
    my $parsed = 0;

    my $short_usage =
        "Usage:\n" .
        "  Fortran::F90Namelist::Group->new()\n" ;
    #    unless($file) {
    #   croak $short_usage;
    #   return;
    #    }

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
    my $debug  = ($args{debug} || 0);

    ##
    ## Finish the object
    ##
    # public data of object
    $self->{DATA}   = \%data;
    $self->{NAMES}  = \@names;
    $self->{NLISTS} = $nlists;

    # internal data
    $self->{PARSED_} = $parsed;
    $self->{DEBUG_}  = $debug;

    bless($self,$class);
    return $self;
}

# ====================================================================== #

##
##  Methods
##

sub parse {
#
#   $obj->parse($text)
#   $obj->parse(file => $filename|$filehandle)
#   $obj->parse({text => $textstring)
#   $obj->parse(... , append => 1)
#
# Parse text or file containing F90 namelist(s)
#
    my $self = shift;
    my @args = @_;              # can't use shift() since we change value
                                # of $text
    my %args;
    my $text = '';

    # Parse arguments (file => <filename>, etc.); may be single string,
    # list, hash or hashref
    if (ref($args[0]) eq 'HASH') { # parse($hashref)
        %args = %{$args[0]};
    } else {
        if (@_ == 1) {          # parse($string)
            $text = $args[0];
        } else {                # parse(%hash) or parse(@list)
            %args = @args;
        }
    }
    my $file   = ($args{file}   || ''    );
    $text      = ($args{text}   || $text );
    my $append = ($args{append} || 0     );

    my $new_nlists = 0;         # counter

    if (! $append) {            # clear all data
        $self->{DATA}    = {};
        $self->{NAMES}   = [];
        $self->{NLISTS}  = 0;
        $self->{PARSED_} = 0;
    }

    # Get text from file if necessary
    if ($text eq '') {
        croak "Fortran::F90Namelist::Group->parse(): need text or file argument\n"
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

    while (length($text) > 0) {
        # Note: shouldn't be difficult to reuse F90Namelist object here
        my $nl = Fortran::F90Namelist->new()
          or croak "Couldn't get F90Namelist object\n";
        $nl->debug(1) if $self->{DEBUG_}; # hand down debugging flag
        $nl->parse($text);
        my $name = $nl->name;
        if (defined($name) and $name ne '') {
            $self->{DATA}->{$name} = $nl;
            push @{$self->{NAMES}}, $name;
            $self->{NLISTS}++;
            $new_nlists++;
        }
    }

    $self->{PARSED_} = 1;

    return $new_nlists;
}

# ---------------------------------------------------------------------- #

sub nlists {
#
# Get number of namelists in group
#
    my $self = shift();
    return $self->{NLISTS};
}

# ---------------------------------------------------------------------- #

sub names {
#
# Return array of namelist names in order
#
    my $self = shift();
    return $self->{NAMES};
}

# ---------------------------------------------------------------------- #

sub insert {
#
# Insert namelist at given position
#
#   $nlgrp->insert($nlobj)
#   $nlgrp->insert($nlobj, $pos)
#
    my $self = shift();
    my ($nlobj,$pos) = @_;

    $pos = $self->nlists() unless defined($pos);

    croak "Usage:  Fortran::F90Namelist::Group::insert(\$nlobj,\$pos)\n"
      unless (ref($nlobj) eq "Fortran::F90Namelist");

    my $nlname = $nlobj->name();

    # Insert namelist data
    $self->{DATA}->{$nlname} = $nlobj;
    # Insert name in order list
    my @names = @{$self->{NAMES}};
    splice @names, $pos, 0, $nlname;
    $self->{NAMES} = \@names;
    # Increment counter
    $self->{NLISTS}++;

    return 1;
}

# ---------------------------------------------------------------------- #

sub delete {                    ##no critic ProhibitBuiltinHomonyms
#
# Delete namelist from group
#
#   $nlgrp->delete($nlobj)
#   $nlgrp->delete($nlname)
#   $nlgrp->delete($pos)
#

    my $self = shift();
    my ($arg1) = @_;

    # Extract nl name from whatever argument we have
    my $nlname;
    if (ref($arg1) eq 'Fortran::F90Namelist') { # F::NL object
        $nlname = $arg1->name();
    } elsif ($arg1 =~ /^[0-9]+$/) { # number (=position)
        $nlname = $self->{NAMES}[$arg1];
    } else {                    # name
        $nlname = $arg1;
    }

    # Delete namelist data
    delete $self->{DATA}->{$nlname} or return 0;

    # Remove name from list
    my @names = grep { !/^${nlname}$/ } @{$self->{NAMES}};
    $self->{NAMES} = \@names;

    # Decrement counter
    $self->{NLISTS}--;

    return 1;
}

# ---------------------------------------------------------------------- #

sub first {
#
# Return the first namelist in the group as Fortran::F90Namelist object.
#
    my $self = shift();

    return $self->nth(0);
}

# ---------------------------------------------------------------------- #

sub nth {
#
# Return the n-th namelist in the group as Fortran::F90Namelist object.
# Counts from 0.
#
    my $self = shift();
    my $n    = shift();

    my $nlists = $self->nlists();
    if (($n < 0) || ($n >= $nlists)) {
        croak "Fortran::F90Namelist::Group->nth(): "
          . "Need 0 <= n < $nlists, but n=$n";
    }
    my $nlname = $self->{NAMES}[$n];
    my $nl   = $self->{DATA}->{$nlname};

    return $nl;
}

# ---------------------------------------------------------------------- #

sub pop {                       ##no critic ProhibitBuiltinHomonyms
#
# Return the first namelist in the group as Fortran::F90Namelist object (or
# undef if there is no first namelist) and remove it from the group,
# pretty much like pop(@arr) does for arrays.
#
    my $self = shift();

    if ($self->nlists() < 1) {
        return;
    } else {
        my $nl = $self->first();
        $self->delete(0);
        return $nl;
    }
}

# ---------------------------------------------------------------------- #

sub flatten {
#
# Merge all namelist data in the group into one Fortran::F90Namelist object.
#
    my $self = shift();
    my @argv = @_;         # will just be handed over to F90Namelist::merge()

    my $nl = $self->first();

    foreach my $i (1..$self->{NLISTS}-1) {
        my $nl2 = $self->nth($i);
        $nl->merge($nl2, @argv);
    }

    return $nl;
}

# ---------------------------------------------------------------------- #

sub hash {
#
# Return hash with parsed namelist contents
#
    my $self = shift;

    # Collect the $nl->hash() values of the individual namelists
    my $hashref = {};
    foreach my $name (@{$self->{NAMES}}) {
        my $nl = $self->{DATA}->{$name};
        my $name = $nl->name();

        $hashref->{$name} = $nl->hash;
    }

    return $hashref;
}

# ---------------------------------------------------------------------- #

sub output {
#
# Write namelist in specified format (defaults to 'f90')
# Optional arguments (see Fortran::F90Namelist.pm):
#   format   => format   ('f90' [default] or 'idl')
#   trim     => 0/1      (trim strings?)
#   double   => 0/1      (mark all floats as double precision?)
#   oneline  => 0/1      (write all in one line? [only for some formats])
#   maxslots => N        (similar to oneline, but split every N slots)
#
    my $self = shift();
    my @argv = @_;        # will just be handed over to F90Namelist::output()

    my $string='';

    foreach my $nlname (@{$self->{NAMES}}) {
        my $nl = $self->{DATA}->{$nlname};
        $string .= $nl->output(@argv);
    }

    return $string;

}

# ---------------------------------------------------------------------- #

sub debug {
#
#   $obj->debug(1)     # debugging on
#   $obj->debug(0)     # debugging off
#
# Purposefully undocumented: Set/get debug flag
#
    my $self = shift();
    if (@_) { $self->{DEBUG_} = shift };

    return $self->{DEBUG_}
}


# ====================================================================== #

## Private utility subroutines:

# ---------------------------------------------------------------------- #

## Done.

1;

# End of file
