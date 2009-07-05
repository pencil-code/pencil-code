#
#                            ConfigParser.pm
#                            ---------------
#
# Description:

#   Parse a Pencil Code .conf file, recursively expanding %include macros
#   and checking that %begin and %end are balanced within each file.
#
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Pencil::ConfigParser;

use warnings;
use strict;
use Carp;
use vars qw($VERSION);

##use critic

$VERSION = '0.1';

# ---------------------------------------------------------------------- #
##
## Object constructor
##

sub new {
#
#   Pencil::ConfigParser->new($filename);
#
    my $proto = shift;          # either classref or object ref or string
    my @argv  = @_;

    my $self = {};
    my $class;
    my $parent = {};
    if (ref($proto)) {
        $class = ref($proto);
        # If this is a call of the form $object->new(), extract properties
        # of original object and use them as default values.
        if ($proto->isa('Pencil::ConfigParser')) {
            $parent = $proto;
        }
    } else {
        $class = $proto;
    }

    croak "Usage: ConfigParser->new(\$filename)" unless (@argv == 1);
    $self->{FILENAME}  = pop(@argv);

    $self->{DEBUG}   = 0;
    $self->{VERBOSE} = 0;
    # debugging implies verbosity:
    $self->{VERBOSE} = $self->{DEBUG} || 0;

    $self->{PARSED} = 0;

    bless($self, $class);
    return($self);
}

# ====================================================================== #

##
## Public methods
##

sub get_makefile_params {
#
# Return a hashref containing the Makefile parameters
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{MAKE_PARAMS_REF};
}

# ---------------------------------------------------------------------- #

sub get_section_hash {
#
# For debugging only -- to be removed
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{SECTIONS};
}

# ====================================================================== #

##
## Private methods
##

sub parse {
#
# Parse the config file
#
    my $self = shift;
    carp "Why am I parsing this twice?" if ($self->{PARSED});

    $self->{SECTIONS} = get_sections($self->{FILENAME}, '__GLOBAL__');
}

# ====================================================================== #


##
## Utility (class) methods
##

# ---------------------------------------------------------------------- #
sub get_sections {
#
# Read file line by line,
# - discarding comment lines
# - connecting continuation lines
# - recursively expanding %include statements
#
# Return the given section (or all sections foun, if $section is the
# global section) as a hash
#   { name1 => [line11, line12, ...],
#     name2 => [line21, line22, ...],
#     ... }
#
    my ($file, $section0) = @_;

    my @sections = ($section0);

    my %section_map;
    my $line_fragment = '';

    open(my $fh, "< $file") or croak "Cannot open $file for reading: $!";
    line: while (defined(my $line = <$fh>)) {

        $line =~ s{#.*}{};       # remove comments

        if ($line =~ /^\s*$/) { # whitespace line
            next line;
        }

        if ($line =~ /(.*)\\\s*$/) { # continuation line
            $line_fragment .= $1;
            next line;
        }

        if ($line =~ /^\s*%section\s+(\S+)\s*$/) { # start section
            push @sections, $1;
            next line;
        }

        if ($line =~ /^\s*%endsection\s+(\S+)\s*$/) { # end section
            my $ending_section = $1;
            unless ($ending_section eq pop(@sections)) {
                croak "Ending section <$ending_section> that is not open:\n"
                  . "  $file:$.: $line\n";
            }
            next line;
        }

        my $complete_line = $line_fragment . $line;
        push @{$section_map{$sections[-1]}}, normalize_line($complete_line);
        $line_fragment = '';
    }
    close $fh;

    if (@sections != 1) {
        croak "Section <" . pop(@sections) . "> not closed"
          . " in file $file\n";
    }

    return \%section_map;
}

# ---------------------------------------------------------------------- #
sub normalize_line {
#
# Strip irrelevant whitespace from $line
#
    my ($line) = @_;

    $line =~ s{^\s*(.*?)\s*$}{$1};

    return $line;
}
# ---------------------------------------------------------------------- #

1;

__END__


=head1 NAME

Pencil::ConfigParser - Parse Pencil Code configuration files

=head1 SYNOPSIS

  use Pencil::ConfigParser;
  my $parser = Pencil::ConfigParser->new('mycomputer.conf');

  my %make_params = $parser->get_makefile_params();
  print $parser->get_makefile_param('FFLAGS');

  my %runtime_params = $parser->get_runtime_params();
  print $parser->get_runtime_param('local_disc');

  my @all_make_keys = $parser->get_makefile_keys();
  my @all_run_keys = $parser->get_runtime_keys();


=head1 DESCRIPTION

Pencil::ConfigParser parses Pencil Code configuration files, recursively
expanding `%include' macros.

=head2 Methods

=over 4


=item B<$doc-E<gt>new>($filename)

Create a new object.
I<$filename> is

=item B<$doc-E<gt>parse>(file)

Bla, bla...

=item B<$doc-E<gt>path(@path)

Set the internal config search path to @path (an array of path name
strings).
If @path is I<undef>, don't change the search path.
Returns the new search path.

[Note: this doesn't really work: we need to set @path before the file
name...]

=head1 BUGS AND LIMITATIONS

=over 4

=item *

None worth mentioning (so far).

=back


=head1 AUTHOR

Wolfgang Dobler <wdobler [at] cpan [dot] org


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2009, Wolfgang Dobler <wdobler [at] cpan [dot] org>.
All rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same conditions as Perl or under the GNU General Public
License, version 3 or later.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=cut


# End of file
