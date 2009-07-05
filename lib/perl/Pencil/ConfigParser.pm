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

sub debug {
#
# Get/set debugging flag
#
    my $self = shift();
    my ($debug) = @_;

    if (defined($debug)) {
        $self->{DEBUG} = $debug;
    }
    return $self->{DEBUG};
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

    my %section_map;
    get_sections(
        $self->{FILENAME}, '__GLOBAL__', \%section_map, $self->{DEBUG}
    );
    $self->{SECTIONS} = \%section_map;

    $self->{PARSED} = 1;
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
# Store the given section (or all sections found, if $section is the
# global section) in the hash referenced by the third argument:
#   { name1 => [line11, line12, ...],
#     name2 => [line21, line22, ...],
#     ... }
#
    my ($file, $enclosing_section, $section_map_ref, $debug) = @_;

    print STDERR "get_sections($file, $enclosing_section)\n" if ($debug);

    my @sections = ($enclosing_section);

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
            print STDERR "  %section $1\n" if ($debug);
            push @sections, $1;
            next line;
        }

        if ($line =~ /^\s*%endsection\s+(\S+)\s*$/) { # end section
            print STDERR "  %endsection $1\n" if ($debug);
            my $ending_section = $1;
            unless ($ending_section eq pop(@sections)) {
                croak "Ending section <$ending_section> that is not open:\n"
                  . "  $file:$.: $line\n";
            }
            next line;
        }

        if ($line =~ /^\s*%include\s+(\S+)\s*$/) { # include config file
            my $include_file = find_include_file($1, $file);
            print STDERR "  %include $1 -> $include_file\n" if ($debug);
            get_sections($include_file, $sections[-1], $section_map_ref, $debug);
            next line;
        }

        my $complete_line = $line_fragment . $line;
        if ($debug) {
            my $output = substr($complete_line, 0, 30);
            chomp($output);
            print STDERR "  $output\n" if ($debug);
        }
        push @{$section_map_ref->{$sections[-1]}},
            normalize_line($complete_line);
        $line_fragment = '';
    }
    close $fh;

    if (@sections != 1) {
        croak "Section <" . pop(@sections) . "> not closed"
          . " in file $file\n";
    }
}

# ---------------------------------------------------------------------- #
sub find_include_file {
#
# Locate a file for $fragment, looking in the following directories:
#
# If $fragment starts with `./', only the directory part of $base_file is
# considered.
# Otherwise, the search path consists of
# 1. ${HOME}/.pencil-config ,
# 2. ${PENCIL_HOME}/config .
#
    my ($fragment, $base_file) = @_;

    my @path = ();
    if ($fragment =~ m{^\./(.*)}) {
        $fragment = $1;
        push @path, ".";
    } else {
        for my $dir ("$ENV{HOME}/.pencil-config", "$ENV{PENCIL_HOME}/config") {
            push @path, $dir if (-d $dir);
        }
    }
    if (@path < 1) {
        croak("No elements in path\n");
    }

    for my $dir (@path) {
        my $file = "$dir/$fragment.conf";
        if (-e $file) {
            return $file;
        }
    }

    croak "Couldn't find file $fragment.conf in path (@path)\n";
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
  print $parser->get_makefile_params->{'FFLAGS'};

  my %runtime_params = $parser->get_runtime_params();
  print $parser->get_runtime_params->{'local_disc'};

  my @all_make_keys = $parser->get_makefile_keys();
  my @all_run_keys = $parser->get_runtime_keys();


=head1 DESCRIPTION

Pencil::ConfigParser parses Pencil Code configuration files, recursively
expanding `%include' macros.

=head2 Methods

=over 4


=item B<new>($filename)

Create a new object. I<$filename> will be parsed lazily, i.e. at the point
where results from parsing are actually needed.

=item B<debug>([$debug])

With argument: set debugging flag.
Returns the debugging flag.

=item B<get_makefile_params>()

Return a hashref of the parmeters defined in the `Makefile' section:

  { KEY1 => value1, KEY2 => value2, ...}

=item B<get_runtime_params>()

Return a hashref of the parmeters defined in the `runtime' section:

  { key1 => value1, key1 => value2, ...}

=item B<get_makefile_keys>()

Return arrayref of all keys in the `Makefile' section, in the order in
which they occured int the config files.

=item B<get_runtime_keys>()

Return arrayref of all keys in the `runtime' section, in the order in
which they occured int the config files.


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
