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
#   Pencil::ConfigParser->new(@filenames);
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

    croak "Usage: ConfigParser->new(\@filenames)" unless (@argv >= 1);
    $self->{FILENAMES}  = \@argv;

    $self->{DEBUG}   = 0;
    $self->{VERBOSE} = 0;
    # debugging implies verbosity:
    $self->{VERBOSE} = $self->{DEBUG} || 0;

    my %empty_map;
    my @empty_list;
    $self->{MAKEFILE_PARAMS} = \%empty_map;
    $self->{MAKEFILE_KEYS}   = \@empty_list;
    $self->{RUNTIME_PARAMS}  = \%empty_map;
    $self->{RUNTIME_KEYS}    = \@empty_list;

    $self->{PARSED} = 0;

    bless($self, $class);
    return($self);
}

# ====================================================================== #

##
## Public methods
##

# ---------------------------------------------------------------------- #

sub get_makefile_params {
#
# Return hash ref of Makefile parameters
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{MAKEFILE_PARAMS};
}

# ---------------------------------------------------------------------- #

sub get_makefile_keys {
#
# Return array ref of Makefile keys
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{MAKEFILE_KEYS};
}

# ---------------------------------------------------------------------- #

sub get_makefile_args {
#
# Return array ref of Makefile arguments
#   [ 'VAR1=val1', 'VAR2='val2', ... ]
# that can be interpolated into a `make' command line.
#
# Note that whitespace in a value is currently _not_ escaped, so
# the array is best used in a system() call with list argument.
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    my @args;
    for my $key (@{$self->{MAKEFILE_KEYS}}) {
        push @args, "$key=$self->{MAKEFILE_PARAMS}->{$key}";
    }

    return \@args;
}

# ---------------------------------------------------------------------- #

sub get_runtime_params {
#
# Return hash ref of runtime parameters
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{RUNTIME_PARAMS};
}

# ---------------------------------------------------------------------- #

sub get_runtime_keys {
#
# Return array ref of Makefile keys
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{RUNTIME_KEYS};
}

# ---------------------------------------------------------------------- #

sub get_environment_params {
#
# Return hash ref of runtime parameters
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{ENVIRONMENT_PARAMS};
}

# ---------------------------------------------------------------------- #

sub get_environment_keys {
#
# Return array ref of Makefile keys
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    return $self->{ENVIRONMENT_KEYS};
}

# ---------------------------------------------------------------------- #

sub get_environment_args {
#
# Return array ref of environment variable settings
#   [ 'VAR1=val1', 'VAR2='val2', ... ]
# that can be interpolated into an `env' or `export' command line.
#
# Note that whitespace in a value is currently _not_ escaped, so
# the array is best used in a system() call with list argument.
#
    my $self = shift();

    $self->parse() unless ($self->{PARSED});

    my @args;
    for my $key (@{$self->{ENVIRONMENT_KEYS}}) {
        push @args, "$key=$self->{ENVIRONMENT_PARAMS}->{$key}";
    }

    return \@args;
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

    # Parse everything into %section_map:
    #   { section1 => [line11, line12, ...],
    #     section2 => [line21, line22, ...], ... }
    my %section_map;
    for my $filename (@{$self->{FILENAMES}}) {
        get_sections(
            $filename, '__GLOBAL__', \%section_map, $self->{DEBUG}
        );
    }
    if ($self->{DEBUG}) {
        eval {
            require Data::Dumper;
            print STDERR '\%section_map = ',
              Data::Dumper::Dumper(\%section_map);
        }
    }

    # Parse %section_map into XYZ_params hashes and keys (for retaining
    # order information)
    foreach my $section (keys %section_map) {
        my ($map_ref, $keys_ref)
          = parse_lines($section_map{$section}, $self->{DEBUG});
        if ($section eq 'Makefile') {
            $self->{MAKEFILE_PARAMS} = $map_ref;
            $self->{MAKEFILE_KEYS} = $keys_ref;
        } elsif ($section eq 'runtime') {
            $self->{RUNTIME_PARAMS} = $map_ref;
            $self->{RUNTIME_KEYS} = $keys_ref;
        } elsif ($section eq 'environment') {
            $self->{ENVIRONMENT_PARAMS} = $map_ref;
            $self->{ENVIRONMENT_KEYS} = $keys_ref;
        } else {
            carp "Warning: Unknown section <$section>\n";
        }

    }

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

        if ($line =~ /(.*) \\ \s* $/x) { # continuation line
            $line_fragment = join_fragments($line_fragment, $1);
            next line;
        }

        if ($line =~ /^\s*%section\s+(\S+)\s*$/) { # start section
            my $sect = $1;
            print STDERR "  %section $sect\n" if ($debug);
            push @sections, $sect;
            next line;
        }

        if ($line =~ /^\s*%endsection\s+(\S+)\s*$/) { # end section
            my $sect = $1;
            print STDERR "  %endsection $sect\n" if ($debug);
            my $ending_section = $sect;
            unless ($ending_section eq pop(@sections)) {
                croak "Ending section <$ending_section> that is not open:\n"
                  . "  $file:$.: $line\n";
            }
            next line;
        }

        # Ignore current line if it is in the wrong section:
        unless ( ($sections[-1] eq $enclosing_section)
             || ($enclosing_section eq '__GLOBAL__')) {
            next line;
        }

        if ($line =~ /^\s*%include\s+(\S+)\s*$/) { # include config file
            my $include_file = find_include_file($1, $file);
            print STDERR "  %include $1 -> $include_file\n" if ($debug);
            get_sections($include_file, $sections[-1], $section_map_ref, $debug);
            next line;
        }

        my $complete_line = join_fragments($line_fragment, $line);
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

sub parse_lines {
#
# Parse an arrayref of lines
#   ['VAR1 = rhs1', 'VAR2 = rhs2a', 'VAR3 = rhs3', 'VAR2 += rhs2b', ...]
# into a hash
#   { 'VAR1' => 'rhs1', 'VAR2' => 'rhs2a rhs2b', 'VAR3' => 'rhs3', ... }
# and the ordered list of keys
#   ['VAR1', 'VAR2', 'VAR3', ...]
#
    my ($lines_ref, $debug) = @_;

    my (%map, @keys);
    foreach my $line (@$lines_ref) {
        ($line =~ /^\s*([^=]*?)\s*(\+?=)\s*(.*?)\s*$/)
          or croak "Cannot parse line <$line>\n";
        my ($key, $op, $val) = ($1, $2, $3);

        push @keys, $key unless defined($map{$key});

        if ($op eq '=') {
            $map{$key} = $val;
        } elsif ($op eq '+=') {
            my $oldval = $map{$key};
            if (defined $oldval && $oldval !~ /^\s*$/) {
                $map{$key} .= " " . $val;
            } else {
                $map{$key} = $val;
            }
        } else {
            die "Unexpected assignment operator `$op'\n";
        }
    }

    return (\%map, \@keys);
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

sub join_fragments {
#
# Join two strings, replacing any whitespace around the contact point with
# exactly one space chacater.
#
    my ($left, $right) = @_;

    $left  =~ s{\s*$}{};
    $right =~ s{^\s*}{};

    return $left . ' ' . $right;

}
# ---------------------------------------------------------------------- #

1;

__END__


=head1 NAME

Pencil::ConfigParser - Parse Pencil Code configuration files

=head1 SYNOPSIS

  use Pencil::ConfigParser;
  my $parser = Pencil::ConfigParser->new('mycomputer.conf');
  my $parser = Pencil::ConfigParser->new('os/GNU_Linux.conf', 'mpi/mpich.conf');

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


=item B<new>(@filenames)

Create a new object that will accumulate all entries from the files listed
in @filenames (and their included files).
Parsing is done lazily, i.e. occurs at the point where results from
parsing are actually needed.

=item B<debug>([$debug])

With argument: set debugging flag.
Returns the debugging flag.

=item B<get_makefile_params>()

Return a hashref of the parmeters defined in the `Makefile' section:

  { KEY1 => value1, KEY2 => value2, ...}

=item B<get_makefile_keys>()

Return arrayref of all keys in the `Makefile' section, in the order in
which they occured int the config files.

=item B<get_makefile_args>()

Return array ref of Makefile arguments

  [ 'VAR1=val1', 'VAR2='val2', ... ]

that can be interpolated into a `make' command line

Note that whitespace in a value is currently _not_ escaped, so
the array is best used in a system() call with list argument.

=item B<get_runtime_params>()

Return a hashref of the parmeters defined in the `runtime' section:

  { key1 => value1, key1 => value2, ...}

=item B<get_runtime_keys>()

Return arrayref of all keys in the `runtime' section, in the order in
which they occured int the config files.

=item B<get_environment_params>()

Return a hashref of the parmeters defined in the `environment' section:

  { KEY1 => value1, KEY1 => value2, ...}

=item B<get_environment_keys>()

Return arrayref of all keys in the `environment' section, in the order in
which they occured int the config files.

=item B<get_environment_args>()

Return array ref of environment variable settings

  [ 'VAR1=val1', 'VAR2='val2', ... ]

that can be interpolated into an `env' or `export' command line.

Note that whitespace in a value is currently _not_ escaped, so
the array is best used in a system() call with list argument.

=back


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

This file is part of the Pencil Code and licensed under the GNU Public
License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=cut


# End of file
