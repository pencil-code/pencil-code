#
#                            ConfigFinder.pm
#                            ---------------
#
# Description:

#   Find a Pencil Code appropriate configuration file for the current machine
#   and directory.
#
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Pencil::ConfigFinder;

use warnings;
use strict;
use Carp;
use vars qw($VERSION);
use constant NO_CONFIG_FILE_FOUND => 15;

##use critic

$VERSION = '0.1';

my $quiet = 0;
my $debug = 1;

# ---------------------------------------------------------------------- #

sub find_config_file {
#
# Try all host IDs listed in ConfigFinder's specs and return the first
# config file found.
# If no config file is found, exit with status 15.
#
    for my $host_id (get_host_ids()) {
        my $config_file = find_config_file_for_host_id($host_id);
        return $config_file if defined($config_file);
    }

    exit NO_CONFIG_FILE_FOUND;
}

# ---------------------------------------------------------------------- #

sub find_config_file_for_host_id {
#
# Return config file for the given host ID, or undef.
#
    my ($host_id) = @_;

    my @path = (
                "$ENV{HOME}/.pencil/config/computers",
                "$ENV{PENCIL_HOME}/config/computers"
               );
    for my $dir (@path) {
        debug("Trying dir <$dir>");

        my $file = "${dir}/${host_id}.conf";
        unless (-e $file) {
            debug("No such file: <$file>\n");
            next;
        }

        debug("Found file: <$file>\n");
        if (-f $file) {
            return $file;
        } else {
            warn "Not a regular file: <$file>\n";
        }
    }

    return undef;               # no file found
}

# ---------------------------------------------------------------------- #

sub get_host_ids {
#
# Return list of host IDs to try
#
    my @ids = ();
    add_host_id_from_file("./host-ID", \@ids);
    add_host_id_from_file("$ENV{HOME}/.pencil/host-ID", \@ids);
    add_host_id_from_fqdn(\@ids);
    add_host_id_from_scraping_system_info(\@ids);
    push @ids, strip_whitespace(`uname -s`);

    debug("get_host_ids: <" . join(">, <", @ids) . ">");
    return @ids;
}

# ---------------------------------------------------------------------- #

sub add_host_id_from_file {
#
# Take the first line from $file and append it to array
#
    my ($file, $ids_ref) = @_;

    unless (-x $file) {
        debug("No such file <$file>");
        return;
    }

    debug("Reading host id from file <$file>");
    my $line = first_line_from_file($file);
    push @$ids_ref, strip_whitespace($line) if defined($line);
}

# ---------------------------------------------------------------------- #

sub add_host_id_from_fqdn {
#
# Try finding the fully-qualified domain name and append it to array
#
    my ($ids_ref) = @_;

    my $fqdname = `hostname --fqdn`;
    chomp($fqdname);

    if ($fqdname =~ /^ [^.]+ \. .* \. [^.]+$/x) {
        debug("Fully-qualified domain name: <$fqdname>");
        push @$ids_ref, strip_whitespace($fqdname);
    } else {
        debug("Not a fully-qualified domain name: <$fqdname>");
    }
}

# ---------------------------------------------------------------------- #

sub strip_whitespace {
#
# Remove leading and trailing whitespace from a host ID
#
    my ($text) = @_;

    return $text unless defined($text);

    chomp($text);
    $text =~ m{^\s*(.*?)\s*$};
    return $1;

}

# ---------------------------------------------------------------------- #

sub first_line_from_file {
#
# Extract the first line from the given file
#
    my ($file) = @_;

    unless (-x $file) {
        log_msg("No such file: <$file>");
        return undef;
    }

    if (-r $file) {
        my $fh;
        unless (open($fh, "< $file")) {
            warn "Cannot open file <$file>\n";
            return undef;
        }
        my $line = <$fh>;
        chomp($line);
        return $line;
    } else {
        log_msg("Not readable: <$file>");
        return undef;
    }
}

# ---------------------------------------------------------------------- #

sub debug {
#
# Print the given line to STDERR if the $debug flag is set.
#
    my ($text) = @_;

    if ($debug) {
        chomp($text);
        print STDERR "DEBUG: $text\n";
    }
}

# ---------------------------------------------------------------------- #

sub log_msg {
#
# Print the given line to STDERR if the $debug flag is set.
#
    my ($text) = @_;

    if ((! $quiet) || $debug) {
        chomp($text);
        print STDERR "$text\n";
    }
}

# ---------------------------------------------------------------------- #

1;

__END__


=head1 NAME

Pencil::ConfigFinder - Find the appropriate Pencil Code configuration file


=head1 SYNOPSIS

  use Pencil::ConfigFinder;

  # Find config file using default algorithm
  my $config_file = Pencil::ConfigFinder::find_config_file();

  # Try host ID "toto" first, then fall back on default algorithm
  my $file = ( Pencil::ConfigFinder::find_config_file_for_host_id('toto')
               || Pencil::ConfigFinder::find_config_file() );


=head1 DESCRIPTION

Pencil::ConfigFinder locates the best configuration file for the given
computer and directory.

=head2 Functions

C<Pencil::ConfigFinder> provides only two functions:

=over 4

=item B<find_config_file()>

Return the full path name of the (first) matching config file for the
curent host.

If no matching file is found, exit with status 15.

=item B<find_config_file_for_host_id($host_ID)>

Return the full path name of the (first) matching config file for the
given host ID.
If not such file is found, return undef.

=back


=head1 ALGORITHM

=head2 The Host ID

A host ID is supposed to uniquely identify a computer.

For a computer with a permanent, global IP address, the host ID is
normally the fully-qualified domain name, like C<workhorse.pencil.org>,
but this can be overridden.

For computers without a fully qualified domain name (compute nodes on
inernal subnets or laptops), other sources of information are evaluated.

C<find_config_file()> tries the following host IDs, in this order:

=over 4

=item 1.
Command line options I<[not yet implemented]>

=item 2.
If the file C<./host-ID> exists, its first line (without
leading/trailing whitespace) is the host ID.
[This should become part of some larger per-run-directory configuration
setting, either in one file C<./pencil-config> with sections and the like,
or in its own file C<./pencil-config/host-ID>]

=item 3.
If the file ~/.pencil/host-ID exists, its first line (without
leading/trailing whitespace)is the host ID.
[Again: should this be one file under C<~/.pencil> or one section/line in
C<~/.pencil-config> or C<~/.pencilrc>?]

=item 4.
If the IP number of the current computer is not in the private range [it
looks like there is no fast, portable way to get hold of the ip number
(and there could be several)] and it is possible to determine its
fully-qualified host name (i.e. the host and domain name), then this is
used as host ID.

=item 5.
[Describe some desperate measures, generating a host ID from the uname,
/etc/issue, etc. uname will be needed anyway for the `default Linux' setup
and the like.]

=item 6.
If no configuration file for that host ID is found, the output
from `C<uname -s>' is tried as host ID.

=item 7.
If still no configuration file for that host ID is found, the host ID
`C<default>' is tried.

=item 8.
If still no configuration is found, Pencil::ConfigFinder aborts with a
specific error code.


=back

For each host ID, Pencil::ConfigFinder looks for a corresponding
configuration file (see L</"Locating the config file"> below).
If such a file is found, C<find_config_file()> exits and returns that
file's name.


=head2 Locating the config file

For a given host ID, C<find_config_file()> looks for a config file.

E.g. if the host ID is workhorse.pencil.org, C<find_config_file()> will
look for a file C<workhorse.pencil.org.conf>. in the following
directories:

=over 4

=item 1.
C<~/.pencil/config/computers>

=item 2.
${PENCIL_HOME}/config/computers

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

This program is free software; you can redistribute it and/or modify it
under the same conditions as Perl or under the GNU General Public
License, version 3 or later.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=cut


# End of file
