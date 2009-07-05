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

package Pencil::ConfigParser;

use warnings;
use strict;
use Carp;
use vars qw($VERSION);

##use critic

$VERSION = '0.1';

# ---------------------------------------------------------------------- #
# ---------------------------------------------------------------------- #

1;

__END__


=head1 NAME

Pencil::ConfigFinder - Find the appropriate Pencil Code configuration file


=head1 SYNOPSIS

  use Pencil::ConfigFinder;
  my $config_file = Pencil::ConfigFinder::find_config_file();


=head1 DESCRIPTION

Pencil::ConfigFinder locates the best configuration file for the given
computer and directory.

=head2 Functions

C<Pencil::ConfigFinder> provides only one function:

=over 4

=item B<find_config_file()>

=item B<find_config_file($host_ID)>

Return the full path name of the (first) matching config file.

If $host_ID is defined, try it as host ID string before falling back on
other methods (see L<"The Host ID"> below).

If no matching file is found, exit with status 15.

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
I<[Command line options]>

=item 2.
If the file C<./.pencil/host-ID> exists, its first line (without
leading/trailing whitespace) is the host ID.

=item 3.
If the file ~/.pencil/host-ID exists, its first line (without
leading/trailing whitespace)is the host ID.

=item 4.
If the IP number of the current computer is not in the private
range and it is possible to determine its fully-qualified host name (i.e.
the host and domain name), then this is used as host ID.

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

For a given host ID, C<find_config_file()> looks for a config file. E.g. if
the host ID is workhorse.pencil.org, C<find_config_file()> will look for a
file C<workhorse.pencil.org.conf>. in the . E.g. if the host ID is
workhorse.pencil.org, configure will look for a file in the following
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
