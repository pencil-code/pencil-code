#
#                            Util.pm
#                            -------
#
# Description:
#   Some utility routines we repeatedly need from Perl
# Author: wd (wdobler#cpan:org =~ tr/:#/.@/)
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Pencil::Util;

use warnings;
use strict;
use Carp;
use vars qw($VERSION);

##use critic

$VERSION = '0.1';

# ---------------------------------------------------------------------- #

sub use_pencil_perl_module {
# Try to use a module from the Pencil Code lib/perl directory and give
# useful advice if the module cannot be found.
#
    my $module = shift;

    my $modulefile = $module;
    $modulefile =~ s{::}{/}g;
    $modulefile .= '.pm';

    eval "use $module";
    my $msg = $@;
    if ($msg) {               # Error loading module -- give useful advice
        my $error = join("\n", $msg);
        if ($error =~ m|Can't locate $modulefile in \@INC|) {
            croak "Perl: Can't locate $modulefile in \@INC\n"
                . "You need to check out the Pencil library directory"
                . " (including lib/perl):\n"
                . "  (cd \$PENCIL_HOME; cvs up -dA lib)\n";
        } else {                # unknown error message
            croak "$@\n";
        }
    } else {
        return 1;               # everything was fine
    }
}

# ---------------------------------------------------------------------- #

# End of file
