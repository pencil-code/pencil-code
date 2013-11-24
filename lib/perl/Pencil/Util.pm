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
# E.g.
#   use Pencil::Util;
#   Pencil::Util::use_pencil_perl_module('Pencil::DocExtractor') or die;
#   my $diag_doc = Pencil::DocExtractor->new([...]);
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
                . "  (cd \$PENCIL_HOME; svn update lib)\n";
        } else {                # unknown error message
            croak "$@\n";
        }
    } else {
        return 1;               # everything was fine
    }
}

# ---------------------------------------------------------------------- #

sub notify {
# Inform the user that we are done.
    my ($phase) = @_;

    notify_audibly($phase);
    notify_visually($phase);
}

# ---------------------------------------------------------------------- #

sub notify_audibly {
# Tell the user that we are done.
    my ($phase) = @_;

    my $player = '';
    foreach my $pl (qw/play aplay mplayer/) {
        if (in_PATH($pl)) {
            $player = $pl;
            last;
        }
    }
    return unless $player;

    my @sounds = (
      "/usr/share/games/rocksndiamonds/levels/BD Dream/sounds/exit_opening.wav",
      "/usr/share/sounds/question.wav",
      "/usr/share/gnome-games/sounds/laughter.ogg"
    );
    foreach my $sound (@sounds) {
        if (-e $sound) {
            run_silently_in_bg($player, "'$sound'");
            return;
        }
    }

    print "\a";
    sleep 1;
    print "\a";
}

# ---------------------------------------------------------------------- #

sub notify_visually {
# Show the user that we are done.
    my ($phase) = @_;

    use Cwd;
    my $cwd = getcwd();
    if (in_PATH('zenity')) {
        run_silently_in_bg(
            'zenity', '--info', '--text', "'Done $phase in $cwd'"
        );
    } elsif (in_PATH('xmessage')) {
        run_silently_in_bg('xmessage', "Done $phase in $cwd");
    }
}

# ---------------------------------------------------------------------- #

sub in_PATH {
# Check whether an executable is available in the execution PATH
    my ($file) = @_;

    foreach my $path (split(/:/,$ENV{PATH})) {
        if (-x "$path/$file") { return 1; }
    }
    return 0;
}

# ---------------------------------------------------------------------- #

sub run_silently_in_bg {
# Run the given command in the background and suppress its output
    my @cmd = @_;

    warn "@cmd\n";
    system("@cmd > /dev/null 2>&1 &");
}

# ---------------------------------------------------------------------- #

# End of file
