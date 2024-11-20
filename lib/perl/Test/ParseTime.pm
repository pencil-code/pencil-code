# 20-Nov-2024: Kishore: outsourced from bin/reaper so that parse_time can be used in pc_auto-test

package Test::ParseTime;

use warnings;
use strict;
use vars qw($VERSION);

$VERSION = '0.1';

sub parse_time {
# Parse a time interval of the form
# `1h[our] 5m[in] 7s[ec]'
# `1:05:07'  (ditto)
# `1:05'     (same as `1h 5m')
#  into seconds.
    my ($string) = (@_);

    my ($hour, $min, $sec);

    if ($string =~ /:/) {
        ($hour, $min, $sec)
          = ($string =~ /^\s*([0-9]+):([0-9]+)(?::([0-9]+))?\s*$/)
            or die "Cannot parse <$string> in hh:mm[:ss] format\n";
    } else {
        ($hour, $min, $sec)
          = ($string =~
             /
                 ^\s*
                 (?:
                     ([0-9]+)    # capture hours
                     h[^0-9]* # h[ours]
                     \s*
                 )?              # hours are optional
                 (?:
                     ([0-9]+)    # capture minutes
                     m[^0-9]* # m[in[utes]]
                     \s*         # minutes are optional
                 )?
                 (?:
                     ([0-9]+)    # capture seconds
                     s[^0-9]* # s[ec[onds]]
                     \s*         # seconds are optional
                 )?
                 \s*
                 $
             /x
            ) or die
              "Cannot parse <$string> in [<hh>h] [<mm>m] [<ss>s] format\n";
    }

    die "Bad time format <string>\n"
      unless (defined($hour) || defined($min) || defined($sec));

    my $time = numeric_value($hour);
    $time = 60*$time + numeric_value($min);
    $time = 60*$time + numeric_value($sec);

    return $time;
}
# ---------------------------------------------------------------------- #
sub numeric_value {
# Return argument as number; if undefined, return 0
    my ($arg) = (@_);
    return ( defined($arg) ? "0" + $arg : "0" );
}

1;
__END__
