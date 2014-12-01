#
#                         ScriptTester.pm
#                         ---------------
#
# Description:
#   Run script tests inside a set of directories.
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Test::ScriptTester;

use warnings;
use strict;
use Cwd qw/getcwd abs_path/;
use Carp;
use File::Basename;
use File::Copy 'move';
use File::Find;
use Test::NumericFileComparator;
use vars qw($VERSION);

use feature 'say';

##use critic

$VERSION = '0.1';


=head1 NAME

Test::ScriptTester - Run script tests in a set of directories

=head1 SYNOPSIS

  use Test::ScriptTests;

  $tester = Test::ScriptTester->new(
      [dir1, dir2, ...],
      {'python' => 'python', 'idl' => '/usr/bin/gdl'}
      );

  $tester->run();

  my @tests = $tester->list_tests();

  $tester->run_tests(@tests[0..2);

  my %default_interpreters = Test::ScriptTests::get_default_interpreters()

  my $idl_interpreter = Test::ScriptTests::find_interpreter_for('idl')


=head1 DESCRIPTION

Scan the given directories for subdirectories named 'tests/'; in each
tests directory, run all script tests.

A I<script test> consists of a test file that

=over 4

=item *

is named <test_name>.<suff>, where <suff> is a known suffix <suff>
(currently supported: py, pro),

=item *

is executable and can be run from the tests directory (thus, many of the
test scripts read data from '../data/'),

=item *

when run, writes a file <test_name>.out in the same directory.

=item *

There exist a file <test_name>.ref in the same directory that defines the
reference data and possibly accuracy.

=back

Each test script is run (if we find an appropriate interpreter)
and the <test_name>.out file is compared to the reference data
<test_name>.ref .


=head2 Methods

=over 4

=cut


our (@default_types, @types, %type_map);


=item B<Test::ScriptTester-E<gt>new>($dirs, $interpreters)
      B<Test::ScriptTester-E<gt>new>($dirs)

Create a new object that searches the given directories

  $dirs = [dir1, dir2, ...]

and uses the given interpreters, e.g.

  {'python' => 'python', 'idl' => '/usr/bin/gdl'} )

If no interpreter map ref is given, use the default map as returned by
I<Test::ScriptTests::get_default_interpreters>().

Only test types listed in the interpreters map are run.

=cut

sub new {
#
#   Test::ScriptTester->new($dirs_ref [, $interpreters_ref]);
#   $tester->new($dirs_ref [, $interpreters_ref]);
#
    my $proto = shift;          # either classref or object ref or string
    my ($dirs_ref, $interpreters_ref)  = @_;

    my $self = {};
    my $class;
    my $parent = {};
    if (ref($proto)) {
        $class = ref($proto);
    } else {
        $class = $proto;
    }

    $self->{DIRS} = $dirs_ref;

    if (defined $interpreters_ref) {
        $self->{INTERPRETERS} = $interpreters_ref;
    } else {
        my %interpreters = get_default_interpreters();
        $self->{INTERPRETERS} = \%interpreters;
    }

    my %full_type_map = (       # all known suffixes and types
        '.py'  => 'python',
        '.pro' => 'idl',
        );
    my %type_map = ();
    while (my ($suffix, $type) = each %full_type_map) {
        if (exists $self->{INTERPRETERS}->{$type}) {
            $type_map{$suffix} = $type;
        }
    }
    $self->{TYPE_MAP} = \%type_map;

    $self->{SUFFIXES} = [ keys %{$self->{TYPE_MAP}} ];

    $self->{DEBUG} = 1;  # @@@

    bless $self, $class;
    return $self;
}


=item B<$tester-E<gt>run>()

Run all tests supported by this $tester.

Return counts ($good_count, $bad_count) of successful and failed tests.

=cut

sub run {
    my ($self) = @_;

    return $self->run_tests($self->list_tests());
}


=item B<$tester-E<gt>list_tests>()

List all tests supported by this $tester.
Returns a list of list refs

  ([$tests_dir1, $script1], ...)

=cut

sub list_tests {
    my ($self) = @_;

    my @test_dirs;
    foreach my $dir (@{$self->{DIRS}}) {
        push @test_dirs, $self->_find_test_dirs($dir);
    }

    my @tests;
    foreach my $testdir (@test_dirs) {
        foreach my $test_script ($self->_find_test_scripts($testdir)) {
            push @tests, [$testdir, $test_script];
        }
    }

    return @tests;
}


=item B<$tester-E<gt>run_tests>(@tests)

Run the given tests.

Each test is a list ref [$tests_dir, $script], where $script is the path
to an executable script file, either relative to $tests_dir, or absolute.

Return counts ($good_count, $bad_count) of successful and failed tests.

=cut

sub run_tests {
    my ($self, @tests) = @_;

    my ($good_count, $bad_count) = (0, 0);
    foreach my $test (@tests) {
        my $good = $self->_run_test($test);
        if ($good) {
            $good_count++;
        } else {
            $bad_count++;
        }
    }
    return ($good_count, $bad_count);
}


sub _run_test {
    my ($self, $test_ref) = @_;

    my ($tests_dir, $test) = @$test_ref;
    my ($file, $type) = @$test;
    my ($base, $dir, $suffix) = fileparse($file, @{$self->{SUFFIXES}});
    my $script = $base . $suffix;
    my $outfile = "${base}.out";
    my $reffile = "${base}.ref";

    $self->debug("Running $type script $script from $tests_dir / $dir");
    my $workdir = getcwd();
    chdir $tests_dir;

    backup($outfile);
    my $ok = 1;

    my $interpreter = $self->{INTERPRETERS}{$type};
    my @cmd = split(/\s+/, $interpreter);
    $ok &= (system(@cmd, $file) == 0);

    if (-e $reffile) {
        if (-e $outfile) {
            $ok &= compare_files($reffile, $outfile);
        } else {
            warn "Script $script did not write expected file $outfile\n";
            $ok = 0;
        }
    }

    chdir $workdir;
    return $ok;
}


sub debug {
    my ($self, @args) = @_;

    if ($self->{DEBUG}) {
        my $string = join(' ', @args);
        chomp($string);
        print "$string\n";
    }
}


=item @default_types

All supported test types.

=cut



=item B<Test::ScriptTester::get_default_interpreters>()

Return a hash

  ( type1 => interpreter1,
    type2 => interpreter1, ...
  )

representing the default interpreters (values) for all known test
types (keys).

=cut

sub get_default_interpreters {
    my %interpreters = (
        'python' => 'python',
        );

    my $idl_interpreter;
    if (_in_PATH('idl')) {
        $idl_interpreter = 'idl';
    } elsif (_in_PATH('gdl')) {
        $idl_interpreter = 'gdl';
    } elsif (_in_PATH('gnudl')) {
        $idl_interpreter = 'gnudl';
    }

    if (defined $idl_interpreter) {
        $interpreters{'idl'} = $idl_interpreter;
    }

    return %interpreters;
}


sub compare_files {
    my ($reference, $actual) = @_;

    my $comparator = Test::NumericFileComparator->new($reference);

    # Compare file to reference data
    my @message = $comparator->compare($actual);
    if (@message) {
        print "File $actual differs: @message\n";
        return 0;
    } else {
        return 1;
    }
}


sub backup {
# Move $file to $file.old if applicable
# An existing backup file will be overwritten without further ado.
    my ($file) = @_;

    if (-e $file) {
        move $file, "${file}.old";
    }
}


sub _in_PATH {
# Check whether an executable is available in the execution PATH
    my $file = shift;

    foreach my $path (split(/:/,$ENV{PATH})) {
        if (-x "$path/$file") { return 1; }
    }
    return 0;
}


=item B<find_interpreter_for>($test_type)

Return the default interpreter (path of an executable) for the given
$test_type, or undef.

=cut

sub find_interpreter_for {
    my ($self, $test_type) = @_;

}


# ---------------------------------------------------------------------- #


sub _find_test_dirs {
# Find all test directories at or below the given @top_dirs.
# We do not recurse further into identified test directories.
    my ($self, @top_dirs) = @_;

    my @dirs;
    for my $dir (@top_dirs) {
        $dir = abs_path($dir);
        if ($self->_is_test_dir($dir)) {
            push @dirs, $dir;
        } else {
            File::Find::find({
                    wanted => sub {
                        my $name = $File::Find::name;
                        if ($self->_is_test_dir($name)) {
                            push @dirs, $name;
                            my $dummy = $File::Find::prune;  # suppress
                                                             # 'used only
                                                             # once'
                                                             # warning
                            $File::Find::prune = 1;
                        }
                    },
                    follow => 1,       # follow symlinks
                    follow_skip => 2,  # ignore duplicates
                },
                $dir
            );
        }
    }

    return @dirs;
}


sub _is_test_dir {
# Does the given $path represent a test directory?
    my ($self, $path) = @_;

    my ($name, undef, undef) = fileparse($path);
    return (-d $path) && ($name =~ '^tests$');
}


sub _find_test_scripts {
# Find all test scripts below directory $dir.
# A test script has a known suffix and is executable.
# Return an array ref [$filename, $type], where
# $filename is the file name of the script relative to the test dir.
    my ($self, $dir) = @_;

    my @scripts;
    File::Find::find({
            wanted => sub {
                my $name = $File::Find::name;
                my $type = $self->_is_test_script($name);
                if ($type) {
                    push @scripts, [_make_relative($name, $dir), $type];
                }
            },
            follow => 1,       # follow symlinks
            follow_skip => 2,  # ignore duplicates
        },
        $dir
    );

    return @scripts;
}


sub _is_test_script {
# Does the given $path represent a test script?
# Return the test type if it does, '' otherwise
    my ($self, $path) = @_;

    my ($name, $dir, $suffix) = fileparse($path, @{$self->{SUFFIXES}});

    if ((-x $path) && $suffix) {
        return $self->{TYPE_MAP}->{$suffix};
    } else {
        return '';
    }
}


sub _get_suffixes {
# Map a list of file types to a list of supported suffices.
    my $self = shift;
    my @types = @_;

    my @suffixes;
    for my $type (@types) {
        push @suffixes, $self->_type_to_suffix($type);
    }
    return @suffixes;
}


sub _type_to_suffix {
# Map test type (e.g. 'idl') to suffix of executable files ('pro' in the
# example)
    my ($self, $type) = @_;

    my $suffix = $self->{TYPE_MAP}->{$type};
    if (defined $suffix) {
        return $suffix;
    } else {
        die "Unknown test type <$type>\n";
    }
}


sub _make_relative {
# Return $path relative to $dir.
# This is a dumb version that assumes (and verifies) that $path is a
# subpath of $dir.
# In cases where a good relative path would require leading '..', this
# routine will just cry foul.
    my ($path, $dir) = @_;

    my $abs_path = abs_path($path);
    my $abs_dir = abs_path($dir);
    my ($undef, $rel_path) = ($abs_path =~ m{^($abs_dir)/*(.*)});
    if (defined $rel_path) {
        return $rel_path;
    } else {
        croak("_make_relative() requires that $path is below \$dir");
    }
}


# ---------------------------------------------------------------------- #

=back

=head1 EXAMPLES

  use Test::ScriptTester;

  my $tester = Test::ScriptTester->new(
      ["$env{PENCIL_HOME}/tests"]
      );
  $tester->run();  # Run all test scripts under $PENCIL_HOME/tests/

=cut


1;

__END__



# End of file
