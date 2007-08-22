#
#                            DocExtractor.pm
#                            ---------------
#
# Description:
#   Extract Documentation strings following variable declarations and
#   Collect in LaTeX longtable environment for inclusion in the manual.
# Author: wd (wdobler [at] cpan.org)
# $Date: 2007-08-22 15:56:46 $
# $Revision: 1.1 $
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Pencil::DocExtractor;

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
# Takes configuration parameters as argument, either in a (flattened) hash
#
#   Pencil::DocExtractor->new(marker   => qr/!\s*DIAG_DOC:/,
#                             whatelse => somethingelse);
#
# or in a hashref
#
#   Pencil::DocExtractor->new({ marker => qr/!\s*BC[XYZ]_DOC:/});
#
    my $proto = shift;          # either classref or object ref or string
    my @argv  = @_;

    my $class = ref($proto) || $proto;
    my $self = {};

    # Parse argument(s) (name => <nlname>); may be list or hashref
    my %config = ( marker => qr/!\s*[A-Z]+_DOC:/,
                   debug  => 0,
                   quiet  => 0,
                 );
    my %args;
    if (@argv) {
        if (ref($argv[0]) eq 'HASH') { # parse($hashref)
            %args = %{$argv[0]};
        } else {                # parse(%hash) or parse(@list)
            %args = @argv;
        }
    }
    # Overwrite default settings with explicitly given config parameters
    while (my ($key, $val) = each %args) {
        $config{$key} = $val;
    }
    #

    my $quiet = ($config{quiet} || 0);
    my $debug = ($config{debug} || 0);

    $self->{MARKER} = $config{marker};
    $self->{DEBUG}  = $config{debug};
    $self->{QUIET}  = $config{quiet};

    bless($self, $class);
    return($self);
}
# ====================================================================== #

##
## Methods
##

sub parse {
#
# Parse a file and collect the docstrings.
#
    my $self = shift();
    my $file = shift();

    carp("Cannot read file <$file>\n") unless (-r $file);

    my @localdoc = get_docs_from_file($file,
                                      $self->{MARKER},
                                      $self->{DEBUG},
                                      $self->{QUIET});
    (my $sfile = $file) =~ s{.*/}{}; # remove path
    $self->{DOC}{$sfile} = \@localdoc if (@localdoc);
    my $count = scalar @localdoc;

    return $count;
}

# ---------------------------------------------------------------------- #

sub longtable {
#
#   $doc->longtable()
#   $doc->longtable(sort_files  => 1/0,
#                   print_empty => 0/1)
#
# Output docstrings in a LaTeX {longtable} environment.
#
    my $self = shift;
    my @args = @_;

    my %args;
    # Parse arguments (sort_files => <true/false>, etc.); may be hash or hashref
    if (ref($args[0]) eq 'HASH') { # longtable($hashref)
        %args = %{$args[0]};
    } else {                    # longtable(%hash)
        %args = @args;
    }

    my $docref = $self->{DOC};
    my @files = keys %$docref;

    my $sort = 1;
    $sort = $args{sort_files} if defined $args{sort_files};

    my $print_empty = $args{print_empty} || 0;

    # Sort file names in pre-defined order
    if ($sort) {
        @files = sort { compare_modulenames_for_sorting($a)
                          <=> compare_modulenames_for_sorting($b)
                      } @files;
    }

    my $text  = header(@files);

    foreach my $module (@files) {
        # Header line for each section of table
        $text .=
            "\\midrule\n"
          . "  \\multicolumn{2}{c}{Module \\file{$module}} \\\\\n"
          . "\\midrule\n";

        # Loop through variables
        my @file_docs = @{$docref->{$module}}; # (['var1', 'doc1'],
                                               #  ['var2', 'doc2'], ...)
        foreach my $vardocref (@file_docs) {
            my ($var,$doc) = @$vardocref;

            next unless ($print_empty || $doc =~ /\S/);

            # Indent continued lines, so LaTeX code is easier to read:
            $doc =~ s{\n}{\n                  }g;

            $text .= sprintf "  %-15s & %s \\\\\n", "\\var{$var}", $doc;
        }

    }

    $text .= footer();

}

# ====================================================================== #

##
## Private utility subroutines:
##

sub get_docs_from_file {
#
# Extract documentation lines from one file and return list of array refs
# ( [var1, doc1], [var2, doc2], ... )
#
    my $module = shift;
    my $marker = shift;
    my $debug  = shift || 0;
    my $quiet  = shift || 0;

    my @localdoc;
    my $file  = $module;
    my $count = 0;

    unless (open(MODULE, "< $file")) {
        carp "Cannot open $file for reading: $!\n";
        return ();
    }
    print STDERR "$module:\n" unless ($quiet);
    LINE: while(defined(my $line = <MODULE>)) {
        next unless $line =~ /$marker/;
        my ($var,$misc,$docstring) = ('', '', '');

        my ($decl, $latex)
          = ($line =~ /^\s*(.*?)\s*$marker\s*(.*?)\s*$/);

        if ($decl ne '') {      # there is a declaration part
            print STDERR "docstring at ${module}:$.\n" if ($debug);
            ($var,$misc) = 
              ($decl =~
               /^integer\s*::\s*idiag_(\S+)(.*?)(?:\s*=\s*[-+0-9]+\s*)/i);
            if ($misc =~ /idiag_/i) {
                carp "In line $. of $file: "
                  . "multiple diagnostic variables in one line:\n";
                carp "  $var, $misc\n";
                next LINE;
            }
            print STDERR "    $var -> $latex\n" if ($debug);
            push @localdoc, [$var, $latex];
            $count++;

        } else {              #  no declaration part --> continuation line
            print STDERR "..continuation line at ${module}:$.\n" if ($debug);
            ## Append latex part to previous entry
            my ($var1,$latex1) = @{ pop @localdoc || [] }
                or next LINE; # nothing to append to
            push @localdoc, [$var1, "$latex1\n  $latex"];
        }
    }

    if ($count) {
        print STDERR "Found documentation for $count diagnostic variables\n"
          if ($debug);
    } else {
        print STDERR "Hmm, no documentation found in $file\n" unless ($quiet);
    }

    return @localdoc;
}

# ---------------------------------------------------------------------- #

sub filter_doc {
# Remove (or not) items with empty documentation lines
    my $docref = shift;
    my $quiet  = shift || 0;

    my %newdoc;
    my $empty    = 0;
    my $nonempty = 0;
    foreach my $module (keys %$docref) {
        my @file_docs = @{$docref->{$module}}; # (['var1', 'doc1'],
                                               #  ['var2', 'doc2'], ...)
        foreach my $vardocref (@file_docs) {
            my ($var,$doc) = @$vardocref;
            if ($doc =~ /\S/) {
                push @{$newdoc{$module}}, [$var,$doc];
                $nonempty++;
            } else {
                $empty++;
            }
        }
    }

    # Give statistical feedback
    unless ($quiet) {
        print STDERR "  (doc, undoc, tot) = ($nonempty, $empty, ",
          $nonempty+$empty, ")\n";
    }

    return \%newdoc;
}

# ---------------------------------------------------------------------- #

sub header {
#
# Print LaTeX longtable header
#
    my @files =  @_;

    my $string =
        "%% This file was automatically generated by Pencil::DocExtractor\n"
      . "%% at " . scalar localtime() . " from\n%%   "
      . join("\n%%   ", @files) . "\n"
      . "%% So think twice before you modify it.\n\n";

    $string .=
        "% ---------------------------------------------------------------- %\n"
      . "\\begin{longtable}{lp{0.7\\textwidth}}\n"
      . "\\toprule\n"
      . "  \\multicolumn{1}{c}{\\emph{Variable}} \& {\\emph{Meaning}} \\\\\n";

    return $string;
}

# ---------------------------------------------------------------------- #

sub footer {
#
# Return LaTeX longtable footer
#
    my $string =
        "%\n"
      . "\\bottomrule\n"
      . "\\end{longtable}\n\n";

    return $string;
}

# ---------------------------------------------------------------------- #

sub compare_modulenames_for_sorting {
# Comparison function that makes sure we get interesting modules (as
# defined by the author of this script) before more boring ones
    my $module = shift;

    (my $short = $module) =~ s/\.f90$//; # remove suffix

    my %mapping =
      (
       'cdata'    => 1,
       'hydro'    => 2,
       'density'  => 3,
       'entropy'  => 4,
       'magnetic' => 5,
       'pscalar'  => 6,
       'equ'      => 7,
      );
    my $infty = 100000;         # or will we get more modules...?

    return $mapping{$short} || $infty;
}

# ---------------------------------------------------------------------- #

1;
__END__


=head1 NAME

Pencil::DocExtractor - Extract doc strings from F90 files and create LaTeX
longtable

=head1 SYNOPSIS

  use Pencil::DocExtractor;
  my $diag = Pencil::DocExtractor->new();
  my $count = $diag->parse('hydro.f90');
  my $count = $diag->parse('entropy.f90');
  print $diag->longtable();


=head1 DESCRIPTION

Pencil::DocExtractor extracts F90 comments of the form

  integer :: idiag_uxmz=0       ! DIAG_DOC: $\left< u_x \right>_{x,y}$
                                ! DIAG_DOC:   \quad(horiz. averaged $x$
                                ! DIAG_DOC:   velocity)
  integer :: idiag_uxmxz=0      ! DIAG_DOC:  $\left< u_x \right>_{y}$ \\

or

  case ('p')
    ! BCX_DOC: periodic
    call bc_per_x(f,topbot,j)
  case ('s')
    ! BCX_DOC: symmetry, $f_{N+i}=f_{N-i}$;
    ! BCX_DOC: implies $f'(x_N)=f'''(x_0)=0$
    call bc_sym_x(f,+1,topbot,j)

and creates a LaTeX {longtable} environment fo inclusion into the manual:

  %% This file was automatically generated by extract-diag-doc
  %% from the src/*.f90 files.
  %% So think twice before you modify it.

  % ---------------------------------------------------------------- %
  \begin{longtable}{lp{0.7\textwidth}}
  \toprule
    \multicolumn{1}{c}{\emph{Variable}} & {\emph{Meaning}} \\
  \midrule
    \multicolumn{2}{c}{Module \file{hydro.f90}} \\
  \midrule
    \var{uxmz}      & $\left< u_x \right>_{x,y}$
                      \quad(horiz. averaged $x$
                      velocity) \\
    \var{uxmxz}     & $\left< u_x \right>_{y}$ \\
  %
  \bottomrule
  \end{longtable}

or

  %% This file was automatically generated by extract-diag-doc
  %% from the src/*.f90 files.
  %% So think twice before you modify it.

  % ---------------------------------------------------------------- %
  \begin{longtable}{lp{0.7\textwidth}}
  \toprule
    \multicolumn{1}{c}{\emph{Variable}} & {\emph{Meaning}} \\
  \midrule
    \multicolumn{2}{c}{Module \file{hydro.f90}} \\
  \midrule
    \var{p}       & periodic \\
    \var{p}       & symmetry, $f_{N+i}=f_{N-i}$;
                    implies $f'(x_N)=f'''(x_0)=0$ \\
  %
  \bottomrule
  \end{longtable}

=head2 Methods

=over 4


=item B<$doc-E<gt>new>()

Create a new object

=item B<$doc-E<gt>parse>(file)

Parse the file.

=item B<$doc-E<gt>longtable>()

=item B<$doc-E<gt>longtable>(%options)

=item B<$doc-E<gt>longtable>(\%options)

Return string containing the extracted documentation strings in a LaTeX
{longtable} environment.
Write this to a file and include it from a LaTeX document.

Additional I<options> are:

=over 8

=item B<sort_files>

If true (the default), sort file names, putting predefined important
modules first and sort remaining file names alphabetically.

=item B<print_empty>

If true, print lines for empty diagnostics, e.g. corresponding to

  integer :: idiag_dtd=0        ! DIAG_DOC:

etc.
If false (the default), skip such lines.

=back


=head1 BUGS AND LIMITATIONS

=over 4

=item *

None worth mentioning (so far).

=back


=head1 AUTHOR

Wolfgang Dobler <wdobler [at] cpan.org>


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2007, Wolfgang Dobler <wdobler [at] cpan.org>.
All rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same conditions as Perl or under the GNU General Public
License, version 3 or later.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=cut


# End of file
