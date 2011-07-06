#
#                            DocExtractor.pm
#                            ---------------
#
# Description:
#   Extract Documentation strings following variable declarations and
#   Collect in LaTeX longtable environment for inclusion in the manual.
# Author: wd (wdobler#cpan:org =~ tr/:#/.@/)
# $Id$
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
#                             prefix   => qr/integer\s*::\s*idiag_/,
#                             whatelse => somethingelse);
#
# or in a hashref
#
#   Pencil::DocExtractor->new({ marker => qr/!\s*BC[XYZ]_DOC:/});
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
        if ($proto->isa('Pencil::DocExtractor')) {
            $parent = $proto;
        }
    } else {
        $class = $proto;
    }

    # Parse argument(s) (name => <nlname>); may be list or hashref
    my %args;
    if (@argv) {
        if (ref($argv[0]) eq 'HASH') { # parse($hashref)
            %args = %{$argv[0]};
        } else {                # parse(%hash) or parse(@list)
            %args = @argv;
        }
    }
    # Set object fields based on
    # 1. explicitly given parameters,
    # 2. existing fields (this could be an object constructor [as in
    #   $bcx_doc = $bc_doc->new(marker => qr/!\s*BCX_DOC:/);] or a class
    #   constructor)
    # 3. a default value

    $self->{MARKER}  = $args{marker}  || $parent->{MARKER} || qr/!\s*[A-Z]+_DOC:/ ;
    $self->{PREFIX}  = $args{prefix}  || $parent->{PREFIX} || qr/\s*/;
    $self->{DEBUG}   = $args{debug}   || $parent->{DEBUG}  || 0;
    $self->{VERBOSE} = $args{verbose} || $parent->{DEBUG}  || 0;
    # debugging implies verbosity:
    $self->{VERBOSE} = $self->{DEBUG} || 0;

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
                                      $self->{PREFIX},
                                      $self->{DEBUG},
                                      $self->{VERBOSE});
    (my $sfile = $file) =~ s{.*/}{}; # remove path
    $self->{DOC}{$sfile} = \@localdoc if (@localdoc);
    my $count = scalar @localdoc;

    return $count;
}

# ---------------------------------------------------------------------- #

sub write_mod_to_file {
#
#   $doc->write_mod_to_file(file)
#   $doc->write_mod_to_file(file          => 'filename',
#                       sort_files    => 1/0,
#                       print_empty   => 0/1,
#                       descr_width   => '0.7\textwidth',
#                       selfcontained => 0/1)
#
# Write LaTeX {longtable} environment of docstrings to given file.
# Just a convenience wrapper around longtable().
#
    my $self = shift();
    my @args = @_;
    my %args;
    # Parse arguments (sort_files => <true/false>, etc.); may be hash or hashref
    if (ref($args[0]) eq 'HASH') { # longtable($hashref)
        %args = %{$args[0]};
    } else {                    # longtable(%hash)
        %args = @args;
    }

    my $file = $args{file} or croak "write_to_file() needs a <file> argument";
    open(my $fh, "> $file") or croak "Cannot open $file for writing: $!";
    print $fh $self->module_table(@args);
    close $fh;
}

# ---------------------------------------------------------------------- #

sub write_to_file {
#
#   $doc->write_to_file(file)
#   $doc->write_to_file(file          => 'filename',
#                       sort_files    => 1/0,
#                       print_empty   => 0/1,
#                       descr_width   => '0.7\textwidth',
#                       selfcontained => 0/1)
#
# Write LaTeX {longtable} environment of docstrings to given file.
# Just a convenience wrapper around longtable().
#
    my $self = shift();
    my @args = @_;
    my %args;
    # Parse arguments (sort_files => <true/false>, etc.); may be hash or hashref
    if (ref($args[0]) eq 'HASH') { # longtable($hashref)
        %args = %{$args[0]};
    } else {                    # longtable(%hash)
        %args = @args;
    }

    my $file = $args{file} or croak "write_to_file() needs a <file> argument";
    open(my $fh, "> $file") or croak "Cannot open $file for writing: $!";
    print $fh $self->longtable(@args);
    close $fh;
}

# ---------------------------------------------------------------------- #

sub longtable {
#
#   $doc->longtable()
#   $doc->longtable(sort_files    => 1/0,
#                   print_empty   => 0/1,
#                   descr_width   => '0.7\textwidth',
#                   selfcontained => 0/1)
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

    my $text  = header(\@files, $args{selfcontained}, $args{descr_width});

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

    $text .= footer($args{selfcontained});

}

# ====================================================================== #
sub module_table {
#
#   $doc->module_table()
#   $doc->module_table(sort_files    => 1/0,
#                   print_empty   => 0/1,
#                   descr_width   => '0.7\textwidth',
#                   selfcontained => 0/1)
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

    my $text  = mod_header(\@files, $args{selfcontained}, $args{descr_width});

    foreach my $module (@files) {
        # Header line for each section of table
        $text .=
            "\\midrule\n";
#          . "  \\multicolumn{2}{c}{Module \\file{$module}} \\\\\n"
#          . "\\midrule\n";

        # Loop through variables
        my @file_docs = @{$docref->{$module}}; # (['var1', 'doc1'],
                                               #  ['var2', 'doc2'], ...)
        foreach my $vardocref (@file_docs) {
            my ($var,$doc) = @$vardocref;

            next unless ($print_empty || $doc =~ /\S/);

            # Indent continued lines, so LaTeX code is easier to read:
            $doc =~ s{\n}{\n                  }g;

            $text .= sprintf "  %-15s & %s \\\\\n", "\\var{$module}", $doc;
            $module=''
        }

    }

    $text .= footer($args{selfcontained});

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
    my $prefix = shift;
    my $debug  = shift || 0;
    my $verbose  = shift || 0;

    my @localdoc;
    my $file  = $module;
    my $count = 0;

    unless (open(MODULE, "< $file")) {
        carp "Cannot open $file for reading: $!\n";
        return ();
    }
    print STDERR "$module:\n" if ($verbose);
    my $no_line_to_continue = 1;

    # Store prefix for processing with the next line to cover cases
    # like
    #   case ('cop')
    #   ! BCZ_DOC: copy value
    my $saved_prefix_line = '';
    LINE: while(defined(my $line = <MODULE>)) {
        print STDERR "----------- input line $. ----------\n" if ($debug);
        if ($line =~ /\s*($prefix)/) {
            $saved_prefix_line = $1;
            print STDERR "Saving prefix line #",
              $., " = <", printable_substring($line), ">\n"
                if ($debug);
        }
        unless ($line =~ /$marker/) {
            # Whatever comes next, cannot be a continuation line:
            $no_line_to_continue = 1;
            next;
        }

        print STDERR "\$line = <", printable_substring($line), ">\n"
          if ($debug);
        my ($var,$misc,$docstring) = ('', '', '');

        my ($decl, $latex)
          = ($line =~ /^\s*(.*?)\s*$marker\s*(.*?)\s*$/);
        warn "Weird... \$line=<$line>\n" unless defined($latex);
        print STDERR "\$decl=<$decl>, \$latex=<$latex>\n" if ($debug);

        # Get stored declaration line if there was no $prefix in the
        # actual line:
        if ($decl eq '') {
            print STDERR "\$decl <-- <",
              printable_substring($saved_prefix_line), ">\n"
                if ($debug);
            $decl = $saved_prefix_line;
        }

        if ($decl =~ /$prefix/) {
            # there is a declaration part
            # -> not a continuation line
            print STDERR "docstring at ${module}:$.\n" if ($debug);

            my ($var,$rest)
              = ($decl =~ /$prefix\s*(.*?)\s*/);
            print STDERR "## \$var=<$var>, \$rest=<$rest>\n" if ($debug);
            $misc
              = ($rest =~ /^(.*?)(?:\s*=\s*[-+0-9]+\s*)/i);
            print STDERR "\$var=<$var>, \$msc=<$misc>\n" if ($debug);
            if (defined($misc) && $misc =~ /idiag_/i) {
                carp "In line $. of $file: "
                  . "multiple diagnostic variables in one line:\n";
                carp "  $var, $misc\n";
                next LINE;
            }
            unless (defined($var)) {
                carp "In line $. of $file: "
                  . "variable name not found:\n";
                next LINE;
            }
            print STDERR "    $var -> $latex\n" if ($debug);
            push @localdoc, [$var, $latex];
            $no_line_to_continue = 0;
            $count++;
        } else {
            #  no declaration part --> continuation line
            print STDERR "..continuation line at ${module}:$.\n" if ($debug);
            ## Append latex part to previous entry
            my ($var1,$latex1) = @{ pop @localdoc || [] }
                or next LINE; # nothing to append to
            push @localdoc, [$var1, "$latex1\n  $latex"];
        }
        $saved_prefix_line = '';
    }

    if ($count) {
        print STDERR "Found documentation for $count diagnostic variables\n"
          if ($debug);
    } else {
        print STDERR "No documentation found in $file\n"
          if ($debug);
    }

    return @localdoc;
}

# ---------------------------------------------------------------------- #

sub filter_doc {
# Remove (or not) items with empty documentation lines
    my $docref  = shift;
    my $verbose = shift || 0;

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
    if ($verbose) {
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
    my @files =  @{shift()};
    my $selfcontained = shift;
    my $descr_width = (shift || '0.7\textwidth');

    my $string =
        '%% $Id$' . "\n"
      . "%% This file was automatically generated by Pencil::DocExtractor,\n"
      . "%% so think twice before you modify it.\n"
      . "%%\n"
      . "%% Source files:\n%%   "
      . join("\n%%   ", @files) . "\n\n";

    if ($selfcontained) {
        $string .=
            "\n\\documentclass[12pt]{article}\n"
          . "\n"
          . "\\usepackage{longtable,booktabs}\n"
          . "\\usepackage{underscore}\n"
          . "\\usepackage{amsmath,newcent,helvet}\n"
          . "\\renewcommand{\\ttdefault}{cmtt} % Courier is too broad\n"
          . "\n"
          . "\\newcommand{\\file}[1]{`\\texttt{#1}'}\n"
          . "\\newcommand{\\var}[1]{\\textsl{#1}}\n"
          . "\n"
          . "\\newcommand{\\vekt}[1] {\\mathbf{#1}}\n"
          . "\n"
          . "\\newcommand{\\Av}{\\vekt{A}}\n"
          . "\\newcommand{\\Bv}{\\vekt{B}}\n"
          . "\\newcommand{\\cs}{c_{\rm s}}\n"
          . "\\newcommand{\\curl}{\\nabla\\times}\n"
          . "\\newcommand{\\jv}{\\vekt{j}}\n"
          . "\\newcommand{\\Strain}{\\boldsymbol{\\mathsf{S}}}\n"
          . "\\newcommand{\\uv}{\\vekt{u}}\n"
          . "\n"
          . "\\begin{document}\n\n"
    }

    $string .=
        "% ---------------------------------------------------------------- %\n"
      . "\\begin{longtable}{lp{$descr_width}}\n"
      . "\\toprule\n"
      . "  \\multicolumn{1}{c}{\\emph{Variable}} \& {\\emph{Meaning}} \\\\\n";

    return $string;
}

# ---------------------------------------------------------------------- #
sub mod_header {
#
# Print LaTeX longtable header
#
    my @files =  @{shift()};
    my $selfcontained = shift;
    my $descr_width = (shift || '0.9\textwidth');

    my $string =
        '%% $Id$' . "\n"
      . "%% This file was automatically generated by Pencil::DocExtractor,\n"
      . "%% so think twice before you modify it.\n"
      . "%%\n"
      . "%% Source files:\n%%   "
      . join("\n%%   ", @files) . "\n\n";

    if ($selfcontained) {
        $string .=
            "\n\\documentclass[12pt]{article}\n"
          . "\n"
          . "\\usepackage{longtable,booktabs}\n"
          . "\\usepackage{underscore}\n"
          . "\\usepackage{amsmath,newcent,helvet}\n"
          . "\\renewcommand{\\ttdefault}{cmtt} % Courier is too broad\n"
          . "\n"
          . "\\newcommand{\\file}[1]{`\\texttt{#1}'}\n"
          . "\\newcommand{\\var}[1]{\\textsl{#1}}\n"
          . "\n"
          . "\\newcommand{\\vekt}[1] {\\mathbf{#1}}\n"
          . "\n"
          . "\\newcommand{\\Av}{\\vekt{A}}\n"
          . "\\newcommand{\\Bv}{\\vekt{B}}\n"
          . "\\newcommand{\\cs}{c_{\rm s}}\n"
          . "\\newcommand{\\curl}{\\nabla\\times}\n"
          . "\\newcommand{\\jv}{\\vekt{j}}\n"
          . "\\newcommand{\\Strain}{\\boldsymbol{\\mathsf{S}}}\n"
          . "\\newcommand{\\uv}{\\vekt{u}}\n"
          . "\n"
          . "\\begin{document}\n\n"
    }

    $string .=
        "% ---------------------------------------------------------------- %\n"
      . "\\begin{longtable}{lp{$descr_width}}\n"
      . "\\toprule\n"
      . "  \\multicolumn{1}{c}{\\emph{Module}} \& {\\emph{Description}} \\\\\n";

    return $string;
}

# ---------------------------------------------------------------------- #

sub footer {
#
# Return LaTeX longtable footer
#
my $selfcontained = shift;
    my $string =
        "%\n"
      . "\\bottomrule\n"
      . "\\end{longtable}\n\n";

    if ($selfcontained) {
        $string .= "\n\\end{document}\n\n";
    }

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

sub printable_substring {
# Extract substring and quote newlines for diagnostic printing
    my $string = shift;
    my $length = shift || 60;

    $string =~ s{\n}{\\n}g;      # quote newlines
    $ string =~ s{\s\s+}{\\s\+}; # compactify whitespace
    my $oldlen = length($string);
    $string = substr($string,0,$length);
    substr($string,-3,3) = '...' if ($length<$oldlen);

    return $string;
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

and creates a LaTeX {longtable} environment for inclusion into the manual:

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
