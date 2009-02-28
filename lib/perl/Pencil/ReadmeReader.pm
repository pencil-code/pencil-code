#
#                            ReadmeReader.pm
#                            ---------------
#
# Description:
#   Parse samples/XYZ/README files
# Author: wd (wdobler#cpan:org =~ tr/:#/.@/)
# $Id$
#
# This file is part of the Pencil Code and licensed under the GNU Public
# License version 3 or later; see $PENCIL_HOME/license/GNU_public_license.txt.
#

package Pencil::ReadmeReader;

use warnings;
use strict;
use Carp;
use vars qw($VERSION);

##use critic

$VERSION = '0.1';


=head1 NAME

Pencil::ReadmeReader - Parse README files from samples/XYZ directories

=head1 SYNOPSIS

  use Pencil::ReadmeReader;
  my $reader = Pencil::DocExtractor->new('samples/conv-slab/README');
  my $section = 'Maintainer';

  # Get list of content lines
  my @lines = $reader->get_content($section);

  # Get join("\n", @lines);
  my $nl_separated = $reader->get_content();

  # Hash ( $section => [content] )
  my %hash = $reader->hash();
  my @lines = @{$hash{$section}};

=head1 DESCRIPTION

Pencil::ReadmeReader parses README files of the form

  Section:
      Content line 1
      Content line 2
  [...]

e.g.

  Directory:
      ${PENCIL_HOME}/samples/conv-slab
  Maintainer:
      Maint Ainer 1 <maintainer1@somewhere.net>
      Maint Ainer 2 <maintainer2@somewhere.net>
  Added:
      04-Aug-2003
  Status:
      has been working since August 2003
  Recommended resolution:
      128x128x128 (nu=eta=1.3e-3)
  Comments:
      This was one of the first sample setups used in the Pencil Code,
      but has not been used in production papers. The setup is motivated
      by the early paper of Hurlburt & Toomre (1984). The present setup
      is described in Brandenburg et al. (1996).

=head2 Methods

=over 4

=cut

# ---------------------------------------------------------------------- #
##
## Object constructor
##

=item B<Pencil::ReadmeReader-E<gt>new>('I<file>')

Create a new object.
Parses I<file> at construction time.

=cut

sub new {
#
#   Pencil::ReadmeReader->new(file);
#   $reader->new(file);
#   $reader->new();
#
    my $proto = shift;          # either classref or object ref or string
    my @argv  = @_;

    my $self = {};
    my $class;
    my $parent = {};
    if (ref($proto)) {
        $class = ref($proto);
    } else {
        $class = $proto;
    }

    if (@argv == 1) {
        $self->{FILE}  = pop(@argv);
    } else {
        croak("Usage: Pencil::ReadmeReader->new(file)");
    }

    $self->{DEBUG} = 0;

    my %hash = parse($self->{FILE});
    $self->{HASH} = \%hash;

    bless($self, $class);
    return($self);
}
# ====================================================================== #

##
## Methods
##

=item B<$reader-E<gt>hash>()

Return hash of the form

   ( 'Section1' => \@content1, 'Section2' => \@content2, [...] )

where @content1 is an array containing the individual lines (without
leading whitespace or trailing newline) bbelonging to that field.

=cut

sub hash {
    my $self = shift();

    return %{$self->{HASH}};
}

# ---------------------------------------------------------------------- #

=item B<$reader-E<gt>get_content>($section)

Return content of the given section, or undef if the section doesn't exist.

In list context, return a list of content lines.
The lines are stripped of leading whitespace and the trailing newline
character.

In scalar context, return a string containing all lines joined with a
newline character.

=cut

sub get_content {
    my $self = shift();
    my ($section) = @_;

    my @lines = @{$self->{HASH}->{$section}};

    if (wantarray()) {
        return @lines;
    } else {
        return join("\n", @lines);
    }
}

# ====================================================================== #

##
## Private utility subroutines:
##

sub parse {
#
# Extract documentation lines from one file and return list of array refs
# ( [var1, doc1], [var2, doc2], ... )
#
    my ($file) = @_;

    unless (-r $file) {
        carp "Cannot open $file\n";
        return undef;
    }

    unless (open(FILE, "< $file")) {
        warn "Cannot open $file: $!\n";
        return undef;
    }

    my %hash;
    my $section = undef;
    while (defined(my $line = <FILE>)) {
        if ($line =~ /^(\S+):/) { # new section
            $section = $1;
            $hash{$section} = [];
            next;
        } else {                # content line
            next unless defined($section);
            chomp($line);
            $line =~ s/^\s*//;
            push @{$hash{$section}}, $line;
        }
    }

    return %hash;
}

# ---------------------------------------------------------------------- #

1;

__END__


=back


=head1 BUGS AND LIMITATIONS

=over 4

=item *

None worth mentioning (so far).

=back


=head1 AUTHOR

Wolfgang Dobler <wdobler#cpan:org =~ tr/:#/.@/>


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2007, Wolfgang Dobler.
All rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same conditions as Perl or under the GNU General Public
License, version 3 or later.


=head1 DISCLAIMER OF WARRANTY

Use completely at your own risk.


=cut


# End of file
