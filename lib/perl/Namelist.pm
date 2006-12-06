#
#                             Namelist.pm
#                             -----------
#
# Description:
#   Parse F90 namelist into a hash and export in different formats.
# Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
# $Date: 2006-12-06 22:46:17 $
# $Revision: 1.7 $

# Current test statistics:
# All tests successful, 1 subtest skipped.
# Files=21, Tests=118,  2 wallclock secs ( 1.75 cusr +  0.14 csys =  1.89 CPU)

package Namelist;

=head1 NAME

Ext::Namelist - Parse F90 namelist into hash and export in different formats

=head1 SYNOPSIS

  use Ext::Namelist;

  my $nl = Ext::Namelist->new() or die "Couldn't get object\n";

  $nl->parse("&runpars\nx=2,y=3\nvec1=1,2,3\nvec2=3*1.3\n/");

  # Dump in arbitrary order
  use Data::Dumper;
  print "Namelist ", $nl->name(), " has ", $nl->nslots(), " slots:\n";
  print Dumper($nl->hash());

  # Retain original order
  print "&",$nl->name(),"\n";
  my $nl_hash = $nl->hash();
  foreach my $var ($nl->order()) {
    print "  $var = $nl_hash{$var}\n";
  }
  print "/\n";

  # Read one namelist from file `one_list.nml'
  $nl->parse(file => "one_list.nml")

  # Read one namelist from file handle HANDLE [¡not implemented yet!]
  $nl->parse(file => HANDLE)

  # Operate on each namelist in file `start.in'
  while ($nl->parse(file => "start.in")) {
      print $nl->name(), "\n";
  }

  # Same thing using file handles
  open(FILE, "< $fname") or die "Couldn't open "$fname for reading\n";
  # or [open file using File::Handle [¡not implemented yet!]]
  while (! eof(FILE)) {
      $nl->parse(file => FILE);
      print $nl->name(), ", ";
  }

  # Read all namelists from file `start.in' and merge into one namelist
  # called `nlist'
  $nl->parse(file => 'start.in', deep => 1, namelist => 'nlist')

  # Write namelist in F90 namelist format
  print $nl->output()

  # Write namelist as IDL structure
  print $nl->output(format => 'idl', name => 'par2')


=head1 DESCRIPTION

Ext::Namelist is a module for parsing Fortran90 namelists into hashs and
re-exporting these hashs in different formats. Currently, the following
data types are supported:

=over 4

=item *

integer

=item *

float/double

=item *

complex numbers

=item *

strings [character(LEN=*)], possibly containing all sorts of quotation
marks

=item *

logical

=back

The following backends exist for re-exporting (or importing into other
languages):

=over 4

=item *

F90 namelist

=item *

IDL struct

=back

=head2 Methods

=over 4

=item B<$fh-E<gt>new>()

Create a new namelist object

=item B<$fh-E<gt>parse>(I<string>)

=item B<$fh-E<gt>parse>(file =>(I<fname>|I<FHANDLE>))

=item B<$fh-E<gt>parse>(file => (I<fname>|I<FHANDLE>) [, I<options> ])

Parse I<string> or the file represented by I<fname> or I<FHANDLE> (a file
handle or File::Handle object [not yet implemeted]); return number of
parsed namelists (often 1 or 0 [e.g. for an empty file]), or undef if
parsing failed.
Additional I<options> are:

=over 8

=item B<deep>

If true, read all namelists from string or file and return hash

    { nl1 => $nl1_obj,
      nl2 => $nl2_obj,
      ...,
      nl_order_ => ['nl1', 'nl2', ...]
    }

[Implement nl_order_ stuff]

=item B<merge>

If true, read all nemelists from string or file and merge into one
namelist

=item B<name>

With B<merge>, set name of resulting namelist (default: name of first
namelist read)

=item B<dups_ok>

With B<merge>, don't warn if new slots have same names as existing slots
[not yet implemented]

=back

=item B<$fh-E<gt>hash>()

Return namelists as Perl hash

=item B<$fh-E<gt>nslots>()

Return number of slots in namelist

=item B<$fh-E<gt>name>()

Return name of namelist

=item B<$fh-E<gt>order>()

Return list of variables in original order

=item B<$fh-E<gt>output>(format => I<format>)

=item B<$fh-E<gt>output>({format => I<format>})

Write namelist in given I<format>. Currently supported formats are
 `f90' (default), and `idl'


=back


=head1 HASH FORMAT

The B<hash> method returns a hash reference of the following structure:

    { 'name of var1' => { 'value' => [ value1, value2, ..],
                          'type'  => numerical_type,
                          'stype' => "type string"
                        },
      'name of var2' => { 'value' => [ value1, value2, ..],
                          'type'  => numerical_type
                          'stype' => "type string"
                        },
      ...
    }

e.g.

    { 'xyz0' => { 'value' => [ 0., -3.141593, 0. ],
                  'type'  => 6,
                  'stype' => 'single precision float'
                },
      'nt'   => { 'value' => [ '1000' ],
                  'type'  => 4,
                  'stype' => 'integer'
                }
    }


Here I<numerical_type> is a number identifying each data type, while
I<stype> is a textual description of the given data type.


=head1 TO DO

=over 4

=item 1.

Method for reading file / file handle

=item 2.

More output methods:

=over 8

Matlab, octave , ...

=back

=back


=head1 BUGS AND LIMITATIONS

=over 4

=item *

No user-defined types (records) are supported, so if you have these LaTeX
comment characters in your namelist data, you are out of luck.

=back

=cut


use strict;
use Carp;
use vars qw($VERSION);

# Cannot use use Perl5.8's constant { x => 1, y=>2 , ..} because 5.6
# is very popular still
#
# Possible states of parser [used at all?]
use constant  UNDEF   => -1;
use constant  START   =>  0;	# initial state of parser
use constant  VAR     =>  1;	# at beginning of variable name
use constant  VALUE   =>  2;	# at beginning of value
use constant  SQUOTE  =>  3;	# in string after opening single quote
use constant  DQUOTE  =>  4;	# in string after opeing double quote
use constant  BRACKET =>  5;	# after opening bracket (e.g. complex number)
use constant  COMMENT =>  6;	# after exclamation mark (F90 comment)
use constant  NL_END  =>  7;	# after closing `/'
#
# F90 data types
use constant  UNKNOWN   => 0;
use constant  SQ_STRING => 1;
use constant  DQ_STRING => 2;
use constant  LOGICAL   => 3;
use constant  INTEGER   => 4;
use constant  FLOAT     => 5;	# a float here can be single or double
use constant  SINGLE    => 6;
use constant  DOUBLE    => 7;
use constant  COMPLEX   => 8;
use constant  DCOMPLEX  => 9;
use constant  MULTIPLE  => 20;
#
use constant  ID        => 100;	# variable name (_not_ a data type)


$VERSION = '0.3';

# general float:
my $numeric = '(?:[+-]?)(?=\d|\.\d)\d*(?:\.\d*)?(?:[EeDd](?:[+-]?\d+))?';
# float:
my $numeric_e = '(?:[+-]?)(?=\d|\.\d)\d*(?:\.\d*)?(?:[Ee](?:[+-]?\d+))?';
# double:
my $numeric_d = '(?:[+-]?)(?=\d|\.\d)\d*(?:\.\d*)?(?:[Dd](?:[+-]?\d+))?';
# float with decimal point, but w/o  exponential part:
my $float = '(?:[-+]?(?:\d+\.\d*|\d*\.\d+))';

## Regexps for the different data type values. Make sure all brackets are
## marked grouping-but-non-capturing (?:...), or else the parsing
## algorithm will fail.
my @regexp;
$regexp[SQ_STRING] = "'(?:[^']|'')*'"; # even covers 'toto''s quote'
$regexp[DQ_STRING] = "\"(?:[^\"]|\"\")*\""; # ditto for double quotes
$regexp[DCOMPLEX]  = "\\(\\s*$numeric_d\\s*,\\s*$numeric_d\\s*\\)";
$regexp[COMPLEX]   = "\\(\\s*$numeric\\s*,\\s*$numeric\\s*\\)";
$regexp[LOGICAL]   = "(?:T|F|\\.(?:true|TRUE|false|FALSE)\\.)";
$regexp[MULTIPLE]  = "[0-9]+\\*"; # also need special treatment...
$regexp[INTEGER]   = "[+-]?[0-9]+";
$regexp[DOUBLE]    = "$numeric_d";
$regexp[SINGLE]    = "$numeric_e";
$regexp[FLOAT]     = "$float";
$regexp[ID]        = "[a-zA-Z](?:[a-zA-Z0-9_])*"; # allowed namelist/var. names

## Corresponding regexp for compatible type class (numeric, complex, ...)
my @regexp2 = @regexp;		# same regexp by default
$regexp2[DCOMPLEX]  = "\\(\\s*$numeric\\s*,\\s*$numeric\\s*\\)";
$regexp2[COMPLEX]   = "\\(\\s*$numeric\\s*,\\s*$numeric\\s*\\)";
$regexp2[INTEGER]   = "$numeric";
$regexp2[DOUBLE]    = "$numeric";
$regexp2[SINGLE]    = "$numeric";
$regexp2[FLOAT]     = "$numeric";

# Hash for looking up symbolic names for type constants. The constants are
# only expanded as numbers if adding 0.
my %stypes = ( UNKNOWN   + 0 => 'unknown',
	       SQ_STRING + 0 => 'single-quote string',
	       DQ_STRING + 0 => 'double-quote string',
	       LOGICAL   + 0 => 'logical',
	       INTEGER   + 0 => 'integer',
	       FLOAT     + 0 => 'unspecified float',
	       SINGLE    + 0 => 'single precision float',
	       DOUBLE    + 0 => 'double precision float',
	       COMPLEX   + 0 => 'complex number',
	       DCOMPLEX  + 0 => 'double precision complex number',
	       MULTIPLE  + 0 => 'multiple data (array)',
	     );

# Global variables related to output() method:
my ($cmplx_pref,$cmplx_suff);

# ---------------------------------------------------------------------- #
##
## Object constructor
##
## Internal structure of Namlist objects (update me):
##   DATA    -- variable names, values, and types (see below)
##   ORDER   -- ordered list of variable names
##   NSLOTS  -- number of slots
##   NAME    -- name of namelist
##   PARSED_ -- flag indicating that argument has been parsed
##   DEBUG_  -- debug flag
##
## Structure of DATA slot: Note: One namelist object holds only one
## namelist -- use {$nl1, $nl2, ..} to group them.
##
##   $self->DATA  = data_hash;
##   data_hash = { 'name of var1' => { 'value' => [ value1, value2, ..],
##                                     'type'  => numerical_type,
##                                     'stype' => "type string"
##                                   }
##                 'name of var2' => { 'value' => [ value1, value2, ..],
##                                     'type'  => numerical_type
##                                     'stype' => "type string"
##                                   }
##                 ...
##               };
sub new {
    my $proto = shift;		# either classref or object ref or string
    my $class = ref($proto) || $proto;
    my $self = {};

    my %data   = ();
    my @order  = ();
    my $nslots = undef;
    my $name   = undef;
    my $parsed = 0;
    my $debug  = 0;

    my $short_usage =
        "Usage:\n" .
        "  Namelist::new()\n" ;
    #    unless($file) {
    #	carp $short_usage;
    #	return undef;
    #    }


    ##
    ## Finish the object
    ##
    # public data of object
    $self->{DATA}   = \%data;
    $self->{ORDER}  = \@order;
    $self->{NSLOTS} = $nslots;
    $self->{NAME}   = $name;

    # internal data
    $self->{PARSED_} = $parsed;
    $self->{DEBUG_}  = $debug;

    bless($self,$class);
    return($self);
}

# ---------------------------------------------------------------------- #

##
##  Methods
##

sub parse {
#
#   $obj->parse($text)
#   $obj->parse(file => $filename)
#   $obj->parse(file => $filehandle)
#   $obj->parse(file => $filename|$filehandle, merge => 1[, name => $name])
#   $obj->parse({file => $filename|$filehandle, merge => 1[, name => $name]})
#   $obj->parse({text => $textstring, merge => 1[, name => $name]})
#   $obj->parse(file => $filename|$filehandle, deep => 1)
#   $obj->parse({file => $filename|$filehandle, deep => 1})
#   $obj->parse({text => $textstring, deep => 1})
#
# Parse text or file containing F90 namelist(s)
#
    my $self = shift;
    my @args = @_;		# can't use shift() since we change value
                                # of $text

    my $state = START;
    my $debug = $self->{DEBUG_};

    my %args;
    my $text = '';

    # Parse arguments (file => <filename>, etc.); may be single string,
    # list, hash or hashref
    if (ref($args[0]) eq 'HASH') { # parse($hashref)
	%args = %{$args[0]};
    } else {
	if (@_ == 1) {		# parse($string)
	    $text = $args[0];
	} else {		# parse(%hash) or parse(@list)
	    %args = @args;
	}
    }
    my $file  = ($args{file}  || ''  );
    my $merge = ($args{merge} || '0' );
    my $name  = ($args{name}  || '' );
    my $deep  = ($args{deep}  || '0' );

    # Get text from file if necessary
    $text     = ($args{text}  || $text );
    if ($text eq '') {
	die "\$self->parse(): need text or file argument\n"
	    unless ($file ne '');
	local $/ = undef;
	open(FH, "< $file") or croak "Cannot open file <$file> for reading";
	$text = <FH>;
	close(FH);
    }

    my ($nslots,%hash,@order);
    my $hashref = \%hash;

    if ($deep) {
	# Not implemented yet; I am not sure this can work here with all
	# the $self stuff, as parse(deep => 1) will not return a namelist
	# object...
	croak "option `deep' not implemented yet\n";
    }

    if ($merge) {		# Read all namelists in string/file
	while ($text) {
	    my ($name1,$nslots1,@order1);
	    my $href = parse_namelist(\$text,\$name1,\$nslots1,\@order1,$debug);
	    if ($name1 ne '') {	# really read a namelist
		$name ||= $name1; # choose first name if not set yet
		# push @order, @order1;
		# %$hashref = (%$hashref, %$href);
		$nslots += $nslots1;
		# Eliminating repeated slots and warn if values don't agree
		slot: foreach my $slot (@order1) {
		    if (defined($$hashref{$slot})) { # slot already known
			my @val1=@{$$hashref{$slot}{'value'}};
			my @val2=@{$$href{$slot}{'value'}};
			while (@val1 and @val2) {
			    my $v1 = pop(@val1);
			    my $v2 = pop(@val2);
			    if ($v1 ne $v2) {
				carp "WARNING: Conflicting slots $slot = [@{$$hashref{$slot}{'value'}}] vs. [@{$$href{$slot}{'value'}}]\n";
				next slot;
			    }
			}
		    } else {	# new slot
			push @order, $slot;
			$$hashref{$slot} = $$href{$slot};
			$nslots++;
		    }
		}
	    }
	}
    } else {
	$hashref = parse_namelist(\$text,\$name,\$nslots,\@order,$debug);
    }

    $self->{NAME}   = $name;
    $self->{NSLOTS} = $nslots;
    $self->{ORDER}  = \@order;
    $self->{DATA}   = $hashref;

    if ($debug) {
	print STDERR
	  "Namelist::parse: Successfully read namelist <$name>\n";
	print STDERR "=================================\n";
    }

    $self->{PARSED_} = 1;

    # Return number of slots read
    $self->{NSLOTS};
}

sub hash {
# Return hash with parsed namelist contents
    my $self = shift;
    $self->{DATA};
}

sub nslots{
# Get number of slots in namelist
    my $self = shift();
    $self->{NSLOTS}
}

sub name{
# Get name of parsed namelists
    my $self = shift();
    $self->{NAME};
}

sub order{ # FIXME
# Return array of variable names in order
    my $self = shift();
    $self->{ORDER}
}

sub output{
# Write namelist in specified format (defaults to 'f90')
    my $self = shift();

    # Optional arguments:
    #   format   => format   ('f90' [default] or 'idl')
    #   name     => nl_name  (name of nlist/struct [default: get from nlist])
    #   trim     => 0/1      (trim strings?)
    #   double   => 0/1      (mark all floats as double precision?)
    #   oneline  => 0/1      (write all in one line? [only for some formats])
    #   maxslots => N        (similar to oneline, but split every N slots)
    my @args = @_;

    # Parse arguments (file => <filename>, etc.); may be list, hash or hashref
    my %args;
    if (ref($args[0]) eq 'HASH') {
	%args = %{$args[0]};
    } else {
	%args = @args;
    }
    my $format   = ($args{format}   || 'f90');
    my $name     = ($args{name}     || '');
    my $trim     = ($args{trim}     || 0);
    my $double   = ($args{double}   || 0);
    my $oneline  = ($args{oneline}  || 0);
    my $maxslots = ($args{maxslots} || 0);
    $oneline = 1 if ($maxslots);

    # Sanity check
    unless ($self->{PARSED_}) {
	croak "Called method output() on unparsed namelist";
	return undef;
    }

    # Get name of namelist(s?)
    my ($name1,$hashref) = (%{$self->{DATA}}); # hash (name=>valhash) ->
                                               # 2-element array; should
                                               # possibly be a loop over all
                                               # name=>hash pairs?
    # Format-dependent settings
    my $header = '';
    my ($footer,$lpref,$lpref2,$lsuff);
    ($cmplx_pref,$cmplx_suff) = ("", ""); # default complex delimiters
    if      (lc($format) eq 'f90') {
	$header = "&$name\n" if ($name);
	$lpref  = "";		# line prefix
	$lpref2 = "  ";		# line prefix for slot lines
	$lsuff  = ",\n";	# line suffix
	$footer = "\n/\n";
    } elsif (lc($format) eq 'idl') {
	$header = "$name = { \$\n" if ($name);
	$lpref  = "  ";		# line prefix
	$lpref2 = "    ";	# line prefix for slot lines
	$lsuff  = ", \$\n";	# line suffix
	$footer = " \$\n$lpref}\n";
	$cmplx_pref = "complex"; # complex number prefix
	if ($oneline) {
	    $header = "{ ";
	    $lpref  = "";
	    $lpref2 = "";
	    $lsuff  =", ";
	    $footer =" }";
	}
    } else                         {
	croak "output(): Format <$format> unknown";
	return undef;
    }

    # Combine slots, one per line
    my @slots;
    my $slot;
    foreach my $var (@{$self->{ORDER}}) {
	my $valhash = $self->{DATA}->{$var};
	my @vals = @{$$valhash{'value'}};
	my $type = $$valhash{'type'};
	if ($trim) {
	    @vals = map { s/\s*$//; $_ } @vals;
	}
	if ($double) { # replace E[0-9]+ by, or append `D0' where necessary
	    if (($type == FLOAT)   ||
		($type == SINGLE)  ||
		($type == DOUBLE)  ||
		($type == COMPLEX) ||
		($type == DCOMPLEX))  {
		@vals = map { s/[eEdD]/D/; $_ } @vals;
		@vals = map { s/(^|\s|,)($float)($|\s|,)/$1$2D0$3/g; $_ } @vals;
		@vals = map { s/(\(\s*)($float)(\s*,\s*)($float)(\s*\))/$1$2D0$3$4D0$5/g; $_ } @vals;
	    }
	}
	$slot = assign_slot_val($var,\@vals,$format,$type);
	push @slots, $slot;
    }

    # Combine data
    my $string='';
    unless ($maxslots) {
	$string = "$lpref$header$lpref2"
	          . join("$lsuff$lpref2", @slots)
		  . "$footer";
    } else {			# chop into chunks of $maxslots
	while (@slots) {
	    my @firstslots;
	    foreach my $i (1..$maxslots) {
		push @firstslots, shift @slots if (@slots);
	    }
	    $string .= "$lpref$header$lpref2"
	               . join("$lsuff$lpref2", @firstslots)
		       . "$footer ";
	}
    }

    $string;
}

sub debug{
#
#   $obj->debug(1)     # debugging on
#   $obj->debug(0)     # debugging off
#
# Undocumented: Set/get debug flag
    my $self = shift();
    if (@_) { $self->{DEBUG_} = shift };
    $self->{DEBUG_}
}


# ====================================================================== #

## Private utility subroutines:

sub parse_namelist {
#
# Parse first F90 namelist from text string; return reference to hash
#
#   parse_namelist(\$text,\$name,\$nslots,\@order,$debug);
#

    my $textref   = shift;
    my $nameref   = shift;
    my $nslotsref = shift;
    my $orderref  = shift;
    my $debug     = shift;

#$debug=1;

    my %hash;
    my $nslots = 0;
    my $state  = START;
    my $id = $regexp[ID];	# allowed namelist/variable names

    my ($status,$var,@values,$type);

    ## Reset to reasonable default values
    $$nslotsref = 0;
    @$orderref  = ();
    $$nameref   = '';

    ## Get name of nl
    $$nameref = extract_nl_name($textref,$debug) or return;

    $status = VAR;
    ## Extract variable slots
    my $text = $$textref;
    while ($text ne '') {
	print STDERR "--------------------\nTop of while loop..\n" if ($debug);
	strip_space_and_comment($text);
	if ($text =~ s/^($id)\s*=\s*//s) { # string starts with <var=...>
	    $var = lc($1);
	    $status = VALUE;
	    if ($debug) {
		print STDERR "parse_namelist 1: \$var=<$var>\n";
		print STDERR "parse_namelist 1: \$text=<",
		  printable_substring($text,50), ">\n";
	    }
	    @values = get_value(\$text,\$type,$var,$debug); # drop $debug here..
	    $nslots++;
	    push @$orderref, $var;
	} elsif ($text =~ s{\s*/\s*}{}) { # string is </>
	    $status = NL_END;
	    last;		# end of namelist
	} else {
	    show_error("Expected var=[...] not found ","",$text,1);
	    return;
	}

	print STDERR "[",join(',',@values), "] -> \$hash{$var}\n" if ($debug);
	my $stype = ($stypes{$type} || 'Type inconsistency!');
	$hash{$var} = { type  => $type,
			stype => $stype,
			value => [@values]
		      };
    }

    carp "Aborted parsing at <",
         printable_substring($text,50),">\n",
         "trying to read slot `$var'\n" unless ($status == NL_END);

    print STDERR
      "parse_namelist: Namelist <$$nameref> parsed succesfully\n"
	if ($debug);

    $$textref   = $text;	# propagate remainder of $text back
    $$nslotsref = $nslots;	# propagate number of slots back
    \%hash;
}

# ---------------------------------------------------------------------- #
sub extract_nl_name {
# Extract namelist name (the part starting with `&')

    my $textref = shift;
    my $debug   = shift;

    my $text = $$textref;
    my $name;
    my $id = $regexp[ID];	# allowed namelist/variable names

    print STDERR "extract_nl_name 1: \$text=<",
                 printable_substring($text,50),">\n" if ($debug);
    strip_space_and_comment($text);

    print STDERR "extract_nl_name 2: \$text=<",
                 printable_substring($text,50), ">\n" if ($debug);
    if ($text =~ s/^&($id)//) {
	$name = lc($1);
    } else {			# empty (comment/whitespace) or erroneous
	if ($text eq '') {
	    print STDERR "Empty text (at most some comments)" if ($debug);
	    $$textref = $text; # propagate remainder of $text back
	    return '';
	} else {
	    show_error("Namelist does not start with &\n","",$text,1);
	    return '';
	}
    }
    strip_space_and_comment($text);

    if ($debug) {
	print STDERR "extract_nl_name 3: \$name=<$name>\n";
	print STDERR "extract_nl_name 3: \$text=<",
		     printable_substring($text,50), ">\n";
    }

    $$textref = $text; # propagate remainder of $text back
    $name;
}

# ---------------------------------------------------------------------- #

sub strip_space_and_comment {
# Strip leading space and anything from possible leading exclamation mark
# til end of line.
    $_[0] =~ s/^(\s*(![^\n]*)?)*//s;
}

# ---------------------------------------------------------------------- #

sub get_value {
# Extract one or several values from string that starts off immediately
# after the equal sign of a slot assignment
    my $txtptr  = shift;
    my $typeptr = shift;
    my $varname = shift;
    my $debug   = shift;    # Need to somewhow get rid of this argument...
    my $text = $$txtptr;
    my @values;

    strip_space_and_comment($text); # (are comments really allowed here?)
    my $type = infer_data_type($text);
    if ($debug) {		# pretty-printing of type
	my @tp;
	$tp[UNKNOWN]   = 'UNKNOWN';
	$tp[SQ_STRING] = 'SQ_STRING';
	$tp[DQ_STRING] = 'DQ_STRING';
	$tp[LOGICAL  ] = 'LOGICAL';
	$tp[INTEGER  ] = 'INTEGER';
	$tp[FLOAT    ] = 'FLOAT';
	$tp[SINGLE   ] = 'SINGLE';
	$tp[DOUBLE   ] = 'DOUBLE';
	$tp[COMPLEX  ] = 'COMPLEX';
	$tp[DCOMPLEX ] = 'DCOMPLEX';
	$tp[MULTIPLE ] = 'MULTIPLE';
	print STDERR
	  "Found data of type $type (", $tp[$type], ") in <",
	  printable_substring($text,40), ">\n";
    }

    if ($type == UNKNOWN) {
	show_error("Cannot identify data type","$varname=","$text",1);
    }

    # Extract data
    my $multiregexp = $regexp[MULTIPLE];
    my $re_type = $regexp2[$type];

    while ($text =~ s/^($multiregexp)?($re_type)\s*(?:(,|\s+|\n|!)|(?=\/))\s*//s) {
	my $mul = ($1 || 1);
	my $val = $2;
	my $rest = ($3 || '');
	$mul =~ s/\*//;
	if ($debug) {
	    print STDERR "\$mul=<$mul>, \$val=<$val>\n";
	    print STDERR "\$rest=<", printable_substring($rest,2),
	      ">, \$text=<", printable_substring($text,40), ">\n";
	}

	# `Widen' data type if necessary (e.g. integer -> double for
	# `var = 1, 2.D0')
	my $valtype = infer_data_type($val);
	$type = $valtype if ($valtype > $type);
	if ($debug) {		# pretty-printing of type
	    my @tp;
	    $tp[UNKNOWN]   = 'UNKNOWN';
	    $tp[SQ_STRING] = 'SQ_STRING';
	    $tp[DQ_STRING] = 'DQ_STRING';
	    $tp[LOGICAL  ] = 'LOGICAL';
	    $tp[INTEGER  ] = 'INTEGER';
	    $tp[FLOAT    ] = 'FLOAT';
	    $tp[SINGLE   ] = 'SINGLE';
	    $tp[DOUBLE   ] = 'DOUBLE';
	    $tp[COMPLEX  ] = 'COMPLEX';
	    $tp[DCOMPLEX ] = 'DCOMPLEX';
	    $tp[MULTIPLE ] = 'MULTIPLE';
	    print STDERR
	      "Data type is now ($valtype) $type (", $tp[$type], ")\n";
	}

	# Remove quotes around (and doubled in) strings
	if ($type == SQ_STRING) {
	    $val =~ s/^'(.*)'$/$1/s;
	    $val =~ s/''/'/gs;
	}
	if ($type == DQ_STRING) {
	    $val =~ s/^"(.*)"$/$1/s;
	    $val =~ s/""/"/gs;
	}

	# Remove embedded newlines from strings (Anders' strange namelist
	# samples from Pencil Code with dust density)
	if (($type == SQ_STRING) || ($type == DQ_STRING)) {
	    $val =~ s/\n//g;
	}

	push @values, ($val) x $mul;
	$text =~ s/.*\n// if ($rest eq '!'); # comment
	print STDERR "<<", ($mul||'1'), "x>><<$val>> <<",
	  printable_substring($text), ">>\n" if ($debug);
    }

    $$txtptr = $text;		# return remaining unparsed string
    $$typeptr = $type;		# return type
    @values;
}

# ---------------------------------------------------------------------- #

sub infer_data_type {
# Determine the F90 data type of first item in string, skipping multiplier
# if present
    my $text = shift;
    $text =~ s/^\s*[0-9]+\*//;	# ignore multiplier for finding data type
    if      ($text =~ /^\s*'/)                     { SQ_STRING;
    } elsif ($text =~ /^\s*"/)                     { DQ_STRING;
    } elsif ($text =~ /^\s*\(\s*$numeric_e\s*,/)   { COMPLEX;
    } elsif ($text =~ /^\s*\(/)                    { DCOMPLEX;
    } elsif ($text =~ /^\s*(T|F|.(true|false).)/i) { LOGICAL;
    } elsif ($text =~ /^\s*[+-]?[0-9]+(\s|,|!|$)/)   { INTEGER;
    } elsif ($text =~ /^\s*$float(\s|,|!|$)/)        { FLOAT;
    } elsif ($text =~ /^\s*$numeric_e(\s|,|!|$)/)    { SINGLE;
    } elsif ($text =~ /^\s*$numeric_d(\s|,|!|$)/)    { DOUBLE;
    } else                                         { UNKNOWN;
    }
}

# ---------------------------------------------------------------------- #

sub show_error {
# Print error message and beginning of string with marker of current
# position [Slightly ridiculous, as the marker will always point to
# beginning of line]
    my $errmsg = shift;
    my $prefix = shift;
    my $text   = shift;
    my $die    = (shift || 0);

    chomp($errmsg);
    # print "$errmsg:\n";
    print "\e[01m$errmsg:\e[00m\n";

    # Escape newlines and only print 75 chars:
    my $subtext = $prefix . $text;
    $subtext =~ s/\n/\\n/g;
    $subtext = substr($subtext,0,75) . "\n";

    # Splice in marker line
    my $marker = " " x length($prefix) . "^------  HERE\n";
    $subtext =~ s/\n/\n$marker/;

    # Now die
    if ($die) {
	croak "$subtext";	# die
    } else {
	carp "$subtext";	# warn
    }
}

# ---------------------------------------------------------------------- #

sub printable_substring {
# Extract substring and quote newlines for diagnostic printing
    my $string = shift;
    my $length = shift || 40;

    $string =~ s/\n/\\n/g;
    my $oldlen = length($string);
    $string = substr($string,0,$length);
    substr($string,-3,3) = '...' if ($length<$oldlen);
    $string;
}

# ---------------------------------------------------------------------- #

sub assign_slot_val {
# Assignment of value to slot for output (format-dependent)
    my $var    = shift;
    my @vals   = @{shift()};
    my $format = shift;
    my $type   = shift;

    my $assmnt = "$var";

    if ($format eq 'f90') {
	$assmnt .= "=";
    } elsif ($format eq 'idl') {
	$assmnt .= ": ";	# structure syntax
    } else {
	die "assign_slot_val: Unknown format <$format>\n";
    }

    encapsulate_values(\@vals,$format,$type); # preprocess values
    if (@vals > 1) {
	$assmnt .= array_bracket(join(",", @vals), $format);
    } else {
	$assmnt .= $vals[0];
    }

    $assmnt;
}

# ---------------------------------------------------------------------- #

sub encapsulate_values {
# Format-specific preprocessing of data values, e.g. quoting strings,
# mapping logicals to integers for IDL, etc.
    my $valref = shift;
    my $format = shift;
    my $type   = shift;
    my @vals = @$valref;

    ## Actions for all formats
    if ($type==COMPLEX or $type==DCOMPLEX) {
	@vals = map { "${cmplx_pref}$_${cmplx_suff}" } @vals;
    }

    ## Actions specific for some formats
    if ($format eq 'f90') {
	#
	#  F90 output format:
	#  - quote strings
	#
	if      ($type==SQ_STRING or $type==DQ_STRING) {
	    @vals = map { quote_string_f90($_) } @vals;
	}
    } elsif ($format eq 'idl') {
	#
	#  IDL output format:
	#  - convert logicals to integers
	#  - quote strings
	#
	if      ($type==LOGICAL) {
	    @vals = map { encaps_logical_idl($_) } @vals;
	} elsif ($type==SQ_STRING or $type==DQ_STRING) {
	    @vals = map { quote_string_f90($_) } @vals;
	}
    } else {
	#
	#  Invalid format
	#
	die "encapsulate_values: Unknown format <$format>\n";
    }

    @$valref = @vals;
}

# ---------------------------------------------------------------------- #

sub array_bracket {
# Add format-specific array delimiters around string
    my $string = shift;
    my $format = shift;

    if     ($format eq 'f90') {
	# No delimiters
    } elsif ($format eq 'idl') {
	$string = "[$string]";
    } else {
	#
	#  Invalid format
	#
	die "array_bracket: Unknown format <$format>\n";
    }

    $string;
}

# ---------------------------------------------------------------------- #

sub encaps_logical_idl {
# Convert logical string to integer for IDL
    my $val = shift;

    $val =~ s/(\.false\.|F)/0L/i;
    $val =~ s/(\.true\.|T)/-1L/i;

    $val;
}

# ---------------------------------------------------------------------- #

sub quote_string_f90 {
# Enclose string by quotation marks, doubling any existing quotation marks
# for Fortran and IDL
    my $val = shift;

    $val =~ s/'/''/g;
    quote_string($val);
}

# ---------------------------------------------------------------------- #

sub quote_string {
# Enclose string by quotation marks
    "'$_[0]'";
}

# ---------------------------------------------------------------------- #

## Done.

1;

# End of file
