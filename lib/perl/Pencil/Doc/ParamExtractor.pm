package Pencil::Doc::ParamExtractor;

use parent 'Pencil::Doc::Extractor';

use warnings;
use strict;
use Carp;
use vars qw($VERSION);

sub new {
    my $class = shift;
    my @argv  = @_;
    my $self = $class->SUPER::new(@_);
    $self->{SRCDIR} = $self->{ARGS}->{src_dir};
    $self->populate_global_docs();
    return $self;
}

sub split_nml {
# Returns an array containing the names of all the variables in the specified namelist. $file is the file to parse, while $nml_name ('init' or 'run') specifies which namelist to look in.
    my $self = shift;
    my $nml_name = shift;
    my $file = shift;
    
    my @var_names;
    my $nml = $self->get_pars_namelist($nml_name, $file);
    if ($nml) {
        my @tmp = split('/', $nml);
        my $vars = $tmp[-1];
        $vars =~ s/!.*$//g; #remove comments
        $vars =~ s/[&\s]+//gs; #remove line continuations
        @var_names = split(',', $vars);
    } else {
        @var_names = ();
    }
    
    return @var_names
}

sub get_named_variables {
# given a filename, return a hash containing all the variables that appear in init_pars or run_pars.
    my $self = shift;
    my $file = shift;

    my %named_vars;

    foreach my $var_name ($self->split_nml('init', $file)) {
        $named_vars{$var_name} = {
            is_start => 1,
        };
    }
    foreach my $var_name ($self->split_nml('run', $file)) {
        if (exists $named_vars{$var_name}) {
            $named_vars{$var_name}{is_run} = 1;
        } else {
            $named_vars{$var_name} = {
                is_run => 1,
            };
        }
    }

    return %named_vars;
}

sub parse {
    my $self = shift();
    my $file = shift();

    carp("Cannot read file <$file>\n") unless (-r $file);

    my %named_vars = $self->get_named_variables($file);
    (my $sfile = $file) =~ s{.*/}{}; # remove path
    
    unless (%named_vars) {
        return 0;
    };

    $self->{DOC}{$sfile} = \%named_vars;

    #Match docstrings to variables
    my @localdoc = $self->get_docs_from_file($file,
                                      $self->{MARKER},
                                      $self->{PREFIX},
                                      $self->{DEBUG},
                                      $self->{VERBOSE}
    );
    my @globaldoc = @{$self->{GLOBAL_DOCS}};
    foreach my $vardocref (@localdoc, @globaldoc) {
        next unless ($vardocref);
        my $var = $vardocref->{var};
        if (exists $named_vars{$var}) {
            foreach my $key (keys %{$vardocref}) {
                $named_vars{$var}{$key} = $vardocref->{$key};
            }
        }
    }
    
    my $count = scalar %named_vars;
    return $count;
}

sub populate_global_docs {
# handles docstrings in cdata.f90 and general.f90. Needs to be done
# separately because the variables defined in these files may appear
# in namelists of other files.
    my $self = shift;
    
    my $path = $self->{SRCDIR};
    
    # TODO: density_init_pars contains cs2top and cs2bot which are defined in the
    # EOS modules. If both eos_idealgas.f90 and eos_idealgas_vapor.f90 document
    # a parameter that appears in density_init_pars, which of them should be
    # prioritized? Note that it is not correct to add eos_idealgas.f90 below,
    # since that would result in the documentation strings from eos_idealgas.f90
    # being indiscriminately applied to all the EOS modules.
    
    my @gd;
    foreach my $file ("$path/cdata.f90", "$path/general.f90") {
        my @doc = $self->get_docs_from_file($file,
                                     $self->{MARKER},
                                     $self->{PREFIX},
                                     $self->{DEBUG},
                                     $self->{VERBOSE}
        );
        push @gd, @doc;
    }
    $self->{GLOBAL_DOCS} = \@gd;
}

# ---------------------------------------------------------------------- #
sub longtable {
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
        @files = sort { $self->smart_compare($a, $b) } @files;
    }

    my $text  = $self->header(\@files, $args{selfcontained}, $args{descr_width});

    foreach my $module (@files) {
        # Header line for each section of table
        my $sectionheader =
            "\\midrule\n"
          . "  \\multicolumn{2}{c}{Module \\file{$module}} \\\\\n"
          . "\\midrule\n";

        my $this_text = "";
        # Loop through variables
        my $file_docs_ref = $docref->{$module}; # (['var1', 'doc1'],
                                               #  ['var2', 'doc2'], ...)
        #sort keys alphabetically, case insensitive
        foreach my $var (sort { "\L$a" cmp "\L$b" } keys %{$file_docs_ref}) {
            my $doc = $file_docs_ref->{$var}->{doc} || '';
            my $is_start = $file_docs_ref->{$var}->{is_start} || 0;
            my $is_run = $file_docs_ref->{$var}->{is_run} || 0;
            my $default_value = $file_docs_ref->{$var}->{default_value} || '';

            next unless ($print_empty || (defined($doc) && $doc =~ /\S/));

            # Indent continued lines, so LaTeX code is easier to read:
            $doc =~ s{\n}{\n                  }g;

            my $nml_label;
            if ($is_start && $is_run) {
                $nml_label = 'init,run';
            } elsif ($is_start) {
                $nml_label = 'init';
            } elsif ($is_run) {
                $nml_label = 'run';
            } else {
                croak "$var was found in neither start nor run";
            }

            if ($default_value && ($default_value ne 'impossible')) {
                $this_text .= sprintf "  %-15s & %s \\\\\n & {[%s; default = %s]} \\\\\n", "\\var{$var}", $doc, $nml_label, "\\verb|$default_value|";
            } else {
                $this_text .= sprintf "  %-15s & %s \\\\\n & {[%s]} \\\\\n", "\\var{$var}", $doc, $nml_label;
            }
        }

        if ($this_text) {
            $text .= $sectionheader . $this_text
        }

    }

    $text .= $self->footer($args{selfcontained});
    return $text;
}

# ---------------------------------------------------------------------- #

sub get_pars_namelist {
# Given a source file, extract the portion of it that specifies the start/run pars
    my $self = shift;
    my $name = shift;
    my $file = shift;

    unless (open(NML_FILE, "< $file")) {
        carp "Cannot open $file for reading: $!\n";
        return ();
    }

    my $processing_namelist = 0;
    my $namelist = '';
    while(defined(my $line = <NML_FILE>)) {
        if ($line =~ /^\s*namelist \/[a-zA-Z_]*_?${name}_pars\//) {
           $processing_namelist = 1;
        }
        if ($processing_namelist) {
            #Assumes each line is terminated by \n; we don't need to manually insert it.
            $namelist = $namelist . $line;

            unless ($line =~ /&\s*$/) {
                $processing_namelist = 0;
                last;
            }
        }
    }

    return $namelist;
}

# ---------------------------------------------------------------------- #


1;
__END__
