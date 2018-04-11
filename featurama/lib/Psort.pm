#!/usr/bin/perl -w

package Psort;

use strict;
use base 'Pred';

sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'PS';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{PATH}, $name);
    bless $self, $class;

    $self->{'exe'}     = "$cfg->{'PSORT'}/runWolfPsortSummary";
    $self->{'k'}       = 0; # Value of k in k-NN method used in run().
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.wolfpsort";
    $self->{'cmd'}     = "$self->{'exe'} animal" .
	                 " < $self->{'path'}/$self->{'md5'}.fsa" .
                         " > $self->{'outfile'}";

    $self->{$self->name()} = {}; # This weird organisation comes from legacy code.

    return $self;
}


sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    # See WoLF-PSORT's runWolfPsortSummary man page & http://wolfpsort.org/aboutWoLF_PSORT.html.en .
    my $type_list = {
                      'chlo' => 'chloroplast',             # GO:0009507, GO:0009543
                      'cysk' => 'cytoskeleton',            # GO:0005856, all children and all grandchildren
                      'cyto' => 'cytosol',                 # GO:0005829 (not documented in the man page due to a bug, as of v0.2).
                      'E.R.' => 'endoplasmic reticulum',   # GO:0005783
                      'extr' => 'extracellular',           # GO:0005576, GO:0005618
                      'golg' => 'Golgi apparatus',         # GO:0005794 and all children
                      'lyso' => 'lysosome',                # GO:0005764
                      'mito' => 'mitochondria',            # GO:0005739
                      'nucl' => 'nuclear',                 # GO:0005634
                      'pero' => 'peroxisome',              # GO:0005777, all children and all grandchildren
                      'plas' => 'plasma membrane',         # GO:0005886
                      'vacu' => 'vacuolar membrane'        # GO:0005774, all children and all grandchildren
                     };
    my $reg_exp = ''; # This will be e.g. 'chlo|cysk|cyto|E\.R\.|extr|golg|lyso|mito|nucl|pero|plas|vacu'.

    foreach my $type (keys %$type_list)
    {
        $RES->{$type} = $CFG->{$type_list->{$type}} = 0;
        if ($type eq 'E.R.') # Handled as exceptional case due to their weird choice of non-word characters in this label.
        {
            $reg_exp .= 'E\.R\.|';
        }
        else
        {
            $reg_exp .= "$type|";
        }
    }
    $reg_exp =~ s/\|$//;

    open(PSORT, "<", $self->outfile()) or die "Cannot open WoLFPSORT file... $!\n";
    my $k_line = <PSORT>;
    my $results_line = <PSORT>;
    close(PSORT);

    if ($k_line =~ /k used for kNN is:\s*(\d+)/)
    {
        $self->{'k'} = $1;
        foreach my $result ($results_line =~ /\b(?:$reg_exp)\s+\d+(?:\.\d+)?/g)
        {
            my ($type, $value) = split(/\s+/, $result);
            $value /= $self->{'k'};
            $RES->{$type} = $CFG->{$type_list->{$type}} = $value;
        }
    }

    $self->{'results'} = $RES;
}


# Old command to run PSORT II.
sub cmd_PSORTII
{
    my ($self) = @_;

    my $cmd = "$self->{exe}/psort " .
	      "$self->{path}/$self->{md5}.fsa > " .
              "$self->{path}/$self->{md5}.psort";

    return $cmd;
}

# Old legacy parser for PSORT II output.
sub parse_PSORTII
{
    my ($self) = @_;
    my $RES = {}; # This will contain 11 features.

    $RES->{'cytoplasmic'}   = $self->{PS}->{'cytoplasmic'}->{val}   = 0;
    $RES->{'nuclear'}       = $self->{PS}->{'nuclear'}->{val}       = 0;
    $RES->{'mitochondrial'} = $self->{PS}->{'mitochondrial'}->{val} = 0;
    $RES->{'peroxisomal'}   = $self->{PS}->{'peroxisomal'}->{val}   = 0;
    $RES->{'extracellular'} = $self->{PS}->{'extracellular'}->{val} = 0;
    $RES->{'golgi'}         = $self->{PS}->{'golgi'}->{val}         = 0;
    $RES->{'plasma'}        = $self->{PS}->{'plasma'}->{val}        = 0;
    $RES->{'vacuolar'}      = $self->{PS}->{'vacuolar'}->{val}      = 0;
    $RES->{'cytoskeletal'}  = $self->{PS}->{'cytoskeletal'}->{val}  = 0;
    $RES->{'vesicles'}      = $self->{PS}->{'vesicles'}->{val}      = 0;
    $RES->{'endoplasmic'}   = $self->{PS}->{'endoplasmic'}->{val}   = 0;

    if (-s "$self->{path}/$self->{md5}.psort")
    {
        open(PSORT, "< $self->{path}/$self->{md5}.psort");

        while(<PSORT>)
        {
            chomp $_;

            if($_ =~ /^\s+(\d+\.\d+)\s+\%\:\s+(\w+)/)
            {
                my ($pcent, $type) = ($1, $2);
                $type = lc($type);
                $type =~ s/\,//g;
                $RES->{$type} = $pcent/100;
                $self->{PS}->{$type}->{val} = $pcent/100;
            }
        }

        close(PSORT);
    }

    $self->{results} = $RES;
}



1;
