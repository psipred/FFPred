#!/usr/bin/perl -w

package PEST;

use strict;
use base 'Pred';

sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'PE';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;

    $self->{'exe'}     = "$cfg->{'PEST'}/epestfind";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.pest";
    $self->{'cmd'}     = "$self->{'exe'}" .
                         " -sequence $self->{'path'}/$self->{'md5'}.fsa" .
                         " -sformat1 fastA" .
                         " -window 10" .
                         " -order 2" .
                         " -outfile $self->{'outfile'}" .
                         " -graph none";

    $self->{'numSegments'} = 8;
    $self->{'segments'}    = $self->createSegments($self->{'numSegments'});
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    return $self;
}


sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_pest_regions'} = log(1+$self->{'results'}{'num_pest_regions'}) / log(20);
    $self->{'results'}{'num_pest_res'} = log(1+$self->{'results'}{'num_pest_res'}) / log(1000);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    my $min_pest_length = 1; # Minimum number of consecutive PEST residues forming a PEST region.

    $RES->{'num_pest_regions'} = 0;
    $RES->{'num_pest_res'}     = 0;
    $RES->{'pest_nterm'}       = 0;
    $RES->{'pest_cterm'}       = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
	$RES->{"pest_S$i"} = 0;
    }

    open(PEST, "<", $self->outfile()) or die "Cannot open PEST file... $!\n";
    while (defined(my $line = <PEST>))
    {
        if ($line =~ /Potential\s+PEST\s+motif\s+with\s+(\d+)\s+amino\s+acids\s+between\s+positions?\s+(\d+)\s+and\s+\d+/)
        {
            my ($length_pest, $start_pest, $stop_pest) = ($1, ($2 + 1), ($2 + $1)); # Peculiar output of epestfind.

            $RES->{'num_pest_res'} += $length_pest;

            foreach my $position ($start_pest..$stop_pest)
            {
                $self->addResidue($RES, $position, 'pest');
            }

            $self->addThresholdedRegion($RES, $CFG, 'pest', [ $min_pest_length ], $start_pest, $length_pest);
        }
    }
    close(PEST);

    $self->{'results'} = $RES;
}


1;
