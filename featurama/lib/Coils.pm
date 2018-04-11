#!/usr/bin/perl -w

package Coils;

use strict;
use base 'Pred';


sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'CC';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;

    $self->{'exe'}     = "$cfg->{'COILS'}/ncoils";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.coils";
    $self->{'cmd'}     = "$self->{'exe'}" .
                         " -f -w < $self->{'path'}/$self->{'md5'}.fsa >" .
                         " $self->{'outfile'}";

    $self->{'numSegments'} = 8;
    $self->{'segments'}    = $self->createSegments($self->{'numSegments'});
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    $ENV{'COILSDIR'} = $cfg->{'COILS'}; # Set up environment for COILS.

    return $self;
}


sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_coil_regions'} = log(1+$self->{'results'}{'num_coil_regions'}) / log(30);
    $self->{'results'}{'num_coil_res'} = log(1+$self->{'results'}{'num_coil_res'}) / log(1000);
}


sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    my $min_coil_length = 5; # Minimum number of consecutive coil residues forming a coiled coil.

    $RES->{'num_coil_regions'} = 0;
    $RES->{'num_coil_res'}     = 0;
    $RES->{'coil_nterm'}       = 0;
    $RES->{'coil_cterm'}       = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"coil_S$i"} = 0;
    }

    my $coils_seq = "";
    open(COILS, "<", $self->outfile()) or die "Cannot open COILS file... $!\n";
    while (defined(my $line = <COILS>))
    {
        $coils_seq .= $line unless ($line =~ /^>/);
    }
    close(COILS);

    $coils_seq =~ s/\s//g;
    die "Bad output in COILS file - has the output format changed ?\n" unless (length($coils_seq) == $self->len());

    my $position = 0;
    my $start_coil = 0;
    my $consecutive_coil_res = 0;
    foreach my $residue (split("", $coils_seq))
    {
        $position++;
        if ($residue eq 'x')
        {
            $start_coil = $position unless ($start_coil);
            $RES->{'num_coil_res'}++;
            $self->addResidue($RES, $position, 'coil');
            $consecutive_coil_res++;
        }
        else
        {
            if ($consecutive_coil_res)
            {
                $self->addThresholdedRegion($RES, $CFG, 'coil', [ $min_coil_length ], $start_coil, $consecutive_coil_res);
                $start_coil = $consecutive_coil_res = 0;
            }
        }
    }

    $self->addThresholdedRegion($RES, $CFG, 'coil', [ $min_coil_length ], $start_coil, $consecutive_coil_res) if ($consecutive_coil_res); # Add last region if it includes the last residue.

    $self->{'results'} = $RES;
}




1;
