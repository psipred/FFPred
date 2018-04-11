#!/usr/bin/perl -w

package LowComplex;

use strict;
use base 'Pred';

sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'LC';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;

    $self->{'exe'} = "$cfg->{'LOWC'}/pfilt";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.lowc";
    $self->{'cmd'} = "$self->{'exe'} -tc" .
                     " $self->{'path'}/$self->{'md5'}.fsa >" .
                     " $self->{'outfile'}";

    $self->{'numSegments'} = 8;
    $self->{'segments'}    = $self->createSegments($self->{'numSegments'});
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    return $self;
}

sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_lowc_regions'} = log(1+$self->{'results'}{'num_lowc_regions'}) / log(30);
    $self->{'results'}{'num_lowc_res'} = log(1+$self->{'results'}{'num_lowc_res'}) / log(1000);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    my $min_lowc_length = 5; # Minimum number of consecutive low complexity residues forming a low complexity region.

    $RES->{'num_lowc_regions'} = 0;
    $RES->{'num_lowc_res'}     = 0;
    $RES->{'lowc_nterm'}       = 0;
    $RES->{'lowc_cterm'}       = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"lowc_S$i"} = 0;
    }

    my $lowc_seq = "";
    open(LOWC, "<", $self->outfile()) or die "Cannot open LOWC file... $!\n";
    while (defined(my $line = <LOWC>))
    {
        $lowc_seq .= $line unless ($line =~ /^>/);
    }
    close(LOWC);

    $lowc_seq =~ s/\s//g;

    die "Bad output in LOWC file - has the output format changed ?\n" unless (length($lowc_seq) == $self->len());

    my $position = 0;
    my $start_lowc = 0;
    my $consecutive_lowc_res = 0;
    foreach my $residue (split("", $lowc_seq))
    {
        $position++;
        if ($residue eq 'X')
        {
            $start_lowc = $position unless ($start_lowc);
            $RES->{'num_lowc_res'}++;
            $self->addResidue($RES, $position, 'lowc');
            $consecutive_lowc_res++;
        }
        else
        {
            if ($consecutive_lowc_res)
            {
                $self->addThresholdedRegion($RES, $CFG, 'lowc', [ $min_lowc_length ], $start_lowc, $consecutive_lowc_res);
                $start_lowc = $consecutive_lowc_res = 0;
            }
        }
    }

    $self->addThresholdedRegion($RES, $CFG, 'lowc', [ $min_lowc_length ], $start_lowc, $consecutive_lowc_res) if ($consecutive_lowc_res); # Add last region if it includes the last residue.

    $self->{'results'} = $RES;
}


1;
