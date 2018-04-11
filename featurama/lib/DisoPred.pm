#!/usr/bin/perl -w

package DisoPred;

use strict;
use PsiBlast;
use MakeMat;
use Seq2Mtx;
use base 'BlastPred';


sub new
{
    my ($class, $aa, $id, $md5, $cfg, $web_control) = @_;
    my $name = 'DI';
    my $mask = 'unmasked';
    my $iter = 3;
    my $evalue = 1e-3;
    my $hvalue = 1e-3;

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg, $name, $mask, $iter, $evalue, $hvalue);

    $self->{'exe'}     = "$cfg->{'DISOPRED'}/bin/disopred";
    $self->{'data'}    = "$cfg->{'DISOPRED'}/data";
    if($web_control)
    {
      $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.diso";
    }
    else
    {
      $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.$self->{'type'}.diso2"; # Manually synchronise this with cmd below! # Edited in the web server version.
    }
    $self->{'cmd'}     = "$self->{'exe'}" .
                         " $self->{'path'}/$self->{'md5'}.$self->{'type'}" .
                         " $self->{'path'}/$self->{'md5'}.$self->{'type'}.mtx" .
                         " $self->{'data'}/" .
                         " 2";
    $self->{'web_control'} = $web_control;
    $self->{'numSegments'} = 8;
    $self->{'segments'}    = $self->createSegments($self->{'numSegments'});
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    return $self;
}


sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_diso_res'} = log(1+$self->{'results'}{'num_diso_res'}) / log(2001);

    $self->{'results'}{'num_diso_5+'}   = log(1+$self->{'results'}{'num_diso_5+'}) / log(40);
    $self->{'results'}{'num_diso_50+'}  = log(1+$self->{'results'}{'num_diso_50+'}) / log(10);
    $self->{'results'}{'num_diso_100+'} = log(1+$self->{'results'}{'num_diso_100+'}) / log(6);
    $self->{'results'}{'num_diso_150+'} = log(1+$self->{'results'}{'num_diso_150+'}) / log(4);
    $self->{'results'}{'num_diso_200+'} = log(1+$self->{'results'}{'num_diso_200+'}) / log(4);
    $self->{'results'}{'num_diso_300+'} = log(1+$self->{'results'}{'num_diso_300+'}) / log(4);
    $self->{'results'}{'num_diso_500+'} = log(1+$self->{'results'}{'num_diso_500+'}) / log(3);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    my $diso_lengths = [ 5, 50, 100, 150, 200, 300, 500 ]; # Threshold values for disordered regions' lengths.

    foreach my $threshold (@$diso_lengths)
    {
        $RES->{"num_diso_${threshold}+"} = 0;
    }

    $RES->{'num_diso_res'} = 0;
    $RES->{'diso_nterm'}   = 0;
    $RES->{'diso_cterm'}   = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"diso_S$i"} = 0;
    }

    my $diso_out = {};
    open(DISO, "<", $self->outfile()) or die "Cannot open DISOPRED .diso2 file... $!\n"; # Edited in the web server version.
    while (defined(my $line = <DISO>))
    {
        $diso_out->{$1} = $2 if ($line =~ /^\s*(\d+)\s+\w\s+(\*|\.)\s/);
    }
    close(DISO);

    die "Bad output in DISOPRED .diso2 file - has the output format changed ?\n" unless (scalar(keys %$diso_out) == $self->len()); # Edited in the web server version.

    my $start_diso = 0;
    my $consecutive_diso_res = 0;
    foreach my $position (sort {$a <=> $b} keys %$diso_out)
    {
        if ($diso_out->{$position} eq '*')
        {
            $start_diso = $position unless ($start_diso);
            $RES->{'num_diso_res'}++;
            $self->addResidue($RES, $position, 'diso');
            $consecutive_diso_res++;
        }
        else
        {
            if ($consecutive_diso_res)
            {
                $self->addThresholdedRegion($RES, $CFG, 'diso', $diso_lengths, $start_diso, $consecutive_diso_res);
                $start_diso = $consecutive_diso_res = 0;
            }
        }
    }

    $self->addThresholdedRegion($RES, $CFG, 'diso', $diso_lengths, $start_diso, $consecutive_diso_res) if ($consecutive_diso_res); # Add last region if it includes the last residue.

    $self->{'results'} = $RES;
}


1;
