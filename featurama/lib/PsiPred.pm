#!/usr/bin/perl -w

package PsiPred;

use strict;
use PsiBlast;
use MakeMat;
use Seq2Mtx;
use base 'BlastPred';


sub new
{
    my ($class, $aa, $id, $md5, $cfg, $web_control) = @_;
    my $name = 'SS';
    my $mask = 'masked';
    my $iter = 3;
    my $evalue = 1e-3;
    my $hvalue = 1e-3;

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg, $name, $mask, $iter, $evalue, $hvalue);

    $self->{'exe_psipred'} = "$cfg->{'PSIPRED'}/bin/psipred";
    $self->{'exe_psipass'} = "$cfg->{'PSIPRED'}/bin/psipass2";
    $self->{'data'}        = "$cfg->{'PSIPRED'}/data";
    if($web_control)
    {
      $self->{'outfile'}     = "$self->{'path'}/$self->{'md5'}.ss2";
    }
    else{
      $self->{'outfile'}     = "$self->{'path'}/$self->{'md5'}.$self->{'type'}.ss2"; # Manually synchronise this with cmd below!
    }
    $self->{'cmd'}         = "$self->{'exe_psipred'}" .
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.mtx" .
                             " $self->{'data'}/weights.dat" .
                             " $self->{'data'}/weights.dat2" .
                             " $self->{'data'}/weights.dat3 >" .
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.ss";
    $self->{'cmd_psipass'} = "$self->{'exe_psipass'}" .
                             " $self->{'data'}/weights_p2.dat" .
                             " 1 1.0 1.0" .
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.ss2" .
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.ss >" .
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.horiz";
    $self->{'web_control'} = $web_control;
    $self->{'numSegments'} = 8;
    $self->{'segments'}    = $self->createSegments($self->{'numSegments'});
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    return $self;
}


sub cmd_psipass
{
    my ($self, $cmd) = @_;

    $self->{'cmd_psipass'} = $cmd if (defined($cmd));
    return $self->{'cmd_psipass'};
}

sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_psipredH_res'} = log(1+$self->{'results'}{'num_psipredH_res'}) / log(1700);
    $self->{'results'}{'num_psipredE_res'} = log(1+$self->{'results'}{'num_psipredE_res'}) / log(1200);
    $self->{'results'}{'num_psipredC_res'} = log(1+$self->{'results'}{'num_psipredC_res'}) / log(2001);

    $self->{'results'}{'num_psipredH_1+'}   = log(1+$self->{'results'}{'num_psipredH_1+'}) / log(80);
    $self->{'results'}{'num_psipredH_10+'}  = log(1+$self->{'results'}{'num_psipredH_10+'}) / log(70);
    $self->{'results'}{'num_psipredH_15+'}  = log(1+$self->{'results'}{'num_psipredH_15+'}) / log(50);
    $self->{'results'}{'num_psipredH_20+'}  = log(1+$self->{'results'}{'num_psipredH_20+'}) / log(40);
    $self->{'results'}{'num_psipredH_30+'}  = log(1+$self->{'results'}{'num_psipredH_30+'}) / log(25);
    $self->{'results'}{'num_psipredH_50+'}  = log(1+$self->{'results'}{'num_psipredH_50+'}) / log(15);
    $self->{'results'}{'num_psipredH_70+'}  = log(1+$self->{'results'}{'num_psipredH_70+'}) / log(6);
    $self->{'results'}{'num_psipredH_100+'} = log(1+$self->{'results'}{'num_psipredH_100+'}) / log(5);

    $self->{'results'}{'num_psipredE_1+'}  = log(1+$self->{'results'}{'num_psipredE_1+'}) / log(200);
    $self->{'results'}{'num_psipredE_10+'} = log(1+$self->{'results'}{'num_psipredE_10+'}) / log(40);
    $self->{'results'}{'num_psipredE_15+'} = log(1+$self->{'results'}{'num_psipredE_15+'}) / log(15);
    $self->{'results'}{'num_psipredE_20+'} = log(1+$self->{'results'}{'num_psipredE_20+'}) / log(15);
    $self->{'results'}{'num_psipredE_25+'} = log(1+$self->{'results'}{'num_psipredE_25+'}) / log(3);
    $self->{'results'}{'num_psipredE_30+'} = log(1+$self->{'results'}{'num_psipredE_30+'}) / log(2);
    $self->{'results'}{'num_psipredE_40+'} = log(1+$self->{'results'}{'num_psipredE_40+'}) / log(2);

    $self->{'results'}{'num_psipredC_regions'} = log(1+$self->{'results'}{'num_psipredC_res'}) / log(200);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    my $psipred_lengths = {
                            'H' => [ 1, 10, 15, 20, 30, 50, 70, 100 ],   # Threshold values for psipred helical regions' lengths.
                            'E' => [ 1, 10, 15, 20, 25, 30, 40 ]         # Threshold values for psipred sheet regions' lengths.
                          };

    foreach my $type_HEC ('H', 'E')
    {
        foreach my $threshold (@{$psipred_lengths->{$type_HEC}})
        {
            $RES->{"num_psipred${type_HEC}_${threshold}+"} = 0;
        }
    }
    $RES->{'num_psipredC_regions'} = 0;

    $RES->{'num_psipredH_res'} = 0;
    $RES->{'num_psipredE_res'} = 0;
    $RES->{'num_psipredC_res'} = 0;
    $RES->{'psipredH_nterm'} = 0;
    $RES->{'psipredE_nterm'} = 0;
    $RES->{'psipredC_nterm'} = 0;
    $RES->{'psipredH_cterm'} = 0;
    $RES->{'psipredE_cterm'} = 0;
    $RES->{'psipredC_cterm'} = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"psipredH_S$i"} = 0;
        $RES->{"psipredE_S$i"} = 0;
        $RES->{"psipredC_S$i"} = 0;
    }

    my $psipred_out = {};
    open(PSIPRED, "<", $self->outfile()) or die "Cannot open PSIPRED .ss2 file... $!\n";
    while (defined(my $line = <PSIPRED>))
    {
        $psipred_out->{$1} = $2 if ($line =~ /^\s*(\d+)\s+\w\s+(H|E|C)\s/);
    }
    close(PSIPRED);

    die "Bad output in PSIPRED .ss2 file - has the output format changed ?\n" unless (scalar(keys %$psipred_out) == $self->len());

    my $type_HEC = "";
    my $start_HEC = 0;
    my $consecutive_HEC_res = 0;
    foreach my $position (sort {$a <=> $b} keys %$psipred_out)
    {
        $type_HEC = $psipred_out->{$position} unless ($type_HEC);

        unless ($psipred_out->{$position} eq $type_HEC)
        {
            if ($type_HEC eq 'C')
            {
                $RES->{'num_psipredC_regions'}++ if ($consecutive_HEC_res > 1); # This is done here because we don't want C regions to appear in $CFG.
            }
            else
            {
                $self->addThresholdedRegion($RES, $CFG, "psipred${type_HEC}", $psipred_lengths->{$type_HEC}, $start_HEC, $consecutive_HEC_res);
            }

            $type_HEC = $psipred_out->{$position};
            $start_HEC = $consecutive_HEC_res = 0;
        }

        $start_HEC = $position unless ($start_HEC);
        $RES->{"num_psipred${type_HEC}_res"}++;
        $self->addResidue($RES, $position, "psipred${type_HEC}");
        $consecutive_HEC_res++;
    }

    # Add last region.
    if ($type_HEC eq 'C')
    {
        $RES->{'num_psipredC_regions'}++; # This is done here because we do not want C regions to appear in $CFG.
    }
    else
    {
        $self->addThresholdedRegion($RES, $CFG, "psipred${type_HEC}", $psipred_lengths->{$type_HEC}, $start_HEC, $consecutive_HEC_res);
    }

    $self->{'results'} = $RES;
}

sub run
{
    my ($self, $error) = @_;
    $error = 0 unless ($error);

    if (-s "$self->{path}/$self->{md5}.$self->{'type'}.chk")
    {
        $self->{'makeMat'}->init();
        $self->{'makeMat'}->run();
        $error += $self->{'makeMat'}->err();
    }
    else
    {
        $self->{'seq2mtx'}->run();
        $error += $self->{'seq2mtx'}->err();
    }

    print STDERR $self->cmd() . "\n";
    $error += system($self->cmd());

    if (-s "$self->{path}/$self->{md5}.$self->{'type'}.ss")
    {
        print STDERR $self->cmd_psipass(), "\n";
        $error += system($self->cmd_psipass());
    }
    else
    {
        print STDERR "Could not find .ss file !\n";
        $error = 1;
    }

    $self->err($error);
}



1;
