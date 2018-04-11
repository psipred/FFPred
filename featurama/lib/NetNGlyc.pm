#!/usr/bin/perl -w

package NetNGlyc;

use strict;
use base 'Pred';


sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'NG';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;

    $self->{'exe'}     = "$cfg->{'NETNGLYC'}/netNglyc";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.netNglyc";
    $self->{'cmd'}     = "$self->{'exe'}" .
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

    $self->{'results'}{'num_nglyc_res'} = log(1+$self->{'results'}{'num_nglyc_res'}) / log(40);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{'num_nglyc_res'} = 0;
    $RES->{'nglyc_nterm'}   = 0;
    $RES->{'nglyc_cterm'}   = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"nglyc_S$i"} = 0;
    }

    if ($self->{'aa'} =~ /N/)
    {
        open(NGLYC, "<", $self->outfile()) or die "Cannot open netNglyc file... $!\n";
        while (defined(my $line = <NGLYC>))
        {
            if ($line =~ /^\s*\S+\s+(\d+)\s+\w+\s+(\S+).+\+/)
            {
                my ($position, $score) = ($1, $2);

                $RES->{'num_nglyc_res'}++;
                $self->addResidue($RES, $position, 'nglyc');

                my $n = scalar @$CFG;
                $CFG->[$n]{'type'} = 'nglyc';
                $CFG->[$n]{'from'} = $position;
                $CFG->[$n]{'to'} = $position;
                $CFG->[$n]{'score'} = $score;
            }
        }
        close(NGLYC);
    }

    $self->{'results'} = $RES;
}

sub run
{
    my ($self, $error) = @_;
    $error = 0 unless ($error);

    $self->SUPER::run($error) if ($self->{'aa'} =~ /N/);
}



1;
