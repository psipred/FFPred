#!/usr/bin/perl -w

package NetOGlyc;

use strict;
use base 'Pred';

sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'OG';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;

    $self->{'exe'}     = "$cfg->{'NETOGLYC'}/netOglyc";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.netOglyc";
    $self->{'cmd'}     = "$self->{'exe'}" .
	                 " $self->{path}/$self->{md5}.fsa >" .
                         " $self->{'outfile'}";

    $self->{'numSegments'} = 8;
    $self->{'segments'}    = $self->createSegments($self->{'numSegments'});
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    return $self;
}

sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_oglycS_res'} = log(1+$self->{'results'}{'num_oglycS_res'}) / log(150);
    $self->{'results'}{'num_oglycT_res'} = log(1+$self->{'results'}{'num_oglycT_res'}) / log(600);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{'num_oglycS_res'} = 0;
    $RES->{'num_oglycT_res'} = 0;
    $RES->{'oglycS_nterm'}   = 0;
    $RES->{'oglycT_nterm'}   = 0;
    $RES->{'oglycS_cterm'}   = 0;
    $RES->{'oglycT_cterm'}   = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"oglycS_S$i"} = 0;
        $RES->{"oglycT_S$i"} = 0;
    }

    open(OGLYC, "<", $self->outfile()) or die "Cannot open netOglyc file... $!\n";
    while (defined(my $line = <OGLYC>))
    {
        if ($line =~ /^\s*\S+\s+\w\s+(\d+)\s+(\S+)\s+\S+\s+(S|T)\s/)
        {
            my ($position, $score, $oglyctype) = ($1, $2, $3);
            my $label = "oglyc${oglyctype}";

            $RES->{"num_${label}_res"}++;
            $self->addResidue($RES, $position, $label);

            my $n = scalar @$CFG;
            $CFG->[$n]{'type'} = $label;
            $CFG->[$n]{'from'} = $position;
            $CFG->[$n]{'to'} = $position;
            $CFG->[$n]{'score'} = $score;
        }
    }
    close(OGLYC);

    $self->{'results'} = $RES;
}



1;
