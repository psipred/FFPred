#!/usr/bin/perl -w

package SeqFeat;

use strict;
use base 'Pred';
use Data::Dumper;

sub new
{
    my ($class, $aa, $id, $md5, $cfg, $web_control) = @_;
    my $name = 'SF';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;

    $self->{'exe'}     = "$cfg->{'SEQFEAT'}/features";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.seqfeat";
    $self->{'cmd'}     = "$self->{'exe'}" .
                         " -i $self->{'path'}/$self->{'md5'}.fsa >" .
                         " $self->{'outfile'}";
    $self->{'web_control'} = $web_control;
    $self->{$self->name()} = {}; # This weird organisation comes from legacy code.
    $self->{'AA'}          = {};
    return $self;
}


sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'len'}       = log(1+$self->{'results'}{'len'}) / log(2001);
    $self->{'results'}{'mwt'}       = log(1+$self->{'results'}{'mwt'}) / log(300000);
    $self->{'results'}{'vol'}       = log(1+$self->{'results'}{'vol'}) / log(300000);
    $self->{'results'}{'surf'}      = log(1+$self->{'results'}{'surf'}) / log(400000);
    $self->{'results'}{'ave_hydro'} = log(5.5+$self->{'results'}{'ave_hydro'}) / log(5.5+70);
    $self->{'results'}{'charge'}    = (300+$self->{'results'}{'charge'}) / (300+200);
    $self->{'results'}{'mol_ext'}   = log(1+$self->{'results'}{'mol_ext'}) / log(500000);
    $self->{'results'}{'iso_pt'}    = log(1+$self->{'results'}{'iso_pt'}) / log(15);
    $self->{'results'}{'aliphatic'} = log(1+$self->{'results'}{'aliphatic'}) / log(200);
    $self->{'results'}{'num_atoms'} = log(1+$self->{'results'}{'num_atoms'}) / log(36000);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    my @AminoAcids   = ('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V');
    my @Atoms = ('atomC', 'atomH', 'atomO', 'atomN', 'atomS');

    open(SEQFEAT, "<", $self->outfile()) or die "Cannot open seqfeat file... $!\n";
    my $seqfeat_out = <SEQFEAT>;
    close(SEQFEAT);

    $seqfeat_out =~ s/^\s+//;
    my @seqfeat = split(/\s+/, $seqfeat_out);

    if($self->{'web_control'}){
      die "Bad output in seqfeat file - has the output format changed ?\n" unless (($seqfeat[2] == $self->len()));
    }
    else{
      die "Bad output in seqfeat file - has the output format changed ?\n" unless (($seqfeat[0] eq $self->md5()) && ($seqfeat[2] == $self->len()));
    }
    # See "$cfg->{'SEQFEAT'}/Sequence.cpp", Sequence::print().
    $CFG->{'length'}                       = $RES->{'len'}       = $seqfeat[2];
    $CFG->{'molecular weight'}             = $RES->{'mwt'}       = $seqfeat[3];
    $CFG->{'volume'}                       = $RES->{'vol'}       = $seqfeat[4];
    $CFG->{'surface area'}                 = $RES->{'surf'}      = $seqfeat[5];
    $CFG->{'hydrophobicity'}               = $RES->{'ave_hydro'} = $seqfeat[7];
    $CFG->{'charge'}                       = $RES->{'charge'}    = $seqfeat[8];
    $CFG->{'molar extinction coefficient'} = $RES->{'mol_ext'}   = $seqfeat[9];
    $CFG->{'isoelectric point'}            = $RES->{'iso_pt'}    = $seqfeat[10];
    $CFG->{'aliphatic index'}              = $RES->{'aliphatic'} = $seqfeat[11];
    $CFG->{'fraction positive residues'}   = $RES->{'ppos'}      = $seqfeat[14];
    $CFG->{'fraction negative residues'}   = $RES->{'pneg'}      = $seqfeat[15];
    $CFG->{'number of atoms'}              = $RES->{'num_atoms'} = $seqfeat[72];

    my $x = 0;
    for (my $i = 17; $i < 57; $i += 2, $x++)
    {
        $self->{'AA'}{$AminoAcids[$x]} = $RES->{$AminoAcids[$x]} = $seqfeat[$i];
    }

    $x = 0;
    for (my $i = 63; $i < 73; $i += 2, $x++)
    {
        $CFG->{$Atoms[$x]} = $RES->{$Atoms[$x]} = $seqfeat[$i];
    }

    $self->{'results'} = $RES;
}



1;
