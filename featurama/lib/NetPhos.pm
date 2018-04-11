#!/usr/bin/perl -w

package NetPhos;

use strict;
use base 'Pred';

sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'PH';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;

    $self->{'exe'}     = "$cfg->{'NETPHOS'}/ape";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.netphos";
    $self->{'cmd'}     = "$self->{'exe'} -m netphos" .
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

    $self->{'results'}{'num_phosS_res'} = log(1+$self->{'results'}{'num_phosS_res'}) / log(1000);
    $self->{'results'}{'num_phosT_res'} = log(1+$self->{'results'}{'num_phosT_res'}) / log(1000);
    $self->{'results'}{'num_phosY_res'} = log(1+$self->{'results'}{'num_phosY_res'}) / log(100);

    $self->{'results'}{'num_ATM_res'}     = log(1+$self->{'results'}{'num_ATM_res'}) / log(30);
    $self->{'results'}{'num_CaM-II_res'}  = log(1+$self->{'results'}{'num_CaM-II_res'}) / log(5);
    $self->{'results'}{'num_CKI_res'}     = log(1+$self->{'results'}{'num_CKI_res'}) / log(350);
    $self->{'results'}{'num_CKII_res'}    = log(1+$self->{'results'}{'num_CKII_res'}) / log(200);
    $self->{'results'}{'num_DNAPK_res'}   = log(1+$self->{'results'}{'num_DNAPK_res'}) / log(40);
    $self->{'results'}{'num_EGFR_res'}    = log(1+$self->{'results'}{'num_EGFR_res'}) / log(15);
    $self->{'results'}{'num_GSK3_res'}    = log(1+$self->{'results'}{'num_GSK3_res'}) / log(120);
    $self->{'results'}{'num_INSR_res'}    = log(1+$self->{'results'}{'num_INSR_res'}) / log(40);
    $self->{'results'}{'num_PKA_res'}     = log(1+$self->{'results'}{'num_PKA_res'}) / log(50);
    $self->{'results'}{'num_PKB_res'}     = log(1+$self->{'results'}{'num_PKB_res'}) / log(30);
    $self->{'results'}{'num_PKC_res'}     = log(1+$self->{'results'}{'num_PKC_res'}) / log(200);
    $self->{'results'}{'num_PKG_res'}     = log(1+$self->{'results'}{'num_PKG_res'}) / log(50);
    $self->{'results'}{'num_RSK_res'}     = log(1+$self->{'results'}{'num_RSK_res'}) / log(50);
    $self->{'results'}{'num_SRC_res'}     = log(1+$self->{'results'}{'num_SRC_res'}) / log(10);
    $self->{'results'}{'num_cdc2_res'}    = log(1+$self->{'results'}{'num_cdc2_res'}) / log(300);
    $self->{'results'}{'num_cdk5_res'}    = log(1+$self->{'results'}{'num_cdk5_res'}) / log(120);
    $self->{'results'}{'num_p38MAPK_res'} = log(1+$self->{'results'}{'num_p38MAPK_res'}) / log(120);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{'num_phosS_res'} = 0;
    $RES->{'num_phosT_res'} = 0;
    $RES->{'num_phosY_res'} = 0;
    $RES->{'phosS_nterm'}   = 0;
    $RES->{'phosT_nterm'}   = 0;
    $RES->{'phosY_nterm'}   = 0;
    $RES->{'phosS_cterm'}   = 0;
    $RES->{'phosT_cterm'}   = 0;
    $RES->{'phosY_cterm'}   = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"phosS_S$i"} = 0;
        $RES->{"phosT_S$i"} = 0;
        $RES->{"phosY_S$i"} = 0;
    }

    $RES->{'num_ATM_res'}     = 0;
    $RES->{'num_CaM-II_res'}  = 0;
    $RES->{'num_CKI_res'}     = 0;
    $RES->{'num_CKII_res'}    = 0;
    $RES->{'num_DNAPK_res'}   = 0;
    $RES->{'num_EGFR_res'}    = 0;
    $RES->{'num_GSK3_res'}    = 0;
    $RES->{'num_INSR_res'}    = 0;
    $RES->{'num_PKA_res'}     = 0;
    $RES->{'num_PKB_res'}     = 0;
    $RES->{'num_PKC_res'}     = 0;
    $RES->{'num_PKG_res'}     = 0;
    $RES->{'num_RSK_res'}     = 0;
    $RES->{'num_SRC_res'}     = 0;
    $RES->{'num_cdc2_res'}    = 0;
    $RES->{'num_cdk5_res'}    = 0;
    $RES->{'num_p38MAPK_res'} = 0;

    my $good_positions = {}; # Needed because the same residue may be positive for several kinases.

    open(PHOS, "<", $self->outfile()) or die "Cannot open netPhos file... $!\n";
    while (defined(my $line = <PHOS>))
    {
        if ($line =~ /^#\s*\S+\s+(\d+)\s+(S|T|Y)\s+\S+\s+(\S+)\s+(\S+).+YES/)
        {
            my ($position, $phosres, $score, $phostype) = ($1, $2, $3, $4);
            my $label = "phos${phosres}";

            unless (exists($good_positions->{$position}))
            {
                $RES->{"num_${label}_res"}++;
                $self->addResidue($RES, $position, $label);
                $good_positions->{$position} = 1;
            }

            $RES->{"num_${phostype}_res"}++ unless ($phostype =~ /unsp/);

            my $n = scalar @$CFG;
            $CFG->[$n]{'type'} = "${label}_${phostype}";
            $CFG->[$n]{'from'} = $position;
            $CFG->[$n]{'to'} = $position;
            $CFG->[$n]{'score'} = $score;
        }
    }
    close(PHOS);

    $self->{'results'} = $RES;
}



1;
