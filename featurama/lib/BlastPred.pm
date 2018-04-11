#!/usr/bin/perl -w

package BlastPred;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg, $name, $mask, $iter, $Evalue) = @_;

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    bless $self, $class;
    $self->{'mask'}     = $mask;
    $self->{'iter'}     = $iter;
    $self->{'type'}     = $self->return_type($self->{'mask'}, $self->{'iter'});
    $self->{'db'}       = $self->{'mask'} eq 'masked' ? $cfg->{'MASK_DB'} : $cfg->{'DB'};
    $self->{'psiblast'} = new PsiBlast($self->{'iter'}, $self->{'mask'}, $cfg->{'BLAST'}, "$self->{'path'}/$self->{'md5'}.fsa", $self->{'db'}, $Evalue, $Evalue);
    $self->{'makeMat'}  = new MakeMat($cfg->{'MAKEMAT'}, "$self->{'path'}/$self->{'md5'}.$self->{'type'}", "$self->{'md5'}.$self->{'type'}.chk", "$self->{'md5'}.fsa");
    $self->{'seq2mtx'}  = new Seq2Mtx($cfg->{'SEQ2MTX'}, $cfg->{'PFILT'}, "$self->{'path'}/$self->{'md5'}.$self->{'type'}", "$self->{'path'}/$self->{'md5'}.fsa");
    return $self;
}

sub checkPsiblast
{
    my ($self) = @_;

    return $self->{'psiblast'}->check();
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

    $self->SUPER::run($error);
}

sub runPsiblast
{
    my ($self) = @_;

    $self->{'psiblast'}->run();
    $self->err($self->err() + $self->{'psiblast'}->err());
}



1;
