#!/usr/bin/perl -w

package DisoPred;

use strict;
use PsiBlast;
use MakeMat;
use Seq2Mtx;
use base 'BlastPred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'DI';
    my $mask = 'unmasked';
    my $iter = 3;
    my $Evalue = 1e-3;

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg, $name, $mask, $iter, $Evalue);

    $self->{'exe'}  = "$cfg->{'DISOPRED'}/bin/disopred2";
    $self->{'data'} = "$cfg->{'DISOPRED'}/data";
    $self->{'cmd'}  = "$self->{'exe'}" .
                      " $self->{'path'}/$self->{'md5'}.$self->{'type'}" .
                      " $self->{'path'}/$self->{'md5'}.$self->{'type'}.mtx" .
                      " $self->{'data'}/" .
                      " 2";

    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    return $self;
}

sub normalise
{
    my ($self) = @_;

    $self->{results}->{num_diso} = log(1+$self->{results}->{num_diso})/log(190);

    $self->{results}->{diso}->{50}    = log(1+$self->{results}->{diso}->{50})/log(177);
    $self->{results}->{diso}->{100}   = log(1+$self->{results}->{diso}->{100})/log(46);
    $self->{results}->{diso}->{150}   = log(1+$self->{results}->{diso}->{150})/log(17);
    $self->{results}->{diso}->{200}   = log(1+$self->{results}->{diso}->{200})/log(4);
    $self->{results}->{diso}->{300}   = log(1+$self->{results}->{diso}->{300})/log(7);
    $self->{results}->{diso}->{500}   = log(1+$self->{results}->{diso}->{500})/log(4);
    $self->{results}->{diso}->{10000} = log(1+$self->{results}->{diso}->{10000})/log(4);

}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{'num_diso'}    = 0;
    $RES->{'diso_res'}    = 0;
    $RES->{'diso_nterm'}  = 0;
    $RES->{'diso_cterm'}  = 0;
    $RES->{'diso'}{50}    = 0;
    $RES->{'diso'}{100}   = 0;
    $RES->{'diso'}{150}   = 0;
    $RES->{'diso'}{200}   = 0;
    $RES->{'diso'}{300}   = 0;
    $RES->{'diso'}{500}   = 0;
    $RES->{'diso'}{10000} = 0;

    for (my $i = 1; $i < 9; $i++)
    {
        $RES->{"diso_S$i"} = 0;
    }

    open(DISO, "<", "$self->{path}/$self->{md5}.$self->{'type'}.diso");

    my $last = "";
    my ($from, $to) = (1, 0);

    my @F = sort {$a <=> $b} keys %{$RES->{'diso'}};

    while (<DISO>)
    {
        $_ =~ s/^\s*//;
        chomp $_;

        next if $_ =~ /^\s*$/;

        my ($idx,$aa,$diso) = split(/\s+/,$_);

        $to = $idx if ($last eq "*");

        if($diso eq "*")
        {
            $RES->{'diso_res'}++;
            $self->addNterm($RES, 'diso') if ($idx <= 50);
            $self->addMidSegment($RES, $idx, 'diso') if (($idx > 50) && ($idx <= ($self->len() - 50)));
            $self->addCterm($RES, 'diso') if ($idx > ($self->len() - 50));
            $from = $idx if ($last ne "*");
        }
        else
        {
            if (($last eq "*") && (($to - $from) >= 5))
            {
                #------add in diso regions by length -----#
                $RES->{'num_diso'}++;
                my $n = exists($CFG->[0]) ? @$CFG : 0;
                $CFG->[$n]{'from'} = $from;
                $CFG->[$n]{'to'}   = $to;

                for (my $i = 0; $i < @F; $i++)
                {
                    if($to-$from < $F[$i] && (!$i || $to-$from >= $F[($i-1)]))
                    {
		        $RES->{'diso'}->{$F[$i]}++;
                    }
                }
            }
        }

        $last = $diso;
    }

    close(DISO);

    if(($last eq "*") && (($to - $from) >= 5))
    {
        #------add in diso regions by length -----#
        my $n = exists($CFG->[0]) ? @$CFG : 0;
        $CFG->[$n]{'from'} = $from;
        $CFG->[$n]{'to'}   = $to;

        for(my $i = 0; $i < @F; $i++)
        {
            if((($to-$from) < $F[$i]) && (!$i || (($to-$from) >= $F[($i-1)])))
            {
                $RES->{'diso'}->{$F[$i]}++;
            }
        }
    }

    $RES->{'diso_res'} /= $self->len();
    $RES->{'diso_nterm'} /= 50;
    $RES->{'diso_cterm'} /= 50;

    for (my $i=1; $i < 9; $i++)
    {
        $RES->{"diso_S$i"} /= $self->seg8();
    }

    $self->{'results'} = $RES;
}



1;
