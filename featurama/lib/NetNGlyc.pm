#!/usr/bin/perl -w

package NetNGlyc;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'NG';
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    
    $self->{'exe'} = "$cfg->{'NETNGLYC'}/netNglyc";
    $self->{'cmd'} = "$self->{'exe'}" . 
                     " $self->{'path'}/$self->{'md5'}.fsa >" . 
                     " $self->{'path'}/$self->{'md5'}.netNglyc";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    
    return $self;
}

sub addMidSegment
{
    my ($self, $RES, $idx, $name) = @_;
  
    return if ($self->len() < 103);
  
    my $seg = 1 + int(($idx-50) / $self->seg3());
    $RES->{"${name}_S$seg"}++;
}

sub normalise
{
    my ($self) = @_;

    $self->{results}->{num_nglyc} = log(1+$self->{results}->{num_nglyc})/log(43);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{'num_nglyc'}   = 0;
    $RES->{'nglyc_nterm'} = 0;
    $RES->{'nglyc_cterm'} = 0;
    
    for(my $i=1; $i < 4; $i++)
    { 
	$RES->{"nglyc_S$i"} = 0;
    }

    $self->{results} = $RES;

    return if $self->{'aa'} !~ /N/;

    if (-s "$self->{path}/$self->{md5}.netNglyc")
    {
	open(NGLYC, "< $self->{path}/$self->{md5}.netNglyc");
        
        while(<NGLYC>)
        {
	    chomp $_;

            #--if no signal peptide then remove ---#

            if ( $_ =~ /\#\s+name/ )
            {
		$_ =<NGLYC>;

                my ($sid,$cmax,$cpos,$cy,$ypos,$yy,
                         $smax,$spos,$sy,$smean,
                         $smy,$d,$dy) = split(/\s+/,$_);

                return if $d < 0.5;
            }

            next if $_ !~ /^[a-zA-Z0-9]+\s+\d+\s+[A-Z]+\s+\d+\.\d+/;

            my ($sid,$idx,$reg,$score,$jury,$res) = split(/\s+/,$_);

            if( $score >= 0.5 )
            {
		$RES->{num_nglyc}++;

                my $n = exists($CFG->[0]) ? @$CFG : 0;

                $CFG->[$n]{'from'}  = $idx;
                $CFG->[$n]{'to'}    = $idx;
                $CFG->[$n]{'score'} = $score;

                $self->addNterm($RES, 'nglyc') if ($idx <= 50);
                $self->addMidSegment($RES, $idx, 'nglyc') if (($idx > 50) && ($idx <= ($self->len() - 50)));
                $self->addCterm($RES, 'nglyc') if ($idx > ($self->len() - 50));
            }
        }

        close(NGLYC);
    }

    $RES->{'nglyc_nterm'} /= 50;
    $RES->{'nglyc_cterm'} /= 50;
    
    for(my $i=1; $i < 4; $i++)
    {
	$RES->{"nglyc_S$i"} /= $self->seg3();
    }
    
    $self->{results} = $RES;
}

sub run
{
    my ($self, $error) = @_;
    $error = 0 unless ($error);

    $self->SUPER::run($error) if ($self->{'aa'} =~ /N/);
}



1;
