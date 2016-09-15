#!/usr/bin/perl -w

package Coils;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'CC';
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    
    $self->{'exe'} = "$cfg->{'COILS'}/ncoils";
    $self->{'cmd'} = "$self->{'exe'}" . 
                     " -f < $self->{'path'}/$self->{'md5'}.fsa >" . 
                     " $self->{'path'}/$self->{'md5'}.coils";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    $ENV{'COILSDIR'} = $cfg->{'COILS'}; # Set up environment for COILS.
    
    return $self;
}

sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_coils'} = log(1+$self->{'results'}{'num_coils'}) / log(48);
}

sub parse
{
    my ($self) = @_;
    
    my $RES = {};
    my $CFG = $self->{$self->name()};
    
    $RES->{'num_coils'} = 0;
    $RES->{'coil_res'} = 0;
    $RES->{'coil_nterm'} = 0;
    $RES->{'coil_cterm'} = 0;
    
    for(my $i = 1; $i < 9; $i++)
    {
        $RES->{"coil_S$i"} = 0;
    }
    
    local $/ = "\n";
    
    if (-s "$self->{'path'}/$self->{'md5'}.coils")
    {
        open(COIL,"< $self->{'path'}/$self->{'md5'}.coils");
        my $coils = "";
        
        while (defined(my $line = <COIL>))
        {
            chomp $line;
            $coils .= $line unless ($line =~ /^>/);
        }
        
        close(COIL);
        
        if ($coils =~ /x/)
        {
            my $i = 1;
            my $last = "";
            my ($from, $to) = (1, 0);
            my @coils = split(//, $coils);
            
            while (my $coil = shift @coils)
            {
                if ($coil =~ /x/)
                {
                    $RES->{'coil_res'}++;
                    $self->addNterm($RES, 'coil') if ($i <= 50);
                    $self->addMidSegment($RES, $i, 'coil') if (($i > 50) && ($i <= ($self->len() - 50)));
                    $self->addCterm($RES, 'coil') if ($i > ($self->len() - 50));
                    $from = $i if ($last !~ /x/);
                }
                else
                {
                    $to = $i if ($last =~ /x/);
                    if ($to-$from >= 5)
                    {
                        my $n = exists($CFG->[0]) ? @$CFG : 0; 
                        $CFG->[$n]{'from'} = $from;
                        $CFG->[$n]{'to'} = $to;
                        $RES->{'num_coils'}++ if ($to-$from);
                    } 
                }
                
                $last = $i;
                $i++;
            }
            
            if ($to-$from >= 5)
            {
                my $n = exists($CFG->[0]) ? @$CFG : 0; 
                $CFG->[$n]{'from'} = $from;
                $CFG->[$n]{'to'} = $to;
                $RES->{'num_coils'}++ if ($to-$from >= 5);
            }
        }
        
        $RES->{'coil_res'} /= $self->len();
        $RES->{'coil_nterm'} /= 50;
        $RES->{'coil_cterm'} /= 50;
        
        for (my $i=1; $i <9; $i++)
        {
            $RES->{"coil_S$i"} /= $self->seg8();
        }
    }
    
    $self->{'results'} = $RES;
}



1;
