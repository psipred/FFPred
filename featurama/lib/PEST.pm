#!/usr/bin/perl -w

package PEST;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'PE';
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    
    $self->{'exe'} = "$cfg->{'PEST'}/epestfind";
    $self->{'cmd'} = "$self->{exe}" . 
                     " -graph NO" . 
                     " -window 10" . 
                     " -order 2" . 
                     " -threshold 5.0" . 
                     " $self->{path}/$self->{md5}.fsa" . 
                     " -outfile $self->{path}/$self->{md5}.pest";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    
    return $self;
}

sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_pest'} = log(1+$self->{'results'}{'num_pest'}) / log(123);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{'num_pest'}   = 0;
    $RES->{'pest_res'}   = 0;
    $RES->{'pest_nterm'} = 0;
    $RES->{'pest_cterm'} = 0;
    
    for(my $i=1; $i<9; $i++)
    {
	$RES->{"pest_S$i"} = 0;
    }

    $/="\n";

    if(-s "$self->{path}/$self->{md5}.pest")
    {
	open(PST, "< $self->{path}/$self->{md5}.pest");

        while(<PST>)
        {
            $_ =~ s/^\s*//;

            $RES->{num_pest} = $+ if $_ =~ /^(\d+)\s+PEST\s+motifs*\s+were\s+identified/;

            if( $_ =~ /PEST\s+motif\s+with\s+\d+\s+amino\s+acids\s+between\s+position\s+(\d+)\s+and\s+(\d+)/ )
            {
              my ($from,$to) = ($1,$2);

              my $n = exists($CFG->[0]) ? @$CFG : 0;
              $CFG->[$n]{'from'} =$from;
              $CFG->[$n]{'to'}   =$to;

              for(my $idx=$from; $idx < $to; $idx++)
              {
	       $RES->{'pest_res'}++;
               $self->addNterm($RES, 'pest') if ($idx <= 50);
               $self->addMidSegment($RES, $idx, 'pest') if (($idx > 50) && ($idx <= ($self->len() - 50)));
               $self->addCterm($RES, 'pest') if ($idx > ($self->len() - 50));
 	      }

	   }
 
        }

        close(PST);
    }
    
    $RES->{pest_nterm} /= 50;
    $RES->{pest_cterm} /= 50;
    
    for(my $i=1; $i<9; $i++)
    {
	$RES->{"pest_S$i"} /= $self->seg8();
    }

    $RES->{pest_res} /= $self->len();

    $self->{results} = $RES;
}



1;
