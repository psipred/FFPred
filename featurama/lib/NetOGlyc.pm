#!/usr/bin/perl -w

package NetOGlyc;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'OG';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
 
    $self->{'exe'} = "$cfg->{'NETOGLYC'}/netOglyc";
    $self->{'cmd'} = "$self->{exe}" . 
	             " $self->{path}/$self->{md5}.fsa >" . 
                     " $self->{path}/$self->{md5}.netOglyc";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    
    return $self;
}

sub normalise
{
    my ($self) = @_;

    $self->{results}->{num_oglyc} = log(1+$self->{results}->{num_oglyc})/log(1704);
    
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{'num_oglyc'}   = 0;
    $RES->{'oglyc_nterm'} = 0;
    $RES->{'oglyc_cterm'} = 0;
    
    for(my $i = 1; $i < 9; $i++)
    {
	$RES->{"oglyc_S$i"} = 0;
    } 
 
    $/="\n";

   if ( -s "$self->{path}/$self->{md5}.netOglyc")
   {
     open(OGLYC, "< $self->{path}/$self->{md5}.netOglyc");

     while(<OGLYC>)
     {
	 chomp $_;

         next if $_ !~ /^[a-z0-9A-Z]+\s+[ST]\s+\d+/;

         my ($sid,$res,$idx,$s1,$s2,$yes) = split(/\s+/,$_);

        if($yes ne'.')
        {
	  my $s = $s1 > $s2 ? $s1 : $s2;
          my $n = exists($CFG->[0]) ? @$CFG : 0;
          
          $CFG->[$n]{'from'}  = $idx;
          $CFG->[$n]{'to'}    = $idx;
          $CFG->[$n]{'score'} = $s;
          
          $RES->{'num_oglyc'}++;
          $self->addNterm($RES, 'oglyc') if ($idx <= 50);
          $self->addMidSegment($RES, $idx, 'oglyc') if (($idx > 50) && ($idx <= ($self->len() - 50)));
          $self->addCterm($RES, 'oglyc') if ($idx > ($self->len() - 50));
        } 
     }

     close(OGLYC);
   }
    
    $RES->{oglyc_nterm} /= 50;
    $RES->{oglyc_cterm} /= 50;
    
    for(my $i=1; $i < 9; $i++)
    {
     $RES->{"oglyc_S$i"} /= $self->seg8();
    }

    $self->{results} = $RES;
}



1;
