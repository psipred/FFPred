#!/usr/bin/perl -w

package LowComplex;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'LC';
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    
    $self->{'exe'} = "$cfg->{'LOWC'}/pfilt";
    $self->{'cmd'} = "$self->{'exe'} -tc" . 
                     " $self->{'path'}/$self->{'md5'}.fsa >" . 
                     " $self->{'path'}/$self->{'md5'}.lowc";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    
    return $self;
}

sub normalise
{
    my ($self)=@_;

    $self->{'results'}{'num_lowc'}=log(1+$self->{'results'}{'num_lowc'})/log(111);
}

sub parse
{

    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{num_lowc}=0;
    $RES->{lowc_res}=0;
    $RES->{lowc_nterm}=0;
    $RES->{lowc_cterm}=0;
    
    for(my $i=1; $i<9; $i++)
    {
	$RES->{"lowc_S$i"} = 0;
    } 

    
    open(LOW, "< $self->{path}/$self->{md5}.lowc");

    $/="\n";
    my $low = "";

    while(<LOW>)
    {
      chomp $_;

      if ( $_ !~ /^\>/ )
      {
	  $low .= $_;
      }
   
    }

    close(LOW);

    if ($low =~ /X/)
    {
      my $i=1;
      my $last="";
      my ($from,$to) = (1,0);

      my @low = split(//,$low);

      foreach my $l (@low)
      {
          if ($l =~ /X/)
          {
	      $RES->{'lowc_res'}++;
              $self->addNterm($RES, 'lowc') if ($i <= 50);
              $self->addMidSegment($RES, $i, 'lowc') if (($i > 50) && ($i <= ($self->len() - 50)));
              $self->addCterm($RES, 'lowc') if ($i > ($self->len() - 50));
              $from = $i if $last !~ /x/;
          }
          else{
           
              $to = $i if $last =~ /x/;

              if($to-$from >=5)
              {
                  my $n = exists($CFG->[0]) ? @$CFG : 0;
                  $CFG->[$n]->{'from'} = $from;
                  $CFG->[$n]->{'to'}   = $to;
		  $RES->{'num_lowc'}++;
              }
 
	  }
  	   $last = $i;
	   $i++;
      }

      if($to-$from >=5)
      { 
       my $n = exists($CFG->[0]) ? @$CFG : 0;
       $CFG->[$n]->{'from'} = $from;
       $CFG->[$n]->{'to'}   = $to;
        
       $RES->{'num_lowc'}++;
      }
    }

    $RES->{lowc_res} /= $self->len();

    $RES->{lowc_nterm} /= 50;
    $RES->{lowc_cterm} /= 50;

    for(my $i=1; $i < 9; $i++)
    {
	$RES->{"lowc_S$i"} /= $self->seg8();
    }

    $self->{results} = $RES;
  
}



1;
