#!/usr/bin/perl -w

package PsiPred;

use strict;
use PsiBlast;
use MakeMat;
use Seq2Mtx;
use base 'BlastPred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'SS';
    my $mask = 'masked';
    my $iter = 3;
    my $Evalue = 1e-3;
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg, $name, $mask, $iter, $Evalue);
    
    $self->{'exe_psipred'} = "$cfg->{'PSIPRED'}/bin/psipred";
    $self->{'exe_psipass'} = "$cfg->{'PSIPRED'}/bin/psipass2";
    $self->{'data'}        = "$cfg->{'PSIPRED'}/data";
    $self->{'cmd'}         = "$self->{'exe_psipred'}" . 
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.mtx" . 
                             " $self->{'data'}/weights.dat" . 
                             " $self->{'data'}/weights.dat2" . 
                             " $self->{'data'}/weights.dat3 >" . 
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.ss";
    $self->{'cmd_psipass'} = "$self->{'exe_psipass'}" . 
                             " $self->{'data'}/weights_p2.dat" . 
                             " 1 1.0 1.0" . 
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.ss2" . 
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.ss >" . 
                             " $self->{'path'}/$self->{'md5'}.$self->{'type'}.horiz";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    
    return $self;
}

sub cmd_psipass
{
    my ($self, $cmd) = @_;
    
    $self->{'cmd_psipass'} = $cmd if (defined($cmd));
    return $self->{'cmd_psipass'};
}

sub addCterm
{
    my ($self,$RES,$ss) = @_;

    if($ss =~ /H/)
    {
        $RES->{helix_cterm}++;
    }
    elsif($ss =~ /E/)
    {
	$RES->{sheet_cterm}++;
    }
    else{

        $RES->{rcoil_cterm}++;
    }
}

sub addMidSegment
{
  my ($self,$RES,$idx,$ss) = @_; 

  return if $self->len() < 108;

  my $seg = int(($idx-50) / $self->seg8())+1;

  if($ss =~ /H/)
  {
      $RES->{"helix_S$seg"}++;
  }
  elsif($ss =~ /E/)
  {
      $RES->{"sheet_S$seg"}++;
  }
}

sub addNterm
{
    my ($self,$RES,$ss) = @_;

    if($ss =~ /H/)
    {
     $RES->{helix_nterm}++;
    }
    elsif($ss =~ /E/)
    {
	$RES->{sheet_nterm}++;
    }
    else{
        $RES->{rcoil_nterm}++;
    }
}

sub addSSRegion
{
    my ($self,$ss,$RES) = @_;

    if ($ss =~ /H/)
    {
	$RES->{num_helix}++;
    }
    elsif($ss =~ /E/)
    {
        $RES->{num_sheet}++;
    }
    else{
        $RES->{num_rcoils}++;
    }
}

sub normalise
{
    my ($self) = @_;

    $self->{results}->{num_helix}  = log(1+$self->{results}->{num_helix})/log(234);
    $self->{results}->{num_sheet}  = log(1+$self->{results}->{num_sheet})/log(276);
    $self->{results}->{num_rcoils} = log(1+$self->{results}->{num_rcoils})/log(300);
 
    $self->{results}->{h}->{10}    = log(1+$self->{results}->{h}->{10})/log(144);
    $self->{results}->{h}->{15}    = log(1+$self->{results}->{h}->{15})/log(74);
    $self->{results}->{h}->{20}    = log(1+$self->{results}->{h}->{20})/log(58);
    $self->{results}->{h}->{30}    = log(1+$self->{results}->{h}->{30})/log(152);
    $self->{results}->{h}->{50}    = log(1+$self->{results}->{h}->{50})/log(70);
    $self->{results}->{h}->{70}    = log(1+$self->{results}->{h}->{70})/log(67);
    $self->{results}->{h}->{100}   = log(1+$self->{results}->{h}->{100})/log(9);
    $self->{results}->{h}->{10000} = log(1+$self->{results}->{h}->{10000})/log(7);

    $self->{results}->{e}->{10}    = log(1+$self->{results}->{e}->{10})/log(355);
    $self->{results}->{e}->{15}    = log(1+$self->{results}->{e}->{15})/log(90);
    $self->{results}->{e}->{20}    = log(1+$self->{results}->{e}->{20})/log(25);
    $self->{results}->{e}->{25}    = log(1+$self->{results}->{e}->{25})/log(6);
    $self->{results}->{e}->{30}    = log(1+$self->{results}->{e}->{30})/log(3);
    $self->{results}->{e}->{40}    = log(1+$self->{results}->{e}->{40})/log(3);
    $self->{results}->{e}->{10000} = log(1+$self->{results}->{e}->{10000})/log(2);
   
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{helix_res}=0;
    $RES->{sheet_res}=0;
    $RES->{rcoil_res}=0;

    $RES->{num_helix} = 0;
    $RES->{num_rcoils} = 0;
    $RES->{num_sheet} = 0;
    $RES->{helix_nterm}  = 0;
    $RES->{helix_cterm}  = 0;
    $RES->{sheet_nterm}  = 0;
    $RES->{sheet_cterm}  = 0;
    $RES->{rcoil_nterm}   = 0;
    $RES->{rcoil_cterm}   = 0;

    for(my $i = 0; $i < 8; $i++)
    {
	$RES->{"helix_S".($i+1)} = 0;
        $RES->{"sheet_S".($i+1)} = 0;
    }

    $RES->{h}->{10}    = 0;
    $RES->{h}->{15}    = 0;
    $RES->{h}->{20}    = 0;
    $RES->{h}->{30}    = 0;
    $RES->{h}->{50}    = 0;
    $RES->{h}->{70}    = 0;
    $RES->{h}->{100}   = 0;
    $RES->{h}->{10000} = 0;     

    $RES->{e}->{10}    = 0;
    $RES->{e}->{15}    = 0;
    $RES->{e}->{20}    = 0;
    $RES->{e}->{25}    = 0;
    $RES->{e}->{30}    = 0;
    $RES->{e}->{40}    = 0;
    $RES->{e}->{10000} = 0;
 
    my $last="";
    my ($from,$to) = (1,0);
    $/="\n";

    if(-s "$self->{path}/$self->{md5}.$self->{'type'}.ss2")
    {
	open(IN, "< $self->{path}/$self->{md5}.$self->{'type'}.ss2");

        while(<IN>)
        {
	    chomp $_;

            $_ =~ s/^\s+//;

            next if $_ !~ /^\s*\d+\s*[A-Z]/;

            my ($idx,$aa,$ss,$h,$c,$e) = split(/\s+/,$_);

            $self->addNterm($RES, $ss) if ($idx <= 50);
            $self->addMidSegment($RES, $idx, $ss) if (($idx > 50) && ($idx <= ($self->len() - 50)));
            $self->addCterm($RES, $ss) if ($idx > ($self->len() - 50));           

            if($ss eq"$last")
	    {
		$to = $idx;
            }
            else{

                my @F = sort {$a <=> $b} keys %{$RES->{lc($last)}};

                if( $to-$from >= 5 )
                {
                  $self->addSSRegion($last,$RES);
                  if($last ne"C")
                  {
                   my $n = exists($CFG->[0]) ? @$CFG : 0;
                   $CFG->[$n]{'from'} = $from;
                   $CFG->[$n]{'to'}   = $to;
                   $CFG->[$n]{'type'} = $last;
                  }
		  for(my $i = 0; $i <@F; $i++)
                  {
                    if($to-$from < $F[$i] && (!$i || $to-$from >= $F[($i-1)]))
                    {
			$RES->{lc($last)}->{$F[$i]}++;
                    }
                  }   
                }
                $from = $idx;
	    }

            $last = $ss;
                 
            if($ss =~ /H/)
            {
             $RES->{helix_res}++;
	    }
            elsif($ss =~ /E/)
            {
             $RES->{sheet_res}++;
	    }
            else{
                  $RES->{rcoil_res}++;
	        }
        }
        close(IN);

        if(length($self->{aa})-$from >=5)
        {
            $self->addSSRegion($last,$RES);
            if( $last ne"C")
            {
             my $n = exists($CFG->[0]) ? @$CFG : 0;
             $CFG->[$n]{'from'} = $from;
             $CFG->[$n]{'to'}   = $to;
             $CFG->[$n]{'type'} = $last;
            }
	    my @F = sort {$a<=>$b} keys %{$RES->{lc($last)}};

            for(my $i = 0; $i <@F; $i++)
            {
		$RES->{lc($last)}->{$F[$i]}++ if $to - $from < $F[$i] && (!$i || $to-$from >= $F[$i]); 
            }
        }

        $RES->{helix_res} /= length($self->{aa});
        $RES->{rcoil_res} /= length($self->{aa});
        $RES->{sheet_res} /= length($self->{aa});

        $RES->{helix_nterm} /= 50;
        $RES->{helix_cterm} /= 50;
        $RES->{sheet_nterm} /= 50;
        $RES->{sheet_cterm} /= 50;
        $RES->{rcoil_nterm}  /= 50;
        $RES->{rcoil_cterm}  /= 50;

        for( my $i = 1; $i < 9; $i++ )
        {
	    $RES->{"helix_S$i"} /= $self->seg8();
            $RES->{"sheet_S$i"} /= $self->seg8();
        }
       
    }
    else
    {
        $self->err(1);
    }

    $self->{'results'} = $RES;
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
    
    print STDERR $self->cmd() . "\n";
    $error += system($self->cmd());
    
    if (-s "$self->{path}/$self->{md5}.$self->{'type'}.ss")
    {
        print STDERR $self->cmd_psipass(), "\n";
        $error += system($self->cmd_psipass());
    }
    else
    {
        print STDERR "Could not find .ss file !\n";
        $error = 1;
    }
    
    $self->err($error);
}



1;
