#!/usr/bin/perl -w

package Memsat;

use strict;
use File::Copy;
use PsiBlast;
use MakeMat;
use Seq2Mtx;
use base 'BlastPred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'TM';
    my $mask = 'unmasked';
    my $iter = 2;
    my $Evalue = 1e-6;
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg, $name, $mask, $iter, $Evalue);
    
    $self->{'folder'}   = $cfg->{'MEMSAT'};
    $self->{'exe'}      = "$self->{'folder'}/run_memsat-svm.pl";
    $self->{'data'}     = "$self->{'folder'}/data";
    $self->{'cmd'}      = "$self->{'exe'}" . 
                          " -mtx 1 -e 1 -j $self->{path}/ -o $self->{md5}.$self->{'type'} -w $self->{'folder'}/" . 
                          " $self->{path}/$self->{md5}.$self->{'type'}.mtx" . 
                          " > $self->{path}/$self->{md5}.memsatSVMlog";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    
    return $self;
}

sub normalise
{
    my ($self) = @_;
    
    $self->{results}->{num_tm} = log(1+$self->{results}->{num_tm})/log(21);
}

sub parse
{
    my ($self) = @_;
    
    my $RES = {};
    my $CFG = $self->{$self->name()};
    
    $RES->{'in/out'}   = 0;
    $RES->{'tm_nterm'} = 0;
    $RES->{'tm_cterm'} = 0;
    $RES->{'num_tm'}   = 0;
    $RES->{'tm_res'}   = 0;    
    
    for(my $i=1; $i < 9; $i++)
    {
        $RES->{"tm_S$i"}=0;
    }
    
    my $hTMData = {};
    open(MEM, "<", "$self->{path}/$self->{md5}.$self->{'type'}.memsat_svm") or die "Cannot open memsat_svm file... $!\n";
    
    my $helix_number = 0;
    my $direction  = '';
    my $best_score = 0;
    while (defined(my $line = <MEM>))
    {
        chomp $line;
        if ($line =~ /Processing\s+(\d+)\s+heli/)
        {
            $helix_number = $1;
            my $first_line = <MEM>;
            if($first_line =~ /Transmembrane\shelix\s1\s+from\s+(\d+)\s+\((.+?)\)\s+to\s+(\d+)\s+.+:\s+(.+)/)
            {
                my $start = $1;
                my $dir = $2;
                my $end = $3;
                my $score = $4;
                if($dir =~ /in/)
                {
                    $direction = "+";
                }
                else
                {
                    $direction = "-";
                }
                $hTMData->{$helix_number}{$direction}{1}{START} = $start;
                $hTMData->{$helix_number}{$direction}{1}{END} = $end;
                $hTMData->{$helix_number}{$direction}{1}{SCORE} = $score;
            }
        }
        
        if($line =~ /Transmembrane\shelix\s(\d+)\s+from\s+(\d+)\s+\((.+?)\)\s+to\s+(\d+)\s+.+:\s+(.+)/)
        {
            my $helix_count = $1;
            my $start = $2;
            my $dir = $3;
            my $end = $4;
            my $score = $5;
            $hTMData->{$helix_number}{$direction}{$helix_count}{START} = $start;
            $hTMData->{$helix_number}{$direction}{$helix_count}{END} = $end;
            $hTMData->{$helix_number}{$direction}{$helix_count}{SCORE} = $score;
        }
        
        if($line =~ /Score\s=\s(.+)/)
        {
            $hTMData->{$helix_number}{$direction}{SET_SCORE} = $1;
        }
        
        if($line =~ /Score:\s+(.+)/)
        {
            $best_score = $1;
        }
    }
    
    my $highest_score= -1000000;
    my $highest_set = 0;
    my $highest_direction = '';
    foreach my $helices (keys %$hTMData)
    {
        foreach my $direction (keys %{$hTMData->{$helices}})
        {
            if($hTMData->{$helices}{$direction}{SET_SCORE} > $highest_score)
            {
                $highest_score = $hTMData->{$helices}{$direction}{SET_SCORE};
                $highest_set = $helices;
                $highest_direction = $direction;
            }
        }
    }
    
    delete $hTMData->{$highest_set}{$highest_direction}{SET_SCORE};
    foreach my $helix (sort {$a <=> $b} keys %{$hTMData->{$highest_set}{$highest_direction}})
    {
        next if $helix =~ /SET_SCORE/;
        my $start =$hTMData->{$highest_set}{$highest_direction}{$helix}{START};
        my $end = $hTMData->{$highest_set}{$highest_direction}{$helix}{END};
        my $score = $hTMData->{$highest_set}{$highest_direction}{$helix}{SCORE};
        
        if ($score != 0 && $helix == 1)
        {
            $RES->{'in/out'} = 1 if $highest_direction eq "-";
            $RES->{'in/out'} = 0.5 if $highest_direction eq "+";
        }
        
        if ($score != 0)
        {
            $RES->{'num_tm'}++;
            my $n = $helix-1;
            $CFG->[$n]{'from'} = $start;
            $CFG->[$n]{'to'}   = $end;         
            
            for (my $i = $start; $i <= $end; $i++)
            {
                $RES->{'tm_res'}++;
                $self->addNterm($RES, 'tm') if ($i <= 50);
                $self->addMidSegment($RES, $i, 'tm') if (($i > 50) && ($i <= ($self->len() - 50)));
                $self->addCterm($RES, 'tm') if ($i > ($self->len() - 50));
            }
        }
    }
    
    close(MEM);
    
    $RES->{tm_res}  /= $self->len();
    $RES->{tm_nterm}/= 50;
    $RES->{tm_cterm}/= 50;
    
    for (my $i=1; $i <9; $i++)
    {
       	$RES->{"tm_S$i"} = $RES->{"tm_S$i"} > $self->seg8() ? 1 : $RES->{"tm_S$i"}/$self->seg8();
    }
    
    $self->{results} = $RES;
}



# Old subroutine used when running memsat3.
sub globCMD
{
    my ($self) = @_;

    my $cmd="$self->{'folder'}/bin/globmem ".
            "$self->{data}/glob_weights.dat ".
            "$self->{path}/$self->{md5}.unmasked2.mtx > ".
            "$self->{path}/$self->{md5}.memsat.globmem";

    return $cmd;
}

# Old subroutine used when running memsat3.
sub memsatCMD
{
    my ($self) = @_;

    my $cmd="$self->{'folder'}/bin/mem_pred ".
	    "$self->{data}/weights.dat ".
            "$self->{path}/$self->{md5}.unmasked2.mtx > ".
            "$self->{path}/$self->{md5}.memsat.nn";

    return $cmd;
}

# Old subroutine used when running memsat3.
sub nnCMD
{
    my ($self) = @_;

    my $cmd="$self->{'folder'}/bin/nnsat ".
            "$self->{path}/$self->{md5}.memsat.nn > ".
            "$self->{path}/$self->{md5}.memsat";

    return $cmd;
}

# Old "runMemsat" subroutine to run memsat3.
sub runMemsat3
{
    my ($self) = @_;
    
    if (-s "$self->{path}/$self->{md5}.unmasked2.chk")
    {
        $self->{'makeMat'}->init();
        $self->{'makeMat'}->run();
    }
    elsif (-s "$self->{path}/$self->{md5}.unmasked3.mtx") # Use DisoPred's files instead.
    {
        copy("$self->{path}/$self->{md5}.unmasked3.mtx", "$self->{path}/$self->{md5}.unmasked2.mtx");   
    }
    
    my $e = system($self->globCMD());
       $e = system($self->memsatCMD()) if -s "$self->{path}/$self->{md5}.memsat.globmem";
       $e = system($self->nnCMD()) if -s "$self->{path}/$self->{md5}.memsat.nn";
}

# Old "parse" subroutine to use with memsat3.
sub parse_memsat3
{
    my ($self) = @_;

    my $RES = {};

    $RES->{'in/out'}   = 0;
    $RES->{'tm_nterm'} = 0;
    $RES->{'tm_cterm'} = 0;
    $RES->{'num_tm'}   = 0;
    $RES->{'tm_res'}   = 0;    

    for(my $i=1; $i < 9; $i++)
    {
	$RES->{"tm_S$i"}=0;
    }

    open(MEM, "< $self->{path}/$self->{md5}.memsat") or die "can't open $!\n";

    while(<MEM>)
    {

      if( $_ =~ /^1\:\s+\((.*?)\)\s+(\d+)\-(\d+)\s+\((\d+\.\d+)\)/)
      {
       
          if( $4 > 0 )
          {
           $self->{TM}->[0]->{from}=$2;
           $self->{TM}->[0]->{to}  =$3;

	   $RES->{'in/out'} = 1 if $+ eq "out";
           $RES->{'in/out'} = 0.5 if $+ eq "in";
          }
      }
      
      if( $_ =~ /^(\d+)\:\s+(\d+)\-(\d+)\s+\((\d+\.\d+)\)/ )
      {

         if( $4 > 0 )
         { 
          $RES->{'num_tm'}++;
          my $n = $1-1;
          $self->{TM}->[$n]->{from}=$2;
          $self->{TM}->[$n]->{to}  =$3;         

          for (my $i = $2; $i <= $3; $i++)
          {
	      $RES->{'tm_res'}++;
              $self->addNterm($RES, 'tm') if ($i <= 50);
              $self->addMidSegment($RES, $i, 'tm') if (($i > 50) && ($i <= ($self->len() - 50)));
              $self->addCterm($RES, 'tm') if ($i > ($self->len() - 50));
 
          }

         }
      }

                 
    }
   
    close(MEM);


    $RES->{tm_res}  /= $self->len();
    $RES->{tm_nterm}/= 50;
    $RES->{tm_cterm}/= 50;
    
    for(my $i=1; $i <9; $i++)
    {
       	$RES->{"tm_S$i"} = $RES->{"tm_S$i"} > $self->seg8() ? 1 : $RES->{"tm_S$i"}/$self->seg8();
        
    }

    $self->{results} = $RES;
}



1;
