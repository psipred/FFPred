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
    my ($class, $aa, $id, $md5, $cfg, $web_control) = @_;
    my $name = 'TM';
    my $mask = 'unmasked';
    my $iter = 2;
    my $evalue = 1e-6;
    my $hvalue = 1e-6;

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg, $name, $mask, $iter, $evalue, $hvalue);

    $self->{'folder'}  = $cfg->{'MEMSAT'};
    $self->{'exe'}     = "$self->{'folder'}/run_memsat-svm.pl";
    $self->{'data'}    = "$self->{'folder'}/data";
    $self->{'globmem'} = "$self->{'path'}/$self->{'md5'}.globmem_svm"; # Manually synchronise this with cmd below!
    if($web_control)
    {
      $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.memsat_svm";
    }
    else{
      $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.$self->{'type'}.memsat_svm"; # Manually synchronise this with cmd below!
    }
    $self->{'cmd'}      = "$self->{'exe'}" .
                          " -mtx 1 -e 1 -j $self->{path}/ -o $self->{md5}.$self->{'type'} -w $self->{'folder'}/" .
                          " $self->{path}/$self->{md5}.$self->{'type'}.mtx" .
                          " > $self->{path}/$self->{md5}.memsatSVMlog";

    $self->{'web_control'} = $web_control;
    $self->{'numSegments'} = 8;
    $self->{'segments'}    = $self->createSegments($self->{'numSegments'});
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.

    return $self;
}

sub globmem
{
    my ($self, $globmem) = @_;

    $self->{'globmem'} = $globmem if (defined($globmem));
    return $self->{'globmem'};
}

sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'num_tm_hel'} = log(1+$self->{'results'}{'num_tm_hel'}) / log(3);
    $self->{'results'}{'num_re_hel'} = log(1+$self->{'results'}{'num_re_hel'}) / log(30);
    $self->{'results'}{'num_tm_res'} = log(1+$self->{'results'}{'num_tm_res'}) / log(600);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    my $min_TM_residues = 7; # Bypass default MEMSAT-SVM1.3 threshold for labelling a protein as transmembrane (= 1) with a more stringent criterion.

    $RES->{'glob/tm'}    = 0.5; # Default value if undefined.
    $RES->{'num_tm_hel'} = 0;
    $RES->{'in/out'}     = 0.5; # Default value if undefined.
    $RES->{'num_re_hel'} = 0;
    $RES->{'num_tm_res'} = 0;
    $RES->{'tm_nterm'}   = 0;
    $RES->{'tm_cterm'}   = 0;

    for (my $i = 1; $i <= $self->{'numSegments'}; $i++)
    {
        $RES->{"tm_S$i"} = 0;
    }

    local $/;
    open(GLOBMEM, "<", $self->globmem()) or die "Cannot open globmem_svm file... $!\n";
    my $globmem_out = <GLOBMEM>;
    close(GLOBMEM);

    my $num_TM_residues = $1 if ($globmem_out =~ /Transmembrane residues found:\s*(\d+)/);

    if (defined($num_TM_residues) && ($num_TM_residues < $min_TM_residues))
    {
        $RES->{'glob/tm'} = 0;
    }
    elsif (defined($num_TM_residues))
    {
        $RES->{'glob/tm'} = 1;

        local $/ = "Results:";
        open(MEMSAT, "<", $self->outfile()) or die "Cannot open memsat_svm file... $!\n";
        my @memsat_out = <MEMSAT>;
        close(MEMSAT);

        my $memsat_result = pop @memsat_out;

        if ($memsat_result =~ /^Topology:\s+(.*)/m)
        {
            my $helices = $1;
            my $n = 0;
            foreach my $helix (split(',', $helices))
            {
                if ($helix =~ /(\d+)-(\d+)/)
                {
                    my ($first, $last) = ($1, $2);

                    $CFG->[$n]{'type'} = 'tm_helix';
                    $CFG->[$n]{'from'} = $first;
                    $CFG->[$n]{'to'}   = $last;

                    foreach my $residue ($first..$last)
                    {
                        $RES->{'num_tm_res'}++;
                        $self->addResidue($RES, $residue, 'tm');
                    }

                    $RES->{'num_tm_hel'}++;
                    $n++;
                }
            }
        }

        if ($memsat_result =~ /^Re-entrant helices:\s+(.*)/m)
        {
            my $reentrant = $1;
            $RES->{'num_re_hel'}++ while ($reentrant =~ /-/g);
        }

        $RES->{'in/out'} = 0 if ($memsat_result =~ /^N-terminal:\s+in/m);
        $RES->{'in/out'} = 1 if ($memsat_result =~ /^N-terminal:\s+out/m);
    }

    $self->{'results'} = $RES;
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
