#!/usr/bin/perl -w

package NetPhos;

use strict;
use base 'Pred';



sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'PH';
    
    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);
    
    $self->{'exe'} = "$cfg->{'NETPHOS'}/ape";
    $self->{'cmd'} = "$self->{'exe'}" . 
                     " -r $self->{'path'}/$self->{'md5'}.fsa >" . 
                     " $self->{'path'}/$self->{'md5'}.netphos";
    
    $self->{$self->name()} = []; # This weird organisation comes from legacy code.
    
    return $self;
}

sub addCterm
{
    my ($self, $RES, $type) = @_;

    if($type eq"S")
   {
    $RES->{cterm_phosS}++;
   }
   elsif($type eq"T")
   {
    $RES->{cterm_phosT}++;
   }
    else
   {
    $RES->{cterm_phosY}++;
   }
}

sub addMidSegment
{
    my ($self, $RES, $idx, $type) = @_;

    return if $self->len() < 108;
 
    my $seg = int(($idx-50) / $self->seg8())+1;

    if($type eq"S")
    {
     $RES->{"phosS_S$seg"}++;
    }
    elsif($type eq"T")
    {
     $RES->{"phosT_S$seg"}++;
    }
    else
    {
     $RES->{"phosY_S$seg"}++;
    }
}

sub addNterm
{
    my ($self, $RES, $type) = @_;

   if($type eq"S")
   {
    $RES->{nterm_phosS}++;
   }
   elsif($type eq"T")
   {
    $RES->{nterm_phosT}++;
   }
    else
   {
    $RES->{nterm_phosY}++;
   }
}

sub addPhosRes
{
    my ($self, $RES, $type) = @_;

    if($type =~ /S/)
    {
	$RES->{num_phosS}++;
    }
    elsif($type =~ /T/)
    {
        $RES->{num_phosT}++;
    }
    else{
        $RES->{num_phosY}++;
    }
}

sub normalise
{
    my ($self) = @_;

    $self->{results}->{num_phosS} = log(1+$self->{results}->{num_phosS})/log(763);
    $self->{results}->{num_phosT} = log(1+$self->{results}->{num_phosT})/log(1057);
    $self->{results}->{num_phosY} = log(1+$self->{results}->{num_phosY})/log(109);
     
    $self->{results}->{num_CKI}    = log(1+$self->{results}->{num_CKI})/log(219);
    $self->{results}->{num_ATM}    = log(1+$self->{results}->{num_ATM})/log(105);
    $self->{results}->{num_CKII}   = log(1+$self->{results}->{num_CKII})/log(173);
    $self->{results}->{num_CaMII}  = log(1+$self->{results}->{num_CaMII})/log(10);
    $self->{results}->{num_DNAPK}  = log(1+$self->{results}->{num_DNAPK})/log(112);
    $self->{results}->{num_EGFR}   = log(1+$self->{results}->{num_EGFR})/log(27); 
    $self->{results}->{num_GSK3}   = log(1+$self->{results}->{num_GSK3})/log(125);
    $self->{results}->{num_INSR}   = log(1+$self->{results}->{num_INSR})/log(67);
    $self->{results}->{num_PKA}    = log(1+$self->{results}->{num_PKA})/log(147);
    $self->{results}->{num_PKB}    = log(1+$self->{results}->{num_PKB})/log(97);
    $self->{results}->{num_PKC}    = log(1+$self->{results}->{num_PKC})/log(428);
    $self->{results}->{num_PKG}    = log(1+$self->{results}->{num_PKG})/log(127);
    $self->{results}->{num_RSK}    = log(1+$self->{results}->{num_RSK})/log(110);
    $self->{results}->{num_SRC}    = log(1+$self->{results}->{num_SRC})/log(40);
    $self->{results}->{num_cdc2}   = log(1+$self->{results}->{num_cdc2})/log(332); 
    $self->{results}->{num_cdk5}   = log(1+$self->{results}->{num_cdk5})/log(606);
    $self->{results}->{num_p38MAPK}= log(1+$self->{results}->{num_p38MAPK})/log(435);
}

sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    $RES->{num_phosS}  =0;
    $RES->{num_phosT}  =0;
    $RES->{num_phosY}  =0;
    $RES->{num_CKI}    =0;
    $RES->{num_ATM}    =0;
    $RES->{num_CKII}   =0;
    $RES->{num_CaMII}  =0;
    $RES->{num_DNAPK}  =0;
    $RES->{num_EGFR}   =0;
    $RES->{num_GSK3}   =0;
    $RES->{num_INSR}   =0;
    $RES->{num_PKA}    =0;
    $RES->{num_PKB}    =0;
    $RES->{num_PKC}    =0;
    $RES->{num_PKG}    =0;
    $RES->{num_RSK}    =0;
    $RES->{num_SRC}    =0;
    $RES->{num_cdc2}   =0;
    $RES->{num_cdk5}   =0;
    $RES->{num_p38MAPK}=0;
    

    for(my $i=1; $i < 9; $i++)
    {
	$RES->{"phosS_S$i"}=0;
        $RES->{"phosT_S$i"}=0;
        $RES->{"phosY_S$i"}=0;
    }

    $RES->{nterm_phosY}=0;
    $RES->{nterm_phosT}=0;
    $RES->{nterm_phosS}=0;
    $RES->{cterm_phosY}=0;
    $RES->{cterm_phosT}=0;
    $RES->{cterm_phosS}=0;

    my $tmp = {};

    if (-s "$self->{path}/$self->{md5}.netphos")
    {
	open(PHOS,"< $self->{path}/$self->{md5}.netphos");

        while(<PHOS>)
        {
	    chomp $_;
            next if $_ !~ /\./;

            my ($dot,$res,$sid,$idx,$s1,$s2,$type) = split(/\s+/,$_);

            if($s1 >= 0.5)
            {
		if(!exists($tmp->{$idx}))
                {
                 my $n = exists($CFG->[0]) ? @$CFG : 0;
                 $CFG->[$n]{'score'} = $s1;
                 $CFG->[$n]{'from'}  = $idx;
                 $CFG->[$n]{'to'}    = $idx;
                 $CFG->[$n]{'type'}  = $res;
	        }

                $type =~ s/\-//g;
		$tmp->{$idx}{'res'} = $res;
                $tmp->{$idx}{'type'}{$type} = 1 if $type !~ /unsp/;
            }

        }

        close(PHOS);

       foreach my $idx(sort {$a<=>$b} keys %$tmp)
       {
    	 $self->addNterm($RES, $tmp->{$idx}{'res'}) if ($idx <= 50);
	 $self->addMidSegment($RES, $idx, $tmp->{$idx}{'res'}) if (($idx > 50) && ($idx <= ($self->len() - 50)));
         $self->addCterm($RES, $tmp->{$idx}{'res'}) if ($idx > ($self->len() - 50));
         $self->addPhosRes($RES, $tmp->{$idx}{'res'});

         foreach my $type (keys %{$tmp->{$idx}{'type'}})
         {
           $RES->{"num_$type"}++;
         }
       }
    }

    for (my $i=1; $i<9; $i++)
    {
	$RES->{"phosS_S$i"} /= $self->seg8();
        $RES->{"phosT_S$i"} /= $self->seg8();
        $RES->{"phosY_S$i"} /= $self->seg8();
    }

    $RES->{nterm_phosS} /= 50;
    $RES->{nterm_phosY} /= 50;
    $RES->{nterm_phosT} /= 50;
    
    $RES->{cterm_phosS} /= 50;
    $RES->{cterm_phosY} /= 50;
    $RES->{cterm_phosT} /= 50;

    $self->{results} = $RES;
}



1;
