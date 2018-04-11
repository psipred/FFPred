#!/usr/bin/perl -w

package SigP;

use strict;
use base 'Pred';


sub new
{
    my ($class, $aa, $id, $md5, $cfg) = @_;
    my $name = 'SP';

    my $self = $class->SUPER::new($aa, $id, $md5, $cfg->{'PATH'}, $name);

    $self->{'exe'}     = "$cfg->{'SIGP'}/signalp";
    $self->{'outfile'} = "$self->{'path'}/$self->{'md5'}.sigp";
    $self->{'cmd'}     = "$self->{'exe'}" .
	                 " -t euk -f summary" .
                         " $self->{'path'}/$self->{'md5'}.sig >" .
                         " $self->{'outfile'}";
    bless $self, $class;
    $self->{$self->name()} = {}; # This weird organisation comes from legacy code.
    $self->printSig();

    return $self;
}

sub normalise
{
    my ($self) = @_;

    $self->{'results'}{'maxCpos'} = log(1+$self->{'results'}{'maxCpos'}) / log(71);
    $self->{'results'}{'maxYpos'} = log(1+$self->{'results'}{'maxYpos'}) / log(71);
    $self->{'results'}{'maxSpos'} = log(1+$self->{'results'}{'maxSpos'}) / log(71);
}

sub printSig
{
    my ($self) = @_;

    open(OUT, ">", "$self->{'path'}/$self->{'md5'}.sig");
    print OUT ">$self->{'md5'}\n" . substr($self->{'aa'}, 0, 70) . "\n";
    close(OUT);
}


# Old command to run SignalP v3.0.
sub cmd_SignalP3
{
    my ($self) = @_;

    my $cmd = "$self->{exe}/signalp ".
	      "-t euk ".
              "$self->{path}/$self->{md5}.sig > ".
              "$self->{path}/$self->{md5}.sigp";

    return $cmd;
}


sub parse
{
    my ($self) = @_;

    my $RES = {};
    my $CFG = $self->{$self->name()};

    # See SignalP manual in "$cfg->{'SIGP'}/signalp.1".
    $RES->{'maxCpos'} = 0;
    $RES->{'maxC'}    = 0;
    $RES->{'maxYpos'} = 0;
    $RES->{'maxY'}    = 0;
    $RES->{'maxSpos'} = 0;
    $RES->{'maxS'}    = 0;
    $RES->{'meanS'}   = 0;
    $RES->{'signal'}  = 0.5; # Default value if undefined.

    open(SIGP, "<", $self->outfile()) or die "Cannot open SignalP file... $!\n";
    while (defined(my $line = <SIGP>))
    {
        if ($line =~ /max\.?\s*C\s*(\d+)\s+(\d+(?:\.\d*)?)/)
        {
            $RES->{'maxCpos'} = $1;
            $RES->{'maxC'} = $2;
            $CFG->{'C'} = "$1\t$2";
        }
        elsif ($line =~ /max\.?\s*Y\s*(\d+)\s+(\d+(?:\.\d*)?)/)
        {
            $RES->{'maxYpos'} = $1;
            $RES->{'maxY'} = $2;
            $CFG->{'Y'} = "$1\t$2";
        }
        elsif ($line =~ /max\.?\s*S\s*(\d+)\s+(\d+(?:\.\d*)?)/)
        {
            $RES->{'maxSpos'} = $1;
            $RES->{'maxS'} = $2;
            $CFG->{'S'} = "$1\t$2";
        }
        elsif ($line =~ /mean\s*S\s*\d+-\d+\s+(\d+(?:\.\d*)?)/)
        {
            $RES->{'meanS'} = $1;
            $CFG->{'S'} .= "\t$1" if (exists($CFG->{'S'}));
        }
        elsif ($line =~ /\s*D.+(YES|NO)/)
        {
            $CFG->{'signal'} = $RES->{'signal'} = ($1 eq 'YES') ? 1 : 0;
        }
    }
    close(SIGP);

    $self->{'results'} = $RES;
}



# Old parser for SignalP v3.0 output.
sub parse_SignalP3
{
    my ($self) = @_;

    my $RES={};

    $RES->{maxy}   =0;
    $RES->{maxc}   =0;
    $RES->{anchor} =0;
    $RES->{means}  =0;
    $RES->{siglen} =0;

    $/="\n";

    if (-s "$self->{path}/$self->{md5}.sigp")
    {
	open(SIG, "< $self->{path}/$self->{md5}.sigp");

        while(<SIG>)
        {
            $_ =~ s/^\s*//;

            $RES->{siglen}  = $+ if $_ =~ /^D\s+\d+\-(\d+)\s+\d+/;
	    $RES->{maxy}    = $+ if $_ =~ /^max\.\s+Y\s+\d+\s+(\d+\.\d+)\s+/;
            $RES->{means}   = $+ if $_ =~ /^mean\s+S\s+\d+\-\d+\s+(\d+\.\d+)\s+/;
            $RES->{sigpave} = $+ if $_ =~ /^D\s+\d+\-\d+\s+(\d+\.\d+)\s+/;
            $RES->{anchor}=$+ if $_ =~ /^Signal\s+anchor\s+probability\:\s+(\d+\.\d+)/;

        }
	$self->{SP}->{SIGNAL}->{val}=1 if $RES->{sigpave} >=0.5;
        $self->{SP}->{ANCHOR}->{val}=1 if $RES->{anchor}  >=0.5;

        close(SIG);
    }

    $self->{results} = $RES;
}

# Old command to normalise SignalP v3.0 output.
sub normalise_SignalP3
{
    my ($self) = @_;

    $self->{results}->{siglen} =  log(1+$self->{results}->{siglen})/log(70);

}



1;
