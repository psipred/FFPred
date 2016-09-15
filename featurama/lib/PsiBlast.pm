#!/usr/bin/perl -w

package PsiBlast;

use strict;
use base 'Tool';



sub new
{
    my ($class, $iter, $mask, $exe_path, $in, $db, $Eval, $multi_Eval) = @_;
    
    my $self = {
                 'iter'       => $iter,
                 'mask'       => $mask,
                 'exe'        => "$exe_path/blastpgp",
                 'infile'     => $in,
                 'db'         => $db,
	         'Eval'       => $Eval,
                 'multi_Eval' => $multi_Eval,
                 'Y'          => 9000000000, # Old value from legacy code.
                 'err'        => 0
               };

    bless $self, $class;
    
    my $prefix = $in;
    $prefix =~ s/\.fsa$//;
    
    $self->{'type'}    = $self->return_type($self->{'mask'}, $self->{'iter'});
    $self->{'outfile'} = "$prefix.$self->{'type'}.psi";
    $self->{'chk'}     = "$prefix.$self->{'type'}.chk";
    $self->{'cmd'}     = "$self->{'exe'}" . 
                         " -d $self->{'db'}" . 
                         " -i $self->{'infile'}" . 
                         " -o $self->{'outfile'}" . 
                         " -e $self->{'Eval'}" . 
                         " -h $self->{'multi_Eval'}" . 
                         " -j $self->{'iter'}" . 
#                         " -Y $self->{Y}" . 
                         " -a 1" . 
                         " -b 0" . 
                         " -C $self->{chk}";
#    $self->{cmd} .= " -F T" if ($self->{'mask'} eq 'masked');
    
    return $self;
}

sub check
{
    my ($self) = @_;
    
    if (-s $self->outfile())
    {
        print STDERR "Checking psiblast output: ";
        
	local $/;
	open(BLAST, "<", $self->outfile());
        my $text = <BLAST>;
        close(BLAST);
        
        my @runs = ($text =~ /Searching\.+done/g);
        my $num_runs = scalar @runs;
        print STDERR "returning $num_runs\n";
        
        my $OK_condition = ( ($num_runs == $self->iter()) || 
                             ($text =~ /CONVERGED!/) || 
                             ($text =~ /No hits found/) );
        return 1 if ($OK_condition);
    }
    
    return 0;
}

sub chk
{
    my ($self, $chk) = @_;
    
    $self->{'chk'} = $chk if (defined($chk));
    return $self->{'chk'};
}

sub infile
{
    my ($self, $infile) = @_;
    
    $self->{'infile'} = $infile if (defined($infile));
    return $self->{'infile'};
}

sub iter
{
    my ($self, $iter) = @_;
    
    $self->{'iter'} = $iter if (defined($iter));
    return $self->{'iter'};
}

sub outfile
{
    my ($self, $outfile) = @_;
    
    $self->{'outfile'} = $outfile if (defined($outfile));
    return $self->{'outfile'};
}

sub run
{
    my ($self, $error) = @_;
    $error = 0 unless ($error);
    
    if (-f $self->infile())
    {
        my $run_condition = 1;
        my $run_attempt = 0;
        while ($run_condition && ($run_attempt < 10))
        {
            $run_attempt++;
            unlink($self->outfile(), $self->chk());
            
            print STDERR "Running PsiBlast\n" . $self->cmd() . "\n";
            my $blast_error = system($self->cmd());
            unlink('error.log');
            
            $blast_error = 999 if (!defined($blast_error) || ($blast_error eq ""));
            print STDERR "blast error state: $blast_error\n";
            
            $self->err($error + $blast_error);
            
            next if ($blast_error);
            $run_condition = 0 if ($self->check());
        }
        print STDERR "Checkpoint - ALREADY CHECKED down to here\n\n";
    }
    else
    {
        $self->err(1);
    }
}



1;
